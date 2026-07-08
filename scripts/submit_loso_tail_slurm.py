#!/usr/bin/env python3
"""
Submit the LOSO-O TAIL experiment (small-serogroup batches) as SLURM jobs.

Same machinery as submit_loso_slurm.py, pre-pointed at the batched-singleton
folds (folds_loso_tail.csv, fold labels batch_1..6) built by
build_loso_singleton_batches.py, writing to loso_tail_results. Nothing to edit --
just build the folds then run this. (Original submit_loso_slurm.py stays pointed
at the completed 19-serogroup run; keep them separate so jobs/logs don't collide.)

Decomposes the shared-clustering nested-cv run into:
  1. ONE shared clustering job (genophi.workflows.run_shared_clustering) over all
     genomes -> shared MMseqs2 DBs + presence/absence matrices.
  2. A job ARRAY, one task per fold (genophi.workflows.run_predefined_fold), each
     reading its fold from the predefined folds CSV and reusing the shared
     artifacts. Depends afterok on (1).
  3. ONE aggregation job (genophi.workflows.aggregate_predefined_folds) pooling
     all completed folds -> global + per-fold metrics, PR/ROC curves. Depends
     afterok on the whole array.

This script is INTENTIONALLY outside the genophi package -- the package stays
scheduler-agnostic and exposes the three functions; this driver wires them into
SLURM. Edit the CONFIG block, then run:  python submit_loso_slurm.py
"""
import os
import subprocess
import pandas as pd

# ============================ CONFIG (edit me) ============================
# Inputs (paths as they exist on the Lawrencium compute nodes)
input_strain_dir = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/strain_AAs_update"
input_phage_dir = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/phage_AAs"
interaction_matrix = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/Gaborieau_interaction_matrix_long_mod.csv"
folds_file = "/global/scratch/users/anoonan/BRaVE/LOSO/tables/folds_loso_tail.csv"

# Output root for this experiment (small-serogroup TAIL batches, batch_1..6)
base_output_dir = "/global/scratch/users/anoonan/BRaVE/LOSO/loso_tail_results"

# Modeling params (passed to each fold)
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"
# Class weights + feature filtering (explicit; not the function defaults)
use_dynamic_weights = True
weights_method = "inverse_frequency"
use_clustering = True
cluster_method = "hierarchical"
max_ram = 100

# SLURM (matches the lab template)
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"
threads = 8
mem_per_job = "120G"
time_limit = "24:00:00"
# Clustering job can be lighter/faster, but keep same for simplicity:
cluster_time_limit = "24:00:00"
# =========================================================================

shared_dir = os.path.join(base_output_dir, "shared_clustering")
folds_out_dir = os.path.join(base_output_dir, "folds")   # one subdir per fold label
logs_dir = os.path.join(base_output_dir, "slurm_logs")


def header(f, job, dep=None, array=None, time_lim=time_limit):
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name={job}\n")
    f.write(f"#SBATCH --account={account}\n#SBATCH --partition={partition}\n#SBATCH --qos={qos}\n")
    f.write("#SBATCH --nodes=1\n#SBATCH --ntasks=1\n")
    f.write(f"#SBATCH --cpus-per-task={threads}\n#SBATCH --mem={mem_per_job}\n#SBATCH --time={time_lim}\n")
    if array:
        f.write(f"#SBATCH --array={array}\n")
        f.write(f"#SBATCH --output={logs_dir}/{job}_%A_%a.out\n#SBATCH --error={logs_dir}/{job}_%A_%a.err\n")
    else:
        f.write(f"#SBATCH --output={logs_dir}/{job}_%j.out\n#SBATCH --error={logs_dir}/{job}_%j.err\n")
    if dep:
        f.write(f"#SBATCH --dependency=afterok:{dep}\n")
    f.write("\nmodule load miniforge3\n")
    f.write('eval "$(conda shell.bash hook)"\n')
    f.write(f"conda activate {environment}\n\n")


def submit(sh_path):
    return subprocess.check_output(["sbatch", "--parsable", sh_path]).decode().strip()


def write_cluster_job():
    """Write + return the path to the shared-clustering job script."""
    w = os.path.join(logs_dir, "LOSOTAIL_cluster.sh")
    with open(w, 'w') as f:
        header(f, "LOSOTAIL_cluster", time_lim=cluster_time_limit)
        f.write(
            "python -c \""
            "from genophi.workflows import run_shared_clustering; "
            f"run_shared_clustering("
            f"input_strain_dir='{input_strain_dir}', "
            f"input_phage_dir='{input_phage_dir}', "
            f"output_dir='{shared_dir}', threads={threads}, "
            f"strain_column='{strain_column}')\"\n"
        )
    return w


def write_single_fold_job(label, dep):
    """Write + return a job script that runs ONE fold by label (afterok dep)."""
    w = os.path.join(logs_dir, f"LOSOTAIL_fold_{label}.sh")
    with open(w, 'w') as f:
        header(f, f"LOSOTAIL_fold_{label}", dep=dep)
        f.write(
            "python -c \""
            "import sys; from genophi.workflows import run_predefined_fold; "
            f"run_predefined_fold("
            "fold_label=sys.argv[1], "
            f"folds_file='{folds_file}', "
            f"shared_dir='{shared_dir}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"input_strain_dir='{input_strain_dir}', "
            f"fold_output_dir='{folds_out_dir}/' + sys.argv[1], "
            f"strain_column='{strain_column}', "
            f"num_runs_fs={num_runs_fs}, num_runs_modeling={num_runs_modeling}, "
            f"use_dynamic_weights={use_dynamic_weights}, "
            f"weights_method='{weights_method}', "
            f"use_clustering={use_clustering}, "
            f"cluster_method='{cluster_method}', "
            f"max_ram={max_ram}, "
            f"threads={threads})\" "
            f"\"{label}\"\n"
        )
    return w


def write_fold_array_job(labels_file, n_folds, dep):
    """Write + return the array job script (one task per fold)."""
    w = os.path.join(logs_dir, "LOSOTAIL_fold.sh")
    with open(w, 'w') as f:
        header(f, "LOSOTAIL_fold", dep=dep, array=f"1-{n_folds}")
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write('echo "Running fold: $LABEL"\n\n')
        f.write(
            "python -c \""
            "import sys; from genophi.workflows import run_predefined_fold; "
            f"run_predefined_fold("
            "fold_label=sys.argv[1], "
            f"folds_file='{folds_file}', "
            f"shared_dir='{shared_dir}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"input_strain_dir='{input_strain_dir}', "
            f"fold_output_dir='{folds_out_dir}/' + sys.argv[1], "
            f"strain_column='{strain_column}', "
            f"num_runs_fs={num_runs_fs}, num_runs_modeling={num_runs_modeling}, "
            f"use_dynamic_weights={use_dynamic_weights}, "
            f"weights_method='{weights_method}', "
            f"use_clustering={use_clustering}, "
            f"cluster_method='{cluster_method}', "
            f"max_ram={max_ram}, "
            f"threads={threads})\" \"$LABEL\"\n"
        )
    return w


def write_aggregate_job(dep):
    """Write + return the aggregation job script (afterok the fold array)."""
    w = os.path.join(logs_dir, "LOSOTAIL_aggregate.sh")
    with open(w, 'w') as f:
        header(f, "LOSOTAIL_aggregate", dep=dep, time_lim="2:00:00")
        f.write(
            "python -c \""
            "from genophi.workflows import aggregate_predefined_folds; "
            f"aggregate_predefined_folds("
            f"folds_dir='{folds_out_dir}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"output_dir='{base_output_dir}', "
            f"strain_column='{strain_column}')\"\n"
        )
    return w


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Submit the LOSO nested-CV experiment to SLURM.")
    ap.add_argument('--test-fold', metavar='LABEL',
                    help="Submit ONLY the clustering job + one fold (by label, e.g. O6) "
                         "as batch jobs. Use to validate the pipeline end-to-end before "
                         "the full array. No array, no aggregation.")
    args = ap.parse_args()

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)

    folds = pd.read_csv(folds_file)
    fold_labels = list(dict.fromkeys(folds['fold_label'].astype(str)))
    n_folds = len(fold_labels)

    # Shared clustering is always step 1 (reused/resumed if already done).
    cl_id = submit(write_cluster_job())
    print(f"  clustering job: {cl_id}")

    if args.test_fold:
        if args.test_fold not in fold_labels:
            raise SystemExit(f"--test-fold '{args.test_fold}' not in folds file. "
                             f"Available: {fold_labels}")
        fold_id = submit(write_single_fold_job(args.test_fold, dep=cl_id))
        print(f"  test fold '{args.test_fold}': {fold_id} (afterok {cl_id})")
        print(f"\n=== Test submitted ===")
        print(f"  cluster {cl_id} -> fold {args.test_fold} {fold_id}")
        print(f"  Result: {folds_out_dir}/{args.test_fold}/model_validation/predict_results/")
        print(f"  Monitor: squeue -u $USER ; "
              f"tail -f {logs_dir}/LOSO_fold_{args.test_fold}_*.out")
        return

    # Full run: array over all folds, then aggregate.
    labels_file = os.path.join(base_output_dir, "fold_labels.txt")
    with open(labels_file, 'w') as fh:
        fh.write("\n".join(fold_labels) + "\n")
    print(f"{n_folds} folds: {fold_labels}")

    arr_id = submit(write_fold_array_job(labels_file, n_folds, dep=cl_id))
    print(f"  fold array: {arr_id} (1-{n_folds})")

    agg_id = submit(write_aggregate_job(dep=arr_id))
    print(f"  aggregate job: {agg_id}")

    print(f"\n=== Submitted ===")
    print(f"  cluster {cl_id} -> fold array {arr_id} (1-{n_folds}) -> aggregate {agg_id}")
    print(f"  Results: {base_output_dir}/performance/  (global + per-fold metrics, curves)")
    print(f"  Per-fold: {folds_out_dir}/<O-group>/")
    print(f"  Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
