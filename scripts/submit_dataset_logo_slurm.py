#!/usr/bin/env python3
"""
Submit a leave-one-dataset-out (genus) nested-CV experiment as SLURM jobs.

Same 3-step pipeline as submit_loso_slurm.py (cluster once -> fold array ->
aggregate), but folds come from a predefined dataset-folds CSV (one fold per
genus: validation = that genus's strains, modeling = all others). Build the
folds file first with build_dataset_folds.py.

  python submit_dataset_logo_slurm.py            # full run (4-fold array)
  python submit_dataset_logo_slurm.py --test-fold "Vibrio"   # cluster + one fold
"""
import os
import subprocess
import pandas as pd

# ============================ CONFIG (edit me) ============================
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"
interaction_matrix = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
folds_file = "/global/home/groups/pc_phiml/embeddings/combined/folds_dataset.csv"

base_output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_logo"

# Modeling params
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"
# Class weights + feature filtering (explicit; not the function defaults)
use_dynamic_weights = True
weights_method = "inverse_frequency"
use_clustering = True
cluster_method = "hierarchical"
max_ram = 100

# SLURM
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"
threads = 32
mem_per_job = "240G"
time_limit = "24:00:00"
cluster_time_limit = "24:00:00"
# =========================================================================

shared_dir = os.path.join(base_output_dir, "shared_clustering")
folds_out_dir = os.path.join(base_output_dir, "folds")
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
    w = os.path.join(logs_dir, "LOGO_cluster.sh")
    with open(w, 'w') as f:
        header(f, "LOGO_cluster", time_lim=cluster_time_limit)
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


def _fold_py(label_expr):
    """The python -c body for one fold; label_expr is a shell expr for the label."""
    return (
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
        f"{label_expr}\n"
    )


def write_single_fold_job(label, dep):
    w = os.path.join(logs_dir, f"LOGO_fold_{label}.sh".replace(' ', '_'))
    with open(w, 'w') as f:
        header(f, f"LOGO_fold_{label}".replace(' ', '_'), dep=dep)
        f.write(_fold_py(f'"{label}"'))
    return w


def write_fold_array_job(labels_file, n_folds, dep):
    w = os.path.join(logs_dir, "LOGO_fold.sh")
    with open(w, 'w') as f:
        header(f, "LOGO_fold", dep=dep, array=f"1-{n_folds}")
        # Dataset labels can contain spaces (e.g. "E. coli"); read the whole line.
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write('echo "Running fold: $LABEL"\n\n')
        f.write(_fold_py('"$LABEL"'))
    return w


def write_aggregate_job(dep):
    w = os.path.join(logs_dir, "LOGO_aggregate.sh")
    with open(w, 'w') as f:
        header(f, "LOGO_aggregate", dep=dep, time_lim="2:00:00")
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
    ap = argparse.ArgumentParser(description="Submit leave-one-dataset-out nested-CV to SLURM.")
    ap.add_argument('--test-fold', metavar='LABEL',
                    help="Submit only clustering + one fold (by dataset label) as batch jobs.")
    args = ap.parse_args()

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)

    folds = pd.read_csv(folds_file)
    fold_labels = list(dict.fromkeys(folds['fold_label'].astype(str)))
    n_folds = len(fold_labels)

    cl_id = submit(write_cluster_job())
    print(f"  clustering job: {cl_id}")

    if args.test_fold:
        if args.test_fold not in fold_labels:
            raise SystemExit(f"--test-fold '{args.test_fold}' not in folds. Available: {fold_labels}")
        fold_id = submit(write_single_fold_job(args.test_fold, dep=cl_id))
        print(f"  test fold '{args.test_fold}': {fold_id} (afterok {cl_id})")
        print(f"\n=== Test submitted ===")
        print(f"  Result: {folds_out_dir}/{args.test_fold}/model_validation/predict_results/")
        return

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
    print(f"  Results: {base_output_dir}/performance/")
    print(f"  Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
