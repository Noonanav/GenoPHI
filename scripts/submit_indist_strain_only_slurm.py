#!/usr/bin/env python3
"""
GenoPHI in-distribution STRAIN-ONLY baseline (shared phages, new strains).

The pLM-vs-GenoPHI comparison: run GenoPHI on the SAME predefined outer folds the
pLM strain_only models use (leave-strains-out; phages shared between train/val).

Pipeline (mirrors submit_loso_slurm.py):
  1. ONE shared clustering job over the combined 4-genera AAs.
  2. A 5-task array, one task per fold (run_predefined_fold), reading the long
     folds CSV produced by build_indist_strain_only_folds.py. afterok (1).
  3. ONE aggregation job (aggregate_predefined_folds) -> global + per-fold
     metrics, PR/ROC curves. afterok the array.

Modeling settings: min_seq_id=0.4, coverage=0.8, dynamic inverse-frequency
weights, hierarchical clustering, max_ram=100 -- the in-distribution comparison
param set (note min_seq_id=0.4 here vs 0.2 in genophi_combined_full_model).

  python build_indist_strain_only_folds.py            # FIRST: build the folds CSV
  python submit_indist_strain_only_slurm.py --test-fold fold_0   # validate one fold
  python submit_indist_strain_only_slurm.py                      # full 5-fold run
"""
import os
import subprocess
import argparse
import pandas as pd

# ============================ CONFIG (edit me) ============================
# Combined 4-genera inputs (identical to genophi_combined_full_model).
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"
interaction_matrix = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
folds_file = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_only/folds_strain_only.csv"

base_output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_only/genophi_results"

# Modeling params (match the combined model).
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"
min_seq_id = 0.4
coverage = 0.8
sensitivity = 7.5
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
mem_per_job = "490G"          # combined-set folds are heavy in RFE/CatBoost (240G OOM'd)
time_limit = "36:00:00"
cluster_time_limit = "24:00:00"
# =========================================================================

shared_dir = os.path.join(base_output_dir, "shared_clustering")
folds_out_dir = os.path.join(base_output_dir, "folds")
logs_dir = os.path.join(base_output_dir, "slurm_logs")


def header(f, job, dep=None, array=None, time_lim=time_limit, mem=mem_per_job):
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name={job}\n")
    f.write(f"#SBATCH --account={account}\n#SBATCH --partition={partition}\n#SBATCH --qos={qos}\n")
    f.write("#SBATCH --nodes=1\n#SBATCH --ntasks=1\n")
    f.write(f"#SBATCH --cpus-per-task={threads}\n#SBATCH --mem={mem}\n#SBATCH --time={time_lim}\n")
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


def submit(sh):
    return subprocess.check_output(["sbatch", "--parsable", sh]).decode().strip()


def write_cluster_job():
    w = os.path.join(logs_dir, "INDIST_SO_cluster.sh")
    with open(w, 'w') as f:
        header(f, "INDIST_SO_cluster", time_lim=cluster_time_limit)
        f.write(
            "python -c \""
            "from genophi.workflows import run_shared_clustering; "
            f"run_shared_clustering("
            f"input_strain_dir='{input_strain_dir}', "
            f"input_phage_dir='{input_phage_dir}', "
            f"output_dir='{shared_dir}', "
            f"min_seq_id={min_seq_id}, coverage={coverage}, sensitivity={sensitivity}, "
            f"threads={threads}, strain_column='{strain_column}')\"\n"
        )
    return w


def _fold_call():
    """The run_predefined_fold python -c body (label comes from sys.argv[1])."""
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
        f"min_seq_id={min_seq_id}, coverage={coverage}, sensitivity={sensitivity}, "
        f"use_dynamic_weights={use_dynamic_weights}, "
        f"weights_method='{weights_method}', "
        f"use_clustering={use_clustering}, "
        f"cluster_method='{cluster_method}', "
        f"max_ram={max_ram}, "
        f"threads={threads})\" "
    )


def write_single_fold_job(label, dep):
    w = os.path.join(logs_dir, f"INDIST_SO_fold_{label}.sh")
    with open(w, 'w') as f:
        header(f, f"INDIST_SO_fold_{label}", dep=dep)
        f.write(_fold_call() + f"\"{label}\"\n")
    return w


def write_fold_array_job(labels_file, n_folds, dep):
    w = os.path.join(logs_dir, "INDIST_SO_fold.sh")
    with open(w, 'w') as f:
        header(f, "INDIST_SO_fold", dep=dep, array=f"1-{n_folds}")
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write('echo "Running fold: $LABEL"\n\n')
        f.write(_fold_call() + '"$LABEL"\n')
    return w


def write_aggregate_job(dep):
    w = os.path.join(logs_dir, "INDIST_SO_aggregate.sh")
    with open(w, 'w') as f:
        header(f, "INDIST_SO_aggregate", dep=dep, time_lim="2:00:00", mem="120G")
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
    ap = argparse.ArgumentParser(description="Submit the in-distribution strain-only GenoPHI baseline.")
    ap.add_argument('--test-fold', metavar='LABEL',
                    help="Submit ONLY clustering + one fold (e.g. fold_0) to validate the "
                         "pipeline before the full array. No array, no aggregation.")
    args = ap.parse_args()

    if not os.path.exists(folds_file):
        raise SystemExit(f"folds_file not found: {folds_file}\n"
                         f"Run build_indist_strain_only_folds.py first.")

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)

    folds = pd.read_csv(folds_file)
    fold_labels = list(dict.fromkeys(folds['fold_label'].astype(str)))
    n_folds = len(fold_labels)

    cl_id = submit(write_cluster_job())
    print(f"  clustering job: {cl_id}")

    if args.test_fold:
        if args.test_fold not in fold_labels:
            raise SystemExit(f"--test-fold '{args.test_fold}' not in {fold_labels}")
        fid = submit(write_single_fold_job(args.test_fold, dep=cl_id))
        print(f"  test fold '{args.test_fold}': {fid} (afterok {cl_id})")
        print(f"\n=== Test submitted ===")
        print(f"  Result: {folds_out_dir}/{args.test_fold}/model_validation/predict_results/")
        print(f"  Monitor: squeue -u $USER ; tail -f {logs_dir}/INDIST_SO_fold_{args.test_fold}_*.out")
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
    print(f"  cluster {cl_id} -> fold array {arr_id} -> aggregate {agg_id}")
    print(f"  Results: {base_output_dir}/performance/  (global + per-fold metrics, curves)")
    print(f"  Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
