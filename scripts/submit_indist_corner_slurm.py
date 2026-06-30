#!/usr/bin/env python3
"""
GenoPHI in-distribution STRAIN_AND_PHAGE (CORNER) baseline (new strains AND new
phages -- both axes held out).

The pLM-vs-GenoPHI comparison: run GenoPHI's corner workflow on the SAME
predefined outer folds the pLM strain_and_phage models use. Each fold trains on
training_strains x training_phages and predicts the held-out CORNER
(validation_strains x validation_phages -- pairs where BOTH are unseen).

The pLM splits already store each fold in genophi's native corner-fold format:

  strain_and_phage/splits/outer_all4_held/fold_N/
    {training,validation}_{strains,phages}.csv

so run_corner_fold_from_shared consumes them DIRECTLY (no builder needed -- the
four *_file args take these CSV paths).

Pipeline:
  1. ONE shared clustering job over the combined 4-genera AAs.
  2. A 5-task array, one task per fold (run_corner_fold_from_shared). afterok (1).
  3. ONE aggregation job (aggregate_predefined_folds, same aggregator -- corner
     folds write fold_group.txt) -> global + per-fold metrics, PR/ROC curves.

Modeling settings: min_seq_id=0.4, coverage=0.8, dynamic inverse-frequency
weights, hierarchical clustering, max_ram=100 (note min_seq_id=0.4 here vs 0.2
in genophi_combined_full_model).

  python submit_indist_corner_slurm.py --test-fold fold_0   # validate one fold
  python submit_indist_corner_slurm.py                      # full 5-fold run
"""
import os
import subprocess
import argparse
import glob

# ============================ CONFIG (edit me) ============================
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"
interaction_matrix = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"

# Per-fold split dirs (genophi native corner format -- used directly).
splits_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_and_phage/splits/outer_all4_held"
base_output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_and_phage/genophi_results"

# Modeling params (match the combined model).
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"
phage_column = "phage"
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
mem_per_job = "490G"
time_limit = "36:00:00"
cluster_time_limit = "24:00:00"
# =========================================================================

shared_dir = os.path.join(base_output_dir, "shared_clustering")
folds_out_dir = os.path.join(base_output_dir, "folds")
logs_dir = os.path.join(base_output_dir, "slurm_logs")


def fold_labels():
    return [os.path.basename(d) for d in sorted(glob.glob(os.path.join(splits_dir, "fold_*")))]


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
    w = os.path.join(logs_dir, "INDIST_CORNER_cluster.sh")
    with open(w, 'w') as f:
        header(f, "INDIST_CORNER_cluster", time_lim=cluster_time_limit)
        f.write(
            "python -c \""
            "from genophi.workflows import run_shared_clustering; "
            f"run_shared_clustering("
            f"input_strain_dir='{input_strain_dir}', "
            f"input_phage_dir='{input_phage_dir}', "
            f"output_dir='{shared_dir}', "
            f"min_seq_id={min_seq_id}, coverage={coverage}, sensitivity={sensitivity}, "
            f"threads={threads}, strain_column='{strain_column}', "
            f"phage_column='{phage_column}')\"\n"
        )
    return w


def _corner_call():
    """run_corner_fold_from_shared python -c body; label from sys.argv[1].

    The four split CSVs are read from splits_dir/<label>/ at runtime.
    """
    return (
        "python -c \""
        "import sys, os; from genophi.workflows import run_corner_fold_from_shared; "
        "lbl=sys.argv[1]; "
        f"fd=os.path.join('{splits_dir}', lbl); "
        f"run_corner_fold_from_shared("
        "fold_label=lbl, "
        "training_strains_file=os.path.join(fd,'training_strains.csv'), "
        "training_phages_file=os.path.join(fd,'training_phages.csv'), "
        "validation_strains_file=os.path.join(fd,'validation_strains.csv'), "
        "validation_phages_file=os.path.join(fd,'validation_phages.csv'), "
        f"shared_dir='{shared_dir}', "
        f"interaction_matrix='{interaction_matrix}', "
        f"input_strain_dir='{input_strain_dir}', "
        f"input_phage_dir='{input_phage_dir}', "
        f"output_dir=os.path.join('{folds_out_dir}', lbl), "
        f"num_runs_fs={num_runs_fs}, num_runs_modeling={num_runs_modeling}, "
        f"min_seq_id={min_seq_id}, coverage={coverage}, sensitivity={sensitivity}, "
        f"use_dynamic_weights={use_dynamic_weights}, "
        f"weights_method='{weights_method}', "
        f"use_clustering={use_clustering}, "
        f"cluster_method='{cluster_method}', "
        f"max_ram={max_ram}, "
        f"strain_column='{strain_column}', phage_column='{phage_column}', "
        f"threads={threads})\" "
    )


def write_single_fold_job(label, dep):
    w = os.path.join(logs_dir, f"INDIST_CORNER_fold_{label}.sh")
    with open(w, 'w') as f:
        header(f, f"INDIST_CORNER_fold_{label}", dep=dep)
        f.write(_corner_call() + f"\"{label}\"\n")
    return w


def write_fold_array_job(labels_file, n_folds, dep):
    w = os.path.join(logs_dir, "INDIST_CORNER_fold.sh")
    with open(w, 'w') as f:
        header(f, "INDIST_CORNER_fold", dep=dep, array=f"1-{n_folds}")
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write('echo "Running corner fold: $LABEL"\n\n')
        f.write(_corner_call() + '"$LABEL"\n')
    return w


def write_aggregate_job(dep):
    w = os.path.join(logs_dir, "INDIST_CORNER_aggregate.sh")
    with open(w, 'w') as f:
        header(f, "INDIST_CORNER_aggregate", dep=dep, time_lim="2:00:00", mem="120G")
        f.write(
            "python -c \""
            "from genophi.workflows import aggregate_predefined_folds; "
            f"aggregate_predefined_folds("
            f"folds_dir='{folds_out_dir}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"output_dir='{base_output_dir}', "
            f"strain_column='{strain_column}', "
            f"phage_column='{phage_column}')\"\n"
        )
    return w


def main():
    ap = argparse.ArgumentParser(description="Submit the in-distribution CORNER GenoPHI baseline.")
    ap.add_argument('--test-fold', metavar='LABEL',
                    help="Submit ONLY clustering + one fold (e.g. fold_0) to validate. "
                         "No array, no aggregation.")
    args = ap.parse_args()

    labels = fold_labels()
    if not labels:
        raise SystemExit(f"No fold_* dirs under {splits_dir}")
    n_folds = len(labels)

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)

    cl_id = submit(write_cluster_job())
    print(f"  clustering job: {cl_id}")

    if args.test_fold:
        if args.test_fold not in labels:
            raise SystemExit(f"--test-fold '{args.test_fold}' not in {labels}")
        fid = submit(write_single_fold_job(args.test_fold, dep=cl_id))
        print(f"  test fold '{args.test_fold}': {fid} (afterok {cl_id})")
        print(f"\n=== Test submitted ===")
        print(f"  Result: {folds_out_dir}/{args.test_fold}/model_validation/predict_results/")
        print(f"  Monitor: squeue -u $USER ; tail -f {logs_dir}/INDIST_CORNER_fold_{args.test_fold}_*.out")
        return

    labels_file = os.path.join(base_output_dir, "fold_labels.txt")
    with open(labels_file, 'w') as fh:
        fh.write("\n".join(labels) + "\n")
    print(f"{n_folds} folds: {labels}")

    arr_id = submit(write_fold_array_job(labels_file, n_folds, dep=cl_id))
    print(f"  fold array: {arr_id} (1-{n_folds})")
    agg_id = submit(write_aggregate_job(dep=arr_id))
    print(f"  aggregate job: {agg_id}")

    print(f"\n=== Submitted ===")
    print(f"  cluster {cl_id} -> fold array {arr_id} -> aggregate {agg_id}")
    print(f"  Results: {base_output_dir}/performance/")
    print(f"  Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
