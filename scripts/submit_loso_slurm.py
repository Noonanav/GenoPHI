#!/usr/bin/env python3
"""
Submit a leave-one-serotype-out (LOSO) nested-CV experiment as SLURM jobs.

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
interaction_matrix = "/global/scratch/users/anoonan/BRaVE/LOSO/tables/Gaborieau_interaction_matrix_long_mod.csv"
folds_file = "/global/scratch/users/anoonan/BRaVE/LOSO/tables/folds_loso_O.csv"

# Output root for this experiment
base_output_dir = "/global/scratch/users/anoonan/BRaVE/LOSO/loso_O_results"

# Modeling params (passed to each fold)
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"

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


def main():
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)

    # Fold labels, in file order.
    folds = pd.read_csv(folds_file)
    fold_labels = list(dict.fromkeys(folds['fold_label'].astype(str)))
    n_folds = len(fold_labels)
    print(f"{n_folds} folds: {fold_labels}")

    # A label->index file so the array task can map SLURM_ARRAY_TASK_ID -> label.
    labels_file = os.path.join(base_output_dir, "fold_labels.txt")
    with open(labels_file, 'w') as fh:
        fh.write("\n".join(fold_labels) + "\n")

    # ---- Job 1: shared clustering (once) ----
    j1 = "LOSO_cluster"
    w1 = os.path.join(logs_dir, f"{j1}.sh")
    with open(w1, 'w') as f:
        header(f, j1, time_lim=cluster_time_limit)
        f.write(
            "python -c \""
            "from genophi.workflows import run_shared_clustering; "
            f"run_shared_clustering("
            f"input_strain_dir='{input_strain_dir}', "
            f"input_phage_dir='{input_phage_dir}', "
            f"output_dir='{shared_dir}', threads={threads}, "
            f"strain_column='{strain_column}')\"\n"
        )
    cl_id = subprocess.check_output(["sbatch", "--parsable", w1]).decode().strip()
    print(f"  clustering job: {cl_id}")

    # ---- Job 2: fold array (afterok clustering) ----
    j2 = "LOSO_fold"
    w2 = os.path.join(logs_dir, f"{j2}.sh")
    with open(w2, 'w') as f:
        header(f, j2, dep=cl_id, array=f"1-{n_folds}")
        # Map array index (1-based) -> fold label, then run that fold.
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
            f"threads={threads})\" \"$LABEL\"\n"
        )
    arr_id = subprocess.check_output(["sbatch", "--parsable", w2]).decode().strip()
    print(f"  fold array: {arr_id} (1-{n_folds})")

    # ---- Job 3: aggregate (afterok whole array) ----
    j3 = "LOSO_aggregate"
    w3 = os.path.join(logs_dir, f"{j3}.sh")
    with open(w3, 'w') as f:
        # afterok on an array job id waits for ALL tasks to succeed.
        header(f, j3, dep=arr_id, time_lim="2:00:00")
        f.write(
            "python -c \""
            "from genophi.workflows import aggregate_predefined_folds; "
            f"aggregate_predefined_folds("
            f"folds_dir='{folds_out_dir}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"output_dir='{base_output_dir}', "
            f"strain_column='{strain_column}')\"\n"
        )
    agg_id = subprocess.check_output(["sbatch", "--parsable", w3]).decode().strip()
    print(f"  aggregate job: {agg_id}")

    print(f"\n=== Submitted ===")
    print(f"  cluster {cl_id} -> fold array {arr_id} (1-{n_folds}) -> aggregate {agg_id}")
    print(f"  Results: {base_output_dir}/performance/  (global + per-fold metrics, curves)")
    print(f"  Per-fold: {folds_out_dir}/<O-group>/")
    print(f"  Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
