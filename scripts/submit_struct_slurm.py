#!/usr/bin/env python3
"""
Submit a structure-stratified (leave-one-structural-cluster-out) nested-CV
experiment as SLURM jobs -- the O-antigen-CHEMISTRY analogue of
submit_phylo_slurm.py / submit_loso_slurm.py, sharing their exact execution path
(shared clustering -> per-fold array -> aggregate) so the struct-CV vs phylo-CV vs
serotype-CV vs random-CV comparison is clean.

Motivation: LOSO-O never tests structural generalization -- every held-out
serogroup still has a structurally-near training serogroup (median min-GED-to-
training 0.167, vs a panel median PAIRWISE GED of 0.500). Holding out whole
structural clusters raises the median to 0.264 with a 0.167-0.385 spread, giving a
real AUROC-vs-chemistry-distance gradient. Folds from scripts/build_struct_folds.py.

REUSING SHARED CLUSTERING. run_shared_clustering() is sequence-only and
label-independent ("can be run a single time and reused across folds" -- its
docstring); the phylo submitter passes no strain_list/phage_list, so it clustered
every FASTA in the input dirs. Its shared_clustering/ output is therefore valid for
ANY fold design over the same inputs. Pass --shared-clustering <dir> to reuse it and
skip the 24 h clustering job entirely:

    python submit_struct_slurm.py --k 20 \\
        --shared-clustering /global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results/shared_clustering

Omit the flag to run clustering fresh (identical to the phylo/LOSO submitters).

    python submit_struct_slurm.py --k 20 --test-fold struct_01   # validate one fold
    python submit_struct_slurm.py --k 20                          # full run
"""
import os
import subprocess
import pandas as pd

# ============================ CONFIG (edit me) ============================
# Inputs (paths as they exist on the Lawrencium compute nodes) -- SAME as LOSO/phylo.
input_strain_dir = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/strain_AAs_update"
input_phage_dir = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/phage_AAs"
interaction_matrix = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/Gaborieau_interaction_matrix_long_mod.csv"

# Struct fold CSVs (scp built-locally outputs of build_struct_folds.py to here).
struct_folds_dir = "/global/scratch/users/anoonan/BRaVE/struct_cv/tables"
struct_base_dir = "/global/scratch/users/anoonan/BRaVE/struct_cv"

# Modeling params (IDENTICAL to LOSO-O / phylo, for comparability)
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"
use_dynamic_weights = True
weights_method = "inverse_frequency"
use_clustering = True
cluster_method = "hierarchical"
max_ram = 100

# SLURM (matches the lab template / LOSO / phylo)
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"
threads = 8
mem_per_job = "120G"
time_limit = "24:00:00"
cluster_time_limit = "24:00:00"
# =========================================================================


def paths_for_k(k):
    folds_file = os.path.join(struct_folds_dir, f"folds_struct_K{k}.csv")
    base_output_dir = os.path.join(struct_base_dir, f"K{k}_results")
    return folds_file, base_output_dir


def header(f, job, logs_dir, dep=None, array=None, time_lim=time_limit):
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


def write_cluster_job(logs_dir, shared_dir, tag):
    w = os.path.join(logs_dir, f"STRUCT{tag}_cluster.sh")
    with open(w, 'w') as f:
        header(f, f"STRUCT{tag}_cluster", logs_dir, time_lim=cluster_time_limit)
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


def _fold_call(folds_file, shared_dir, folds_out_dir):
    """The python -c body shared by the single-fold and array jobs."""
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
        f"threads={threads})\""
    )


def write_single_fold_job(label, dep, logs_dir, folds_file, shared_dir, folds_out_dir, tag):
    w = os.path.join(logs_dir, f"STRUCT{tag}_fold_{label}.sh")
    with open(w, 'w') as f:
        header(f, f"STRUCT{tag}_fold_{label}", logs_dir, dep=dep)
        f.write(_fold_call(folds_file, shared_dir, folds_out_dir) + f' "{label}"\n')
    return w


def write_fold_array_job(labels_file, n_folds, dep, logs_dir, folds_file, shared_dir, folds_out_dir, tag):
    w = os.path.join(logs_dir, f"STRUCT{tag}_fold.sh")
    with open(w, 'w') as f:
        header(f, f"STRUCT{tag}_fold", logs_dir, dep=dep, array=f"1-{n_folds}")
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write('echo "Running fold: $LABEL"\n\n')
        f.write(_fold_call(folds_file, shared_dir, folds_out_dir) + ' "$LABEL"\n')
    return w


def write_aggregate_job(dep, logs_dir, folds_out_dir, base_output_dir, tag):
    w = os.path.join(logs_dir, f"STRUCT{tag}_aggregate.sh")
    with open(w, 'w') as f:
        header(f, f"STRUCT{tag}_aggregate", logs_dir, dep=dep, time_lim="2:00:00")
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
    ap = argparse.ArgumentParser(
        description="Submit the structure-stratified (leave-one-structural-cluster-out) nested-CV experiment.")
    ap.add_argument('--k', type=int, required=True, choices=[20],
                    help="Which fold set to run (20 = average-linkage GED clusters, min_size=3).")
    ap.add_argument('--shared-clustering', metavar='DIR',
                    help="Reuse an EXISTING shared_clustering dir (e.g. from the phylo run) "
                         "and skip the 24h clustering job. run_shared_clustering is "
                         "sequence-only and label-independent, so this is safe across CV designs.")
    ap.add_argument('--test-fold', metavar='LABEL',
                    help="Submit ONLY one fold (e.g. struct_01) to validate the pipeline "
                         "end-to-end before the full array.")
    args = ap.parse_args()

    tag = str(args.k)
    folds_file, base_output_dir = paths_for_k(args.k)
    folds_out_dir = os.path.join(base_output_dir, "folds")
    logs_dir = os.path.join(base_output_dir, "slurm_logs")

    if not os.path.exists(folds_file):
        raise SystemExit(f"folds file not found: {folds_file}\n"
                         f"scp the build_struct_folds.py output (folds_struct_K{args.k}.csv) here first.")

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)

    folds = pd.read_csv(folds_file)
    fold_labels = list(dict.fromkeys(folds['fold_label'].astype(str)))
    n_folds = len(fold_labels)

    # --- shared clustering: reuse an existing dir, or run it fresh ---
    if args.shared_clustering:
        shared_dir = os.path.abspath(args.shared_clustering)
        if not os.path.isdir(shared_dir):
            raise SystemExit(f"--shared-clustering dir not found: {shared_dir}")
        for sub in ("strain", "phage"):
            if not os.path.isdir(os.path.join(shared_dir, sub)):
                raise SystemExit(f"shared clustering dir missing '{sub}/' subdir: {shared_dir}")
        cl_id = None
        print(f"  REUSING shared clustering: {shared_dir}  (clustering job skipped)")
    else:
        shared_dir = os.path.join(base_output_dir, "shared_clustering")
        cl_id = submit(write_cluster_job(logs_dir, shared_dir, tag))
        print(f"  clustering job: {cl_id}")

    if args.test_fold:
        if args.test_fold not in fold_labels:
            raise SystemExit(f"--test-fold '{args.test_fold}' not in folds file. "
                             f"Available: {fold_labels}")
        fold_id = submit(write_single_fold_job(args.test_fold, cl_id, logs_dir,
                                               folds_file, shared_dir, folds_out_dir, tag))
        dep = f"(afterok {cl_id})" if cl_id else "(no dependency; reused clustering)"
        print(f"  test fold '{args.test_fold}': {fold_id} {dep}")
        print(f"\n=== Test submitted (K={args.k}) ===")
        print(f"  Result: {folds_out_dir}/{args.test_fold}/model_validation/predict_results/")
        print(f"  Monitor: squeue -u $USER ; tail -f {logs_dir}/STRUCT{tag}_fold_{args.test_fold}_*.out")
        return

    labels_file = os.path.join(base_output_dir, "fold_labels.txt")
    with open(labels_file, 'w') as fh:
        fh.write("\n".join(fold_labels) + "\n")
    print(f"K={args.k}: {n_folds} folds: {fold_labels}")

    arr_id = submit(write_fold_array_job(labels_file, n_folds, cl_id, logs_dir,
                                         folds_file, shared_dir, folds_out_dir, tag))
    print(f"  fold array: {arr_id} (1-{n_folds})")

    agg_id = submit(write_aggregate_job(arr_id, logs_dir, folds_out_dir, base_output_dir, tag))
    print(f"  aggregate job: {agg_id}")

    print(f"\n=== Submitted (K={args.k}) ===")
    if cl_id:
        print(f"  cluster {cl_id} -> fold array {arr_id} (1-{n_folds}) -> aggregate {agg_id}")
    else:
        print(f"  [reused clustering] -> fold array {arr_id} (1-{n_folds}) -> aggregate {agg_id}")
    print(f"  Results: {base_output_dir}/performance/")
    print(f"  Per-fold: {folds_out_dir}/<struct_NN>/")


if __name__ == "__main__":
    main()
