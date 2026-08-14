#!/usr/bin/env python3
"""
Submit leave-one-genus-out CORNER nested-CV using PURE k-mer features.

Corner = hold out BOTH a genus's strains AND its phages; train on the training
quadrant; predict the held-out corner (validation_strains x validation_phages).
Uses genophi.workflows.run_corner_kmer_fold (k-mers from sequences; held-out
genomes assigned by substring matching -- NO MMseqs2, NO shared clustering).

Because there is no clustering step, each fold is fully independent: this submits
just a fold array + an aggregate job (no cluster dependency). Unlike the
protein-family corner submitter, there is no (min_seq_id, coverage) grid -- the
k-mer sweep axis is k (set K below).

Fold inputs follow the cross_genus_splits/logo convention (fold_N/ with
training/validation x strains/phages CSVs + held_out_group.txt).

  python submit_corner_kmer_logo_slurm.py                    # all folds
  python submit_corner_kmer_logo_slurm.py --test-fold fold_0 # one fold
"""
import os
import argparse
import subprocess

# ============================ CONFIG (edit me) ============================
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"
interaction_matrix = "/global/home/groups/pc_phiml/set_transformer/manuscript/cross_genus_splits/logo/representative_subset.csv"
folds_root = "/global/home/groups/pc_phiml/set_transformer/manuscript/cross_genus_splits/logo"

base_output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_corner_kmer_logo"

# k-mer length (k=4 per the experiment plan).
K = 4
K_RANGE = False

# Modeling params (match the lab's intended settings).
num_runs_fs = 25
num_runs_modeling = 50
num_features = 100
strain_column = "strain"
phage_column = "phage"
use_dynamic_weights = True
weights_method = "inverse_frequency"
use_clustering = True
cluster_method = "hierarchical"
max_ram = 40

# SLURM
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"
threads = 32
mem_per_job = "320G"
time_limit = "14:00:00"
# =========================================================================

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


def list_fold_dirs():
    folds = [d for d in os.listdir(folds_root)
             if d.startswith("fold_") and os.path.isdir(os.path.join(folds_root, d))]
    return sorted(folds, key=lambda d: int(d.split("_")[1]))


def fold_label(fold_name):
    f = os.path.join(folds_root, fold_name, "held_out_group.txt")
    if os.path.exists(f):
        with open(f) as fh:
            t = fh.read().strip()
            if t:
                return t
    return fold_name


def _kmer_corner_py(label_expr, fold_dir_expr):
    return (
        "python -c \""
        "import sys; from genophi.workflows import run_corner_kmer_fold; "
        "fl, fd = sys.argv[1], sys.argv[2]; "
        "run_corner_kmer_fold("
        "fold_label=fl, "
        "training_strains_file=fd + '/training_strains.csv', "
        "training_phages_file=fd + '/training_phages.csv', "
        "validation_strains_file=fd + '/validation_strains.csv', "
        "validation_phages_file=fd + '/validation_phages.csv', "
        f"interaction_matrix='{interaction_matrix}', "
        f"input_strain_dir='{input_strain_dir}', "
        f"input_phage_dir='{input_phage_dir}', "
        f"output_dir='{folds_out_dir}/' + fl.replace(' ', '_'), "
        f"k={K}, k_range={K_RANGE}, "
        f"num_runs_fs={num_runs_fs}, num_runs_modeling={num_runs_modeling}, num_features={num_features}, "
        f"use_dynamic_weights={use_dynamic_weights}, weights_method='{weights_method}', "
        f"use_clustering={use_clustering}, cluster_method='{cluster_method}', "
        f"max_ram={max_ram}, threads={threads}, "
        f"strain_column='{strain_column}', phage_column='{phage_column}')\" "
        f"{label_expr} {fold_dir_expr}\n"
    )


def write_single_fold_job(fold_name):
    label = fold_label(fold_name)
    fdir = os.path.join(folds_root, fold_name)
    w = os.path.join(logs_dir, f"KMER_fold_{fold_name}.sh")
    with open(w, 'w') as f:
        header(f, f"KMER_fold_{fold_name}")
        f.write(_kmer_corner_py(f'"{label}"', f'"{fdir}"'))
    return w


def write_fold_array_job(fold_names):
    labels_file = os.path.join(base_output_dir, "fold_labels.txt")
    dirs_file = os.path.join(base_output_dir, "fold_dirs.txt")
    with open(labels_file, 'w') as fh:
        fh.write("\n".join(fold_label(fn) for fn in fold_names) + "\n")
    with open(dirs_file, 'w') as fh:
        fh.write("\n".join(os.path.join(folds_root, fn) for fn in fold_names) + "\n")

    w = os.path.join(logs_dir, "KMER_fold.sh")
    with open(w, 'w') as f:
        header(f, "KMER_fold", array=f"1-{len(fold_names)}")
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write(f'FOLDDIR=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {dirs_file})\n')
        f.write('echo "k-mer corner fold: $LABEL ($FOLDDIR)"\n\n')
        f.write(_kmer_corner_py('"$LABEL"', '"$FOLDDIR"'))
    return w


def write_aggregate_job(dep):
    w = os.path.join(logs_dir, "KMER_aggregate.sh")
    with open(w, 'w') as f:
        header(f, "KMER_aggregate", dep=dep, time_lim="2:00:00")
        f.write(
            "python -c \""
            "from genophi.workflows import aggregate_predefined_folds; "
            f"aggregate_predefined_folds("
            f"folds_dir='{folds_out_dir}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"output_dir='{base_output_dir}', "
            f"strain_column='{strain_column}', phage_column='{phage_column}')\"\n"
        )
    return w


def main():
    ap = argparse.ArgumentParser(description="Submit pure-k-mer corner LOGO nested-CV.")
    ap.add_argument('--test-fold', metavar='FOLD', help="Submit only one fold (e.g. fold_0).")
    args = ap.parse_args()

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(folds_out_dir, exist_ok=True)
    fold_names = list_fold_dirs()
    if not fold_names:
        raise SystemExit(f"No fold_* dirs under {folds_root}")
    print(f"folds: {fold_names}  (k={K})")

    if args.test_fold:
        if args.test_fold not in fold_names:
            raise SystemExit(f"--test-fold '{args.test_fold}' not found. Available: {fold_names}")
        fid = submit(write_single_fold_job(args.test_fold))
        print(f"test k-mer corner fold '{args.test_fold}': {fid}")
        print(f"result: {folds_out_dir}/<label>/model_validation/predict_results/")
        return

    arr_id = submit(write_fold_array_job(fold_names))
    print(f"k-mer corner array: {arr_id} (1-{len(fold_names)})")
    agg_id = submit(write_aggregate_job(dep=arr_id))
    print(f"aggregate: {agg_id}")
    print(f"\n=== Submitted. results -> {base_output_dir}/performance/ ===")
    print("Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
