#!/usr/bin/env python3
"""
Submit leave-one-genus-out CORNER nested-CV as SLURM jobs, swept over a grid of
MMseqs2 (min_seq_id, coverage) thresholds.

Corner fold = hold out BOTH a genus's strains AND its phages; train on the
training quadrant; predict the held-out corner (validation_strains x
validation_phages, both unseen). Uses genophi.workflows.run_corner_fold_from_shared.

Each (min_seq_id, coverage) combo is its OWN pipeline (clustering depends on the
thresholds, so it cannot be shared across the grid):
    per combo:  1 shared-clustering job  ->  fold corner array (afterok)  ->  aggregate (afterok)

Fold inputs follow the cross_genus_splits/logo convention: a directory of
fold_N/ subdirs, each containing training_strains.csv, training_phages.csv,
validation_strains.csv, validation_phages.csv (+ held_out_group.txt for a label).
The wrapper reads that convention and passes the four explicit paths to the
(format-agnostic) package function.

  python submit_corner_logo_slurm.py                       # full grid
  python submit_corner_logo_slurm.py --test-fold fold_0    # one combo, one fold
  python submit_corner_logo_slurm.py --min_seq_id 0.4 --coverage 0.8   # single combo
"""
import os
import argparse
import subprocess

# ============================ CONFIG (edit me) ============================
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"
interaction_matrix = "/global/home/groups/pc_phiml/set_transformer/manuscript/cross_genus_splits/logo/representative_subset.csv"
folds_root = "/global/home/groups/pc_phiml/set_transformer/manuscript/cross_genus_splits/logo"

base_output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_corner_logo"

# Threshold grid: 3 min_seq_id x 2 coverage = 6 combos.
MIN_SEQ_IDS = [0.2, 0.3, 0.4]
COVERAGES = [0.5, 0.8]

# Modeling params (match the lab's intended settings).
num_runs_fs = 25
num_runs_modeling = 50
strain_column = "strain"
phage_column = "phage"
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


def combo_tag(msi, cov):
    return f"msi{msi}_c{cov}".replace('.', '')


def combo_dirs(msi, cov):
    root = os.path.join(base_output_dir, combo_tag(msi, cov))
    return {
        'root': root,
        'shared': os.path.join(root, 'shared_clustering'),
        'folds': os.path.join(root, 'folds'),
        'logs': os.path.join(root, 'slurm_logs'),
    }


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


def list_fold_dirs():
    """Fold subdirs (fold_*) under folds_root, sorted by index."""
    folds = [d for d in os.listdir(folds_root)
             if d.startswith("fold_") and os.path.isdir(os.path.join(folds_root, d))]
    return sorted(folds, key=lambda d: int(d.split("_")[1]))


def fold_label(fold_name):
    """Human label = held_out_group.txt contents, else the fold dir name."""
    f = os.path.join(folds_root, fold_name, "held_out_group.txt")
    if os.path.exists(f):
        with open(f) as fh:
            t = fh.read().strip()
            if t:
                return t
    return fold_name


def write_cluster_job(msi, cov, dirs):
    w = os.path.join(dirs['logs'], "CORNER_cluster.sh")
    with open(w, 'w') as f:
        header(f, "CORNER_cluster", dirs['logs'], time_lim=cluster_time_limit)
        f.write(
            "python -c \""
            "from genophi.workflows import run_shared_clustering; "
            f"run_shared_clustering("
            f"input_strain_dir='{input_strain_dir}', "
            f"input_phage_dir='{input_phage_dir}', "
            f"output_dir='{dirs['shared']}', "
            f"min_seq_id={msi}, coverage={cov}, "
            f"threads={threads}, strain_column='{strain_column}')\"\n"
        )
    return w


def _corner_py(msi, cov, dirs, fold_label_expr, fold_dir_expr):
    """python -c body for one corner fold. fold_label_expr / fold_dir_expr are
    shell expressions (already quoted) for the label and the fold's source dir."""
    return (
        "python -c \""
        "import sys; from genophi.workflows import run_corner_fold_from_shared; "
        "fl, fd = sys.argv[1], sys.argv[2]; "
        f"run_corner_fold_from_shared("
        "fold_label=fl, "
        "training_strains_file=fd + '/training_strains.csv', "
        "training_phages_file=fd + '/training_phages.csv', "
        "validation_strains_file=fd + '/validation_strains.csv', "
        "validation_phages_file=fd + '/validation_phages.csv', "
        f"shared_dir='{dirs['shared']}', "
        f"interaction_matrix='{interaction_matrix}', "
        f"input_strain_dir='{input_strain_dir}', "
        f"input_phage_dir='{input_phage_dir}', "
        f"output_dir='{dirs['folds']}/' + fl.replace(' ', '_'), "
        f"min_seq_id={msi}, coverage={cov}, "
        f"num_runs_fs={num_runs_fs}, num_runs_modeling={num_runs_modeling}, "
        f"use_dynamic_weights={use_dynamic_weights}, weights_method='{weights_method}', "
        f"use_clustering={use_clustering}, cluster_method='{cluster_method}', "
        f"max_ram={max_ram}, threads={threads}, "
        f"strain_column='{strain_column}', phage_column='{phage_column}')\" "
        f"{fold_label_expr} {fold_dir_expr}\n"
    )


def write_single_corner_job(msi, cov, dirs, fold_name, dep):
    label = fold_label(fold_name)
    fdir = os.path.join(folds_root, fold_name)
    w = os.path.join(dirs['logs'], f"CORNER_fold_{fold_name}.sh")
    with open(w, 'w') as f:
        header(f, f"CORNER_fold_{fold_name}", dirs['logs'], dep=dep)
        f.write(_corner_py(msi, cov, dirs, f'"{label}"', f'"{fdir}"'))
    return w


def write_corner_array_job(msi, cov, dirs, fold_names, dep):
    # Map array index -> fold name -> (label, source dir) via lookup files.
    labels_file = os.path.join(dirs['root'], "fold_labels.txt")
    dirs_file = os.path.join(dirs['root'], "fold_dirs.txt")
    with open(labels_file, 'w') as fh:
        fh.write("\n".join(fold_label(fn) for fn in fold_names) + "\n")
    with open(dirs_file, 'w') as fh:
        fh.write("\n".join(os.path.join(folds_root, fn) for fn in fold_names) + "\n")

    w = os.path.join(dirs['logs'], "CORNER_fold.sh")
    with open(w, 'w') as f:
        header(f, "CORNER_fold", dirs['logs'], dep=dep, array=f"1-{len(fold_names)}")
        f.write(f'LABEL=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {labels_file})\n')
        f.write(f'FOLDDIR=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {dirs_file})\n')
        f.write('echo "Corner fold: $LABEL ($FOLDDIR)"\n\n')
        f.write(_corner_py(msi, cov, dirs, '"$LABEL"', '"$FOLDDIR"'))
    return w


def write_aggregate_job(dirs, dep):
    w = os.path.join(dirs['logs'], "CORNER_aggregate.sh")
    with open(w, 'w') as f:
        header(f, "CORNER_aggregate", dirs['logs'], dep=dep, time_lim="2:00:00")
        f.write(
            "python -c \""
            "from genophi.workflows import aggregate_predefined_folds; "
            f"aggregate_predefined_folds("
            f"folds_dir='{dirs['folds']}', "
            f"interaction_matrix='{interaction_matrix}', "
            f"output_dir='{dirs['root']}', "
            f"strain_column='{strain_column}', phage_column='{phage_column}')\"\n"
        )
    return w


def submit_combo(msi, cov, fold_names, test_fold=None):
    dirs = combo_dirs(msi, cov)
    os.makedirs(dirs['logs'], exist_ok=True)
    os.makedirs(dirs['folds'], exist_ok=True)
    tag = combo_tag(msi, cov)

    cl_id = submit(write_cluster_job(msi, cov, dirs))
    print(f"[{tag}] cluster: {cl_id}")

    if test_fold:
        fold_id = submit(write_single_corner_job(msi, cov, dirs, test_fold, dep=cl_id))
        print(f"[{tag}] test corner fold '{test_fold}': {fold_id} (afterok {cl_id})")
        print(f"[{tag}] result: {dirs['folds']}/<label>/model_validation/predict_results/")
        return

    arr_id = submit(write_corner_array_job(msi, cov, dirs, fold_names, dep=cl_id))
    print(f"[{tag}] corner array: {arr_id} (1-{len(fold_names)})")
    agg_id = submit(write_aggregate_job(dirs, dep=arr_id))
    print(f"[{tag}] aggregate: {agg_id}")
    print(f"[{tag}] results -> {dirs['root']}/performance/")


def main():
    ap = argparse.ArgumentParser(description="Submit corner LOGO nested-CV threshold sweep.")
    ap.add_argument('--test-fold', metavar='FOLD',
                    help="Submit only clustering + one corner fold (e.g. fold_0) for ONE combo.")
    ap.add_argument('--min_seq_id', type=float, help="Run a single combo at this min_seq_id.")
    ap.add_argument('--coverage', type=float, help="Run a single combo at this coverage.")
    args = ap.parse_args()

    os.makedirs(base_output_dir, exist_ok=True)
    fold_names = list_fold_dirs()
    if not fold_names:
        raise SystemExit(f"No fold_* dirs under {folds_root}")

    # Which threshold combos?
    if args.min_seq_id is not None and args.coverage is not None:
        combos = [(args.min_seq_id, args.coverage)]
    elif args.test_fold:
        # Default the test to the first grid combo unless overridden above.
        combos = [(MIN_SEQ_IDS[-1], COVERAGES[-1])]  # 0.4 / 0.8
    else:
        combos = [(m, c) for m in MIN_SEQ_IDS for c in COVERAGES]

    print(f"folds: {fold_names}")
    print(f"combos: {combos}")
    for msi, cov in combos:
        submit_combo(msi, cov, fold_names, test_fold=args.test_fold)

    print("\n=== Submitted. Monitor: squeue -u $USER ===")
    print(f"Per-combo results under: {base_output_dir}/<combo>/performance/")


if __name__ == "__main__":
    main()
