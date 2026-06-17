#!/usr/bin/env python3
"""
SLURM submission for INDEPENDENT multi-output CV, parallelized as one job per
(fold x target) so it fits cluster walltimes even with many targets.

Three stages (each a SLURM job/array, chained by afterok):
  Stage A  fold tables   --array=1-n_folds          build_fold_table(fold)
  Stage B  per-target    --array=1-(n_folds*targets) train_fold_target(fold,target)
           (afterok Stage A)
  Stage C  aggregate     1 job (afterok Stage B)     aggregate_independent_cv(...)

Stage B task index maps to (fold, target):
    fold   = (idx-1) // n_targets + 1
    target = targets[(idx-1) % n_targets]

Splits are the deterministic _kfold_splits used by the joint run, so folds match
exactly (diff cv_splits.csv to confirm). Use this for --strategy independent
instead of the single-job-per-fold wrapper, which serializes all targets per
fold and will time out.
"""
import os
import sys
import argparse
import subprocess
import time


def submit_job(script_path):
    try:
        r = subprocess.run(['sbatch', '--parsable', script_path],
                           capture_output=True, text=True, check=True)
        return r.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error submitting {script_path}: {e}\n{e.stderr}")
        return None


def _slurm_header(args, jobname, extra_lines):
    """Build the SBATCH header. extra_lines is a list of '#SBATCH ...' strings."""
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={jobname}",
        f"#SBATCH --account={args.account}",
        f"#SBATCH --partition={args.partition}",
        f"#SBATCH --qos={args.qos}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={args.threads}",
    ] + list(extra_lines)
    return "\n".join(lines)


def _conda_preamble(args):
    return (
        "module load miniforge3\n"
        'eval "$(conda shell.bash hook)"\n'
        f"conda activate {args.environment}\n"
    )


def _common_fold_kwargs(args, targets):
    # Used by Stage A (build_fold_table) and Stage C (aggregate_independent_cv),
    # both of which resolve splits -> they get the group args. Stage B
    # (train_fold_target) builds its own call and is split-agnostic.
    return (
        f"input_strain_dir={args.input_strain_dir!r}, "
        f"phenotype_matrix={args.phenotype_matrix!r}, "
        f"output_dir={args.output_dir!r}, "
        f"phenotype_column={targets!r}, "
        f"task_type={args.task_type!r}, "
        f"n_folds={args.n_folds}, cv_rounds={args.cv_rounds}, "
        f"sample_column={args.sample_column!r}, suffix={args.suffix!r}, "
        f"group_metadata={args.group_metadata!r}, "
        f"group_strain_column={args.group_strain_column!r}, "
        f"group_column={args.group_column!r}"
    )


def _n_folds(args):
    """Fold count: group count for phylo folds, else n_folds * cv_rounds."""
    if args.group_metadata:
        import pandas as pd
        g = pd.read_csv(args.group_metadata)
        groups = g[args.group_column].dropna()
        groups = groups[groups.astype(str).str.strip() != '']
        return int(groups.nunique())
    return args.n_folds * args.cv_rounds


def create_stage_a(args, run_dir, targets):
    n = _n_folds(args)
    common = _common_fold_kwargs(args, targets)
    header = _slurm_header(args, 'mo_cv_table', [
        f"#SBATCH --mem={args.mem_table}G",
        f"#SBATCH --time={args.time_table}",
        f"#SBATCH --array=1-{n}",
        "#SBATCH --output=logs/table_%A_%a.out",
        "#SBATCH --error=logs/table_%A_%a.err",
    ])
    body = f"""{header}
echo "=== Stage A: fold table $SLURM_ARRAY_TASK_ID ==="; echo "$(date)"
{_conda_preamble(args)}
python3 << 'PYTHON_SCRIPT'
import logging, os
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
from genophi.workflows.multioutput_cv_parallel import build_fold_table
fold = int(os.environ['SLURM_ARRAY_TASK_ID'])
build_fold_table(fold_idx=fold, {common}, k={args.k}, threads={args.threads}, max_ram={args.max_ram},
    use_feature_clustering={args.use_feature_clustering},
    feature_n_clusters={args.feature_n_clusters},
    feature_min_cluster_presence={args.feature_min_cluster_presence})
print(f"Fold {{fold}} table built.")
PYTHON_SCRIPT
echo "$(date)"
"""
    path = os.path.join(run_dir, "stageA_tables.sh")
    open(path, 'w').write(body)
    os.chmod(path, 0o755)
    return path


def create_stage_b(args, run_dir, targets, dependency):
    n_targets = len(targets)
    n = _n_folds(args) * n_targets
    header = _slurm_header(args, 'mo_cv_train', [
        f"#SBATCH --mem={args.mem_train}G",
        f"#SBATCH --time={args.time_train}",
        f"#SBATCH --array=1-{n}",
        f"#SBATCH --dependency=afterok:{dependency}",
        "#SBATCH --output=logs/train_%A_%a.out",
        "#SBATCH --error=logs/train_%A_%a.err",
    ])
    body = f"""{header}
echo "=== Stage B: train job $SLURM_ARRAY_TASK_ID ==="; echo "$(date)"
{_conda_preamble(args)}
python3 << 'PYTHON_SCRIPT'
import logging, os
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
from genophi.workflows.multioutput_cv_parallel import train_fold_target

targets = {targets!r}
n_targets = len(targets)
idx = int(os.environ['SLURM_ARRAY_TASK_ID'])
fold = (idx - 1) // n_targets + 1
target = targets[(idx - 1) % n_targets]
print(f"Job {{idx}}: fold {{fold}}, target {{target}}")

train_fold_target(
    fold_idx=fold, target=target,
    output_dir={args.output_dir!r}, phenotype_column=targets,
    task_type={args.task_type!r}, n_folds={args.n_folds}, cv_rounds={args.cv_rounds},
    sample_column={args.sample_column!r},
    num_runs_fs={args.num_runs_fs}, num_runs_modeling={args.num_runs_modeling},
    method={args.method!r}, max_features={args.max_features!r},
    threads={args.threads}, max_ram={args.max_ram},
    use_clustering={args.use_clustering},
    filter_by_cluster_presence={args.filter_by_cluster_presence},
    min_cluster_presence={args.min_cluster_presence},
)
print(f"fold {{fold}} target {{target}} done.")
PYTHON_SCRIPT
echo "$(date)"
"""
    path = os.path.join(run_dir, "stageB_train.sh")
    open(path, 'w').write(body)
    os.chmod(path, 0o755)
    return path


def create_stage_c(args, run_dir, targets, dependency):
    common = _common_fold_kwargs(args, targets)
    header = _slurm_header(args, 'mo_cv_agg', [
        "#SBATCH --mem=16G",
        "#SBATCH --time=2:00:00",
        f"#SBATCH --dependency=afterok:{dependency}",
        "#SBATCH --output=logs/aggregate_%j.out",
        "#SBATCH --error=logs/aggregate_%j.err",
    ])
    body = f"""{header}
echo "=== Stage C: aggregate ==="; echo "$(date)"
{_conda_preamble(args)}
python3 << 'PYTHON_SCRIPT'
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
from genophi.workflows.multioutput_cv_parallel import aggregate_independent_cv
res = aggregate_independent_cv({common}, threads={args.threads}, strong_top_frac={args.strong_top_frac})
print("CV predictions:", res['predictions'])
print("CV metrics:", res['metrics'])
PYTHON_SCRIPT
echo "$(date)"
"""
    path = os.path.join(run_dir, "stageC_aggregate.sh")
    open(path, 'w').write(body)
    os.chmod(path, 0o755)
    return path


def main():
    p = argparse.ArgumentParser(description="Submit parallelized independent multi-output CV (job per fold x target)")
    p.add_argument('--input_strain_dir', required=True)
    p.add_argument('--phenotype_matrix', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--phenotype_column', required=True, help='Comma-separated target names')
    p.add_argument('--task_type', default='classification', choices=['classification', 'regression'])
    p.add_argument('--n_folds', type=int, default=10)
    p.add_argument('--cv_rounds', type=int, default=1)
    p.add_argument('--sample_column', default='phage')
    p.add_argument('--suffix', default='faa')
    p.add_argument('--k', type=int, default=5)
    p.add_argument('--num_runs_fs', type=int, default=25)
    p.add_argument('--num_runs_modeling', type=int, default=50)
    p.add_argument('--method', default='rfe')
    p.add_argument('--max_features', default='none')
    p.add_argument('--strong_top_frac', type=float, default=0.2)
    # Anti-overfitting: cluster-aware split + phylogenetic feature filtering
    p.add_argument('--use_clustering', action='store_true',
                   help='Cluster-aware train/test split (group related phages)')
    p.add_argument('--use_feature_clustering', action='store_true',
                   help='Remove phylogenetically-linked features (clade markers) at table build')
    p.add_argument('--feature_n_clusters', type=int, default=20)
    p.add_argument('--feature_min_cluster_presence', type=int, default=2)
    p.add_argument('--filter_by_cluster_presence', action='store_true',
                   help='Filter features by how many sample-clusters they appear in')
    p.add_argument('--min_cluster_presence', type=int, default=2)
    # Phylogenetic / PEQ leave-one-group-out folds (overrides random k-fold).
    p.add_argument('--group_metadata', default=None,
                   help='strain->group CSV (e.g. PEQ clades). When set, folds are '
                        'leave-one-group-out; --n_folds/--cv_rounds are ignored and the '
                        'fold count = number of groups (Stage A/B arrays size from it).')
    p.add_argument('--group_strain_column', default='strain')
    p.add_argument('--group_column', default='group')
    # SLURM
    p.add_argument('--account', default='pc_crispriart')
    p.add_argument('--partition', default='lr7')
    p.add_argument('--qos', default='lr_normal')
    p.add_argument('--environment', default='genophi')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--max_ram', type=int, default=100)
    p.add_argument('--mem_table', type=int, default=120, help='Mem per fold-table job (GB)')
    p.add_argument('--time_table', default='8:00:00')
    p.add_argument('--mem_train', type=int, default=120, help='Mem per (fold,target) train job (GB)')
    p.add_argument('--time_train', default='12:00:00')
    p.add_argument('--run_dir', default=None,
                   help='Where to write the run dir (scripts + logs). Default: '
                        'under --output_dir (on scratch), NOT the cwd/home.')
    p.add_argument('--dry_run', action='store_true')
    args = p.parse_args()

    targets = [c.strip() for c in args.phenotype_column.split(',') if c.strip()]
    if not targets:
        print("Error: no targets parsed"); return 1
    for pth, name in [(args.input_strain_dir, 'FASTA dir'), (args.phenotype_matrix, 'matrix')]:
        if not os.path.exists(pth):
            print(f"Error: {name} not found: {pth}"); return 1
    os.makedirs(args.output_dir, exist_ok=True)

    # Run dir (scripts + logs) defaults UNDER output_dir (scratch), so logs never
    # land in home/cwd by accident.
    base = args.run_dir if args.run_dir else args.output_dir
    run_dir = os.path.join(base, f"mo_cv_indep_run_{time.strftime('%Y%m%d_%H%M%S')}")
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)

    n_folds_total = _n_folds(args)
    n_train = n_folds_total * len(targets)
    print("=== Independent multi-output CV (parallel) ===")
    print(f"targets ({len(targets)}): {targets}")
    if args.group_metadata:
        print(f"folds: leave-one-group-out from {args.group_metadata} ({n_folds_total} groups)")
    print(f"Stage A: {n_folds_total} fold-table jobs")
    print(f"Stage B: {n_train} (fold x target) train jobs")
    print(f"Stage C: 1 aggregate job")
    print(f"output: {args.output_dir}")

    a = create_stage_a(args, run_dir, targets)
    b = create_stage_b(args, run_dir, targets, "STAGE_A_ID")
    c = create_stage_c(args, run_dir, targets, "STAGE_B_ID")

    if args.dry_run:
        print(f"\nDry run -- scripts in {run_dir}/ (not submitted)")
        return 0

    orig = os.getcwd(); os.chdir(run_dir)
    print("\nSubmitting Stage A (fold tables)...")
    a_id = submit_job("stageA_tables.sh"); print(f"  A job: {a_id}")
    if a_id:
        s = open("stageB_train.sh").read().replace("STAGE_A_ID", a_id)
        open("stageB_train.sh", 'w').write(s)
        print("Submitting Stage B (per-target, afterok A)...")
        b_id = submit_job("stageB_train.sh"); print(f"  B job: {b_id}")
        if b_id:
            s = open("stageC_aggregate.sh").read().replace("STAGE_B_ID", b_id)
            open("stageC_aggregate.sh", 'w').write(s)
            print("Submitting Stage C (aggregate, afterok B)...")
            c_id = submit_job("stageC_aggregate.sh"); print(f"  C job: {c_id}")
    os.chdir(orig)

    print(f"\nMonitor: squeue -u $USER ; tail -f {run_dir}/logs/*.out")
    print(f"Results: {args.output_dir}/cv_predictions.csv + {args.output_dir}/cv_performance/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
