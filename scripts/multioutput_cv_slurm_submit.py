#!/usr/bin/env python3
"""
SLURM submission for multi-output cross-validation (one job per fold).

Two stages:
  Stage 1: fold array (1 job per fold) — each runs genophi's run_one_cv_fold:
           builds a fold-specific k-mer feature table on the MODELING genomes,
           runs ensemble FS + modeling, then assigns + ensemble-predicts the
           held-out genomes. Folds are independent (deterministic splits).
  Stage 2: aggregation (afterok Stage 1) — pools per-fold held-out predictions
           and scores per target. ALL-OR-NOTHING: if any fold job failed, the
           aggregation errors; re-run the missing fold(s) then resubmit Stage 2.

Mirrors rbp_slurm_submit.py conventions (account/partition/qos, miniforge3 +
conda, job array, afterok dependency, sbatch --parsable).
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


# Shared kwargs string passed to both run_one_cv_fold and aggregate_cv_results,
# so the deterministic splits + targets match exactly across stages.
def _common_kwargs_literal(args, targets):
    return (
        f"input_strain_dir={args.input_strain_dir!r}, "
        f"phenotype_matrix={args.phenotype_matrix!r}, "
        f"output_dir={args.output_dir!r}, "
        f"phenotype_column={targets!r}, "
        f"task_type={args.task_type!r}, "
        f"target_mode={args.target_mode!r}, "
        f"strategy={args.strategy!r}, "
        f"n_folds={args.n_folds}, cv_rounds={args.cv_rounds}, "
        f"sample_column={args.sample_column!r}, suffix={args.suffix!r}, k={args.k}, "
        f"num_runs_fs={args.num_runs_fs}, num_runs_modeling={args.num_runs_modeling}, "
        f"method={args.method!r}, max_features={args.max_features!r}, "
        f"threads={args.threads}, max_ram={args.max_ram}, "
        f"use_clustering={args.use_clustering}, "
        f"use_feature_clustering={args.use_feature_clustering}, "
        f"feature_n_clusters={args.feature_n_clusters}, "
        f"feature_min_cluster_presence={args.feature_min_cluster_presence}, "
        f"filter_by_cluster_presence={args.filter_by_cluster_presence}, "
        f"min_cluster_presence={args.min_cluster_presence}, "
        f"group_metadata={args.group_metadata!r}, "
        f"group_strain_column={args.group_strain_column!r}, "
        f"group_column={args.group_column!r}"
    )


def _n_fold_jobs(args):
    """Fold-array size: group count for phylo folds, else n_folds * cv_rounds."""
    if args.group_metadata:
        import pandas as pd
        g = pd.read_csv(args.group_metadata)
        col = args.group_column
        groups = g[col].dropna()
        groups = groups[groups.astype(str).str.strip() != '']
        return int(groups.nunique())
    return args.n_folds * args.cv_rounds


def create_fold_stage(args, run_dir, targets):
    n_jobs = _n_fold_jobs(args)
    common = _common_kwargs_literal(args, targets)
    script = f"""#!/bin/bash
#SBATCH --job-name=mo_cv_fold
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_fold}G
#SBATCH --time={args.time_fold}
#SBATCH --array=1-{n_jobs}
#SBATCH --output=logs/fold_%A_%a.out
#SBATCH --error=logs/fold_%A_%a.err

echo "=== Multi-output CV fold $SLURM_ARRAY_TASK_ID ==="
echo "Started: $(date)"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate {args.environment}

python3 << 'PYTHON_SCRIPT'
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
import os
from genophi.workflows.multioutput_cv_workflow import run_one_cv_fold

fold_idx = int(os.environ['SLURM_ARRAY_TASK_ID'])
run_one_cv_fold(fold_idx=fold_idx, {common})
print(f"Fold {{fold_idx}} complete.")
PYTHON_SCRIPT

echo "Completed: $(date)"
"""
    path = os.path.join(run_dir, "stage1_cv_folds.sh")
    with open(path, 'w') as f:
        f.write(script)
    os.chmod(path, 0o755)
    return path


def create_aggregate_stage(args, run_dir, targets, dependency):
    # aggregate only needs output_dir + the split-defining args + scoring args.
    agg = (
        f"output_dir={args.output_dir!r}, "
        f"phenotype_column={targets!r}, task_type={args.task_type!r}, "
        f"n_folds={args.n_folds}, cv_rounds={args.cv_rounds}, "
        f"sample_column={args.sample_column!r}, strong_top_frac={args.strong_top_frac}, "
        f"group_metadata={args.group_metadata!r}, "
        f"input_strain_dir={args.input_strain_dir!r}, "
        f"phenotype_matrix={args.phenotype_matrix!r}, suffix={args.suffix!r}, "
        f"group_strain_column={args.group_strain_column!r}, "
        f"group_column={args.group_column!r}"
    )
    script = f"""#!/bin/bash
#SBATCH --job-name=mo_cv_agg
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --dependency=afterok:{dependency}
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err

echo "=== Multi-output CV aggregation ==="; echo "Started: $(date)"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate {args.environment}

python3 << 'PYTHON_SCRIPT'
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
from genophi.workflows.multioutput_cv_workflow import aggregate_cv_results
res = aggregate_cv_results({agg})
print("CV predictions:", res['predictions'])
print("CV metrics:", res['metrics'])
PYTHON_SCRIPT

echo "Completed: $(date)"
"""
    path = os.path.join(run_dir, "stage2_aggregate.sh")
    with open(path, 'w') as f:
        f.write(script)
    os.chmod(path, 0o755)
    return path


def main():
    p = argparse.ArgumentParser(description="Submit multi-output CV as SLURM jobs (one job per fold)")
    # Required inputs
    p.add_argument('--input_strain_dir', required=True, help='Per-genome FASTA dir')
    p.add_argument('--phenotype_matrix', required=True, help='Phenotype matrix CSV (sample col + targets)')
    p.add_argument('--output_dir', required=True)
    p.add_argument('--phenotype_column', required=True, help='Comma-separated target column names')
    # Task / strategy
    p.add_argument('--task_type', default='classification', choices=['classification', 'regression'])
    p.add_argument('--target_mode', default='multilabel',
                   choices=['auto', 'binary', 'multiclass', 'multilabel', 'single', 'multitarget'])
    p.add_argument('--strategy', default='joint', choices=['joint', 'independent'])
    # CV
    p.add_argument('--n_folds', type=int, default=10)
    p.add_argument('--cv_rounds', type=int, default=1)
    p.add_argument('--sample_column', default='phage')
    p.add_argument('--suffix', default='faa')
    # Modeling pass-through
    p.add_argument('--k', type=int, default=5)
    p.add_argument('--num_runs_fs', type=int, default=25)
    p.add_argument('--num_runs_modeling', type=int, default=50)
    p.add_argument('--method', default='rfe')
    p.add_argument('--max_features', default='none')
    p.add_argument('--strong_top_frac', type=float, default=0.2)
    # SLURM
    p.add_argument('--account', default='pc_phiml')
    p.add_argument('--partition', default='lr7')
    p.add_argument('--qos', default='lr_normal')
    p.add_argument('--environment', default='genophi')
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--max_ram', type=int, default=60)
    p.add_argument('--mem_fold', type=int, default=64, help='Mem per fold job (GB)')
    p.add_argument('--time_fold', default='24:00:00')
    # Anti-overfitting: cluster-aware split + phylogenetic feature filtering
    p.add_argument('--use_clustering', action='store_true',
                   help='Cluster-aware train/test split (group related phages)')
    p.add_argument('--use_feature_clustering', action='store_true',
                   help='Remove phylogenetically-linked features (clade markers) at table build')
    p.add_argument('--feature_n_clusters', type=int, default=20)
    p.add_argument('--feature_min_cluster_presence', type=int, default=2,
                   help='Min feature-clusters a feature must span to be kept')
    p.add_argument('--filter_by_cluster_presence', action='store_true',
                   help='Filter features by how many sample-clusters they appear in')
    p.add_argument('--min_cluster_presence', type=int, default=2)
    # Phylogenetic / PEQ leave-one-group-out folds (overrides random k-fold).
    p.add_argument('--group_metadata', default=None,
                   help='strain->group CSV (e.g. PEQ clades). When set, folds are '
                        'leave-one-group-out (one fold per group); --n_folds/--cv_rounds '
                        'are ignored and the fold count = number of groups.')
    p.add_argument('--group_strain_column', default='strain',
                   help='Sample-name column in --group_metadata (default: strain)')
    p.add_argument('--group_column', default='group',
                   help='Group-label column in --group_metadata (default: group)')
    p.add_argument('--run_dir', default=None,
                   help='Where to write the run dir (scripts + logs). Default: '
                        'under --output_dir (on scratch), NOT the cwd/home.')
    p.add_argument('--dry_run', action='store_true')
    args = p.parse_args()

    targets = [c.strip() for c in args.phenotype_column.split(',') if c.strip()]
    if not targets:
        print("Error: no targets parsed from --phenotype_column"); return 1

    for pth, name in [(args.input_strain_dir, 'FASTA dir'), (args.phenotype_matrix, 'phenotype matrix')]:
        if not os.path.exists(pth):
            print(f"Error: {name} not found: {pth}"); return 1
    os.makedirs(args.output_dir, exist_ok=True)

    # Run dir (scripts + logs) defaults UNDER output_dir (scratch), so logs never
    # land in home/cwd by accident. SLURM --output paths are relative to the
    # submit cwd, so we chdir into run_dir before sbatch (done below).
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    base = args.run_dir if args.run_dir else args.output_dir
    run_dir = os.path.join(base, f"mo_cv_run_{timestamp}")
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)

    n_jobs = _n_fold_jobs(args)
    print(f"=== Multi-output CV SLURM submission ===")
    print(f"targets ({len(targets)}): {targets}")
    print(f"strategy={args.strategy} target_mode={args.target_mode} task={args.task_type}")
    if args.group_metadata:
        print(f"folds: leave-one-group-out from {args.group_metadata} "
              f"= {n_jobs} fold jobs + 1 aggregate")
    else:
        print(f"folds: {args.n_folds} x {args.cv_rounds} rounds = {n_jobs} fold jobs + 1 aggregate")
    print(f"output: {args.output_dir}")

    fold_script = create_fold_stage(args, run_dir, targets)
    agg_script = create_aggregate_stage(args, run_dir, targets, "FOLD_JOB_ID")

    if args.dry_run:
        print(f"\nDry run — scripts in {run_dir}/ (not submitted)")
        return 0

    original = os.getcwd()
    os.chdir(run_dir)
    print("\nSubmitting Stage 1 (fold array)...")
    fold_id = submit_job("stage1_cv_folds.sh")
    print(f"  Job ID: {fold_id}")
    if fold_id:
        with open("stage2_aggregate.sh") as f:
            content = f.read()
        with open("stage2_aggregate.sh", 'w') as f:
            f.write(content.replace("FOLD_JOB_ID", fold_id))
        print("Submitting Stage 2 (aggregate, afterok)...")
        agg_id = submit_job("stage2_aggregate.sh")
        print(f"  Job ID: {agg_id}")
    os.chdir(original)

    print(f"\nMonitor: squeue -u $USER ; tail -f {run_dir}/logs/*.out")
    print(f"Results: {args.output_dir}/cv_predictions.csv and {args.output_dir}/cv_performance/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
