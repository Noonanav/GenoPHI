#!/usr/bin/env python3
"""Steps 2 and 3 of the GenoPHI multi-output workflow: cross-validate, then train.

One entry point over two very different jobs:

  --mode cv     Held-out performance. Dispatches to the existing SLURM submitters
                (scripts/multioutput_cv_independent_slurm_submit.py for
                --strategy independent, scripts/multioutput_cv_slurm_submit.py
                for --strategy joint). Produces per-target metrics you can report.

  --mode final  One model per target trained on ALL samples, for predicting new
                genomes later. Emits a single SLURM job running `genophi
                kmer-workflow`. Produces models, not metrics -- there is no
                held-out set, so never quote numbers from this run.

Both modes read the directory written by prepare_inputs.py and refuse to run if
that step did not finish clean.

Why a wrapper rather than calling the submitters directly:
  * the two CV submitters take different flags (the joint one accepts
    --target_mode and rejects --min_features; the independent one is the
    reverse), which is easy to get wrong;
  * `genophi kmer-workflow` defaults to use_clustering=True while every CV path
    defaults it to False -- training a final model with CLI defaults silently
    produces a model that does not match the CV you reported. This wrapper keeps
    the two consistent by default.

Example
-------
    python submit_run.py --prepared_dir prepared/msp_binary \
        --mode cv --output_dir /global/scratch/users/me/msp_cv --dry_run
"""

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SUBMITTERS = {
    'independent': os.path.join(REPO, 'scripts',
                                'multioutput_cv_independent_slurm_submit.py'),
    'joint': os.path.join(REPO, 'scripts', 'multioutput_cv_slurm_submit.py'),
}


def load_prepared(prepared_dir, skip_check):
    """Read prepare_inputs.py outputs; refuse inputs that did not validate.

    Paths are absolutized: they get baked into SLURM scripts that run from a
    different working directory than the submitting shell.
    """
    prepared_dir = os.path.abspath(prepared_dir)
    pheno = os.path.join(prepared_dir, 'phenotype.csv')
    genomes = os.path.join(prepared_dir, 'genomes')
    targets_file = os.path.join(prepared_dir, 'targets.txt')
    report = os.path.join(prepared_dir, 'prepare_report.json')

    for p in (pheno, genomes, targets_file):
        if not os.path.exists(p):
            sys.exit(f"ERROR: {p} not found. Run prepare_inputs.py first.")

    if os.path.exists(report):
        with open(report) as fh:
            meta = json.load(fh)
        if meta.get('status') != 'ready' and not skip_check:
            sys.exit(f"ERROR: {report} says status={meta.get('status')!r}.\n"
                     f"       Errors: {meta.get('errors')}\n"
                     '       Fix them, or pass --skip_prepare_check to override.')
        for w in meta.get('warnings', []):
            print(f'  (prepare warning) {w}')
    elif not skip_check:
        sys.exit(f'ERROR: {report} not found. Re-run prepare_inputs.py.')

    with open(targets_file) as fh:
        targets = [ln.strip() for ln in fh if ln.strip()]

    return pheno, genomes, targets


def build_cv_argv(args, pheno, genomes, targets):
    """Assemble the argv for whichever CV submitter matches --strategy."""
    submitter = SUBMITTERS[args.strategy]
    if not os.path.exists(submitter):
        sys.exit(f'ERROR: submitter not found: {submitter}')

    argv = [
        sys.executable, submitter,
        '--input_strain_dir', genomes,
        '--phenotype_matrix', pheno,
        '--output_dir', args.output_dir,
        '--phenotype_column', ','.join(targets),
        '--task_type', args.task_type,
        '--n_folds', str(args.n_folds),
        '--cv_rounds', str(args.cv_rounds),
        '--sample_column', args.sample_column,
        '--suffix', args.suffix,
        '--k', str(args.k),
        '--num_runs_fs', str(args.num_runs_fs),
        '--num_runs_modeling', str(args.num_runs_modeling),
        '--method', args.method,
        '--max_features', args.max_features,
        '--account', args.account,
        '--partition', args.partition,
        '--qos', args.qos,
        '--environment', args.environment,
        '--threads', str(args.threads),
        '--max_ram', str(args.max_ram),
    ]

    # The two submitters diverge here -- this asymmetry is the main reason this
    # wrapper exists.
    if args.strategy == 'independent':
        argv += ['--min_features', args.min_features]
        argv += ['--mem_table', str(args.mem), '--time_table', args.time_table]
        argv += ['--mem_train', str(args.mem), '--time_train', args.time_train]
    else:
        argv += ['--strategy', 'joint', '--target_mode', args.target_mode]
        argv += ['--mem_fold', str(args.mem), '--time_fold', args.time_train]

    if args.group_metadata:
        argv += ['--group_metadata', args.group_metadata,
                 '--group_strain_column', args.group_strain_column,
                 '--group_column', args.group_column]
    if args.use_clustering:
        argv += ['--use_clustering']
    if args.use_feature_clustering:
        argv += ['--use_feature_clustering',
                 '--feature_n_clusters', str(args.feature_n_clusters),
                 '--feature_min_cluster_presence', str(args.feature_min_cluster_presence)]
    if args.dry_run:
        argv += ['--dry_run']
    return argv


FINAL_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=mo_final
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={mem}G
#SBATCH --time={time}
#SBATCH --output={run_dir}/logs/final_%j.out
#SBATCH --error={run_dir}/logs/final_%j.err

echo "=== Final model training (all samples, no held-out set) ==="
date
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate {environment}

genophi kmer-workflow \\
  --input_strain_dir {genomes} \\
  --phenotype_matrix {pheno} \\
  --output {output_dir} \\
  --sample_column {sample_column} \\
  --phenotype_columns_file {targets_file} \\
  --strategy {strategy} \\
  --task_type {task_type} \\
  --suffix {suffix} \\
  --k {k} \\
  --num_runs_fs {num_runs_fs} \\
  --num_runs_modeling {num_runs_modeling} \\
  --method {method} \\
  --max_features {max_features} \\
  --min_features {min_features} \\
  --threads {threads} \\
  --max_ram {max_ram}{extra}

rc=$?
echo "exit=$rc"
date
exit $rc
"""


def build_final_script(args, pheno, genomes, targets, run_dir):
    """Write the single-job SLURM script that trains the deployable models."""
    extra = ''
    # kmer-workflow enables clustering by default; every CV path disables it.
    # Keep them consistent unless the user explicitly opts in.
    if not args.use_clustering:
        extra += ' \\\n  --no-clustering'
    if args.use_feature_clustering:
        extra += (f' \\\n  --use_feature_clustering'
                  f' \\\n  --feature_n_clusters {args.feature_n_clusters}'
                  f' \\\n  --feature_min_cluster_presence '
                  f'{args.feature_min_cluster_presence}')
    if args.strategy == 'joint':
        extra += f' \\\n  --target_mode {args.target_mode}'

    body = FINAL_TEMPLATE.format(
        account=args.account, partition=args.partition, qos=args.qos,
        threads=args.threads, mem=args.mem, time=args.time_train,
        environment=args.environment, run_dir=run_dir,
        genomes=genomes, pheno=pheno, output_dir=args.output_dir,
        sample_column=args.sample_column,
        targets_file=os.path.join(args.prepared_dir, 'targets.txt'),
        strategy=args.strategy, task_type=args.task_type, suffix=args.suffix,
        k=args.k, num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling, method=args.method,
        max_features=args.max_features, min_features=args.min_features,
        max_ram=args.max_ram, extra=extra)

    path = os.path.join(run_dir, 'final_train.sh')
    with open(path, 'w') as fh:
        fh.write(body)
    os.chmod(path, 0o755)
    return path


def main():
    ap = argparse.ArgumentParser(
        description='Submit multi-output cross-validation or final model training.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--prepared_dir', required=True,
                    help='Directory written by prepare_inputs.py')
    ap.add_argument('--mode', required=True, choices=['cv', 'final'],
                    help="'cv' for held-out performance, 'final' for deployable models")
    ap.add_argument('--output_dir', required=True,
                    help='Run output directory (use cluster scratch, not $HOME)')
    ap.add_argument('--strategy', default='independent',
                    choices=['independent', 'joint'],
                    help='One model per target, or one shared multi-output model')
    ap.add_argument('--task_type', default='classification',
                    choices=['classification', 'regression'])
    ap.add_argument('--target_mode', default='multilabel',
                    choices=['auto', 'binary', 'multiclass', 'multilabel',
                             'single', 'multitarget'],
                    help='Joint strategy only; ignored for independent')
    ap.add_argument('--sample_column', default='phage')
    ap.add_argument('--suffix', default='faa')
    ap.add_argument('--skip_prepare_check', action='store_true',
                    help='Run even if prepare_inputs.py reported errors')

    g = ap.add_argument_group('cross-validation')
    g.add_argument('--n_folds', type=int, default=10)
    g.add_argument('--cv_rounds', type=int, default=1)
    g.add_argument('--group_metadata', default=None,
                   help='CSV of strain,group for leave-one-group-out (phylo) folds')
    g.add_argument('--group_strain_column', default='strain')
    g.add_argument('--group_column', default='group')

    g = ap.add_argument_group('modeling')
    g.add_argument('--k', type=int, default=5)
    g.add_argument('--num_runs_fs', type=int, default=25)
    g.add_argument('--num_runs_modeling', type=int, default=50)
    g.add_argument('--method', default='rfe')
    g.add_argument('--max_features', default='100')
    g.add_argument('--min_features', default='3',
                   help='Minimum features a cutoff table needs to be modeled. '
                        'The GenoPHI default of 5 silently drops small-n folds; '
                        '3 is the safe value for datasets under ~150 samples.')
    g.add_argument('--use_clustering', action='store_true',
                   help='Enable sample clustering. OFF by default so CV and final '
                        'training match; the published runs used no clustering.')
    g.add_argument('--use_feature_clustering', action='store_true',
                   help='Collapse phylogenetically linked features')
    g.add_argument('--feature_n_clusters', type=int, default=20)
    g.add_argument('--feature_min_cluster_presence', type=int, default=2)

    g = ap.add_argument_group('cluster')
    g.add_argument('--account', default='pc_phiml')
    g.add_argument('--partition', default='lr7')
    g.add_argument('--qos', default='lr_normal')
    g.add_argument('--environment', default='genophi-mo')
    g.add_argument('--threads', type=int, default=16)
    g.add_argument('--max_ram', type=int, default=100)
    g.add_argument('--mem', type=int, default=120, help='Memory per job (GB)')
    g.add_argument('--time_table', default='8:00:00')
    g.add_argument('--time_train', default='24:00:00')
    g.add_argument('--dry_run', action='store_true',
                   help='Write the SLURM scripts but do not sbatch them')
    args = ap.parse_args()

    args.prepared_dir = os.path.abspath(args.prepared_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    if args.group_metadata:
        args.group_metadata = os.path.abspath(args.group_metadata)
    pheno, genomes, targets = load_prepared(args.prepared_dir,
                                            args.skip_prepare_check)

    print(f'prepared_dir : {args.prepared_dir}')
    print(f'  phenotype  : {pheno}')
    print(f'  genomes    : {genomes} ({len(os.listdir(genomes))} files)')
    print(f'  targets    : {targets}')
    print(f'mode         : {args.mode}   strategy: {args.strategy}')
    print(f'output_dir   : {args.output_dir}')
    print()

    if args.mode == 'cv':
        argv = build_cv_argv(args, pheno, genomes, targets)
        print('dispatching to:', os.path.basename(argv[1]))
        print('  ' + ' '.join(argv[2:]))
        print()
        sys.stdout.flush()
        return subprocess.call(argv)

    run_dir = os.path.join(args.output_dir, 'run_scripts')
    os.makedirs(os.path.join(run_dir, 'logs'), exist_ok=True)
    script = build_final_script(args, pheno, genomes, targets, run_dir)
    print(f'wrote {script}')
    if args.dry_run:
        print('(dry run -- not submitted)')
        return 0
    r = subprocess.run(['sbatch', '--parsable', script],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f'sbatch failed: {r.stderr}')
        return r.returncode
    print(f'submitted job {r.stdout.strip()}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
