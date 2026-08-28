#!/usr/bin/env python3
"""Steps 2 and 3 of the GenoPHI multi-output workflow: cross-validate, then train.

Two independent choices:

  --mode      cv     Held-out performance -- the numbers you report.
              final  One model per target trained on ALL samples, for predicting
                     new genomes. No held-out set, so it yields no honest metric.

  --executor  slurm  Submit to a cluster (default). Needs --account.
              local  Run here, in this process, one fold at a time.

Both executors call the SAME functions, so a local run and a cluster run of the
same configuration produce the same results -- local is just serial.

Local is fine for small datasets and for smoke tests. Cross-validation is
n_folds x n_targets independent model fits: a 10-fold, 2-target run is 20, and a
10-fold, 5-target run on 69 genomes previously took about 17 hours serially.
Above roughly 50 samples or 3 targets, use a cluster.

Both modes read the directory written by prepare_inputs.py and refuse to run if
that step did not finish clean.

Why a wrapper rather than calling GenoPHI directly:
  * the two CV submitters take different flags (the joint one accepts
    --target_mode and rejects --min_features; the independent one is the
    reverse), and `genophi multioutput-cv` accepts neither --min_features nor
    --group_metadata, so the CLI cannot express a small-n or phylo run at all;
  * `genophi kmer-workflow` defaults to use_clustering=True while every CV path
    defaults it to False -- training a final model with CLI defaults silently
    produces a model that does not match the CV you reported.

Example
-------
    python submit_run.py --prepared_dir prepared/msp_binary \
        --mode cv --executor local --output_dir runs/msp_cv
"""

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SUBMITTERS = {
    'independent': os.path.join(REPO, 'scripts',
                                'multioutput_cv_independent_slurm_submit.py'),
    'joint': os.path.join(REPO, 'scripts', 'multioutput_cv_slurm_submit.py'),
}

ACCOUNT_HELP = """
--account is required for --executor slurm and has no default: SLURM
allocations are per-project and you can only charge one you belong to.

To find yours:

    sacctmgr show associations user=$USER format=Account,Partition,QOS

or on Lawrencium:

    lrc-associations -u $USER

Then pass, for example:

    --account <your_project> --partition lr7 --qos lr_normal

If you have no allocation, ask the project PI to add you, or use
--executor local.
"""


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


def n_fold_jobs(args):
    """Fold count: number of groups for phylo folds, else n_folds * cv_rounds."""
    if args.group_metadata:
        import pandas as pd
        g = pd.read_csv(args.group_metadata)
        groups = g[args.group_column].dropna()
        groups = groups[groups.astype(str).str.strip() != '']
        return int(groups.nunique())
    return args.n_folds * args.cv_rounds


# --------------------------------------------------------------------------
# SLURM
# --------------------------------------------------------------------------

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


FINAL_HEADER = """#!/bin/bash
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

"""


def build_final_command(args, pheno, genomes):
    """The `genophi kmer-workflow` invocation, as an argv list."""
    cmd = [
        'genophi', 'kmer-workflow',
        '--input_strain_dir', genomes,
        '--phenotype_matrix', pheno,
        '--output', args.output_dir,
        '--sample_column', args.sample_column,
        '--phenotype_columns_file', os.path.join(args.prepared_dir, 'targets.txt'),
        '--strategy', args.strategy,
        '--task_type', args.task_type,
        '--suffix', args.suffix,
        '--k', str(args.k),
        '--num_runs_fs', str(args.num_runs_fs),
        '--num_runs_modeling', str(args.num_runs_modeling),
        '--method', args.method,
        '--max_features', args.max_features,
        '--min_features', args.min_features,
        '--threads', str(args.threads),
        '--max_ram', str(args.max_ram),
    ]
    # kmer-workflow enables clustering by default; every CV path disables it.
    # Keep them consistent unless the user explicitly opts in.
    if not args.use_clustering:
        cmd += ['--no-clustering']
    if args.use_feature_clustering:
        cmd += ['--use_feature_clustering',
                '--feature_n_clusters', str(args.feature_n_clusters),
                '--feature_min_cluster_presence',
                str(args.feature_min_cluster_presence)]
    if args.strategy == 'joint':
        cmd += ['--target_mode', args.target_mode]
    return cmd


def build_final_script(args, cmd, run_dir):
    """Wrap the training command in a single SLURM job script."""
    body = FINAL_HEADER.format(
        account=args.account, partition=args.partition, qos=args.qos,
        threads=args.threads, mem=args.mem, time=args.time_train,
        environment=args.environment, run_dir=run_dir)
    body += ' \\\n  '.join([cmd[0] + ' ' + cmd[1]] + _pairs(cmd[2:]))
    body += '\nrc=$?\necho "exit=$rc"\ndate\nexit $rc\n'

    path = os.path.join(run_dir, 'final_train.sh')
    with open(path, 'w') as fh:
        fh.write(body)
    os.chmod(path, 0o755)
    return path


def _pairs(rest):
    """Render argv tail as '--flag value' strings for a readable shell script."""
    out, i = [], 0
    while i < len(rest):
        if i + 1 < len(rest) and not rest[i + 1].startswith('--'):
            out.append(f'{rest[i]} {rest[i + 1]}')
            i += 2
        else:
            out.append(rest[i])
            i += 1
    return out


# --------------------------------------------------------------------------
# Local execution -- same functions the SLURM stages call, run serially
# --------------------------------------------------------------------------

def run_cv_local(args, pheno, genomes, targets):
    """Run the whole CV here, one fold (and target) at a time."""
    n = n_fold_jobs(args)
    units = n * (len(targets) if args.strategy == 'independent' else 1)
    print(f'Local CV: {n} folds, {len(targets)} target(s) '
          f'-> {units} model fit(s), serial.')
    print('This is the same code the cluster runs; it is just not parallel.')
    print('Interrupting is safe -- completed folds are skipped on re-run.\n')

    common = dict(
        input_strain_dir=genomes, phenotype_matrix=pheno,
        output_dir=args.output_dir, phenotype_column=targets,
        task_type=args.task_type, n_folds=args.n_folds,
        cv_rounds=args.cv_rounds, sample_column=args.sample_column,
        suffix=args.suffix, group_metadata=args.group_metadata,
        group_strain_column=args.group_strain_column,
        group_column=args.group_column)

    t0 = time.time()
    if args.strategy == 'independent':
        from genophi.workflows.multioutput_cv_parallel import (
            build_fold_table, train_fold_target, aggregate_independent_cv)
        for fold in range(1, n + 1):
            print(f'--- fold {fold}/{n}: building feature table '
                  f'({time.time() - t0:.0f}s elapsed)')
            build_fold_table(
                fold_idx=fold, k=args.k, threads=args.threads,
                max_ram=args.max_ram,
                use_feature_clustering=args.use_feature_clustering,
                feature_n_clusters=args.feature_n_clusters,
                feature_min_cluster_presence=args.feature_min_cluster_presence,
                **common)
            for target in targets:
                print(f'--- fold {fold}/{n}, target {target} '
                      f'({time.time() - t0:.0f}s elapsed)')
                train_fold_target(
                    fold_idx=fold, target=target,
                    output_dir=args.output_dir, phenotype_column=targets,
                    task_type=args.task_type, n_folds=args.n_folds,
                    cv_rounds=args.cv_rounds, sample_column=args.sample_column,
                    num_runs_fs=args.num_runs_fs,
                    num_runs_modeling=args.num_runs_modeling,
                    method=args.method, max_features=args.max_features,
                    min_features=args.min_features, threads=args.threads,
                    max_ram=args.max_ram, use_clustering=args.use_clustering)
        print('--- aggregating')
        aggregate_independent_cv(threads=args.threads,
                                 strong_top_frac=args.strong_top_frac, **common)
    else:
        from genophi.workflows.multioutput_cv_workflow import (
            run_one_cv_fold, aggregate_cv_results)
        for fold in range(1, n + 1):
            print(f'--- fold {fold}/{n} ({time.time() - t0:.0f}s elapsed)')
            run_one_cv_fold(
                fold_idx=fold, target_mode=args.target_mode, strategy='joint',
                k=args.k, num_runs_fs=args.num_runs_fs,
                num_runs_modeling=args.num_runs_modeling, method=args.method,
                max_features=args.max_features, threads=args.threads,
                max_ram=args.max_ram, use_clustering=args.use_clustering,
                use_feature_clustering=args.use_feature_clustering,
                feature_n_clusters=args.feature_n_clusters,
                feature_min_cluster_presence=args.feature_min_cluster_presence,
                **common)
        print('--- aggregating')
        aggregate_cv_results(
            output_dir=args.output_dir, phenotype_column=targets,
            task_type=args.task_type, n_folds=args.n_folds,
            cv_rounds=args.cv_rounds, sample_column=args.sample_column,
            strong_top_frac=args.strong_top_frac,
            group_metadata=args.group_metadata, input_strain_dir=genomes,
            phenotype_matrix=pheno, suffix=args.suffix,
            group_strain_column=args.group_strain_column,
            group_column=args.group_column)

    print(f'\nDone in {(time.time() - t0) / 60:.1f} min.')
    print(f'Metrics: {os.path.join(args.output_dir, "cv_performance", "model_performance_metrics.csv")}')
    return 0


def main():
    ap = argparse.ArgumentParser(
        description='Run multi-output cross-validation or final model training, '
                    'locally or on a cluster.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--prepared_dir', required=True,
                    help='Directory written by prepare_inputs.py')
    ap.add_argument('--mode', required=True, choices=['cv', 'final'],
                    help="'cv' for held-out performance, 'final' for deployable models")
    ap.add_argument('--executor', default='slurm', choices=['slurm', 'local'],
                    help="'local' runs here serially; 'slurm' submits jobs")
    ap.add_argument('--output_dir', required=True,
                    help='Run output directory (on a cluster, use scratch not $HOME)')
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
    g.add_argument('--strong_top_frac', type=float, default=0.2,
                   help='Regression only: top fraction treated as strong '
                        'responders for detection metrics')
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
                        'GenoPHI defaults to 5, which silently drops small-n '
                        'folds; 3 is safe under ~150 samples. Independent '
                        'strategy only -- the joint fold runner has no such knob.')
    g.add_argument('--use_clustering', action='store_true',
                   help='Enable sample clustering. OFF by default so CV and final '
                        'training match; the published runs used no clustering.')
    g.add_argument('--use_feature_clustering', action='store_true',
                   help='Collapse phylogenetically linked features')
    g.add_argument('--feature_n_clusters', type=int, default=20)
    g.add_argument('--feature_min_cluster_presence', type=int, default=2)
    g.add_argument('--threads', type=int, default=16)
    g.add_argument('--max_ram', type=int, default=100, help='Max RAM (GB)')

    g = ap.add_argument_group('cluster (--executor slurm)')
    g.add_argument('--account', default=None,
                   help='SLURM allocation to charge. REQUIRED -- no default, '
                        'since allocations are per-project.')
    g.add_argument('--partition', default='lr7')
    g.add_argument('--qos', default='lr_normal')
    g.add_argument('--environment', default='genophi-mo',
                   help='Conda environment (name or full path) to activate in jobs')
    g.add_argument('--mem', type=int, default=120, help='Memory per job (GB)')
    g.add_argument('--time_table', default='8:00:00')
    g.add_argument('--time_train', default='24:00:00')
    g.add_argument('--dry_run', action='store_true',
                   help='Write the SLURM scripts but do not sbatch them')
    args = ap.parse_args()

    if args.executor == 'slurm' and not args.account:
        sys.exit(ACCOUNT_HELP)

    # --target_mode defaults to multilabel, which is a classification mode. For
    # regression the joint path needs multitarget (MultiRMSE); silently running
    # a regression as multilabel would be a confusing failure.
    if args.task_type == 'regression' and args.target_mode == 'multilabel':
        args.target_mode = 'multitarget'
        print("note: --task_type regression -> --target_mode multitarget\n")

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
    print(f'mode         : {args.mode}   strategy: {args.strategy}   '
          f'executor: {args.executor}')
    print(f'output_dir   : {args.output_dir}')
    print()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.mode == 'cv':
        if args.executor == 'local':
            return run_cv_local(args, pheno, genomes, targets)
        argv = build_cv_argv(args, pheno, genomes, targets)
        print('dispatching to:', os.path.basename(argv[1]))
        print('  ' + ' '.join(argv[2:]))
        print()
        sys.stdout.flush()
        return subprocess.call(argv)

    cmd = build_final_command(args, pheno, genomes)
    if args.executor == 'local':
        print('Final training (all samples, no held-out set):')
        print('  ' + ' '.join(cmd) + '\n')
        sys.stdout.flush()
        return subprocess.call(cmd)

    run_dir = os.path.join(args.output_dir, 'run_scripts')
    os.makedirs(os.path.join(run_dir, 'logs'), exist_ok=True)
    script = build_final_script(args, cmd, run_dir)
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
