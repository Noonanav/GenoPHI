#!/usr/bin/env python3
"""
Submit the OOD arm of the k-mer comparison: train once, then predict.

Two stages, submitted separately on purpose:

  train    ONE all-4-genera k-mer model (kmer_train_all4.py). This is the
           largest job in the k-mer set -- every LOGO fold trains on three
           genera, this one trains on four -- so it may hit the 72 h wall
           clock. It resumes: run_kmer_table_workflow checkpoints its feature
           tables, and feature selection skips any run_N that already has
           feature_importances.csv. Just submit `train` again.

  predict  One light job per OOD dataset (kmer_predict_ood.py), applying that
           model to a dataset's hosts x phages. No training, so these are
           quick.

They are NOT chained with afterok. A dependency would be satisfied only by a
clean first run, and the training job is likely to time out at least once --
the dependency would then go DependencyNeverSatisfied and the predictions would
need resubmitting by hand anyway. Run `predict` once `train` genuinely
finishes; it refuses to submit if the model is incomplete.

Scoring: predictions cover the full host x phage grid. Merge
strain_median_predictions.csv against each dataset's <dataset>_pairs.csv on
(strain, phage), as with the protein-family OOD run.

  python submit_kmer_ood_slurm.py train
  python submit_kmer_ood_slurm.py predict
  python submit_kmer_ood_slurm.py predict --dataset cellulophaga
  python submit_kmer_ood_slurm.py status
"""
import os
import argparse
import subprocess

# ============================ CONFIG (edit me) ============================
# Where the all-4-genera k-mer model is written.
ALL4_DIR = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_kmer_all4"

# Training inputs (defaults of kmer_train_all4.py; repeated here so the
# submitted command is explicit in the log).
MATRIX = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
TRAIN_STRAIN_DIR = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
TRAIN_PHAGE_DIR = "/global/home/groups/pc_phiml/data/combined/phage_AAs"

# OOD datasets: <OOD_ROOT>/<dataset>/{host_faa,phage_faa}/
OOD_ROOT = "/global/scratch/users/anoonan/set_transformer/manuscript/OOD"
DATASETS = ["cellulophaga", "enterococcus", "salmonella_adler"]
HOST_SUBDIR = "host_faa"
PHAGE_SUBDIR = "phage_faa"
PRED_ROOT = os.path.join(OOD_ROOT, "kmer_predictions")

# k-mer + modeling params (match the LOGO and in-distribution k-mer arms).
K = 4
num_runs_fs = 25
num_runs_modeling = 50
num_features = 100
max_ram = 40

# Where the scripts live.
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))

# SLURM
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"

train_threads = 32
train_mem = "490G"
train_time = "72:00:00"

predict_threads = 8
predict_mem = "120G"
predict_time = "8:00:00"
# =========================================================================


def header(f, job, logs_dir, threads, mem, time_lim):
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name={job}\n")
    f.write(f"#SBATCH --account={account}\n#SBATCH --partition={partition}\n#SBATCH --qos={qos}\n")
    f.write("#SBATCH --nodes=1\n#SBATCH --ntasks=1\n")
    f.write(f"#SBATCH --cpus-per-task={threads}\n#SBATCH --mem={mem}\n#SBATCH --time={time_lim}\n")
    f.write(f"#SBATCH --output={logs_dir}/{job}_%j.out\n#SBATCH --error={logs_dir}/{job}_%j.err\n")
    f.write("\nmodule load miniforge3\n")
    f.write('eval "$(conda shell.bash hook)"\n')
    f.write(f"conda activate {environment}\n\n")


def submit(sh_path):
    return subprocess.check_output(["sbatch", "--parsable", sh_path]).decode().strip()


def model_is_complete(model_dir):
    """True once training has written its modeling metrics.

    That file is the last thing run_kmer_table_workflow produces, so its
    presence means feature selection AND modeling finished -- not merely that
    the directory exists.
    """
    return os.path.exists(os.path.join(
        model_dir, 'modeling', 'modeling_results', 'model_performance',
        'model_performance_metrics.csv'))


def training_progress(model_dir):
    """(completed_fs_runs, num_runs_fs) for a partially-trained model."""
    fs_dir = os.path.join(model_dir, 'modeling', 'feature_selection')
    if not os.path.isdir(fs_dir):
        return 0, num_runs_fs
    done = sum(
        1 for d in os.listdir(fs_dir)
        if d.startswith('run_') and os.path.exists(
            os.path.join(fs_dir, d, 'feature_importances.csv'))
    )
    return done, num_runs_fs


def cmd_train(args):
    logs_dir = os.path.join(ALL4_DIR, "slurm_logs")
    os.makedirs(logs_dir, exist_ok=True)

    if model_is_complete(ALL4_DIR) and not args.force:
        raise SystemExit(
            f"Model already complete at {ALL4_DIR}. Use --force to retrain, "
            f"or run: {os.path.basename(__file__)} predict")

    done, total = training_progress(ALL4_DIR)
    if done:
        print(f"resuming: {done}/{total} feature-selection runs already done")

    w = os.path.join(logs_dir, "KMERALL4_train.sh")
    with open(w, 'w') as f:
        header(f, "KMERALL4", logs_dir, train_threads, train_mem, train_time)
        f.write(
            f"python {os.path.join(SCRIPTS_DIR, 'kmer_train_all4.py')} \\\n"
            f"    --matrix {MATRIX} \\\n"
            f"    --strain-dir {TRAIN_STRAIN_DIR} \\\n"
            f"    --phage-dir {TRAIN_PHAGE_DIR} \\\n"
            f"    --output {ALL4_DIR} \\\n"
            f"    --k {K} \\\n"
            f"    --num-runs-fs {num_runs_fs} \\\n"
            f"    --num-runs-modeling {num_runs_modeling} \\\n"
            f"    --num-features {num_features} \\\n"
            f"    --max-ram {max_ram} \\\n"
            f"    --threads {train_threads}\n"
        )
    jid = submit(w)
    print(f"training: {jid}  -> {ALL4_DIR}")
    print("This job is expected to need more than one submission. When it "
          "times out, run `train` again -- it resumes.")


def cmd_predict(args):
    if not model_is_complete(ALL4_DIR):
        done, total = training_progress(ALL4_DIR)
        raise SystemExit(
            f"Model at {ALL4_DIR} is not finished ({done}/{total} "
            f"feature-selection runs done; modeling metrics absent). "
            f"Run `train` to completion first.")

    datasets = [args.dataset] if args.dataset else DATASETS
    unknown = [d for d in datasets if d not in DATASETS]
    if unknown:
        raise SystemExit(f"Unknown dataset(s): {unknown}. Known: {DATASETS}")

    logs_dir = os.path.join(PRED_ROOT, "slurm_logs")
    os.makedirs(logs_dir, exist_ok=True)

    for dataset in datasets:
        host_dir = os.path.join(OOD_ROOT, dataset, HOST_SUBDIR)
        phage_dir = os.path.join(OOD_ROOT, dataset, PHAGE_SUBDIR)
        for d in (host_dir, phage_dir):
            if not os.path.isdir(d):
                raise SystemExit(f"Missing genome directory: {d}")

        out_dir = os.path.join(PRED_ROOT, dataset)
        os.makedirs(out_dir, exist_ok=True)

        w = os.path.join(logs_dir, f"KMERPRED_{dataset}.sh")
        with open(w, 'w') as f:
            header(f, f"KMERPRED_{dataset}", logs_dir,
                   predict_threads, predict_mem, predict_time)
            f.write(
                f"python {os.path.join(SCRIPTS_DIR, 'kmer_predict_ood.py')} \\\n"
                f"    --model {ALL4_DIR} \\\n"
                f"    --host-dir {host_dir} \\\n"
                f"    --phage-dir {phage_dir} \\\n"
                f"    --output {out_dir} \\\n"
                f"    --threads {predict_threads}\n"
            )
        jid = submit(w)
        print(f"predict {dataset}: {jid}  -> {out_dir}")

    print("\nScore by merging predict_results/strain_median_predictions.csv "
          "against each <dataset>_pairs.csv on (strain, phage).")


def cmd_status(args):
    print(f"model: {ALL4_DIR}")
    if model_is_complete(ALL4_DIR):
        print("  training COMPLETE -- ready to predict")
    else:
        done, total = training_progress(ALL4_DIR)
        print(f"  training incomplete: {done}/{total} feature-selection runs")
    print(f"\npredictions: {PRED_ROOT}")
    for dataset in DATASETS:
        pred = os.path.join(PRED_ROOT, dataset, 'predict_results',
                            'strain_median_predictions.csv')
        print(f"  {dataset:20s} {'done' if os.path.exists(pred) else '-'}")


def main():
    ap = argparse.ArgumentParser(
        description="Submit the OOD k-mer arm: train one all-4 model, then predict.")
    # dest= without required= so this also runs under python 3.6, which the
    # login shell has; `required` on add_subparsers landed in 3.7.
    sub = ap.add_subparsers(dest='command')

    p_train = sub.add_parser('train', help="Submit the all-4-genera training job.")
    p_train.add_argument('--force', action='store_true',
                         help="Submit even if the model looks complete.")
    p_train.set_defaults(func=cmd_train)

    p_pred = sub.add_parser('predict', help="Submit prediction jobs (needs a trained model).")
    p_pred.add_argument('--dataset', help=f"Just one of: {', '.join(DATASETS)}")
    p_pred.set_defaults(func=cmd_predict)

    p_stat = sub.add_parser('status', help="Show training and prediction progress.")
    p_stat.set_defaults(func=cmd_status)

    args = ap.parse_args()
    if not getattr(args, 'command', None):
        ap.print_help()
        raise SystemExit(1)
    args.func(args)


if __name__ == "__main__":
    main()
