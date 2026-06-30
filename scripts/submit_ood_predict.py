#!/usr/bin/env python3
"""
Predict OOD datasets with the combined all-4-genera GenoPHI model.

Applies the trained combined protein-family model (genophi_combined_full_model)
to out-of-distribution datasets (new host/phage genera) via the corner E1+E2
prediction pattern:

  E1. Assign the dataset's PHAGES against the combined model's phage cluster DB
      (run_assign_features_workflow, genome_type='phage') -> a phage feature
      table in the model's phage feature vocabulary.
  E2. Assign the dataset's HOSTS against the strain cluster DB and predict them
      paired with that phage table (assign_predict_workflow, genome_type='strain')
      -> predictions for host x phage.

The OOD genomes are QUERIES searched against the trained cluster DB, so the DB
need not contain them. duplicate_all=True prefixes protein IDs (genome::protein)
to match the model's namespace. The best cutoff is selected by MCC (genophi's
_select_best_cutoff: highest MCC, tie-break higher cut_off) -- NOT hardcoded.

Per dataset (one SLURM job each): predicts the full host x phage grid. Score
afterward by merging predictions against <dataset>_pairs.csv on (strain, phage).

  python submit_ood_predict.py                       # all configured datasets
  python submit_ood_predict.py --dataset salmonella_adler   # one dataset
"""
import os
import argparse
import subprocess

# ============================ CONFIG (edit me) ============================
# The trained combined model.
MODEL = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_combined_full_model"

# OOD root + the datasets that have BOTH host_faa and phage_faa.
OOD_ROOT = "/global/scratch/users/anoonan/set_transformer/manuscript/OOD"
DATASETS = ["cellulophaga", "enterococcus", "salmonella_adler"]

# Per-dataset layout (relative to OOD_ROOT/<dataset>/).
host_faa_subdir = "host_faa"
phage_faa_subdir = "phage_faa"

# MMseqs assignment params -- MUST match how the combined model was clustered.
min_seq_id = 0.2
coverage = 0.8
sensitivity = 7.5
duplicate_all = True
threads = 32

# SLURM
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"
mem_per_job = "120G"      # prediction is light vs training
time_limit = "8:00:00"
# =========================================================================

base_output_dir = os.path.join(OOD_ROOT, "genophi_combined_predictions")
logs_dir = os.path.join(base_output_dir, "slurm_logs")


def header(f, job):
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name={job}\n")
    f.write(f"#SBATCH --account={account}\n#SBATCH --partition={partition}\n#SBATCH --qos={qos}\n")
    f.write("#SBATCH --nodes=1\n#SBATCH --ntasks=1\n")
    f.write(f"#SBATCH --cpus-per-task={threads}\n#SBATCH --mem={mem_per_job}\n#SBATCH --time={time_limit}\n")
    f.write(f"#SBATCH --output={logs_dir}/{job}_%j.out\n#SBATCH --error={logs_dir}/{job}_%j.err\n")
    f.write("\nmodule load miniforge3\n")
    f.write('eval "$(conda shell.bash hook)"\n')
    f.write(f"conda activate {environment}\n\n")


def predict_py(dataset):
    host_dir = os.path.join(OOD_ROOT, dataset, host_faa_subdir)
    phage_dir = os.path.join(OOD_ROOT, dataset, phage_faa_subdir)
    od = os.path.join(base_output_dir, dataset)
    return (
        "python -c \""
        "import os, pandas as pd; "
        "from genophi.workflows.assign_features_workflow import run_assign_features_workflow; "
        "from genophi.workflows.assign_predict_workflow import assign_predict_workflow; "
        f"M='{MODEL}'; od='{od}'; os.makedirs(od, exist_ok=True); "
        # best cutoff by MCC (genophi _select_best_cutoff logic)
        "mdf=pd.read_csv(os.path.join(M,'modeling_results','model_performance','model_performance_metrics.csv')); "
        "mdf=mdf.sort_values(['MCC','cut_off'],ascending=[False,False]); "
        "cut=mdf['cut_off'].values[0]; print('best cutoff (MCC):', cut); "
        "model_dir=os.path.join(M,'modeling_results',str(cut)); "
        "select_ft=os.path.join(M,'feature_selection','filtered_feature_tables','select_feature_table_'+str(cut)+'.csv'); "
        # E1: assign OOD phages -> phage feature table in the model vocabulary
        "vpdir=os.path.join(od,'phage_features'); os.makedirs(vpdir,exist_ok=True); "
        "run_assign_features_workflow("
        f"input_dir='{phage_dir}', "
        f"mmseqs_db=os.path.join(M,'tmp','phage','mmseqs_db'), "
        "tmp_dir=os.path.join(vpdir,'tmp'), output_dir=vpdir, "
        "feature_map=os.path.join(M,'phage','features','selected_features.csv'), "
        "clusters_tsv=os.path.join(M,'phage','clusters.tsv'), "
        "genome_type='phage', "
        f"sensitivity={sensitivity}, coverage={coverage}, min_seq_id={min_seq_id}, "
        f"threads={threads}, duplicate_all={duplicate_all}); "
        "vpt=os.path.join(vpdir,'phage_combined_feature_table.csv'); "
        # E2: assign OOD hosts, predict host x phage
        "vodir=os.path.join(od,'model_validation'); os.makedirs(vodir,exist_ok=True); "
        "assign_predict_workflow("
        f"input_dir='{host_dir}', "
        f"mmseqs_db=os.path.join(M,'tmp','strain','mmseqs_db'), "
        "clusters_tsv=os.path.join(M,'strain','clusters.tsv'), "
        "feature_map=os.path.join(M,'strain','features','selected_features.csv'), "
        "tmp_dir=os.path.join(vodir,'tmp'), suffix='faa', model_dir=model_dir, "
        "feature_table=select_ft, phage_feature_table_path=vpt, output_dir=vodir, "
        f"threads={threads}, genome_type='strain', sensitivity={sensitivity}, "
        f"coverage={coverage}, min_seq_id={min_seq_id}, duplicate_all={duplicate_all}); "
        "print('DONE', od)"
        "\"\n"
    )


def submit(sh):
    return subprocess.check_output(["sbatch", "--parsable", sh]).decode().strip()


def write_job(dataset):
    j = f"OOD_{dataset}"
    w = os.path.join(logs_dir, f"{j}.sh")
    with open(w, 'w') as f:
        header(f, j)
        f.write(f'echo "OOD predict: {dataset}"\n\n')
        f.write(predict_py(dataset))
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', help="Run one dataset. Default: all.")
    args = ap.parse_args()

    os.makedirs(logs_dir, exist_ok=True)
    targets = [args.dataset] if args.dataset else DATASETS
    for ds in targets:
        if ds not in DATASETS:
            raise SystemExit(f"--dataset '{ds}' not in {DATASETS}")
        # Sanity: required inputs + model artifacts exist.
        host_dir = os.path.join(OOD_ROOT, ds, host_faa_subdir)
        phage_dir = os.path.join(OOD_ROOT, ds, phage_faa_subdir)
        for p in (host_dir, phage_dir,
                  os.path.join(MODEL, 'strain', 'features', 'selected_features.csv'),
                  os.path.join(MODEL, 'phage', 'features', 'selected_features.csv'),
                  os.path.join(MODEL, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')):
            if not os.path.exists(p):
                print(f"  WARNING [{ds}]: missing {p}")
        jid = submit(write_job(ds))
        print(f"  {ds}: {jid} -> {os.path.join(base_output_dir, ds)}")

    print(f"\nPredictions: {base_output_dir}/<dataset>/model_validation/predict_results/")
    print("Score: merge strain_median_predictions.csv against <dataset>_pairs.csv on (strain,phage).")
    print("Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
