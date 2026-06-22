#!/usr/bin/env python3
"""
Re-predict the CORNER on the FULL held-out genus using SUBSET-TRAINED models.

The subset corner scan (genophi_corner_logo/msi02_c08) trained each genus fold on
the subset and predicted on the SUBSET held-out genus. This script reuses those
trained models + selected features, but re-runs ONLY the prediction step (corner
E1+E2) against the COMPLETE held-out genus's strains x phages.

It does NOT retrain and does NOT re-cluster. assign_predict_workflow searches the
new (full) held-out genomes' proteins against the existing shared cluster DB as a
target, so the DB does not need to contain them.

Per genus (one SLURM job each):
  E1. run_assign_features_workflow(genome_type='phage', genome_list=FULL val phages)
      -> a full-held-out-phage feature table in the model's phage vocabulary
  E2. assign_predict_workflow(genome_type='strain', genome_list=FULL val strains,
      phage_feature_table_path=E1 output) -> predict FULL val strains x val phages

Reused per-fold artifacts (from the SUBSET run):
  <SUBSET_FOLDS>/<genus>/strain/features/selected_features.csv
  <SUBSET_FOLDS>/<genus>/phage/features/selected_features.csv
  <SUBSET_FOLDS>/<genus>/modeling_results/<best_cutoff>/        (the model)
  <SUBSET_FOLDS>/<genus>/feature_selection/filtered_feature_tables/select_feature_table_<best_cutoff>.csv
  <SUBSET_SHARED>/{strain,phage}/clusters.tsv, /tmp/{strain,phage}/mmseqs_db   (assignment target)

FULL held-out validation lists come from the full corner folds:
  <FULL_FOLDS>/fold_<genus>/validation_{strains,phages}.csv

  python submit_corner_repredict_full_heldout.py                 # all genera
  python submit_corner_repredict_full_heldout.py --genus Vibrio  # one genus
"""
import os
import argparse
import subprocess

# ============================ CONFIG (edit me) ============================
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"

# SUBSET run: trained models + selected features + shared clustering (msi02_c08).
SUBSET_RUN = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_corner_logo/msi02_c08"
SUBSET_FOLDS = os.path.join(SUBSET_RUN, "folds")          # <genus>/ subdirs (note: dir names use '_' for spaces)
SUBSET_SHARED = os.path.join(SUBSET_RUN, "shared_clustering")

# FULL corner folds: validation lists = the COMPLETE held-out genus.
FULL_FOLDS = "/global/home/groups/pc_phiml/set_transformer/manuscript/genophi_corner_folds_FULL"

# Output (fresh dir; does NOT touch the subset predictions).
base_output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_corner_repredict_FULL"

# Held-out genera. The label (with spaces) must match how the subset fold dirs
# and full fold dirs were named: subset folds = label.replace(' ', '_');
# full folds = 'fold_' + safe_name(label). Map both here.
GENERA = ["E. coli", "Kleb - Hemaa", "Pseudomonas", "Vibrio"]

# MMseqs assignment params -- MUST match what the models were trained at (msi02_c08).
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
mem_per_job = "240G"
time_limit = "12:00:00"   # prediction only -- much faster than training
# =========================================================================


def safe(s):
    import re
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s.strip())


def subset_fold_dir(label):
    # subset run named fold dirs by label with spaces -> underscores.
    return os.path.join(SUBSET_FOLDS, label.replace(' ', '_'))


def full_fold_dir(label):
    # build_corner_folds.py named dirs 'fold_<safe_name>'.
    return os.path.join(FULL_FOLDS, f"fold_{safe(label)}")


def out_dir(label):
    return os.path.join(base_output_dir, "folds", safe(label))


logs_dir = os.path.join(base_output_dir, "slurm_logs")


def header(f, job, time_lim=time_limit):
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name={job}\n")
    f.write(f"#SBATCH --account={account}\n#SBATCH --partition={partition}\n#SBATCH --qos={qos}\n")
    f.write("#SBATCH --nodes=1\n#SBATCH --ntasks=1\n")
    f.write(f"#SBATCH --cpus-per-task={threads}\n#SBATCH --mem={mem_per_job}\n#SBATCH --time={time_lim}\n")
    f.write(f"#SBATCH --output={logs_dir}/{job}_%j.out\n#SBATCH --error={logs_dir}/{job}_%j.err\n")
    f.write("\nmodule load miniforge3\n")
    f.write('eval "$(conda shell.bash hook)"\n')
    f.write(f"conda activate {environment}\n\n")


def repredict_py(label):
    sub = subset_fold_dir(label)
    full = full_fold_dir(label)
    od = out_dir(label)
    # The python -c reproduces corner E1+E2 against existing artifacts.
    return (
        "python -c \""
        "import os, pandas as pd; "
        "from genophi.workflows.assign_features_workflow import run_assign_features_workflow; "
        "from genophi.workflows.assign_predict_workflow import assign_predict_workflow; "
        f"sub='{sub}'; full='{full}'; od='{od}'; "
        "os.makedirs(od, exist_ok=True); "
        # best cutoff = highest-MCC row from the trained model's metrics (matches _select_best_cutoff)
        "mdf=pd.read_csv(os.path.join(sub,'modeling_results','model_performance','model_performance_metrics.csv')); "
        "mdf=mdf.sort_values(['MCC','cut_off'],ascending=[False,False]); "
        "cut=mdf['cut_off'].values[0]; "
        "model_dir=os.path.join(sub,'modeling_results',str(cut)); "
        "select_ft=os.path.join(sub,'feature_selection','filtered_feature_tables','select_feature_table_'+str(cut)+'.csv'); "
        # E1: assign FULL held-out phages -> phage feature table in model vocabulary
        "vpdir=os.path.join(od,'validation_phage_features'); os.makedirs(vpdir,exist_ok=True); "
        "run_assign_features_workflow("
        f"input_dir='{input_phage_dir}', "
        f"mmseqs_db='{os.path.join(SUBSET_SHARED,'tmp','phage','mmseqs_db')}', "
        "tmp_dir=os.path.join(vpdir,'tmp'), output_dir=vpdir, "
        "feature_map=os.path.join(sub,'phage','features','selected_features.csv'), "
        f"clusters_tsv='{os.path.join(SUBSET_SHARED,'phage','clusters.tsv')}', "
        "genome_type='phage', genome_list=os.path.join(full,'validation_phages.csv'), "
        f"sensitivity={sensitivity}, coverage={coverage}, min_seq_id={min_seq_id}, "
        f"threads={threads}, duplicate_all={duplicate_all}); "
        "vpt=os.path.join(vpdir,'phage_combined_feature_table.csv'); "
        # E2: assign FULL held-out strains, predict paired with the phage table
        "vodir=os.path.join(od,'model_validation'); os.makedirs(vodir,exist_ok=True); "
        "assign_predict_workflow("
        f"input_dir='{input_strain_dir}', "
        "genome_list=os.path.join(full,'validation_strains.csv'), "
        f"mmseqs_db='{os.path.join(SUBSET_SHARED,'tmp','strain','mmseqs_db')}', "
        f"clusters_tsv='{os.path.join(SUBSET_SHARED,'strain','clusters.tsv')}', "
        "feature_map=os.path.join(sub,'strain','features','selected_features.csv'), "
        "tmp_dir=os.path.join(vodir,'tmp'), suffix='faa', model_dir=model_dir, "
        "feature_table=select_ft, phage_feature_table_path=vpt, output_dir=vodir, "
        f"threads={threads}, genome_type='strain', sensitivity={sensitivity}, "
        f"coverage={coverage}, min_seq_id={min_seq_id}, duplicate_all={duplicate_all}); "
        "print('DONE', od)"
        "\"\n"
    )


def submit(sh):
    return subprocess.check_output(["sbatch", "--parsable", sh]).decode().strip()


def write_job(label):
    j = f"REPRED_{safe(label)}"
    w = os.path.join(logs_dir, f"{j}.sh")
    with open(w, 'w') as f:
        header(f, j)
        f.write(f'echo "Re-predicting full held-out: {label}"\n\n')
        f.write(repredict_py(label))
    return w


def main():
    ap = argparse.ArgumentParser(description="Re-predict corner on FULL held-out genus (subset-trained models).")
    ap.add_argument('--genus', help="Run one genus (e.g. 'Vibrio'). Default: all.")
    args = ap.parse_args()

    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(os.path.join(base_output_dir, "folds"), exist_ok=True)

    targets = [args.genus] if args.genus else GENERA
    for label in targets:
        if label not in GENERA:
            raise SystemExit(f"--genus '{label}' not in {GENERA}")
        # Sanity: required artifacts exist before submitting.
        sub = subset_fold_dir(label); full = full_fold_dir(label)
        for p in (os.path.join(sub, 'modeling_results', 'model_performance', 'model_performance_metrics.csv'),
                  os.path.join(sub, 'strain', 'features', 'selected_features.csv'),
                  os.path.join(sub, 'phage', 'features', 'selected_features.csv'),
                  os.path.join(full, 'validation_strains.csv'),
                  os.path.join(full, 'validation_phages.csv')):
            if not os.path.exists(p):
                print(f"  WARNING [{label}]: missing {p}")
        jid = submit(write_job(label))
        print(f"  {label}: {jid} -> {out_dir(label)}")

    print(f"\nResults: {base_output_dir}/folds/<genus>/model_validation/predict_results/")
    print("Monitor: squeue -u $USER")


if __name__ == "__main__":
    main()
