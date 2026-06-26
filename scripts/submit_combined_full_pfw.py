#!/usr/bin/env python3
"""
Train ONE GenoPHI protein-family model on all 4 genera combined.

No cross-validation, no held-out prediction, no folds -- a single
`genophi protein-family-workflow` run over the full combined interaction matrix:
MMseqs2 clustering -> feature tables -> feature selection -> modeling.

Settings match the msi02_c08 corner experiment (min_seq_id=0.2, coverage=0.8,
dynamic inverse-frequency weights, hierarchical sample clustering, filter_type
strain, max_ram=100), so this combined model is comparable to the corner folds.

  python submit_combined_full_pfw.py             # submit
  python submit_combined_full_pfw.py --dry_run   # write the .sh, don't sbatch
"""
import os
import argparse
import subprocess

# ============================ CONFIG (edit me) ============================
input_strain_dir = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
input_phage_dir = "/global/home/groups/pc_phiml/data/combined/phage_AAs"
phenotype_matrix = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"

output_dir = "/global/scratch/users/anoonan/set_transformer/manuscript/genophi_combined_full_model"

# Modeling settings (match msi02_c08).
min_seq_id = 0.2
coverage = 0.8
sensitivity = 7.5
num_runs_fs = 25
num_runs_modeling = 50
filter_type = "strain"
weights_method = "inverse_frequency"
cluster_method = "hierarchical"
max_ram = 100
threads = 32

# SLURM
account = "pc_crispriart"
partition = "lr7"
qos = "lr_normal"
environment = "genophi"
mem_per_job = "240G"
time_limit = "48:00:00"   # clustering + FS + modeling on the full combined set is heavy
# =========================================================================

logs_dir = os.path.join(output_dir, "slurm_logs")
tmp_dir = os.path.join(output_dir, "tmp")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry_run', action='store_true', help="Write the .sh but do not submit.")
    args = ap.parse_args()

    os.makedirs(logs_dir, exist_ok=True)

    job = "COMBINED_pfw"
    sh = os.path.join(logs_dir, f"{job}.sh")
    with open(sh, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name={job}\n")
        f.write(f"#SBATCH --account={account}\n#SBATCH --partition={partition}\n#SBATCH --qos={qos}\n")
        f.write("#SBATCH --nodes=1\n#SBATCH --ntasks=1\n")
        f.write(f"#SBATCH --cpus-per-task={threads}\n#SBATCH --mem={mem_per_job}\n#SBATCH --time={time_limit}\n")
        f.write(f"#SBATCH --output={logs_dir}/{job}_%j.out\n#SBATCH --error={logs_dir}/{job}_%j.err\n")
        f.write("\nmodule load miniforge3\n")
        f.write('eval "$(conda shell.bash hook)"\n')
        f.write(f"conda activate {environment}\n\n")
        f.write(
            "genophi protein-family-workflow \\\n"
            f"  --input_path_strain {input_strain_dir} \\\n"
            f"  --input_path_phage {input_phage_dir} \\\n"
            f"  --phenotype_matrix {phenotype_matrix} \\\n"
            f"  --output_dir {output_dir} \\\n"
            f"  --tmp_dir {tmp_dir} \\\n"
            f"  --min_seq_id {min_seq_id} \\\n"
            f"  --coverage {coverage} \\\n"
            f"  --sensitivity {sensitivity} \\\n"
            f"  --filter_type {filter_type} \\\n"
            f"  --num_runs_fs {num_runs_fs} \\\n"
            f"  --num_runs_modeling {num_runs_modeling} \\\n"
            f"  --task_type classification \\\n"
            f"  --use_dynamic_weights \\\n"
            f"  --weights_method {weights_method} \\\n"
            f"  --use_clustering \\\n"
            f"  --cluster_method {cluster_method} \\\n"
            f"  --max_ram {max_ram} \\\n"
            f"  --threads {threads}\n"
        )
    print(f"wrote {sh}")

    if args.dry_run:
        print("dry_run: not submitting. Inspect the .sh above.")
        return
    jid = subprocess.check_output(["sbatch", "--parsable", sh]).decode().strip()
    print(f"submitted: {jid}")
    print(f"output -> {output_dir}/  (models under modeling_results/)")
    print(f"monitor: squeue -u $USER ; tail -f {logs_dir}/{job}_{jid}.out")


if __name__ == "__main__":
    main()
