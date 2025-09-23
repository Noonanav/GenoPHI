#!/usr/bin/env python
"""
Prepare Slurm submission scripts for running the Set-Transformer Optuna workflow
on Lawrencium.

Edit the CONFIG block below (or pass overrides on the CLI) and run:

    python prepare_optuna_jobs.py --array

to make a single job-array submission file, *or*

    python prepare_optuna_jobs.py --split

to emit one .slurm file per trial.

Then run ./submit_optuna_jobs.sh
"""

from __future__ import annotations
import argparse, json, textwrap
from pathlib import Path
from datetime import datetime

# ────────────────────────────────────────────────────────────────── CONFIG ──
CONFIG = dict(
    # ---------- data / workflow ----------
    interactions="/global/scratch/users/anoonan/BRaVE/ecoli/Gaborieau_interaction_matrix_long_mod.csv",
    strain_emb="/global/scratch/users/anoonan/set_transformer/pLMs/ESM2/ecoli/Gaborieau/strain",
    phage_emb="/global/scratch/users/anoonan/set_transformer/pLMs/ESM2/ecoli/Gaborieau/phage",
    output_dir="/global/scratch/users/anoonan/set_transformer/pLMs/modeling/ecoli/ESM2",
    trial_total=500,
    folds=5,
    final_seeds=0,
    search_config='/global/scratch/users/anoonan/set_transformer/scripts/ecoli/update/ecoli_ESM2_params.yaml',
    # ---------- Slurm resources ----------
    account="ac_mak",               # e.g. "es1" or "lr5"
    partition="es1",            # or lr5, lr_bigmem, gpu
    qos="es_normal",
    time="48:00:00",            # HH:MM:SS  (per trial)
    nodes=1,
    ntasks=1,
    cpus_per_task=2,
    mem="80G",
    gres="gpu:H100:1",                    # e.g. "gpu:1"  if needed
    conda_env="phage-set-transformer",  # name OR full conda activate path
)
# ────────────────────────────────────────────────────────────────────────────

TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=pst_trial_{job_id}
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}{maybe_gres}
#SBATCH --output={out_dir}/slurm-%A_%a.out
#SBATCH --error={out_dir}/slurm-%A_%a.err
{array_directive}

TRIAL_ID=${{SLURM_ARRAY_TASK_ID:-{job_id}}}

module load anaconda3
conda activate {conda_env} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {conda_env}
}}

set -euo pipefail
echo "Running trial $TRIAL_ID on $HOSTNAME at $(date)"

python -m phage_set_transformer.cli optimize \\
    --interactions {interactions} \\
    --strain-embeddings {strain_emb} \\
    --phage-embeddings {phage_emb} \\
    --trials 1 \\
    --folds {folds} \\
    --final-seeds {final_seeds} \\
    --output {out_dir} \\
    --study-name pst_lr_{timestamp} \\
    --seed $(( {seed_base} + $TRIAL_ID )){maybe_search_config}

echo "Finished at $(date)"
"""

def make_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def render_slurm(idx: int, cfg: dict, array: bool) -> str:
    # Conditional search config line
    maybe_search_config = ""
    if cfg["search_config"]:
        maybe_search_config = f" \\\n    --search-config {cfg['search_config']}"
    
    return TEMPLATE.format(
        job_id=idx,
        account=cfg["account"],
        partition=cfg["partition"],
        qos=cfg["qos"],
        time=cfg["time"],
        nodes=cfg["nodes"],
        ntasks=cfg["ntasks"],
        cpus_per_task=cfg["cpus_per_task"],
        mem=cfg["mem"],
        maybe_gres=f"\n#SBATCH --gres={cfg['gres']}" if cfg["gres"] else "",
        out_dir=cfg["output_dir"],
        array_directive="" if not array else f"#SBATCH --array=0-{cfg['trial_total']-1}%50",
        conda_env=cfg["conda_env"],
        interactions=cfg["interactions"],
        strain_emb=cfg["strain_emb"],
        phage_emb=cfg["phage_emb"],
        folds=cfg["folds"],
        final_seeds=cfg["final_seeds"],
        timestamp=cfg["timestamp"],
        seed_base=cfg["seed_base"],
        maybe_search_config=maybe_search_config,
    )

def create_optuna_study(output_dir: str, study_name: str, random_seed: int = 42):
    """Pre-create the Optuna study to avoid concurrency issues."""
    try:
        import optuna
        from optuna.storages import RDBStorage
        from optuna.samplers import TPESampler
        from optuna.pruners import MedianPruner
        
        storage = RDBStorage(url=f"sqlite:///{output_dir}/study.db")
        study = optuna.create_study(
            study_name=study_name,
            direction="maximize",
            sampler=TPESampler(seed=random_seed),
            pruner=MedianPruner(n_startup_trials=3),
            storage=storage,
            load_if_exists=True  # Won't error if already exists
        )
        print(f"Created Optuna study: {study_name}")
        print(f"Database: {output_dir}/study.db")
        return True
    except ImportError:
        print("Warning: Optuna not available. Jobs will create study individually (may cause conflicts).")
        return False
    except Exception as e:
        print(f"Warning: Could not create study: {e}")
        print("Jobs will attempt to create study individually.")
        return False

RETRAIN_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=pst_retrain_{timestamp}
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G{maybe_gres}
#SBATCH --output={out_dir}/retrain-%j.out
#SBATCH --error={out_dir}/retrain-%j.err
#SBATCH --dependency=afterany:{array_job_id}

module load anaconda3
conda activate {conda_env} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {conda_env}
}}

set -euo pipefail
echo "Starting final retraining at $(date)"

# Run final retraining with best parameters
python -m phage_set_transformer.cli optimize \\
    --interactions {interactions} \\
    --strain-embeddings {strain_emb} \\
    --phage-embeddings {phage_emb} \\
    --trials 0 \\
    --final-seeds 10 \\
    --output {out_dir} \\
    --study-name pst_lr_{timestamp}

echo "Final retraining completed at $(date)"
"""

# Update the main() function:
def main():
    p = argparse.ArgumentParser(description="Generate Slurm scripts for Optuna CV search")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--array", action="store_true", help="Write array job + retraining job")
    g.add_argument("--split", action="store_true", help="Write individual trial jobs")
    
    # Search config argument
    p.add_argument("--search-config", default=None, 
                   help="Path to YAML search space configuration file")
    
    args = p.parse_args()

    cfg = CONFIG.copy()
    cfg["timestamp"] = datetime.now().strftime("%Y%m%d_%H%M%S")
    cfg["seed_base"] = 42
    cfg["output_dir"] = str(make_dir(Path(cfg["output_dir"]) / cfg["timestamp"]))

    # Override search config from command line
    if args.search_config:
        cfg["search_config"] = args.search_config

    # Force final_seeds=0 for optimization jobs
    cfg["final_seeds"] = 0

    # Pre-create the Optuna study
    study_name = f"pst_lr_{cfg['timestamp']}"
    create_optuna_study(cfg["output_dir"], study_name, cfg["seed_base"])

    scripts_dir = make_dir(Path("slurm_scripts") / cfg["timestamp"])

    if args.array:
        # Generate optimization array job
        script = render_slurm(0, cfg, array=True)
        opt_script_path = scripts_dir / "optuna_array.slurm"
        opt_script_path.write_text(script)
        
        # Generate retraining job (will be submitted with dependency)
        retrain_script = RETRAIN_TEMPLATE.format(
            timestamp=cfg["timestamp"],
            account=cfg["account"],
            partition=cfg["partition"],
            qos=cfg["qos"],
            maybe_gres=f"\n#SBATCH --gres={cfg['gres']}" if cfg["gres"] else "",
            out_dir=cfg["output_dir"],
            array_job_id="PLACEHOLDER",  # Will be replaced in submission script
            conda_env=cfg["conda_env"],
            interactions=cfg["interactions"],
            strain_emb=cfg["strain_emb"],
            phage_emb=cfg["phage_emb"],
        )
        retrain_script_path = scripts_dir / "retrain.slurm"
        retrain_script_path.write_text(retrain_script)
        
        # Generate submission script
        submit_script = f"""#!/bin/bash
# Submit optimization array job
ARRAY_JOB=$(sbatch --parsable {opt_script_path})
echo "Submitted optimization array job: $ARRAY_JOB"

# Submit retraining job with dependency
sed "s/PLACEHOLDER/$ARRAY_JOB/" {retrain_script_path} > {scripts_dir}/retrain_final.slurm
RETRAIN_JOB=$(sbatch --parsable {scripts_dir}/retrain_final.slurm)
echo "Submitted retraining job: $RETRAIN_JOB (depends on $ARRAY_JOB)"

echo "Workflow submitted!"
echo "Monitor: watch 'squeue -u $USER'"
"""
        submit_path = scripts_dir / "submit_all.sh"
        submit_path.write_text(submit_script)
        submit_path.chmod(0o755)
        
        print(f"Generated optimization and retraining jobs in {scripts_dir}")
        print(f"To submit: {submit_path}")

    else:  # split mode unchanged
        for i in range(cfg["trial_total"]):
            s = render_slurm(i, cfg, array=False)
            (scripts_dir / f"trial_{i}.slurm").write_text(s)
        print(f"Wrote {cfg['trial_total']} individual .slurm files to {scripts_dir}")

if __name__ == "__main__":
    main()
