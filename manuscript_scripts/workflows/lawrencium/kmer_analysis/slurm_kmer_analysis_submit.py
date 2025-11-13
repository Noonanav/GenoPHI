#!/usr/bin/env python3
"""
Prepare SLURM submission scripts for k-mer analysis workflow.

Usage:
    # For strain features:
    python prepare_kmer_analysis_jobs.py \
        --modeling-dir k4/modeling \
        --feature-selection-dir k4/modeling/feature_selection \
        --feature2cluster-path k4/feature_tables/selected_features.csv \
        --protein-mapping-csv k4/strain_proteins.csv \
        --aa-sequence-file k4/strain_combined.faa \
        --output-dir k4/kmer_analysis/strain \
        --feature-type strain
    
    # For phage features:
    python prepare_kmer_analysis_jobs.py \
        --modeling-dir k4/modeling \
        --feature-selection-dir k4/modeling/feature_selection \
        --feature2cluster-path k4/feature_tables/phage_selected_features.csv \
        --protein-mapping-csv k4/phage_proteins.csv \
        --aa-sequence-file k4/phage_combined.faa \
        --output-dir k4/kmer_analysis/phage \
        --feature-type phage
"""

import argparse
from pathlib import Path
from datetime import datetime
import os

# ────────────────────────────────────────────────────────────────── CONFIG ──
SLURM_CONFIG = {
    "account": "pc_phiml",
    "partition": "lr7",
    "qos": "lr_normal",
    "time": "12:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 8,
    "mem": "120G",
    "conda_env": "phage_modeling",
}

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=kmer_{feature_type}_{timestamp}
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}
#SBATCH --output={log_dir}/slurm-%j.out
#SBATCH --error={log_dir}/slurm-%j.err

# Load conda (adjust module name if needed for your system)
# Common options: anaconda3, python, anaconda, or remove if conda already in PATH
module load python/3.11 2>/dev/null || module load anaconda3 2>/dev/null || true

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate {conda_env} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {conda_env}
}}

# Verify conda environment is active
if [[ "$CONDA_DEFAULT_ENV" != "{conda_env}" ]]; then
    echo "ERROR: Failed to activate conda environment: {conda_env}"
    echo "Current environment: $CONDA_DEFAULT_ENV"
    echo "Available environments:"
    conda env list
    exit 1
fi

echo "Conda environment activated: $CONDA_DEFAULT_ENV"
echo ""

set -euo pipefail

echo "=================================================="
echo "K-mer Analysis Job - {feature_type}"
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $HOSTNAME"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "=================================================="
echo ""

# Run k-mer analysis
python {script_path} \\
    --modeling_dir {modeling_dir} \\
    --feature_selection_dir {feature_selection_dir} \\
    --feature2cluster_path {feature2cluster_path} \\
    --protein_mapping_csv {protein_mapping_csv} \\
    --aa_sequence_file {aa_sequence_file} \\
    --output_dir {output_dir} \\
    --feature_type {feature_type} \\
    --threads {cpus_per_task}{maybe_cutoff}{maybe_clustering}{maybe_plot}

exit_code=$?

echo ""
echo "=================================================="
echo "Job completed with exit code: $exit_code"
echo "Finished: $(date)"
echo "=================================================="

exit $exit_code
"""

SUBMIT_TEMPLATE = """\
#!/bin/bash
# Submission script for k-mer analysis
# Generated: {timestamp}

echo "=================================================="
echo "K-mer Analysis Submission"
echo "=================================================="
echo "Feature type: {feature_type}"
echo "Output directory: {output_dir}"
echo "Log directory: {log_dir}"
echo "Script: {slurm_script}"
echo "=================================================="
echo ""

# Submit the job
JOB_ID=$(sbatch {slurm_script} | awk '{{print $4}}')

echo "Submitted job: $JOB_ID"
echo ""
echo "Monitor job status:"
echo "  squeue -j $JOB_ID"
echo "  squeue -u $USER"
echo ""
echo "View logs:"
echo "  tail -f {log_dir}/slurm-$JOB_ID.out"
echo "  tail -f {log_dir}/slurm-$JOB_ID.err"
echo ""
echo "Cancel job:"
echo "  scancel $JOB_ID"
echo ""
"""


def resolve_path(path_str):
    """Resolve relative paths to absolute paths."""
    path = Path(path_str).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return str(path.resolve())


def main():
    parser = argparse.ArgumentParser(
        description="Generate SLURM submission scripts for k-mer analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Required arguments
    parser.add_argument(
        "--modeling-dir",
        required=True,
        help="Path to modeling directory (contains modeling_results/)"
    )
    parser.add_argument(
        "--feature-selection-dir",
        required=True,
        help="Path to feature_selection directory (contains filtered_feature_tables/)"
    )
    parser.add_argument(
        "--feature2cluster-path",
        required=True,
        help="Path to selected_features.csv (or strain_selected_features.csv)"
    )
    parser.add_argument(
        "--protein-mapping-csv",
        required=True,
        help="Path to protein mapping CSV (strain_proteins.csv or phage_proteins.csv)"
    )
    parser.add_argument(
        "--aa-sequence-file",
        required=True,
        help="Path to combined FASTA file (strain_combined.faa or phage_combined.faa)"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--feature-type",
        required=True,
        choices=["strain", "phage"],
        help="Type of features to analyze"
    )
    
    # Optional analysis parameters
    parser.add_argument(
        "--cutoff",
        help="Specific cutoff to use (default: auto-detect best)"
    )
    parser.add_argument(
        "--use-clustering",
        action="store_true",
        help="Use MMseqs2 to cluster proteins before alignment"
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=0.5,
        help="MMseqs2 minimum sequence identity (default: 0.5)"
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=0.8,
        help="MMseqs2 minimum coverage (default: 0.8)"
    )
    parser.add_argument(
        "--sensitivity",
        type=float,
        default=7.5,
        help="MMseqs2 sensitivity (default: 7.5)"
    )
    parser.add_argument(
        "--create-plots",
        action="store_true",
        help="Create visualization plots"
    )
    
    # SLURM configuration overrides
    parser.add_argument("--account", help="SLURM account (default: pc_phiml)")
    parser.add_argument("--partition", help="SLURM partition (default: lr7)")
    parser.add_argument("--qos", help="SLURM QOS (default: lr_normal)")
    parser.add_argument("--time", help="Time limit (default: 12:00:00)")
    parser.add_argument("--cpus", type=int, help="CPUs per task (default: 8)")
    parser.add_argument("--mem", help="Memory per job (default: 120G)")
    parser.add_argument("--conda-env", help="Conda environment name (default: phage_modeling)")
    
    # Script location
    parser.add_argument(
        "--script-path",
        default="kmer_analysis.py",
        help="Path to kmer_analysis.py script (default: ./kmer_analysis.py)"
    )
    
    args = parser.parse_args()
    
    # Merge with defaults
    cfg = SLURM_CONFIG.copy()
    if args.account:
        cfg["account"] = args.account
    if args.partition:
        cfg["partition"] = args.partition
    if args.qos:
        cfg["qos"] = args.qos
    if args.time:
        cfg["time"] = args.time
    if args.cpus:
        cfg["cpus_per_task"] = args.cpus
    if args.mem:
        cfg["mem"] = args.mem
    if args.conda_env:
        cfg["conda_env"] = args.conda_env
    
    # Generate timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Resolve all paths to absolute
    modeling_dir = resolve_path(args.modeling_dir)
    feature_selection_dir = resolve_path(args.feature_selection_dir)
    feature2cluster_path = resolve_path(args.feature2cluster_path)
    protein_mapping_csv = resolve_path(args.protein_mapping_csv)
    aa_sequence_file = resolve_path(args.aa_sequence_file)
    output_dir = resolve_path(args.output_dir)
    script_path = resolve_path(args.script_path)
    
    # Validate script exists
    if not Path(script_path).exists():
        print(f"ERROR: Script not found: {script_path}")
        print(f"Please provide correct path with --script-path")
        return 1
    
    # Create output directories
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    log_dir = Path.cwd() / "slurm_logs" / f"{args.feature_type}_{timestamp}"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    scripts_dir = Path.cwd() / "slurm_scripts" / f"{args.feature_type}_{timestamp}"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    
    # Build optional arguments
    maybe_cutoff = f" \\\n    --cutoff {args.cutoff}" if args.cutoff else ""
    
    maybe_clustering = ""
    if args.use_clustering:
        maybe_clustering = f" \\\n    --use_clustering"
        maybe_clustering += f" \\\n    --min_seq_id {args.min_seq_id}"
        maybe_clustering += f" \\\n    --coverage {args.coverage}"
        maybe_clustering += f" \\\n    --sensitivity {args.sensitivity}"
    
    maybe_plot = " \\\n    --create_plots" if args.create_plots else ""
    
    # Generate SLURM script
    slurm_script = SLURM_TEMPLATE.format(
        feature_type=args.feature_type,
        timestamp=timestamp,
        account=cfg["account"],
        partition=cfg["partition"],
        qos=cfg["qos"],
        time=cfg["time"],
        nodes=cfg["nodes"],
        ntasks=cfg["ntasks"],
        cpus_per_task=cfg["cpus_per_task"],
        mem=cfg["mem"],
        conda_env=cfg["conda_env"],
        log_dir=str(log_dir),
        script_path=script_path,
        modeling_dir=modeling_dir,
        feature_selection_dir=feature_selection_dir,
        feature2cluster_path=feature2cluster_path,
        protein_mapping_csv=protein_mapping_csv,
        aa_sequence_file=aa_sequence_file,
        output_dir=output_dir,
        maybe_cutoff=maybe_cutoff,
        maybe_clustering=maybe_clustering,
        maybe_plot=maybe_plot
    )
    
    slurm_script_path = scripts_dir / f"kmer_analysis_{args.feature_type}.slurm"
    slurm_script_path.write_text(slurm_script)
    slurm_script_path.chmod(0o755)
    
    # Generate submission script
    submit_script = SUBMIT_TEMPLATE.format(
        timestamp=timestamp,
        feature_type=args.feature_type,
        output_dir=output_dir,
        log_dir=str(log_dir),
        slurm_script=str(slurm_script_path)
    )
    
    submit_script_path = scripts_dir / f"submit_{args.feature_type}.sh"
    submit_script_path.write_text(submit_script)
    submit_script_path.chmod(0o755)
    
    # Print summary
    print("=" * 70)
    print("K-mer Analysis SLURM Scripts Generated")
    print("=" * 70)
    print(f"Feature type:      {args.feature_type}")
    print(f"Output directory:  {output_dir}")
    print(f"Log directory:     {log_dir}")
    print(f"")
    print(f"SLURM script:      {slurm_script_path}")
    print(f"Submit script:     {submit_script_path}")
    print(f"")
    print("SLURM Configuration:")
    print(f"  Account:         {cfg['account']}")
    print(f"  Partition:       {cfg['partition']}")
    print(f"  QOS:             {cfg['qos']}")
    print(f"  Time limit:      {cfg['time']}")
    print(f"  CPUs:            {cfg['cpus_per_task']}")
    print(f"  Memory:          {cfg['mem']}")
    print(f"  Conda env:       {cfg['conda_env']}")
    print(f"")
    print("To submit:")
    print(f"  {submit_script_path}")
    print(f"")
    print("Or directly:")
    print(f"  sbatch {slurm_script_path}")
    print("=" * 70)
    
    return 0


if __name__ == "__main__":
    exit(main())