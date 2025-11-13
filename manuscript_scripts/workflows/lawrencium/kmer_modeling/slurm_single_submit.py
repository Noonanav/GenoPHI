#!/usr/bin/env python3
"""
SLURM workflow submission for single k-mer-based phage-host interaction modeling run.
Creates a single SLURM job that runs the full k-mer workflow on the complete dataset.
"""

import os
import sys
import argparse
import subprocess
import time

def submit_job(script_path):
    """Submit a SLURM job and return job ID."""
    try:
        result = subprocess.run(['sbatch', '--parsable', script_path], 
                              capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error submitting {script_path}: {e}")
        print(f"Error output: {e.stderr}")
        return None

def create_single_job_script(args, run_dir):
    """Create SLURM job script for single k-mer workflow run"""
    
    # Get the absolute path to the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Build the command-line arguments for kmer_full_workflow.py
    workflow_args = [
        f"-is {args.input_strain_dir}",
        f"-pm {args.interaction_matrix}",
        f"-o {args.output_dir}",
        f"--k {args.k}",
        f"--suffix faa",
        f"--threads {args.threads}",
        f"--strain_list {args.strain_list}",
        f"--phage_list {args.phage_list}",
        f"--strain_column {args.strain_column}",
        f"--phage_column phage",
        f"--sample_column {args.strain_column}",
        f"--phenotype_column interaction",
        f"--num_features 100",
        f"--filter_type {args.filter_type}",
        f"--num_runs_fs {args.num_runs_fs}",
        f"--num_runs_modeling {args.num_runs_modeling}",
        f"--method {args.method}",
        f"--task_type classification",
        f"--max_features none",
        f"--max_ram {args.max_ram}",
        f"--weights_method {args.weights_method}",
        f"--cluster_method {args.cluster_method}",
        f"--n_clusters {args.n_clusters}",
        f"--min_cluster_size {args.min_cluster_size}",
        f"--cluster_selection_epsilon {args.cluster_selection_epsilon}",
        f"--min_cluster_presence {args.min_cluster_presence}",
        f"--feature_cluster_method {args.feature_cluster_method}",
        f"--feature_n_clusters {args.feature_n_clusters}",
        f"--feature_min_cluster_presence {args.feature_min_cluster_presence}",
    ]
    
    # Add optional phage directory
    if args.input_phage_dir:
        workflow_args.append(f"-ip {args.input_phage_dir}")
    
    # Add boolean flags
    if args.k_range:
        workflow_args.append("--k_range")
    if args.one_gene:
        workflow_args.append("--one_gene")
    if args.use_dynamic_weights:
        workflow_args.append("--use_dynamic_weights")
    if args.use_clustering:
        workflow_args.append("--use_clustering")
    if args.min_samples:
        workflow_args.append(f"--min_samples {args.min_samples}")
    if args.check_feature_presence:
        workflow_args.append("--check_feature_presence")
    if args.filter_by_cluster_presence:
        workflow_args.append("--filter_by_cluster_presence")
    if args.use_feature_clustering:
        workflow_args.append("--use_feature_clustering")
    if args.remove_suffix:
        workflow_args.append("--remove_suffix")
    if args.use_shap:
        workflow_args.append("--use_shap")
    if args.run_predictive_analysis:
        workflow_args.append("--run_predictive_analysis")
    
    workflow_cmd = " ".join(workflow_args)
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=kmer_single_run
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_per_job}G
#SBATCH --time={args.time_limit}
#SBATCH --output=logs/kmer_single_%j.out
#SBATCH --error=logs/kmer_single_%j.err

echo "=== K-mer Single Workflow Run ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

# Run the k-mer full workflow
python3 -m phage_modeling.workflows.kmer_full_workflow {workflow_cmd}

echo "K-mer workflow completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "kmer_single_job.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Submit k-mer single workflow as SLURM job")
    
    # Required input arguments
    parser.add_argument('--input_strain_dir', type=str, required=True, help="Directory containing strain FASTA files.")
    parser.add_argument('--input_phage_dir', type=str, help="Directory containing phage FASTA files.")
    parser.add_argument('--interaction_matrix', type=str, required=True, help="Path to the interaction matrix.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--strain_list', type=str, help="Path to strain list file (default: uses interaction_matrix).")
    parser.add_argument('--phage_list', type=str, help="Path to phage list file (default: uses interaction_matrix).")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads to use.")
    parser.add_argument('--strain_column', type=str, default='strain', help="Column in the interaction matrix containing strain names.")
    parser.add_argument('--num_runs_fs', type=int, default=25, help="Number of runs for feature selection.")
    parser.add_argument('--num_runs_modeling', type=int, default=50, help="Number of runs for modeling.")
    parser.add_argument('--method', type=str, default='rfe', choices=['rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'], help="Feature selection method.")
    
    # K-mer specific parameters
    parser.add_argument('--k', type=int, default=5, help="K-mer length.")
    parser.add_argument('--k_range', action='store_true', help="Generate k-mers from 3 to k.")
    parser.add_argument('--one_gene', action='store_true', help="Include features with only one gene.")
    
    # Feature selection parameters
    parser.add_argument('--filter_type', type=str, default='strain', help="Filter type for train-test split ('none', 'strain', 'phage').")
    parser.add_argument('--use_dynamic_weights', action='store_true', help="Use dynamic weights for feature selection.")
    parser.add_argument('--weights_method', type=str, default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help="Method for calculating dynamic weights.")
    parser.add_argument('--use_clustering', action='store_true', help="Use clustering results for feature selection.")
    parser.add_argument('--cluster_method', type=str, default='hdbscan', choices=['hdbscan', 'hierarchical'], help="Clustering method to use.")
    parser.add_argument('--n_clusters', type=int, default=20, help="Number of clusters for clustering feature selection.")
    parser.add_argument('--min_cluster_size', type=int, default=2, help="Minimum cluster size for clustering feature selection.")
    parser.add_argument('--min_samples', type=int, help="Minimum number of samples for clustering feature selection.")
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help="Epsilon value for clustering feature selection.")
    parser.add_argument('--check_feature_presence', action='store_true', help="Check for feature presence for train-test splits.")
    parser.add_argument('--filter_by_cluster_presence', action='store_true', help="Filter features by cluster presence across train/test splits.")
    parser.add_argument('--min_cluster_presence', type=int, default=2, help="Minimum number of clusters that must contain the feature.")
    parser.add_argument('--max_ram', type=float, default=40, help="Maximum RAM usage in GB.")
    parser.add_argument('--use_shap', action='store_true', help="Use SHAP values for analysis.")
    parser.add_argument('--run_predictive_analysis', action='store_true', help="Run predictive analysis workflow.")
    
    # Feature clustering parameters
    parser.add_argument('--use_feature_clustering', action='store_true', help="Enable pre-processing cluster-based feature filtering.")
    parser.add_argument('--feature_cluster_method', type=str, default='hierarchical', choices=['hierarchical'], help="Pre-processing clustering method ('hierarchical' only for now).")
    parser.add_argument('--feature_n_clusters', type=int, default=20, help="Number of clusters for pre-processing feature clustering.")
    parser.add_argument('--feature_min_cluster_presence', type=int, default=2, help="Min clusters a feature must appear in during pre-processing.")
    
    # Optional parameters
    parser.add_argument('--remove_suffix', action='store_true', help="Remove suffix from genome names when merging.")
    
    # SLURM-specific arguments
    parser.add_argument('--account', default='ac_mak', help='SLURM account (default: ac_mak).')
    parser.add_argument('--partition', default='lr7', help='SLURM partition (default: lr7).')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS (default: lr_normal).')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment name (default: phage_modeling).')
    parser.add_argument('--mem_per_job', type=int, default=60, help='Memory per job in GB (default: 60).')
    parser.add_argument('--time_limit', default='12:00:00', help='Time limit per job (default: 12:00:00).')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts but do not submit jobs.')
    
    args = parser.parse_args()
    
    # Set default strain_list and phage_list to interaction_matrix if not provided
    if not args.strain_list:
        args.strain_list = args.interaction_matrix
    if not args.phage_list:
        args.phage_list = args.interaction_matrix
    
    # Validate critical inputs exist
    if not os.path.exists(args.input_strain_dir):
        print(f"Error: Input strain directory not found: {args.input_strain_dir}")
        return 1
    if args.input_phage_dir and not os.path.exists(args.input_phage_dir):
        print(f"Error: Input phage directory not found: {args.input_phage_dir}")
        return 1
    if not os.path.exists(args.interaction_matrix):
        print(f"Error: Interaction matrix not found: {args.interaction_matrix}")
        return 1
    
    # Check if output directory parent exists
    output_parent = os.path.dirname(args.output_dir)
    if output_parent and not os.path.exists(output_parent):
        print(f"Error: Output parent directory not found: {output_parent}")
        return 1
    
    # Validate interaction matrix has expected column
    try:
        import pandas as pd
        df = pd.read_csv(args.interaction_matrix)
        if args.strain_column not in df.columns:
            print(f"Error: Column '{args.strain_column}' not found in interaction matrix")
            print(f"Available columns: {list(df.columns)}")
            return 1
        print(f"Interaction matrix validation: ✓ Found {len(df)} rows with column '{args.strain_column}'")
    except Exception as e:
        print(f"Error reading interaction matrix: {e}")
        return 1
    
    # Check strain files exist
    strain_files = [f for f in os.listdir(args.input_strain_dir) if f.endswith('.faa')]
    if len(strain_files) == 0:
        print(f"Error: No .faa files found in strain directory: {args.input_strain_dir}")
        return 1
    print(f"Strain directory validation: ✓ Found {len(strain_files)} .faa files")
    
    # Check phage files exist if directory provided
    if args.input_phage_dir:
        phage_files = [f for f in os.listdir(args.input_phage_dir) if f.endswith('.faa')]
        if len(phage_files) == 0:
            print(f"Error: No .faa files found in phage directory: {args.input_phage_dir}")
            return 1
        print(f"Phage directory validation: ✓ Found {len(phage_files)} .faa files")
    
    # Create timestamped run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"kmer_single_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== K-mer Single Workflow SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Strain list: {args.strain_list}")
    print(f"Phage list: {args.phage_list}")
    print(f"Threads per job: {args.threads}")
    print(f"K-mer length: {args.k}")
    print(f"K-mer range: {args.k_range}")
    print(f"Filter type: {args.filter_type}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    print(f"Use SHAP: {args.use_shap}")
    print(f"Use feature clustering: {args.use_feature_clustering}")
    if args.use_feature_clustering:
        print(f"Feature cluster method: {args.feature_cluster_method}")
        print(f"Feature n_clusters: {args.feature_n_clusters}")
        print(f"Feature min presence: {args.feature_min_cluster_presence}")
    
    # Resource validation
    if args.mem_per_job < 32:
        print(f"WARNING: {args.mem_per_job}GB memory may be insufficient for k-mer generation on complex datasets")
    
    print()
    
    # Create script
    print("Creating SLURM job script...")
    job_script = create_single_job_script(args, run_dir)
    
    if args.dry_run:
        print("Dry run - script created but not submitted")
        print("Script:")
        print(f"  Single workflow job: {job_script}")
        return
    
    # Change to run directory
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    # Submit job
    print("Submitting k-mer single workflow job...")
    job_id = submit_job("kmer_single_job.sh")
    print(f"Single workflow job: {job_id}")
    
    # Change back to original directory
    os.chdir(original_dir)
    
    print(f"\n=== K-mer Single Workflow Job Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Expected total runtime: {args.time_limit}")
    print("\nMonitor with:")
    print("  squeue -u $USER")
    print("  squeue -u $USER --name=kmer_single_run")
    print("  tail -f logs/kmer_single_*.out")
    print(f"\nResults will be in:")
    print(f"  Output directory: {args.output_dir}/")
    print(f"  Feature tables: {args.output_dir}/feature_tables/")
    print(f"  Modeling results: {args.output_dir}/modeling/")
    print(f"  Reports: {args.output_dir}/kmer_workflow_report.txt")

if __name__ == "__main__":
    sys.exit(main())
