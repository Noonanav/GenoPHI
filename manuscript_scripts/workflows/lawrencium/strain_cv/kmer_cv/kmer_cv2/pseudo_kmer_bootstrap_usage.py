#!/usr/bin/env python3
"""
K-mer bootstrap cross-validation workflow submission script.
Edit the paths and parameters below, then run: python3 kmer_bootstrap_usage_example.py
"""

import subprocess
import sys
import os

def main():
    # =============================================
    # YOUR PATHS - EDIT THESE
    # =============================================
    input_strain_dir = "/global/scratch/users/anoonan/BRaVE/pseudomonas/strain_AAs"
    input_phage_dir = "/global/scratch/users/anoonan/BRaVE/pseudomonas/phage_AAs"
    interaction_matrix = "/global/scratch/users/anoonan/BRaVE/pseudomonas/pseudo_anarita_matrix_long.csv"
    output_dir = "/global/scratch/users/anoonan/BRaVE/pseudomonas/bootstrapping/pseudomonas_modeling_hierarchical"
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"                    # SLURM partition 
    qos = "lr_normal"                    # SLURM QOS
    environment = "phage_modeling"       # Conda environment name
    mem_per_job = "120"                   # Memory per iteration job in GB
    time_limit = "24:00:00"              # Time limit per iteration job
    
    # =============================================
    # K-MER BOOTSTRAP VALIDATION PARAMETERS
    # =============================================
    
    # Core parameters
    n_iterations = "20"                  # Number of cross-validation iterations (20 for 20-fold CV)
    threads = "16"                       # Number of threads per job
    strain_column = "strain"             # Column in interaction matrix with strain names
    max_ram = "60"                       # Max RAM usage in GB
    
    # Feature selection and modeling runs
    num_runs_fs = "25"                   # Number of feature selection runs per iteration
    num_runs_modeling = "50"             # Number of modeling runs per iteration
    
    # Dynamic weights and clustering for feature selection
    use_dynamic_weights = True           # Use dynamic class weights
    weights_method = "inverse_frequency" # Weight calculation method ('log10', 'inverse_frequency', 'balanced')
    use_clustering = True                # Use clustering for feature selection
    cluster_method = "hierarchical"      # Clustering method ('hdbscan' or 'hierarchical')
    n_clusters = "20"                    # Number of clusters for hierarchical clustering
    min_cluster_size = "3"               # Minimum cluster size for HDBSCAN
    min_samples = None                   # Min samples for HDBSCAN (None = auto)
    cluster_selection_epsilon = "0.0"    # Cluster selection epsilon for HDBSCAN
    check_feature_presence = False       # Check feature presence in train-test splits
    filter_by_cluster_presence = True    # Filter features by cluster presence across train/test
    min_cluster_presence = "2"           # Min clusters a feature must appear in
    
    # Pre-processing feature clustering parameters
    use_feature_clustering = False       # Enable pre-processing cluster-based feature filtering
    feature_cluster_method = "hierarchical"  # Pre-processing clustering method
    feature_n_clusters = "20"           # Number of clusters for pre-processing feature clustering
    feature_min_cluster_presence = "2"  # Min clusters a feature must appear in during pre-processing
    
    # K-mer specific parameters
    duplicate_all = True                 # Duplicate all genomes when duplicates found
    
    # Debug options
    dry_run = False                      # Create scripts but don't submit jobs
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "kmer_bootstrap_submit.py",
        
        # Required arguments
        "--input_strain_dir", input_strain_dir,
        "--input_phage_dir", input_phage_dir,
        "--interaction_matrix", interaction_matrix,
        "--output_dir", output_dir,
        
        # SLURM configuration  
        "--account", account,
        "--partition", partition,
        "--qos", qos,
        "--environment", environment,
        "--mem_per_job", mem_per_job,
        "--time_limit", time_limit,
        
        # Bootstrap parameters
        "--n_iterations", n_iterations,
        "--threads", threads,
        "--strain_column", strain_column,
        "--max_ram", max_ram,
        
        # Feature selection and modeling
        "--num_runs_fs", num_runs_fs,
        "--num_runs_modeling", num_runs_modeling,
        "--weights_method", weights_method,
        
        # Clustering for feature selection
        "--cluster_method", cluster_method,
        "--n_clusters", n_clusters,
        "--min_cluster_size", min_cluster_size,
        "--cluster_selection_epsilon", cluster_selection_epsilon,
        "--min_cluster_presence", min_cluster_presence,
        
        # Pre-processing feature clustering
        "--feature_cluster_method", feature_cluster_method,
        "--feature_n_clusters", feature_n_clusters,
        "--feature_min_cluster_presence", feature_min_cluster_presence,
    ]
    
    # Add optional arguments
    if min_samples:
        cmd.extend(["--min_samples", str(min_samples)])
    
    # Add boolean flags
    if use_dynamic_weights:
        cmd.append("--use_dynamic_weights")
    if use_clustering:
        cmd.append("--use_clustering")
    if check_feature_presence:
        cmd.append("--check_feature_presence")
    if filter_by_cluster_presence:
        cmd.append("--filter_by_cluster_presence")
    if use_feature_clustering:
        cmd.append("--use_feature_clustering")
    if duplicate_all:
        cmd.append("--duplicate_all")
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # SUBMIT WORKFLOW
    # =============================================
    print("=" * 80)
    print("K-mer Bootstrap Cross-Validation Workflow Submission")
    print("=" * 80)
    print(f"Input strain dir:    {input_strain_dir}")
    print(f"Input phage dir:     {input_phage_dir}")
    print(f"Interaction matrix:  {interaction_matrix}")
    print(f"Output directory:    {output_dir}")
    print()
    print(f"Number of iterations: {n_iterations}")
    print(f"Feature sel runs:     {num_runs_fs}")
    print(f"Modeling runs:        {num_runs_modeling}")
    print(f"Use dynamic weights:  {use_dynamic_weights}")
    print(f"Use clustering:       {use_clustering}")
    print(f"Cluster method:       {cluster_method}")
    print(f"Filter by cluster:    {filter_by_cluster_presence}")
    print(f"Use feature clustering: {use_feature_clustering}")
    print(f"Duplicate all:        {duplicate_all}")
    print()
    print(f"SLURM account:       {account}")
    print(f"Environment:         {environment}")
    print(f"Memory per job:      {mem_per_job}GB")
    print(f"Time limit per job:  {time_limit}")
    print()
    
    if dry_run:
        print("🧪 DRY RUN MODE - Scripts will be created but not submitted")
        print()
    
    # Estimate costs and runtime
    estimated_hours = float(time_limit.split(":")[0])
    estimated_cost = int(n_iterations) * estimated_hours * 2.0  # Rough estimate at $2/hour per job
    
    print(f"📊 ESTIMATED RESOURCE USAGE:")
    print(f"   Parallel jobs: {n_iterations} (one per iteration)")
    print(f"   Max runtime: {time_limit} per job (all run in parallel)")
    print(f"   Total core-hours: ~{int(n_iterations) * estimated_hours * int(threads)}")
    print(f"   Estimated cost: ~${estimated_cost:.0f}")
    print()
    
    print("Submitting workflow with command:")
    print(" ".join(cmd))
    print()
    
    # Validate paths before submission
    validation_errors = []
    if not os.path.exists(input_strain_dir):
        validation_errors.append(f"Strain directory not found: {input_strain_dir}")
    if not os.path.exists(input_phage_dir):
        validation_errors.append(f"Phage directory not found: {input_phage_dir}")
    if not os.path.exists(interaction_matrix):
        validation_errors.append(f"Interaction matrix not found: {interaction_matrix}")
    
    if validation_errors:
        print("❌ VALIDATION ERRORS:")
        for error in validation_errors:
            print(f"   {error}")
        return 1
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created in kmer_bootstrap_run_* directory")
        else:
            print("\n✅ K-mer bootstrap validation workflow submitted successfully!")
            print("\n📋 Monitor progress with:")
            print("   squeue -u $USER")
            print("   tail -f kmer_bootstrap_run_*/logs/kmer_bootstrap_*.out")
            print("\n⏱️  Expected workflow completion:")
            print(f"   All {n_iterations} iterations will run in parallel")
            print(f"   Total wall time: ~{time_limit}")
            print(f"   Final results: {output_dir}/final_predictions.csv")
            print("\n💡 Tips:")
            print("   - Each iteration processes ~10% of data as validation set")
            print("   - Uses k-mer analysis for more detailed feature characterization")
            print("   - Includes duplicate protein ID detection and resolution")
            print("   - Results from all iterations are automatically aggregated")
            print("   - Individual iteration results saved in iteration_N/ subdirectories")
            print("   - K-mer results in kmer_modeling/ subdirectories")
            print("   - Use 'scancel <job_id>' to cancel if needed")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())