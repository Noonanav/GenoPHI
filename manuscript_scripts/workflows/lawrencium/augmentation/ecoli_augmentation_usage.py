#!/usr/bin/env python3
"""
Augmentation testing workflow submission script.
Edit the paths and parameters below, then run: python3 augmentation_usage_example.py
"""

import subprocess
import sys
import os

def main():
    # =============================================
    # YOUR PATHS - EDIT THESE
    # =============================================
    cv_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/augmentation/ecoli_modeling_hierarchical"
    input_strain_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/strain_AAs"
    input_phage_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/phage_AAs"
    
   # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"                    # SLURM partition 
    qos = "lr_normal"                    # SLURM QOS
    environment = "phage_modeling"       # Conda environment name
    mem_per_job = "240"                   # Memory per iteration job in GB
    time_limit = "36:00:00"              # Time limit per iteration job
    
    # =============================================
    # AUGMENTATION TEST PARAMETERS
    # =============================================
    
    # Augmentation fractions to test
    fractions = [0.005, 0.01, 0.02, 0.05, 0.1]  # Different fractions to test
    fold_increase = 4                            # Fold increase for augmentation
    
    # CV parameters (auto-detected if not specified)
    n_iterations = None                          # None = auto-detect from CV directory
    feature_table_path = "merged/full_feature_table.csv"  # Path within each iteration directory
    
    # Core parameters
    threads = "24"                       # Number of threads per job
    max_ram = "60"                       # Max RAM usage in GB
    
    # Feature selection and modeling runs
    num_runs_fs = "25"                   # Number of feature selection runs per test
    num_runs_modeling = "50"             # Number of modeling runs per test
    num_features = "100"                 # Number of features to select
    
    # Dynamic weights and clustering for feature selection
    use_dynamic_weights = True           # Use dynamic class weights
    weights_method = "inverse_frequency" # Weight calculation method
    use_clustering = True                # Use clustering for feature selection
    cluster_method = "hierarchical"      # Clustering method
    n_clusters = "20"                    # Number of clusters for hierarchical clustering
    min_cluster_size = "3"               # Minimum cluster size for HDBSCAN
    min_samples = None                   # Min samples for HDBSCAN (None = auto)
    cluster_selection_epsilon = "0.0"    # Cluster selection epsilon for HDBSCAN
    check_feature_presence = False        # Check feature presence in train-test splits
    filter_by_cluster_presence = True    # Filter features by cluster presence
    min_cluster_presence = "2"           # Min clusters feature must be present in
    bootstrapping = True                 # Enable bootstrapping for feature selection
    
    # Prediction parameters
    duplicate_all = True                 # Duplicate all genomes in feature table for predictions
    
    # Debug options
    dry_run = False                      # Create scripts but don't submit jobs
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "augmentation_submit.py",
        
        # Required arguments
        "--cv_dir", cv_dir,
        "--input_strain_dir", input_strain_dir,
        "--input_phage_dir", input_phage_dir,
        
        # SLURM configuration  
        "--account", account,
        "--partition", partition,
        "--qos", qos,
        "--environment", environment,
        "--mem_per_job", mem_per_job,
        "--time_limit", time_limit,
        
        # Augmentation parameters
        "--fractions"] + [str(f) for f in fractions] + [
        "--fold_increase", str(fold_increase),
        "--feature_table_path", feature_table_path,
        
        # Modeling parameters
        "--threads", threads,
        "--max_ram", max_ram,
        "--num_runs_fs", num_runs_fs,
        "--num_runs_modeling", num_runs_modeling,
        "--num_features", num_features,
        
        # Feature selection configuration
        "--weights_method", weights_method,
        "--cluster_method", cluster_method,
        "--n_clusters", n_clusters,
        "--min_cluster_size", min_cluster_size,
        "--cluster_selection_epsilon", cluster_selection_epsilon,
        "--min_cluster_presence", min_cluster_presence,
    ]
    
    # Add optional arguments
    if n_iterations:
        cmd.extend(["--n_iterations", str(n_iterations)])
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
    if bootstrapping:
        cmd.append("--bootstrapping")
    if duplicate_all:
        cmd.append("--duplicate_all")
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # SUBMIT WORKFLOW
    # =============================================
    print("=" * 80)
    print("Augmentation Testing Workflow Submission")
    print("=" * 80)
    print(f"CV directory:        {cv_dir}")
    print(f"Input strain dir:    {input_strain_dir}")
    print(f"Input phage dir:     {input_phage_dir}")
    print()
    print(f"Augmentation fractions: {fractions}")
    print(f"Fold increase:          {fold_increase}")
    print(f"Feature sel runs:       {num_runs_fs}")
    print(f"Modeling runs:          {num_runs_modeling}")
    print(f"Use dynamic weights:    {use_dynamic_weights}")
    print(f"Use clustering:         {use_clustering}")
    print(f"Cluster method:         {cluster_method}")
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
    n_fractions = len(fractions)
    
    print(f"📊 ESTIMATED RESOURCE USAGE:")
    print(f"   Augmentation fractions: {n_fractions}")
    print(f"   CV iterations: Auto-detected from {cv_dir}")
    print(f"   Max runtime: {time_limit} per fraction-iteration combo")
    print(f"   Total jobs: {n_fractions} × n_iterations (parallel execution)")
    print()
    
    print("Submitting workflow with command:")
    print(" ".join(cmd))
    print()
    
    # Validate paths before submission
    validation_errors = []
    if not os.path.exists(cv_dir):
        validation_errors.append(f"CV directory not found: {cv_dir}")
    if not os.path.exists(input_strain_dir):
        validation_errors.append(f"Strain directory not found: {input_strain_dir}")
    if not os.path.exists(input_phage_dir):
        validation_errors.append(f"Phage directory not found: {input_phage_dir}")
    
    # Check for at least one iteration directory
    iteration_dirs = [d for d in os.listdir(cv_dir) if d.startswith('iteration_')]
    if not iteration_dirs:
        validation_errors.append(f"No iteration directories found in: {cv_dir}")
    
    if validation_errors:
        print("❌ VALIDATION ERRORS:")
        for error in validation_errors:
            print(f"   {error}")
        return 1
    
    print(f"✅ Found {len(iteration_dirs)} iteration directories in CV output")
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created in augmentation_test_run_* directory")
        else:
            print("\n✅ Augmentation testing workflow submitted successfully!")
            print("\n📋 Monitor progress with:")
            print("   squeue -u $USER")
            print("   tail -f augmentation_test_run_*/logs/aug_*.out")
            print("\n⏱️  Expected workflow completion:")
            print(f"   All fraction-iteration combinations will run in parallel")
            print(f"   Total wall time: ~{time_limit}")
            print(f"   Final results: {cv_dir}/augmentation_test_summary.csv")
            print("\n💡 Tips:")
            print("   - Each job tests one fraction on one CV iteration")
            print("   - Results automatically aggregated across all tests")
            print("   - Individual results saved in iteration_N/augmentation/frac_*pct/")
            print("   - Use 'scancel <job_id>' to cancel if needed")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())