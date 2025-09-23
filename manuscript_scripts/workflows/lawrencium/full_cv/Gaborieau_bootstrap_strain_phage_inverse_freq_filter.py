#!/usr/bin/env python3
"""
Bootstrap cross-validation workflow submission script for STRAIN+PHAGE-based cross-validation.
This version leaves out 10% of strains AND 10% of phages for validation, testing generalization
to completely novel strain-phage combinations.
Edit the paths and parameters below, then run: python3 bootstrap_strain_phage_usage.py
"""

import subprocess
import sys
import os

def main():
    # =============================================
    # YOUR PATHS - EDIT THESE
    # =============================================
    input_strain_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/strain_AAs_update"
    input_phage_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/phage_AAs"
    interaction_matrix = "/global/scratch/users/anoonan/BRaVE/ecoli/Gaborieau_interaction_matrix_long_mod.csv"
    output_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/full_bootstrapping/inverse_frequency_hierarchical_feature_filter"
    
    # Optional: Use existing clustering results to speed up workflow
    clustering_dir = None  # Set to path if you have pre-computed clustering results
    # clustering_dir = "/global/scratch/users/anoonan/EDGE/full_workflow_bootstrap"
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"                    # SLURM partition 
    qos = "lr_normal"                    # SLURM QOS
    environment = "phage_modeling"       # Conda environment name
    mem_per_job = "180"                  # Memory per iteration job in GB (higher for strain+phage processing)
    time_limit = "48:00:00"              # Time limit per iteration job (longer for more complex workflow)
    
    # =============================================
    # BOOTSTRAP VALIDATION PARAMETERS
    # =============================================
    
    # Core parameters
    n_iterations = "20"                  # Number of cross-validation iterations (20 for 20-fold CV)
    threads = "32"                       # Number of threads per job
    strain_column = "strain"             # Column in interaction matrix with strain names
    phage_column = "phage"               # Column in interaction matrix with phage names
    max_ram = "180"                       # Max RAM usage in GB
    
    # Feature selection and modeling runs
    num_runs_fs = "25"                   # Number of feature selection runs per iteration
    num_runs_modeling = "50"             # Number of modeling runs per iteration
    
    # Filter type can be either 'strain' or 'phage' - both are used in this workflow
    filter_type = "strain"               # Primary filter type for feature selection
    
    # Dynamic weights and clustering for feature selection
    use_dynamic_weights = True           # Use dynamic class weights
    weights_method = "inverse_frequency" # Weight calculation method ('log10', 'inverse_frequency', 'balanced')
    use_clustering = True                # Use clustering for feature selection
    cluster_method = "hierarchical"      # Clustering method ('hdbscan' or 'hierarchical')
    n_clusters = "20"                    # Number of clusters for hierarchical clustering
    min_cluster_size = "3"               # Minimum cluster size for HDBSCAN
    min_samples = None                   # Min samples for HDBSCAN (None = auto)
    cluster_selection_epsilon = "0.0"    # Cluster selection epsilon for HDBSCAN
    check_feature_presence = False        # Check feature presence in train-test splits
    filter_by_cluster_presence = False    # Filter features by cluster presence across train/test
    min_cluster_presence = "2"     
    bootstrapping = True                  # Enable bootstrapping for feature selection
    
    # Pre-processing feature clustering parameters
    use_feature_clustering = True        # Enable pre-processing cluster-based feature filtering
    feature_cluster_method = "hierarchical"  # Pre-processing clustering method ('hierarchical' only for now)
    feature_n_clusters = "20"            # Number of clusters for pre-processing feature clustering
    feature_min_cluster_presence = "2"   # Min clusters a feature must appear in during pre-processing
    
    # Prediction parameters
    duplicate_all = True                 # Duplicate all genomes in feature table for predictions
    
    # Debug options
    dry_run = False                      # Create scripts but don't submit jobs
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "slurm_bootstrap_submit_strain_phage_cv.py",
        
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
        "--phage_column", phage_column,
        "--max_ram", max_ram,
        
        # Feature selection and modeling
        "--num_runs_fs", num_runs_fs,
        "--num_runs_modeling", num_runs_modeling,
        "--weights_method", weights_method,
        "--filter_type", filter_type,
        
        # Clustering for feature selection
        "--cluster_method", cluster_method,
        "--n_clusters", n_clusters,
        "--min_cluster_size", min_cluster_size,
        "--cluster_selection_epsilon", cluster_selection_epsilon,
        "--min_cluster_presence", min_cluster_presence,
        
        # Pre-processing feature clustering parameters
        "--feature_cluster_method", feature_cluster_method,
        "--feature_n_clusters", feature_n_clusters,
        "--feature_min_cluster_presence", feature_min_cluster_presence,
    ]
    
    # Add optional arguments
    if clustering_dir:
        cmd.extend(["--clustering_dir", clustering_dir])
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
    if use_feature_clustering:
        cmd.append("--use_feature_clustering")
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # SUBMIT WORKFLOW
    # =============================================
    print("=" * 80)
    print("Bootstrap Cross-Validation Workflow Submission - STRAIN+PHAGE-BASED CV")
    print("=" * 80)
    print(f"Input strain dir:    {input_strain_dir}")
    print(f"Input phage dir:     {input_phage_dir}")
    print(f"Interaction matrix:  {interaction_matrix}")
    print(f"Output directory:    {output_dir}")
    print(f"Clustering dir:      {clustering_dir if clustering_dir else 'None (will compute)'}")
    print()
    print(f"CV TYPE:             STRAIN+PHAGE-based (leaving out BOTH strains AND phages)")
    print(f"Filter type:         {filter_type}")
    print(f"Number of iterations: {n_iterations}")
    print(f"Feature sel runs:     {num_runs_fs}")
    print(f"Modeling runs:        {num_runs_modeling}")
    print(f"Use dynamic weights:  {use_dynamic_weights}")
    print(f"Use clustering:       {use_clustering}")
    print(f"Cluster method:       {cluster_method}")
    print(f"Use feature clustering: {use_feature_clustering}")
    if use_feature_clustering:
        print(f"Feature cluster method: {feature_cluster_method}")
        print(f"Feature n_clusters:     {feature_n_clusters}")
        print(f"Feature min presence:   {feature_min_cluster_presence}")
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
    print(f"   Expected data retention: ~81% per iteration (90% strains × 90% phages)")
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
    if clustering_dir and not os.path.exists(clustering_dir):
        validation_errors.append(f"Clustering directory not found: {clustering_dir}")
    
    if validation_errors:
        print("❌ VALIDATION ERRORS:")
        for error in validation_errors:
            print(f"   {error}")
        return 1
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created in bootstrap_run_* directory")
        else:
            print("\n✅ Bootstrap validation workflow submitted successfully!")
            print("\n📋 Monitor progress with:")
            print("   squeue -u $USER")
            print("   tail -f bootstrap_run_*/logs/bootstrap_*.out")
            print("\n⏱️  Expected workflow completion:")
            print(f"   All {n_iterations} iterations will run in parallel")
            print(f"   Total wall time: ~{time_limit}")
            print(f"   Final results: {output_dir}/final_predictions.csv")
            print("\n💡 Tips:")
            print("   - Each iteration processes ~10% of STRAINS + ~10% of PHAGES as validation set")
            print("   - Training uses 90% strains × 90% phages interactions only")  
            print("   - Validation strains and phages get features assigned from training clustering")
            print("   - Models predict interactions between validation strains × validation phages")
            print("   - This tests generalization to COMPLETELY NOVEL strain-phage combinations")
            print("   - Expected ~81% data retention per iteration (some interactions excluded)")
            print("   - Results from all iterations are automatically aggregated")
            print("   - Individual iteration results saved in iteration_N/ subdirectories")
            print("   - Data retention statistics saved in data_retention_stats.csv")
            print("   - Use 'scancel <job_id>' to cancel if needed")
            print("\n🎯 What this tests:")
            print("   - Can the model generalize to strain-phage pairs it has NEVER seen?")
            print("   - Neither the strain NOR the phage was in the training data")
            print("   - Most stringent test of model generalization capability")
            print("   - Complementary to strain-only and phage-only cross-validation")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())