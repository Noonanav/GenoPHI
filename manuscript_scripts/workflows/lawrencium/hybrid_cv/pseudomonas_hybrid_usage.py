#!/usr/bin/env python3
"""
Hybrid bootstrap cross-validation workflow submission script.
Uses protein families for STRAINS and k-mer features for PHAGES.
Edit the paths and parameters below, then run: python3 bootstrap_usage_hybrid.py
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
    output_dir = "/global/scratch/users/anoonan/BRaVE/pseudomonas/hybrid_bootstrapping/inverse_frequency_hierarchical_filter"

    # Optional: Use existing clustering results to speed up workflow
    clustering_dir = None  # Set to path if you have pre-computed clustering results
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"
    qos = "lr_normal"
    environment = "phage_modeling"
    mem_per_job = "80"
    time_limit = "24:00:00"
    
    # =============================================
    # BOOTSTRAP VALIDATION PARAMETERS
    # =============================================
    
    # Core parameters
    n_iterations = "20"
    threads = "32"
    strain_column = "strain"
    max_ram = "60"
    
    # Feature selection and modeling runs
    num_runs_fs = "25"
    num_runs_modeling = "50"
    
    # Protein family parameters (for strains)
    min_seq_id = "0.4"
    coverage = "0.8"
    sensitivity = "7.5"
    
    # K-mer parameters (for phages)
    k = "5"
    k_range = True
    one_gene = True
    ignore_families = True
    
    # Dynamic weights and clustering for feature selection
    use_dynamic_weights = True
    weights_method = "inverse_frequency"
    use_clustering = True
    cluster_method = "hierarchical"
    n_clusters = "20"
    min_cluster_size = "3"
    min_samples = None
    cluster_selection_epsilon = "0.0"
    check_feature_presence = False
    filter_by_cluster_presence = True
    min_cluster_presence = "2"
    bootstrapping = True
    
    # Feature clustering (pre-processing)
    use_feature_clustering = True
    feature_cluster_method = "hierarchical"
    feature_n_clusters = "20"
    feature_min_cluster_presence = "2"
    
    # Prediction parameters
    duplicate_all = True
    
    # Debug options
    dry_run = False
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "hybrid_slurm.py",
        
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
        
        # Protein family parameters
        "--min_seq_id", min_seq_id,
        "--coverage", coverage,
        "--sensitivity", sensitivity,
        
        # K-mer parameters
        "--k", k,
        
        # Clustering for feature selection
        "--cluster_method", cluster_method,
        "--n_clusters", n_clusters,
        "--min_cluster_size", min_cluster_size,
        "--cluster_selection_epsilon", cluster_selection_epsilon,
        "--min_cluster_presence", min_cluster_presence,
        
        # Feature clustering parameters
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
    if k_range:
        cmd.append("--k_range")
    if one_gene:
        cmd.append("--one_gene")
    if ignore_families:
        cmd.append("--ignore_families")
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
    print("Hybrid Bootstrap Cross-Validation Workflow Submission")
    print("=" * 80)
    print(f"Input strain dir:    {input_strain_dir}")
    print(f"Input phage dir:     {input_phage_dir}")
    print(f"Interaction matrix:  {interaction_matrix}")
    print(f"Output directory:    {output_dir}")
    print(f"Clustering dir:      {clustering_dir if clustering_dir else 'None (will compute)'}")
    print()
    print(f"HYBRID FEATURE APPROACH:")
    print(f"  Strain features:   Protein families (sc_*)")
    print(f"  Phage features:    K-mers (pc_*), k={k}")
    print()
    print(f"Number of iterations: {n_iterations}")
    print(f"Feature sel runs:     {num_runs_fs}")
    print(f"Modeling runs:        {num_runs_modeling}")
    print(f"Use dynamic weights:  {use_dynamic_weights}")
    print(f"Use clustering:       {use_clustering}")
    print(f"Cluster method:       {cluster_method}")
    print(f"Use feature clustering: {use_feature_clustering}")
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
    estimated_cost = int(n_iterations) * estimated_hours * 2.0
    
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
            print("Scripts created in hybrid_bootstrap_run_* directory")
        else:
            print("\n✅ Hybrid bootstrap validation workflow submitted successfully!")
            print("\n📋 Monitor progress with:")
            print("   squeue -u $USER")
            print("   tail -f hybrid_bootstrap_run_*/logs/bootstrap_*.out")
            print("\n⏱️  Expected workflow completion:")
            print(f"   All {n_iterations} iterations will run in parallel")
            print(f"   Total wall time: ~{time_limit}")
            print(f"   Final results: {output_dir}/final_predictions.csv")
            print("\n💡 Tips:")
            print("   - Strain features: Protein families via MMseqs2 clustering")
            print("   - Phage features: K-mer based for all phages")
            print("   - Each iteration validates on ~10% of strains")
            print("   - All phages used in both training and validation")
            print("   - Results automatically aggregated from all iterations")
            print("   - Use 'scancel <job_id>' to cancel if needed")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
