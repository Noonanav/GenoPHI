#!/usr/bin/env python3
"""
K-mer bootstrap submission script for ALL FAMILIES mode.
Run on existing protein family workflow results to add k-mer modeling with all families.

Edit the paths and parameters below, then run: python3 kmer_all_families_submit.py
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
    
    # IMPORTANT: This should point to EXISTING output directory with completed protein family workflows
    output_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/ecoli_hierarchical_inv_freq_filter_update_mod"
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "pc_phiml"
    partition = "lr7"
    qos = "lr_normal"
    environment = "phage_modeling"
    mem_per_job = "240"      # Increased for all families mode
    time_limit = "36:00:00"  # Longer time for all families
    
    # =============================================
    # BOOTSTRAP PARAMETERS
    # =============================================
    n_iterations = "20"      # Should match existing iterations in output_dir
    threads = "32"
    strain_column = "strain"
    max_ram = "240"
    
    # Feature selection and modeling runs
    num_runs_fs = "25"
    num_runs_modeling = "50"
    
    # Dynamic weights and clustering (should match original workflow)
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
    
    # Pre-processing feature clustering parameters
    use_feature_clustering = False
    feature_cluster_method = "hierarchical"
    feature_n_clusters = "20"
    feature_min_cluster_presence = "2"
    
    # Prediction parameters
    duplicate_all = True
    
    # =============================================
    # KEY PARAMETER: ALL FAMILIES MODE
    # =============================================
    # Set to True to use ALL protein families for k-mer modeling
    # Set to False to use only predictive families (standard mode)
    use_all_families_for_kmers = True
    
    # Debug options
    dry_run = False  # Set to True to create scripts without submitting
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "kmer_bootstrap_all_families.py",
        
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
    
    # 🆕 KEY FLAG: Use all families for k-mer modeling
    if use_all_families_for_kmers:
        cmd.append("--use_all_families_for_kmers")
    
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # VALIDATION AND SUBMISSION
    # =============================================
    print("=" * 80)
    print("K-mer Bootstrap Validation - ALL FAMILIES Mode")
    print("=" * 80)
    print()
    print("📁 PATHS:")
    print(f"   Input strain dir:    {input_strain_dir}")
    print(f"   Input phage dir:     {input_phage_dir}")
    print(f"   Interaction matrix:  {interaction_matrix}")
    print(f"   Output directory:    {output_dir}")
    print()
    print("🔧 WORKFLOW PARAMETERS:")
    print(f"   Number of iterations: {n_iterations}")
    print(f"   Feature sel runs:     {num_runs_fs}")
    print(f"   Modeling runs:        {num_runs_modeling}")
    print(f"   Use dynamic weights:  {use_dynamic_weights}")
    print(f"   Use clustering:       {use_clustering}")
    print(f"   Cluster method:       {cluster_method}")
    print(f"   Filter by cluster:    {filter_by_cluster_presence}")
    print(f"   Use feature clustering: {use_feature_clustering}")
    print()
    print("🆕 K-MER MODE:")
    if use_all_families_for_kmers:
        print("   ✅ ALL FAMILIES MODE ENABLED")
        print("   Will use ALL protein families for k-mer modeling")
        print("   Output: iteration_N/kmer_modeling_all_families/")
        print("   Final: all_families_final_predictions.csv")
    else:
        print("   📊 PREDICTIVE FAMILIES ONLY (standard)")
        print("   Will use only predictive protein families")
        print("   Output: iteration_N/kmer_modeling/")
        print("   Final: kmer_final_predictions.csv")
    print()
    print("💻 SLURM CONFIGURATION:")
    print(f"   Account:             {account}")
    print(f"   Partition:           {partition}")
    print(f"   QOS:                 {qos}")
    print(f"   Environment:         {environment}")
    print(f"   Memory per job:      {mem_per_job}GB")
    print(f"   Time limit per job:  {time_limit}")
    print()
    
    if dry_run:
        print("🧪 DRY RUN MODE - Scripts will be created but not submitted")
        print()
    
    # Estimate resources
    estimated_hours = float(time_limit.split(":")[0])
    estimated_cost = int(n_iterations) * estimated_hours * 2.0
    
    print("📊 ESTIMATED RESOURCE USAGE:")
    print(f"   Parallel jobs: {n_iterations} (one per iteration)")
    print(f"   Max runtime: {time_limit} per job (all run in parallel)")
    print(f"   Total core-hours: ~{int(n_iterations) * estimated_hours * int(threads)}")
    print(f"   Estimated cost: ~${estimated_cost:.0f}")
    print()
    
    # Validate paths before submission
    validation_errors = []
    
    if not os.path.exists(input_strain_dir):
        validation_errors.append(f"Strain directory not found: {input_strain_dir}")
    if not os.path.exists(input_phage_dir):
        validation_errors.append(f"Phage directory not found: {input_phage_dir}")
    if not os.path.exists(interaction_matrix):
        validation_errors.append(f"Interaction matrix not found: {interaction_matrix}")
    if not os.path.exists(output_dir):
        validation_errors.append(f"Output directory not found: {output_dir}")
        validation_errors.append("  → This script expects EXISTING output with completed protein family workflows")
    
    # Check for iteration directories
    if os.path.exists(output_dir):
        iteration_dirs = [d for d in os.listdir(output_dir) if d.startswith('iteration_')]
        if not iteration_dirs:
            validation_errors.append(f"No iteration directories found in {output_dir}")
            validation_errors.append("  → Expected directories like iteration_1/, iteration_2/, etc.")
        else:
            print(f"✅ VALIDATION: Found {len(iteration_dirs)} iteration directories")
            
            # Check a sample iteration for required files
            sample_iter = os.path.join(output_dir, 'iteration_1')
            if os.path.exists(sample_iter):
                required_files = [
                    'modeling_strains.csv',
                    'validation_strains.csv',
                    'strain/clusters.tsv',
                    'phage/clusters.tsv',
                    'modeling_results/model_performance/model_performance_metrics.csv'
                ]
                missing = []
                for req_file in required_files:
                    full_path = os.path.join(sample_iter, req_file)
                    if not os.path.exists(full_path):
                        missing.append(req_file)
                
                if missing:
                    validation_errors.append(f"Missing required files in {sample_iter}:")
                    for mf in missing:
                        validation_errors.append(f"  → {mf}")
                else:
                    print("✅ VALIDATION: Required files present in iteration_1")
    
    if validation_errors:
        print()
        print("❌ VALIDATION ERRORS:")
        for error in validation_errors:
            print(f"   {error}")
        print()
        print("Please ensure:")
        print("  1. Output directory exists with completed protein family workflows")
        print("  2. All iteration directories (iteration_1/, etc.) are present")
        print("  3. Each iteration has completed modeling results")
        return 1
    
    print()
    print("✅ All validation checks passed!")
    print()
    print("Submitting workflow with command:")
    print(" ".join(cmd))
    print()
    
    # Confirm if all families mode
    if use_all_families_for_kmers and not dry_run:
        print("⚠️  CONFIRMATION REQUIRED")
        print("You are about to run k-mer modeling with ALL protein families.")
        print("This will:")
        print(f"  - Create new directories: iteration_N/kmer_modeling_all_families/")
        print(f"  - Process 10-100x more features than predictive-only mode")
        print(f"  - Take longer and use more memory")
        print(f"  - NOT modify any existing data")
        print()
        response = input("Continue? (yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("Submission cancelled.")
            return 0
        print()
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created in kmer_bootstrap_run_* directory")
            print("Review the scripts and run without --dry_run to submit")
        else:
            output_prefix = "all_families_" if use_all_families_for_kmers else "kmer_"
            
            print("\n✅ K-mer bootstrap workflow submitted successfully!")
            print("\n📋 MONITORING:")
            print("   Monitor all jobs:     squeue -u $USER")
            print("   Monitor this workflow: squeue -u $USER --name=kmer_bootstrap_all")
            print("   View logs:            tail -f kmer_bootstrap_run_*/logs/kmer_bootstrap_*.out")
            print()
            print("📁 OUTPUT LOCATIONS:")
            print(f"   Individual iterations: {output_dir}/iteration_N/kmer_modeling_all_families/")
            print(f"   Final predictions:     {output_dir}/{output_prefix}final_predictions.csv")
            print(f"   Summary statistics:    {output_dir}/{output_prefix}prediction_summary.csv")
            print()
            print("💡 TIPS:")
            print("   - Each iteration will skip completed protein family workflow steps")
            print("   - Will generate all_AA_seqs.faa from existing clusters.tsv")
            print("   - K-mer modeling will use ALL protein families (not just predictive)")
            print("   - Results automatically aggregated across all iterations")
            print("   - Use 'scancel <job_id>' to cancel if needed")
            print()
            print("📖 For detailed safety information, see: SAFETY_REVIEW.md")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
