#!/usr/bin/env python3
"""
Complete SLURM workflow submission script with all parameters.
Edit the paths and parameters below, then run: python3 submit_my_workflow.py
"""

import subprocess
import sys

def main():
    # =============================================
    # YOUR PATHS - EDIT THESE
    # =============================================
    input_strain = "/global/scratch/users/anoonan/BRaVE/ecoli/strain_AAs_update"
    input_phage = "/global/scratch/users/anoonan/BRaVE/ecoli/phage_AAs"
    phenotype_matrix = "/global/scratch/users/anoonan/BRaVE/ecoli/Gaborieau_interaction_matrix_long_mod.csv"
    output_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/ecoli_phage_hierarchical_inv_freq_filter_update_mod"
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"                    # SLURM partition 
    qos = "lr_normal"                    # SLURM QOS
    environment = "phage_modeling"       # Conda environment name
    
    # =============================================
    # WORKFLOW PARAMETERS (with defaults)
    # =============================================
    
    # Clustering parameters
    clustering_dir = None                # Path to existing clustering directory (None = run new clustering)
    min_seq_id = "0.4"                 # Minimum sequence identity for clustering
    coverage = "0.8"                   # Minimum coverage for clustering  
    sensitivity = "7.5"                # Sensitivity for clustering
    suffix = "faa"                     # File suffix for FASTA files
    compare = False                     # Compare clustering results
    
    # Input filtering
    strain_list = phenotype_matrix               # List of strains for filtering ('none' = use all)
    phage_list = phenotype_matrix                 # List of phages for filtering ('none' = use all)
    strain_column = "strain"            # Column name for strain identifiers
    phage_column = "phage"              # Column name for phage identifiers  
    source_strain = "strain"            # Source prefix for strain features
    source_phage = "phage"              # Source prefix for phage features
    sample_column = "strain"            # Sample column name in phenotype matrix
    phenotype_column = "interaction"    # Phenotype column name
    
    # Feature selection parameters
    filter_type = "phage"                # Filter type ('none', 'strain', 'phage')
    method = "rfe"                      # Feature selection method ('rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap')
    num_features = "none"               # Number of features to select ('none' = auto-determine)
    num_runs_fs = "25"                  # Number of feature selection runs
    max_features = "none"               # Maximum features for modeling ('none' = no limit)
    
    # Modeling parameters  
    num_runs_modeling = "50"            # Number of modeling runs
    task_type = "classification"        # Task type ('classification' or 'regression')
    use_dynamic_weights = True         # Use dynamic class weights
    weights_method = "inverse_frequency"            # Weight calculation method ('log10', 'inverse_frequency', 'balanced')
    
    # Clustering for feature selection
    use_clustering = True               # Use clustering for feature selection
    cluster_method = "hierarchical"          # Clustering method ('hdbscan' or 'hierarchical')  
    n_clusters = "20"                   # Number of clusters for hierarchical clustering
    min_cluster_size = "5"              # Minimum cluster size for HDBSCAN
    min_samples = None                  # Min samples for HDBSCAN (None = auto)
    cluster_selection_epsilon = "0.0"   # Cluster selection epsilon for HDBSCAN
    check_feature_presence = False      # Check feature presence in train-test splits
    filter_by_cluster_presence = True    # Filter features by cluster presence across train/test
    min_cluster_presence = "2"          # NEW: Minimum number of clusters/groups a feature must be present in
    
    # NEW: Pre-processing feature clustering parameters
    use_feature_clustering = False       # Enable pre-processing cluster-based feature filtering
    feature_cluster_method = "hierarchical"  # Pre-processing clustering method ('hierarchical' only for now)
    feature_n_clusters = "20"           # Number of clusters for pre-processing feature clustering
    feature_min_cluster_presence = "2"  # Min clusters a feature must appear in during pre-processing
    
    # K-mer parameters
    k = "6"                             # K-mer length
    k_range = False                     # Use range of k-mer lengths (3 to k)
    one_gene = False                    # Include features with only one gene
    ignore_families = False              # Ignore protein families in k-mer analysis
    remove_suffix = False               # Remove suffix from genome names
    
    # Analysis options
    use_shap = True                    # Calculate SHAP values (increases runtime)
    annotation_table_path = None        # Path to annotation table (None = not used)
    protein_id_col = "protein_ID"       # Protein ID column name
    
    # Resource parameters
    threads = "32"                       # Number of threads per job (will be overridden by SLURM allocation)
    max_ram = "120"                       # Maximum RAM usage in GB
    
    # Cleanup
    clear_tmp = False                   # Clear temporary files after workflow
    modeling = True                     # Run modeling workflow (should be True)
    
    # Debug options
    dry_run = False                     # Create scripts but don't submit jobs
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "slurm_full_modeling_submit_big.py",
        
        # Required arguments
        "--input_strain", input_strain,
        "--input_phage", input_phage,
        "--phenotype_matrix", phenotype_matrix,
        "--output", output_dir,
        
        # SLURM configuration  
        "--account", account,
        "--partition", partition,
        "--qos", qos,
        "--environment", environment,
        
        # Clustering parameters
        "--min_seq_id", min_seq_id,
        "--coverage", coverage,
        "--sensitivity", sensitivity,
        "--suffix", suffix,
        "--strain_list", strain_list,
        "--phage_list", phage_list,
        "--strain_column", strain_column,
        "--phage_column", phage_column,
        "--source_strain", source_strain,
        "--source_phage", source_phage,
        
        # Feature selection parameters
        "--filter_type", filter_type,
        "--method", method,
        "--num_features", num_features,
        "--num_runs_fs", num_runs_fs,
        "--max_features", max_features,
        
        # Modeling parameters
        "--num_runs_modeling", num_runs_modeling,
        "--sample_column", sample_column,
        "--phenotype_column", phenotype_column,
        "--task_type", task_type,
        "--weights_method", weights_method,
        
        # Clustering for feature selection
        "--cluster_method", cluster_method,
        "--n_clusters", n_clusters,
        "--min_cluster_size", min_cluster_size,
        "--cluster_selection_epsilon", cluster_selection_epsilon,
        "--min_cluster_presence", min_cluster_presence,  # NEW: Added this parameter
        
        # NEW: Pre-processing feature clustering parameters
        "--feature_cluster_method", feature_cluster_method,
        "--feature_n_clusters", feature_n_clusters,
        "--feature_min_cluster_presence", feature_min_cluster_presence,
        
        # K-mer parameters
        "--k", k,
        "--protein_id_col", protein_id_col,
        
        # Resource parameters
        "--threads", threads,
        "--max_ram", max_ram,
    ]
    
    # Add optional arguments that have None values
    if clustering_dir:
        cmd.extend(["--clustering_dir", clustering_dir])
    if annotation_table_path:
        cmd.extend(["--annotation_table_path", annotation_table_path])
    if min_samples:
        cmd.extend(["--min_samples", str(min_samples)])
    
    # Add boolean flags
    if compare:
        cmd.append("--compare")
    if use_dynamic_weights:
        cmd.append("--use_dynamic_weights")
    if use_clustering:
        cmd.append("--use_clustering")
    if use_feature_clustering:  # NEW
        cmd.append("--use_feature_clustering")
    if filter_by_cluster_presence:  # NEW: Added this boolean flag
        cmd.append("--filter_by_cluster_presence")
    if k_range:
        cmd.append("--k_range")
    if one_gene:
        cmd.append("--one_gene")
    if ignore_families:
        cmd.append("--ignore_families")
    if remove_suffix:
        cmd.append("--remove_suffix")
    if use_shap:
        cmd.append("--use_shap")
    if check_feature_presence:
        cmd.append("--check_feature_presence")
    if clear_tmp:
        cmd.append("--clear_tmp")
    if modeling:
        cmd.append("--modeling")
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # SUBMIT WORKFLOW
    # =============================================
    print("=" * 60)
    print("SLURM Workflow Submission")
    print("=" * 60)
    print(f"Input strain:      {input_strain}")
    print(f"Input phage:       {input_phage}")
    print(f"Phenotype matrix:  {phenotype_matrix}")
    print(f"Output directory:  {output_dir}")
    print(f"SLURM account:     {account}")
    print(f"Environment:       {environment}")
    print(f"Feature sel runs:  {num_runs_fs}")
    print(f"Modeling runs:     {num_runs_modeling}")
    print(f"K-mer length:      {k}")
    print(f"Use clustering:    {use_clustering}")
    print(f"Use feature clustering: {use_feature_clustering}")  # NEW
    if use_feature_clustering:  # NEW
        print(f"Feature cluster method: {feature_cluster_method}")
        print(f"Feature n_clusters:     {feature_n_clusters}")
        print(f"Feature min presence:   {feature_min_cluster_presence}")
    print(f"Filter by cluster presence: {filter_by_cluster_presence}")  # NEW
    print(f"Min cluster presence:       {min_cluster_presence}")        # NEW
    print(f"Ignore families:   {ignore_families}")
    print()
    
    if dry_run:
        print("🧪 DRY RUN MODE - Scripts will be created but not submitted")
        print()
    
    print("Submitting workflow with command:")
    print(" ".join(cmd))
    print()
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created in slurm_run_* directory")
        else:
            print("\n✅ Workflow submitted successfully!")
            print("\nMonitor progress with:")
            print("  squeue -u $USER")
            print("  tail -f slurm_run_*/logs/stage*_*.out")
            print("\nExpected runtime: 8-12 hours total (6 sequential stages)")
            print("Estimated cost: ~$200-300")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())