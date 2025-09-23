#!/usr/bin/env python3
"""
Hybrid workflow that combines protein family clustering for strains with k-mer features for phages.
Based on the strain-only CV framework from slurm_bootstrap_submit.py.
"""

import os
import logging
import pandas as pd
import shutil
from Bio import SeqIO
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables
from phage_modeling.kmer_table_workflow import construct_feature_table, get_genome_assignments_tables, feature_selection_optimized, feature_assignment
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments

def run_hybrid_protein_family_kmer_workflow(
    input_path_strain, 
    input_path_phage,
    phenotype_matrix,
    output_dir, 
    tmp_dir="tmp",
    # K-mer parameters for phages
    k=5,
    one_gene=True,
    k_range=False,
    ignore_families=False,
    # MMseqs2 parameters for strains  
    min_seq_id=0.6, 
    coverage=0.8, 
    sensitivity=7.5, 
    suffix='faa', 
    threads=4,
    # Strain/phage lists
    strain_list='none', 
    phage_list='none',
    strain_column='strain', 
    phage_column='phage', 
    # Feature selection and modeling
    num_features='none', 
    filter_type='strain', 
    num_runs_fs=10, 
    num_runs_modeling=10,
    sample_column='strain', 
    phenotype_column='interaction', 
    method='rfe',
    task_type='classification', 
    max_features='none', 
    max_ram=8, 
    use_dynamic_weights=False, 
    weights_method='log10',
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    filter_by_cluster_presence=False,
    min_cluster_presence=2, 
    use_shap=False, 
    bootstrapping=False,
    clear_tmp=False,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    use_augmentation=False,
    augmentation_strain_fraction=0.01,
    augmentation_phage_fraction=0.01,
    augmentation_fold_increase=3
):
    """
    Hybrid workflow combining protein family clustering for strains with k-mer features for phages.
    
    Args:
        input_path_strain (str): Path to strain FASTA files
        input_path_phage (str): Path to phage FASTA files  
        phenotype_matrix (str): Path to interaction matrix
        output_dir (str): Output directory
        k (int): K-mer length for phages
        one_gene (bool): Include single-gene features for k-mers
        k_range (bool): Generate k-mer range from 3 to k
        ignore_families (bool): Ignore protein families in k-mer generation
        ... (other parameters same as protein_family_workflow)
    
    Returns:
        str: Path to final merged feature table
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    logging.info("Starting hybrid protein family (strains) + k-mer (phages) workflow...")
    
    # Step 1: Process strains with MMseqs2 protein family clustering
    logging.info("Step 1: Processing strains with protein family clustering...")
    strain_output_dir = os.path.join(output_dir, 'strain')
    strain_tmp_dir = os.path.join(tmp_dir, 'strain')
    
    # Run MMseqs2 clustering for strains
    run_clustering_workflow(
        input_path=input_path_strain,
        output_dir=strain_output_dir,
        tmp_dir=strain_tmp_dir,
        min_seq_id=min_seq_id,
        coverage=coverage,
        sensitivity=sensitivity,
        suffix=suffix,
        threads=threads,
        strain_list=strain_list,
        strain_column=strain_column,
        compare=False,
        bootstrapping=bootstrapping,
        clear_tmp=clear_tmp
    )
    
    # Generate strain features from clustering results
    run_feature_assignment(
        input_file=os.path.join(strain_output_dir, "presence_absence_matrix.csv"),
        output_dir=os.path.join(strain_output_dir, 'features'),
        source='strain',
        select='none',
        select_column=strain_column,
        input_type='directory',
        max_ram=max_ram,
        threads=threads
    )
    
    # Step 2: Process phages with k-mer generation
    logging.info("Step 2: Processing phages with k-mer features...")
    phage_output_dir = os.path.join(output_dir, 'phage')
    phage_tmp_dir = os.path.join(tmp_dir, 'phage')
    phage_feature_dir = os.path.join(phage_output_dir, 'features')
    os.makedirs(phage_feature_dir, exist_ok=True)
    
    # Load phage list if provided
    phage_genome_list = None
    if phage_list != 'none' and phage_list:
        if os.path.exists(phage_list):
            phage_df = pd.read_csv(phage_list)
            if phage_column in phage_df.columns:
                phage_genome_list = list(phage_df[phage_column].unique())
                logging.info(f"Loaded {len(phage_genome_list)} phages from phage list")
    
    # Create a minimal protein CSV for k-mer workflow compatibility
    # This is needed because construct_feature_table expects protein feature data
    protein_csv_phage = os.path.join(phage_tmp_dir, 'phage_proteins.csv')
    os.makedirs(phage_tmp_dir, exist_ok=True)
    
    # Generate protein CSV from FASTA files
    create_phage_protein_csv(input_path_phage, protein_csv_phage, suffix, phage_genome_list)
    
    # Generate k-mer feature table for phages
    phage_feature_table_path = construct_feature_table(
        fasta_file=input_path_phage,
        protein_csv=protein_csv_phage,
        k=k,
        id_col='phage',
        one_gene=one_gene,
        output_dir=phage_feature_dir,
        output_name="phage",
        k_range=k_range,
        ignore_families=ignore_families,
        genome_list=phage_genome_list
    )
    
    # Process k-mer features similar to protein families
    phage_feature_table = pd.read_csv(phage_feature_table_path)
    
    # Get genome assignments for phages
    phage_genome_assignments = get_genome_assignments_tables(
        phage_feature_table, 'phage', phage_feature_dir, prefix='phage'
    )
    
    # Optimize feature selection for phages
    phage_selected_features = feature_selection_optimized(
        phage_feature_table, "phage", 'phage', phage_feature_dir, prefix='phage'
    )
    
    # Assign features to phage genomes  
    phage_assignment_df, phage_final_feature_table, phage_final_feature_table_output = feature_assignment(
        phage_genome_assignments, phage_selected_features, 'phage', 
        phage_feature_dir, prefix='phage', all_genomes=phage_genome_list
    )
    
    # Step 3: Merge strain and phage features with phenotype matrix
    logging.info("Step 3: Merging strain and phage features...")
    strain_feature_table_path = os.path.join(strain_output_dir, 'features', 'feature_table.csv')
    
    merged_table_path = merge_feature_tables(
        strain_features=strain_feature_table_path,
        phenotype_matrix=phenotype_matrix,
        output_dir=output_dir,
        sample_column=sample_column,
        phage_features=phage_final_feature_table_output,
        remove_suffix=False,
        use_feature_clustering=use_feature_clustering,
        feature_cluster_method=feature_cluster_method,
        feature_n_clusters=feature_n_clusters,
        feature_min_cluster_presence=feature_min_cluster_presence
    )
    
    logging.info(f"Hybrid feature table created: {merged_table_path}")
    
    # Step 4: Feature selection
    logging.info("Step 4: Running feature selection...")
    base_fs_output_dir = os.path.join(output_dir, 'feature_selection')
    
    # Determine number of features if not specified
    if num_features == 'none':
        merged_df = pd.read_csv(merged_table_path)
        num_interactions = len(merged_df)
        if num_interactions < 500:
            num_features = int(num_interactions / 10)
        else:
            num_features = int(num_interactions / 20)
        num_features = max(10, min(num_features, 1000))  # Reasonable bounds
        logging.info(f"Auto-determined num_features: {num_features}")
    else:
        num_features = int(num_features)
    
    run_feature_selection_iterations(
        input_path=merged_table_path,
        base_output_dir=base_fs_output_dir,
        threads=threads,
        num_features=num_features,
        filter_type=filter_type,
        num_runs=num_runs_fs,
        method=method,
        sample_column=sample_column,
        phenotype_column=phenotype_column,
        phage_column=phage_column,
        task_type=task_type,
        max_ram=max_ram,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        check_feature_presence=check_feature_presence,
        filter_by_cluster_presence=filter_by_cluster_presence,
        min_cluster_presence=min_cluster_presence,
        use_augmentation=use_augmentation,
        augmentation_strain_fraction=augmentation_strain_fraction,
        augmentation_phage_fraction=augmentation_phage_fraction,
        augmentation_fold_increase=augmentation_fold_increase
    )
    
    # Step 5: Generate feature tables
    logging.info("Step 5: Generating filtered feature tables...")
    max_features_int = None if max_features == 'none' else int(max_features)
    filter_table_dir = os.path.join(base_fs_output_dir, 'filtered_feature_tables')
    
    generate_feature_tables(
        model_testing_dir=base_fs_output_dir,
        full_feature_table_file=merged_table_path,
        filter_table_dir=filter_table_dir,
        phenotype_column=phenotype_column,
        sample_column=sample_column,
        cut_offs=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50],
        binary_data=True,  # Use binary for hybrid features
        max_features=max_features_int,
        filter_type=filter_type
    )
    
    # Step 6: Modeling
    logging.info("Step 6: Running modeling experiments...")
    modeling_output_dir = os.path.join(output_dir, 'modeling_results')
    
    run_experiments(
        input_dir=filter_table_dir,
        base_output_dir=modeling_output_dir,
        threads=threads,
        num_runs=num_runs_modeling,
        set_filter=filter_type,
        sample_column=sample_column,
        phage_column=phage_column,
        phenotype_column=phenotype_column,
        task_type=task_type,
        binary_data=True,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        max_ram=max_ram,
        use_shap=use_shap,
        use_augmentation=use_augmentation,
        augmentation_strain_fraction=augmentation_strain_fraction,
        augmentation_phage_fraction=augmentation_phage_fraction,
        augmentation_fold_increase=augmentation_fold_increase
    )
    
    logging.info("Hybrid workflow completed successfully!")
    return merged_table_path

def create_phage_protein_csv(input_path_phage, output_csv, suffix='faa', phage_genome_list=None):
    """
    Create a protein CSV file from phage FASTA files for k-mer workflow compatibility.
    
    Args:
        input_path_phage (str): Path to phage FASTA files or directory
        output_csv (str): Output CSV file path
        suffix (str): FASTA file suffix
        phage_genome_list (list): List of phages to include
    """
    protein_data = []
    
    if os.path.isfile(input_path_phage):
        # Single file - extract phage name from filename
        phage_name = os.path.basename(input_path_phage).replace(f'.{suffix}', '')
        fasta_files = [input_path_phage]
        phage_names = [phage_name]
    else:
        # Directory of files
        all_files = [f for f in os.listdir(input_path_phage) if f.endswith(f'.{suffix}')]
        phage_names = [f.replace(f'.{suffix}', '') for f in all_files]
        fasta_files = [os.path.join(input_path_phage, f) for f in all_files]
    
    # Filter phages if list provided
    if phage_genome_list:
        filtered_files = []
        filtered_names = []
        for fasta_file, phage_name in zip(fasta_files, phage_names):
            if phage_name in phage_genome_list:
                filtered_files.append(fasta_file)
                filtered_names.append(phage_name)
        fasta_files = filtered_files
        phage_names = filtered_names
        logging.info(f"Filtered to {len(fasta_files)} phage files from genome list")
    
    cluster_counter = 0
    for fasta_file, phage_name in zip(fasta_files, phage_names):
        try:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                protein_data.append({
                    'protein_ID': record.id,
                    'phage': phage_name,
                    'Feature': f'phage_feature_{cluster_counter}',  # Dummy feature name
                    'cluster': f'cluster_{cluster_counter}'  # Dummy cluster
                })
                cluster_counter += 1
        except Exception as e:
            logging.warning(f"Error processing {fasta_file}: {e}")
    
    if not protein_data:
        logging.error("No protein data found - check FASTA files and phage list")
        raise ValueError("No proteins found in phage FASTA files")
    
    # Create DataFrame and save
    protein_df = pd.DataFrame(protein_data)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    protein_df.to_csv(output_csv, index=False)
    logging.info(f"Created protein CSV with {len(protein_df)} proteins for {len(set(protein_df['phage']))} phages")

def hybrid_assign_predict_workflow(
    input_dir, genome_list, mmseqs_db, clusters_tsv, feature_map,
    tmp_dir, suffix, model_dir, feature_table, 
    phage_feature_table_path, output_dir, threads, genome_type='strain',
    sensitivity=7.5, coverage=0.8, min_seq_id=0.4, duplicate_all=True,
    # K-mer parameters for phage validation
    k=5, filtered_kmers=None, aa_sequence_file=None, threshold=0.5,
    # Additional hybrid parameters
    strain_feature_table_path=None, reuse_existing=True
):
    """
    Hybrid assignment and prediction workflow for validation.
    
    This function coordinates between protein family features (strains) and k-mer features (phages)
    to perform hybrid feature assignment and prediction for cross-validation.
    
    For strain-based CV with hybrid features:
    - Validation strains: assign protein family features via MMseqs2
    - Phages: use existing k-mer feature table (all phages used in training)
    - Prediction: combine both feature types
    
    For phage-based CV with hybrid features:
    - Strains: use existing protein family feature table
    - Validation phages: assign k-mer features
    - Prediction: combine both feature types
    
    Args:
        input_dir (str): Directory with validation genome FASTA files
        genome_list (str): Path to validation genome list CSV
        mmseqs_db (str): Path to training MMseqs2 database
        clusters_tsv (str): Path to training clusters file
        feature_map (str): Path to training feature map
        tmp_dir (str): Temporary directory for processing
        suffix (str): FASTA file suffix (default: 'faa')
        model_dir (str): Directory containing trained models
        feature_table (str): Path to filtered feature table for prediction
        phage_feature_table_path (str): Path to phage k-mer feature table
        output_dir (str): Output directory for results
        threads (int): Number of threads for processing
        genome_type (str): 'strain' or 'phage' - type of validation genomes
        sensitivity (float): MMseqs2 sensitivity parameter
        coverage (float): MMseqs2 coverage parameter  
        min_seq_id (float): MMseqs2 minimum sequence identity
        duplicate_all (bool): Duplicate all genomes in feature table
        k (int): K-mer length for phage feature assignment
        filtered_kmers (str): Path to filtered k-mers file (for phage assignment)
        aa_sequence_file (str): Path to amino acid sequences (for phage assignment)
        threshold (float): K-mer matching threshold (for phage assignment)
        strain_feature_table_path (str): Path to existing strain features (for phage CV)
        reuse_existing (bool): Whether to reuse existing assignment results
    
    Returns:
        str: Path to final prediction results
    """
    import os
    import logging
    import shutil
    from phage_modeling.workflows.assign_predict_workflow import assign_predict_workflow
    from phage_modeling.workflows.kmer_assign_predict_workflow import kmer_assign_predict_workflow
    
    logging.info(f"Starting hybrid assignment and prediction workflow for {genome_type} validation...")
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    hybrid_assign_dir = os.path.join(output_dir, "hybrid_assign_results")
    hybrid_predict_dir = os.path.join(output_dir, "predict_results")
    os.makedirs(hybrid_assign_dir, exist_ok=True)
    os.makedirs(hybrid_predict_dir, exist_ok=True)
    
    # Validate required inputs based on genome type and CV strategy
    if genome_type == 'strain':
        # Strain-based CV: need phage features, will assign strain features
        if phage_feature_table_path is None:
            raise ValueError("phage_feature_table_path is required for strain-based hybrid validation")
        if not os.path.exists(phage_feature_table_path):
            raise FileNotFoundError(f"Phage feature table not found: {phage_feature_table_path}")
            
    elif genome_type == 'phage':
        # Phage-based CV: need strain features and k-mer assignment parameters
        if strain_feature_table_path is None:
            raise ValueError("strain_feature_table_path is required for phage-based hybrid validation")
        if not os.path.exists(strain_feature_table_path):
            raise FileNotFoundError(f"Strain feature table not found: {strain_feature_table_path}")
        if filtered_kmers is None or aa_sequence_file is None:
            raise ValueError("filtered_kmers and aa_sequence_file are required for phage k-mer assignment")
    
    # Paths for assigned features
    assigned_feature_table_path = os.path.join(hybrid_assign_dir, f'{genome_type}_combined_feature_table.csv')
    final_prediction_path = os.path.join(hybrid_predict_dir, f"{genome_type}_median_predictions.csv")
    
    # Skip assignment if results already exist and reuse is enabled
    if reuse_existing and os.path.exists(assigned_feature_table_path):
        logging.info(f"Found existing assigned features: {assigned_feature_table_path}")
    else:
        # Perform feature assignment based on genome type
        if genome_type == 'strain':
            logging.info("Assigning protein family features to validation strains...")
            
            # Use standard protein family assignment workflow for strains
            run_assign_features_workflow(
                input_dir=input_dir,
                mmseqs_db=mmseqs_db,
                tmp_dir=os.path.join(tmp_dir, 'strain_assign'),
                output_dir=hybrid_assign_dir,
                feature_map=feature_map,
                clusters_tsv=clusters_tsv,
                genome_type='strain',
                genome_list=genome_list,
                sensitivity=sensitivity,
                coverage=coverage,
                min_seq_id=min_seq_id,
                threads=threads,
                suffix=suffix,
                duplicate_all=duplicate_all
            )
            
        elif genome_type == 'phage':
            logging.info("Assigning k-mer features to validation phages...")
            
            # Use k-mer assignment workflow for phages
            from phage_modeling.workflows.kmer_assign_features_workflow import run_kmer_assign_features_workflow
            
            run_kmer_assign_features_workflow(
                input_dir=input_dir,
                mmseqs_db=mmseqs_db,
                tmp_dir=os.path.join(tmp_dir, 'phage_assign'),
                output_dir=hybrid_assign_dir,
                feature_map=feature_map,
                filtered_kmers=filtered_kmers,
                aa_sequence_file=aa_sequence_file,
                clusters_tsv=clusters_tsv,
                genome_type='phage',
                genome_list=genome_list,
                sensitivity=sensitivity,
                coverage=coverage,
                min_seq_id=min_seq_id,
                threads=threads,
                suffix=suffix,
                threshold=threshold,
                reuse_existing=reuse_existing
            )
    
    # Verify assignment completed successfully
    if not os.path.exists(assigned_feature_table_path):
        raise FileNotFoundError(f"Feature assignment failed - no output found: {assigned_feature_table_path}")
    
    logging.info(f"Feature assignment completed: {assigned_feature_table_path}")
    
    # Skip prediction if results already exist and reuse is enabled
    if reuse_existing and os.path.exists(final_prediction_path):
        logging.info(f"Found existing prediction results: {final_prediction_path}")
        return final_prediction_path
    
    # Set up prediction with hybrid features
    logging.info("Starting hybrid prediction workflow...")
    
    # Create prediction input directory with proper feature table structure
    predict_input_dir = os.path.join(output_dir, "prediction_input")
    os.makedirs(predict_input_dir, exist_ok=True)
    
    if genome_type == 'strain':
        # Strain-based CV: newly assigned strain features + existing phage features
        prediction_strain_features = os.path.join(predict_input_dir, "strain_feature_table.csv") 
        prediction_phage_features = phage_feature_table_path
        
        # Copy newly assigned strain features to prediction input
        if not os.path.exists(prediction_strain_features):
            shutil.copy2(assigned_feature_table_path, prediction_strain_features)
            logging.info("Copied newly assigned strain features for hybrid prediction")
            
    elif genome_type == 'phage':
        # Phage-based CV: existing strain features + newly assigned phage features
        prediction_strain_features = os.path.join(predict_input_dir, "strain_feature_table.csv")
        prediction_phage_features = assigned_feature_table_path
        
        # Copy existing strain features to prediction input
        if not os.path.exists(prediction_strain_features):
            shutil.copy2(strain_feature_table_path, prediction_strain_features)
            logging.info("Copied existing strain features for hybrid prediction")
    
    # Run prediction workflow with hybrid features
    from phage_modeling.workflows.prediction_workflow import run_prediction_workflow
    
    try:
        run_prediction_workflow(
            input_dir=predict_input_dir,
            phage_feature_table_path=prediction_phage_features,
            model_dir=model_dir,
            feature_table=feature_table,
            output_dir=hybrid_predict_dir,
            strain_source='strain',
            phage_source='phage', 
            threads=threads
        )
        
        logging.info(f"Hybrid prediction completed: {final_prediction_path}")
        
    except Exception as e:
        logging.error(f"Prediction workflow failed: {e}")
        # Provide more detailed error information
        logging.error(f"Prediction input dir: {predict_input_dir}")
        logging.error(f"Strain features: {prediction_strain_features}")
        logging.error(f"Phage features: {prediction_phage_features}")
        logging.error(f"Model dir: {model_dir}")
        logging.error(f"Feature table: {feature_table}")
        raise
    
    # Verify prediction completed successfully
    if not os.path.exists(final_prediction_path):
        raise FileNotFoundError(f"Prediction failed - no output found: {final_prediction_path}")
    
    logging.info(f"Hybrid assignment and prediction workflow completed successfully for {genome_type} validation")
    logging.info(f"Results saved to: {final_prediction_path}")
    
    return final_prediction_path


# Helper function to import the assignment workflow
def run_assign_features_workflow(input_dir, mmseqs_db, tmp_dir, output_dir, feature_map, 
                                clusters_tsv, genome_type, genome_list=None, sensitivity=7.5, 
                                coverage=0.8, min_seq_id=0.4, threads=4, suffix='faa', 
                                duplicate_all=True):
    """
    Helper function to run standard protein family feature assignment.
    This is a wrapper around the existing assign_features_workflow.
    """
    from phage_modeling.workflows.assign_features_workflow import run_assign_features_workflow as run_assign
    
    return run_assign(
        input_dir=input_dir,
        mmseqs_db=mmseqs_db,
        tmp_dir=tmp_dir,
        output_dir=output_dir,
        feature_map=feature_map,
        clusters_tsv=clusters_tsv,
        genome_type=genome_type,
        genome_list=genome_list,
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        threads=threads,
        suffix=suffix,
        duplicate_all=duplicate_all
    )
