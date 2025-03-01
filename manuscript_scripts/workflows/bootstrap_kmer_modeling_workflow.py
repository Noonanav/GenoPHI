import os
import pandas as pd
import argparse
import logging
import time
import psutil
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments
from phage_modeling.workflows.feature_annotations_workflow import run_predictive_proteins_workflow
from phage_modeling.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table
from phage_modeling.workflows.kmer_table_workflow import run_kmer_table_workflow 

# Set up logging
def setup_logging(output_dir):
    """
    Set up logging to both console and a log file.
    Ensures duplicate handlers are not added.
    """
    log_file = os.path.join(output_dir, "combined_workflow.log")
    os.makedirs(output_dir, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Prevent duplicate handlers
    if not logger.hasHandlers():
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(console_handler)

        # File handler
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

    logger.info("Logging has been set up. Writing logs to: %s", log_file)

# Main combined workflow
def run_kmer_iteration(
    phenotype_matrix,
    output,
    num_features='none',
    filter_type='none',
    num_runs_fs=10,
    num_runs_modeling=20,
    sample_column='strain',
    phenotype_column='interaction',
    method='rfe',
    protein_id_col='protein_ID',
    task_type='classification',
    max_features='none',
    max_ram=8,
    threads=4,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=True,
    cluster_method='hdbscan',
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    use_shap=False,
    clear_tmp=False,
    k=5,
    k_range=False,
    remove_suffix=False,
    one_gene=False,
    ignore_families=False,
    modeling=False
):
    os.makedirs(output, exist_ok=True)
    setup_logging(output)

    # Start overall workflow timing
    logging.info("Starting combined workflow...")
    start_time = time.time()

    # Step 2: Run K-mer Workflow
    logging.info("Running k-mer feature table workflow...")
    kmer_output_dir = os.path.join(output, "kmer_modeling")
    os.makedirs(kmer_output_dir, exist_ok=True)

    # Hardcoded paths derived from protein family workflow outputs
    metrics_file = os.path.join(output, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
    performance_df = pd.read_csv(metrics_file)
    top_cutoff = performance_df.iloc[0]['cut_off'].split('_')[-1]
    feature_file_path = os.path.join(output, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_cutoff_{top_cutoff}.csv')

    strain_fasta = os.path.join(output, "modeling_results", "model_performance", "predictive_proteins", "strain_predictive_AA_seqs.faa")
    protein_csv = os.path.join(output, "modeling_results", "model_performance", "predictive_proteins", "strain_predictive_feature_overview.csv")
    phage_fasta = os.path.join(output, "modeling_results", "model_performance", "predictive_proteins", "phage_predictive_AA_seqs.faa")
    protein_csv_phage = os.path.join(output, "modeling_results", "model_performance", "predictive_proteins", "phage_predictive_feature_overview.csv")

    # Validate required files
    required_files = {
        'strain_fasta': strain_fasta,
        'protein_csv': protein_csv,
        'phage_fasta': phage_fasta,
        'protein_csv_phage': protein_csv_phage
    }

    for file_desc, file_path in required_files.items():
        if not os.path.exists(file_path):
            logging.error(f"Required file {file_desc} not found at: {file_path}")
            raise FileNotFoundError(f"{file_desc} missing: {file_path}")

    run_kmer_table_workflow(
        strain_fasta=strain_fasta,
        protein_csv=protein_csv,
        k=k,
        id_col='strain',
        one_gene=one_gene,
        output_dir=kmer_output_dir,
        k_range=k_range,
        phenotype_matrix=phenotype_matrix,
        phage_fasta=phage_fasta,
        protein_csv_phage=protein_csv_phage,
        remove_suffix=remove_suffix,
        sample_column=sample_column,
        phenotype_column=phenotype_column,
        modeling=modeling,
        filter_type=filter_type,
        num_features=num_features,
        num_runs_fs=num_runs_fs,
        num_runs_modeling=num_runs_modeling,
        method=method,
        strain_list=feature_file_path,
        phage_list=feature_file_path,
        threads=threads,
        task_type=task_type,
        max_features=max_features,
        ignore_families=ignore_families,
        max_ram=max_ram,
        use_shap=use_shap,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon
    )

    # Final combined report
    total_runtime = time.time() - start_time
    logging.info(f"Combined workflow completed in {total_runtime:.2f} seconds.")

def run_kmer_iterations(
    phenotype_matrix,
    modeling_dir,
    num_features='none',
    filter_type='none',
    num_runs_fs=10,
    num_runs_modeling=20,
    sample_column='strain',
    phenotype_column='interaction',
    method='rfe',
    protein_id_col='protein_ID',
    task_type='classification',
    max_features='none',
    max_ram=8,
    threads=4,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=True,
    cluster_method='hdbscan',
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    use_shap=False,
    clear_tmp=False,
    k=5,
    k_range=False,
    remove_suffix=False,
    one_gene=False,
    ignore_families=False,
    modeling=False
):
    iterations = os.listdir(modeling_dir)
    iterations = [x for x in iterations if 'iteration' in x]

    for iteration in iterations:
        iteration_dir = os.path.join(modeling_dir, iteration)

        try:
            run_kmer_iteration(
                phenotype_matrix=phenotype_matrix,
                output=iteration_dir,
                num_features=num_features,
                filter_type=filter_type,
                num_runs_fs=num_runs_fs,
                num_runs_modeling=num_runs_modeling,
                sample_column=sample_column,
                phenotype_column=phenotype_column,
                method=method,
                protein_id_col=protein_id_col,
                task_type=task_type,
                max_features=max_features,
                max_ram=max_ram,
                threads=threads,
                use_dynamic_weights=use_dynamic_weights,
                weights_method=weights_method,
                use_clustering=use_clustering,
                cluster_method=cluster_method,
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                cluster_selection_epsilon=cluster_selection_epsilon,
                use_shap=use_shap,
                clear_tmp=clear_tmp,
                k=k,
                k_range=k_range,
                remove_suffix=remove_suffix,
                one_gene=one_gene,
                ignore_families=ignore_families,
                modeling=modeling
            )
        except FileNotFoundError as e:
            logging.error(f"Skipping {iteration} due to missing files: {str(e)}")
            continue

# CLI
def main():
    parser = argparse.ArgumentParser(description="Run combined protein family and k-mer workflows with reporting.")
    
    # Input data
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('--phenotype_matrix', required=True, help='Phenotype matrix file path.')

    # Optional input arguments
    optional_input_group = parser.add_argument_group('Optional input arguments')
    optional_input_group.add_argument('--sample_column', default='strain', help='Sample column name.')
    optional_input_group.add_argument('--phenotype_column', default='interaction', help='Phenotype column name.')
    optional_input_group.add_argument('--protein_id_col', default='protein_ID', help='Protein ID column name.')
    optional_input_group.add_argument('--use_shap', action='store_true', help='Use SHAP values for analysis (default: False).')
    optional_input_group.add_argument('--clear_tmp', action='store_true', help='Clear temporary files after workflow.')

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('--modeling_dir', required=True, help='Output directory.')

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--filter_type', default='none', help='Filter type for feature selection.')
    fs_modeling_group.add_argument('--method', default='rfe', choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'], help='Feature selection method.')
    fs_modeling_group.add_argument('--num_features', default='none', help='Number of features for selection.')
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection runs.')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=20, help='Number of modeling runs.')
    fs_modeling_group.add_argument('--task_type', default='classification', choices=['classification', 'regression'], help='Task type for modeling.')
    fs_modeling_group.add_argument('--max_features', default='none', help='Maximum number of features.')
    fs_modeling_group.add_argument('--use_dynamic_weights', action='store_true', help='Use dynamic weights for feature selection and modeling.')
    fs_modeling_group.add_argument('--weights_method', default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help='Method to calculate class weights (default: log10)')
    fs_modeling_group.add_argument('--use_clustering', action='store_true', help='Use clustering for feature selection.')
    fs_modeling_group.add_argument('--cluster_method', default='hdbscan', choices=['hdbscan', 'hierarchical'], help='Clustering method for feature selection.')
    fs_modeling_group.add_argument('--min_cluster_size', type=int, default=5, help='Minimum cluster size for HDBSCAN clustering (default: 5)')
    fs_modeling_group.add_argument('--min_samples', type=int, default=None, help='Min samples parameter for HDBSCAN')
    fs_modeling_group.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help='Cluster selection epsilon for HDBSCAN (default: 0.0)')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads to use.')
    general_group.add_argument('--max_ram', type=float, default=8, help='Maximum RAM usage in GB.')
    general_group.add_argument('--k', type=int, default=5, help='K-mer length.')
    general_group.add_argument('--k_range', action='store_true', help='Use range of k-mer lengths.')
    general_group.add_argument('--remove_suffix', action='store_true', help='Remove suffix from genome names.')
    general_group.add_argument('--one_gene', action='store_true', help='Include features with one gene.')
    general_group.add_argument('--ignore_families', action='store_true', help='Ignore protein families.')
    general_group.add_argument('--modeling', action='store_true', help='Run modeling workflow.')

    args = parser.parse_args()

    run_kmer_iterations(
        phenotype_matrix=args.phenotype_matrix,
        modeling_dir=args.modeling_dir,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        method=args.method,
        protein_id_col=args.protein_id_col,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram,
        threads=args.threads,
        use_shap=args.use_shap,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        clear_tmp=args.clear_tmp,
        k=args.k,
        k_range=args.k_range,
        remove_suffix=args.remove_suffix,
        one_gene=args.one_gene,
        ignore_families=args.ignore_families,
        modeling=args.modeling
    )

if __name__ == "__main__":
    main()
