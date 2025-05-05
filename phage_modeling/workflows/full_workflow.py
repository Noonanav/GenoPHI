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
from phage_modeling.workflows.protein_family_workflow import run_protein_family_workflow
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


def monitor_resources(start_time, func, *args, **kwargs):
    """
    Wrapper to monitor resource usage (CPU, RAM) and execution time for a function.
    """
    process = psutil.Process()
    cpu_usage = []
    ram_usage = 0

    def monitor_cpu():
        while True:
            cpu_usage.append(psutil.cpu_percent(interval=1))
            if stop_monitoring.is_set():
                break

    from threading import Event, Thread
    stop_monitoring = Event()
    monitor_thread = Thread(target=monitor_cpu)
    monitor_thread.start()

    try:
        func_start = time.time()
        result = func(*args, **kwargs)
        ram_usage = process.memory_info().rss  # Max RAM in bytes
        return result, {
            'runtime': time.time() - func_start,
            'avg_cpu': sum(cpu_usage) / len(cpu_usage) if cpu_usage else 0,
            'max_cpu': max(cpu_usage) if cpu_usage else 0,
            'max_ram': ram_usage / (1024 ** 3)  # Convert to GB
        }
    finally:
        stop_monitoring.set()
        monitor_thread.join()

def write_section_report(section, metrics, output_dir):
    """Append section metrics to the workflow report."""
    report_file = os.path.join(output_dir, "workflow_section_metrics.csv")
    if not os.path.exists(report_file):
        with open(report_file, "w") as f:
            f.write("Section,Runtime (s),Avg CPU (%),Max CPU (%),Max RAM (GB)\n")
    with open(report_file, "a") as f:
        f.write(f"{section},{metrics['runtime']:.2f},{metrics['avg_cpu']:.2f},{metrics['max_cpu']:.2f},{metrics['max_ram']:.2f}\n")

# Main combined workflow
def run_full_workflow(
    input_strain,
    input_phage,
    phenotype_matrix,
    output,
    clustering_dir=None,
    min_seq_id=0.4,
    coverage=0.8,
    sensitivity=7.5,
    suffix='faa',
    strain_list='none',
    phage_list='none',
    strain_column='strain',
    phage_column='phage',
    source_strain='strain',
    source_phage='phage',
    compare=False,
    num_features='none',
    filter_type='none',
    num_runs_fs=10,
    num_runs_modeling=20,
    sample_column='strain',
    phenotype_column='interaction',
    method='rfe',
    annotation_table_path=None,
    protein_id_col='protein_ID',
    task_type='classification',
    max_features='none',
    max_ram=8,
    threads=4,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
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

    metrics = {}
    
    # Step 1: Run Full Workflow
    logging.info("Running full protein family workflow...")
    _, metrics['protein_family'] = monitor_resources(start_time, run_protein_family_workflow,
                                                     input_path_strain=input_strain,
                                                     input_path_phage=input_phage,
                                                     phenotype_matrix=phenotype_matrix,
                                                     output_dir=output,
                                                     clustering_dir=clustering_dir,
                                                     min_seq_id=min_seq_id,
                                                     coverage=coverage,
                                                     sensitivity=sensitivity,
                                                     suffix=suffix,
                                                     strain_list=strain_list,
                                                     phage_list=phage_list,
                                                     strain_column=strain_column,
                                                     phage_column=phage_column,
                                                     source_strain=source_strain,
                                                     source_phage=source_phage,
                                                     compare=compare,
                                                     num_features=num_features,
                                                     filter_type=filter_type,
                                                     num_runs_fs=num_runs_fs,
                                                     num_runs_modeling=num_runs_modeling,
                                                     sample_column=sample_column,
                                                     phenotype_column=phenotype_column,
                                                     method=method,
                                                     annotation_table_path=annotation_table_path,
                                                     protein_id_col=protein_id_col,
                                                     task_type=task_type,
                                                     max_features=max_features,
                                                     max_ram=max_ram,
                                                     threads=threads,
                                                     use_shap=use_shap,
                                                     use_dynamic_weights=use_dynamic_weights,
                                                     weights_method=weights_method,
                                                     use_clustering=use_clustering,
                                                     cluster_method=cluster_method,
                                                     n_clusters=n_clusters,
                                                     min_cluster_size=min_cluster_size,
                                                     min_samples=min_samples,
                                                     cluster_selection_epsilon=cluster_selection_epsilon,
                                                     check_feature_presence=check_feature_presence,
                                                     clear_tmp=clear_tmp)
    write_section_report("Protein_Family_Workflow", metrics['protein_family'], output)

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

    _, metrics['kmer_workflow'] = monitor_resources(start_time, run_kmer_table_workflow,
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
                                                    n_clusters=n_clusters,
                                                    min_cluster_size=min_cluster_size,
                                                    min_samples=min_samples,
                                                    cluster_selection_epsilon=cluster_selection_epsilon,
                                                    check_feature_presence=check_feature_presence)
    write_section_report("Kmer_Workflow", metrics['kmer_workflow'], output)

    # Final combined report
    total_runtime = time.time() - start_time
    logging.info(f"Combined workflow completed in {total_runtime:.2f} seconds.")
    with open(os.path.join(output, "combined_workflow_summary.txt"), "w") as report:
        report.write("Combined Workflow Summary\n")
        report.write("=" * 40 + "\n")
        for section, m in metrics.items():
            report.write(f"{section}: Runtime = {m['runtime']:.2f}s, Avg CPU = {m['avg_cpu']:.2f}%, Max CPU = {m['max_cpu']:.2f}%, Max RAM = {m['max_ram']:.2f}GB\n")
        report.write(f"Total Workflow Runtime: {total_runtime:.2f}s\n")

# CLI
def main():
    parser = argparse.ArgumentParser(description="Run combined protein family and k-mer workflows with reporting.")
    
    # Input data
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('--input_strain', required=True, help='Input strain FASTA path.')
    input_group.add_argument('--input_phage', required=True, help='Input phage FASTA path.')
    input_group.add_argument('--phenotype_matrix', required=True, help='Phenotype matrix file path.')

    # Optional input arguments
    optional_input_group = parser.add_argument_group('Optional input arguments')
    optional_input_group.add_argument('--suffix', default='faa', help='Suffix for input FASTA files (default: faa).')
    optional_input_group.add_argument('--strain_list', default='none', help='List of strains for filtering (default: none).')
    optional_input_group.add_argument('--phage_list', default='none', help='List of phages for filtering (default: none).')
    optional_input_group.add_argument('--strain_column', default='strain', help='Column name for strain data (default: strain).')
    optional_input_group.add_argument('--phage_column', default='phage', help='Column name for phage data (default: phage).')
    optional_input_group.add_argument('--source_strain', default='strain', help='Source prefix for strain (default: strain).')
    optional_input_group.add_argument('--source_phage', default='phage', help='Source prefix for phage (default: phage).')
    optional_input_group.add_argument('--sample_column', default='strain', help='Sample column name (default: strain).')
    optional_input_group.add_argument('--phenotype_column', default='interaction', help='Phenotype column name (default: interaction).')
    optional_input_group.add_argument('--annotation_table_path', help='Path to annotation table.')
    optional_input_group.add_argument('--protein_id_col', default='protein_ID', help='Protein ID column name (default: protein_ID).')
    optional_input_group.add_argument('--use_shap', action='store_true', help='Use SHAP values for analysis (default: False).')
    optional_input_group.add_argument('--clear_tmp', action='store_true', help='Clear temporary files after workflow.')

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('--output', required=True, help='Output directory.')

    # Clustering parameters
    clustering_group = parser.add_argument_group('Clustering')
    clustering_group.add_argument('--clustering_dir', help='Path to an existing strain clustering directory.')
    clustering_group.add_argument('--min_seq_id', type=float, default=0.4, help='Minimum sequence identity for clustering (default: 0.4).')
    clustering_group.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering (default: 0.8).')
    clustering_group.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering (default: 7.5).')
    clustering_group.add_argument('--compare', action='store_true', help='Compare clustering results.')

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--filter_type', default='none', choices=['none', 'strain', 'phage'], help='Filter type for feature selection. (default: none)')
    fs_modeling_group.add_argument('--method', default='rfe', choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                                   help='Feature selection method (default: rfe)')
    fs_modeling_group.add_argument('--num_features', default='none', help='Number of features for selection (default: none).')
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection runs (default: 10).')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=20, help='Number of modeling runs (default: 20).')
    fs_modeling_group.add_argument('--task_type', default='classification', choices=['classification', 'regression'], help='Task type for modeling (default: classification)')
    fs_modeling_group.add_argument('--max_features', default='none', help='Maximum number of features for modeling (default: none).')
    fs_modeling_group.add_argument('--use_dynamic_weights', action='store_true', help='Use dynamic weights for feature selection and modeling.')
    fs_modeling_group.add_argument('--weights_method', default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help='Method to calculate class weights (default: log10)')
    fs_modeling_group.add_argument('--use_clustering', action='store_true', help='Use clustering for feature selection.')
    fs_modeling_group.add_argument('--cluster_method', default='hierarchical', choices=['hdbscan', 'hierarchical'], help='Clustering method for feature selection (default: hierarchical)')
    fs_modeling_group.add_argument('--n_clusters', type=int, default=20, help='Number of clusters for hierarchical clustering (default: 20)')
    fs_modeling_group.add_argument('--min_cluster_size', type=int, default=5, help='Minimum cluster size for HDBSCAN clustering (default: 5)')
    fs_modeling_group.add_argument('--min_samples', type=int, default=None, help='Min samples parameter for HDBSCAN')
    fs_modeling_group.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help='Cluster selection epsilon for HDBSCAN (default: 0.0)')
    fs_modeling_group.add_argument('--check_feature_presence', action='store_true', help='Check for presence of features during train-test split.')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')
    general_group.add_argument('--max_ram', type=float, default=8, help='Maximum RAM usage in GB (default: 8).')
    general_group.add_argument('--k', type=int, default=5, help='K-mer length (default: 5).')
    general_group.add_argument('--k_range', action='store_true', help='Use range of k-mer lengths.')
    general_group.add_argument('--remove_suffix', action='store_true', help='Remove suffix from genome names.')
    general_group.add_argument('--one_gene', action='store_true', help='Include features with one gene.')
    general_group.add_argument('--ignore_families', action='store_true', help='Ignore protein families.')
    general_group.add_argument('--modeling', action='store_true', help='Run modeling workflow.')

    args = parser.parse_args()

    run_full_workflow(
        input_strain=args.input_strain,
        input_phage=args.input_phage,
        phenotype_matrix=args.phenotype_matrix,
        output=args.output,
        clustering_dir=args.clustering_dir,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        compare=args.compare,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        method=args.method,
        annotation_table_path=args.annotation_table_path,
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
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
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
