import os
import argparse
import sys
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables

def check_method_task_type_compatibility(method, task_type):
    """
    Checks if the selected feature selection method is compatible with the task type.

    Args:
        method (str): Feature selection method.
        task_type (str): Task type ('classification' or 'regression').

    Raises:
        ValueError: If the method and task_type are incompatible.
    """
    incompatible_methods = {
        'classification': [],
        'regression': ['chi_squared']
    }
    if method in incompatible_methods.get(task_type, []):
        raise ValueError(f"The feature selection method '{method}' is not compatible with task type '{task_type}'.")

def run_feature_selection_workflow(
    input_path, 
    base_output_dir, 
    threads=4, 
    num_features=500, 
    filter_type='none',
    num_runs=50, method='rfe', 
    task_type='classification', 
    phenotype_column=None,
    sample_column='strain', 
    phage_column='phage',
    binary_data=False, 
    max_features='none', 
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    max_ram=8):
    """
    Workflow for running feature selection iterations and generating feature tables.

    Args:
        input_path (str): Path to the input feature table.
        base_output_dir (str): Directory to save results.
        threads (int): Number of threads to use.
        num_features (int): Number of features to select.
        filter_type (str): Filter type for the input data ('strain', 'phage', 'none').
        num_runs (int): Number of runs to perform.
        method (str): Feature selection method ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').
        task_type (str): Task type for modeling ('classification' or 'regression').
        phenotype_column (str or None): Column name for the phenotype or target variable.
        sample_column (str or None): Column name for sample identifiers.
        phage_column (str): Name of the phage column (default: 'phage').
        binary_data (bool): If True, converts feature values to binary (0/1).
        max_features (str): Maximum number of features to include in feature tables.
        use_dynamic_weights (bool): Whether to use dynamic weights for feature selection.
        weights_method (str): Method for calculating weights ('log10', 'inverse_frequency', 'balanced').
        use_clustering (bool): Whether to use clustering for filtering.
        cluster_method (str): Clustering method to use ('hdbscan' or 'hierarchical').
        n_clusters (int): Number of clusters for hierarchical clustering (default: 20).
        min_cluster_size (int): Minimum cluster size for HDBSCAN clustering.
        min_samples (int): Minimum number of samples for HDBSCAN (default: None for same as min_cluster_size).
        cluster_selection_epsilon (float): Epsilon value for HDBSCAN clustering.
        check_feature_presence (bool): If True, only include features present in both train and test sets.
        max_ram (int): Maximum RAM to use in GB.
    """
    # Check compatibility of method and task type
    try:
        check_method_task_type_compatibility(method, task_type)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Run multiple iterations of feature selection
    print("Running feature selection iterations...")
    run_feature_selection_iterations(
        input_path=input_path,
        base_output_dir=base_output_dir,
        threads=threads,
        num_features=num_features,
        filter_type=filter_type,
        num_runs=num_runs,
        method=method,
        phenotype_column=phenotype_column,
        sample_column=sample_column,
        phage_column=phage_column,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        task_type=task_type,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        check_feature_presence=check_feature_presence,
        max_ram=max_ram,
    )
    
    # Generate feature tables based on the results
    print("Generating feature tables...")
    max_features = None if max_features == 'none' else int(max_features)
    filter_table_dir = os.path.join(base_output_dir, 'filtered_feature_tables')
    generate_feature_tables(
        model_testing_dir=base_output_dir,
        full_feature_table_file=input_path,
        filter_table_dir=filter_table_dir,
        phenotype_column=phenotype_column,
        sample_column=sample_column,
        cut_offs=[3, 4, 5, 6, 7, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50],
        binary_data=binary_data,
        max_features=max_features,
        filter_type=filter_type
    )

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run feature selection workflow.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input path for the full feature table.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Base output directory for the results.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use.')
    parser.add_argument('--num_features', default='none', help='Number of features to select during feature selection.')
    parser.add_argument('--filter_type', type=str, default='none', help="Type of filtering to use ('none', 'strain', 'phage').")
    parser.add_argument('--num_runs', type=int, default=50, help='Number of feature selection iterations to run.')
    parser.add_argument('--method', type=str, default='rfe', choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                        help="Feature selection method to use ('rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').")
    parser.add_argument('--task_type', type=str, required=True, choices=['classification', 'regression'],
                        help="Task type for model ('classification' or 'regression').")
    parser.add_argument('--phenotype_column', type=str, help='Optional column name for phenotype/target variable.')
    parser.add_argument('--sample_column', type=str, default='strain', help='Optional column name for sample identifiers.')
    parser.add_argument('--phage_column', type=str, default='phage', help='Optional column name for phage identifiers.')
    parser.add_argument('--binary_data', action='store_true', help='If set, converts feature values to binary (1/0); otherwise, continuous values are kept.')
    parser.add_argument('--max_features', default='none', help='Maximum number of features to include in the feature tables.')
    parser.add_argument('--use_dynamic_weights', action='store_true', help='If set, uses dynamic weights for feature selection.')
    parser.add_argument('--weights_method', type=str, default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help='Method for calculating dynamic weights.')
    parser.add_argument('--use_clustering', action='store_true', help='If set, uses clustering for feature selection.')
    parser.add_argument('--cluster_method', type=str, default='hdbscan', choices=['hdbscan', 'hierarchical'], help='Clustering method to use.')
    parser.add_argument('--n_clusters', type=int, default=20, help='Number of clusters for clustering.')
    parser.add_argument('--min_cluster_size', type=int, default=5, help='Minimum cluster size for clustering')
    parser.add_argument('--min_samples', help='Minimum number of samples for filtering.')
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help='Epsilon value for clustering.')
    parser.add_argument('--check_feature_presence', action='store_true', help='If set, checks for presence of features during train-test split.')
    parser.add_argument('--max_ram', type=int, default=8, help='Maximum amount of RAM to use in GB.')

    args = parser.parse_args()

    # Run the full feature selection workflow
    run_feature_selection_workflow(
        input_path=args.input,
        base_output_dir=args.output,
        threads=args.threads,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs=args.num_runs,
        method=args.method,
        task_type=args.task_type,
        phenotype_column=args.phenotype_column,
        sample_column=args.sample_column,
        phage_column=args.phage_column,
        binary_data=args.binary_data,
        max_features=args.max_features,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        max_ram=args.max_ram
    )

if __name__ == "__main__":
    main()
