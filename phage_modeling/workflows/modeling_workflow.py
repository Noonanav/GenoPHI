import os
import argparse
from phage_modeling.select_feature_modeling import run_experiments

def run_modeling_workflow(
    input_dir, 
    base_output_dir, 
    threads=4, 
    num_runs=100, 
    set_filter='none', 
    sample_column=None, 
    phenotype_column='interaction', 
    phage_column='phage',
    task_type='classification', 
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=True,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    binary_data=False
):
    """
    Workflow to run experiments on selected feature tables using grid search and MCC/R2 optimization.
    
    Args:
        input_dir (str): Directory containing the feature tables to model.
        base_output_dir (str): Directory to save the results of each experiment.
        threads (int): Number of threads to use.
        num_runs (int): Number of runs to perform per feature table.
        set_filter (str): Filter type for dataset ('none', 'strain', 'phage', 'dataset').
        sample_column (str): Column for sample identification.
        phenotype_column (str): Column for phenotype data.
        task_type (str): Specifies if the task is 'classification' or 'regression'.
        binary_data (bool): If True, plot SHAP jitter plot with binary data.
    """
    # Run the experiments on each feature table in the input directory
    print(f"Running {task_type} modeling experiments on feature tables...")
    run_experiments(
        input_dir=input_dir,
        base_output_dir=base_output_dir,
        threads=threads,
        num_runs=num_runs,
        set_filter=set_filter,
        sample_column=sample_column,
        phenotype_column=phenotype_column,
        phage_column=phage_column,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        task_type=task_type,
        use_clustering=use_clustering,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        binary_data=binary_data
    )

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run modeling workflow on selected feature tables.')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Directory containing selected feature tables.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save results of the experiments.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use.')
    parser.add_argument('--num_runs', type=int, default=100, help='Number of runs per feature table.')
    parser.add_argument('--set_filter', type=str, default='none', help="Filter for dataset ('none', 'strain', 'phage', 'dataset').")
    parser.add_argument('--sample_column', type=str, default='strain', help='Column name for the sample identifier (optional).')
    parser.add_argument('--phenotype_column', type=str, default='interaction', help='Column name for the phenotype (optional).')
    parser.add_argument('--phage_column', type=str, default='phage', help='Column name for the phage identifier (optional).')
    parser.add_argument('--task_type', type=str, default='classification', choices=['classification', 'regression'], help="Specify 'classification' or 'regression' task.")
    parser.add_argument('--use_dynamic_weights', action='store_true', help='If True, use dynamic weights for feature selection.')
    parser.add_argument('--weights_method', type=str, default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help='Method for calculating dynamic weights.')
    parser.add_argument('--use_clustering', action='store_true', help='If True, use clustering for feature selection')
    parser.add_argument('--min_cluster_size', type=int, default=5, help='Minimum cluster size for clustering feature selection.')
    parser.add_argument('--min_samples', type=int, help='Minimum number of samples for clustering feature selection.')
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help='Epsilon value for clustering feature selection.')
    parser.add_argument('--binary_data', action='store_true', help='If True, plot SHAP jitter plot with binary data.')

    args = parser.parse_args()

    # Run the full modeling workflow
    run_modeling_workflow(
        input_dir=args.input_dir,
        base_output_dir=args.output_dir,
        threads=args.threads,
        num_runs=args.num_runs,
        set_filter=args.set_filter,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        phage_column=args.phage_column,
        task_type=args.task_type,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        binary_data=args.binary_data
    )

if __name__ == "__main__":
    main()
