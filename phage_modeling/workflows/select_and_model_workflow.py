import os
import argparse
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments

def run_modeling_workflow_from_feature_table(full_feature_table, output_dir, threads=4, num_features=100, filter_type='none', num_runs_fs=10, num_runs_modeling=10, sample_column='strain', phenotype_column=None, method='rfe'):
    """
    Workflow for feature selection and modeling starting from a previously generated full feature table.

    Args:
        full_feature_table (str): Path to the full feature table.
        output_dir (str): Directory to save results.
        threads (int): Number of threads to use.
        num_features (int): Number of features to select.
        filter_type (str): Filter type for the input data ('strain', 'phage', 'none').
        num_runs_fs (int): Number of feature selection iterations.
        num_runs_modeling (int): Number of runs per feature table for modeling.
        sample_column (str): Column name for the sample identifier.
        phenotype_column (str): Column name for the phenotype.
        method (str): Feature selection method ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').
    """
    
    # Step 1: Feature Selection
    print("Step 1: Running feature selection iterations...")
    base_fs_output_dir = os.path.join(output_dir, 'feature_selection')
    run_feature_selection_iterations(
        input_path=full_feature_table,
        base_output_dir=base_fs_output_dir,
        threads=threads,
        num_features=num_features,
        filter_type=filter_type,
        num_runs=num_runs_fs,
        method=method,
        sample_column=sample_column,
        phenotype_column=phenotype_column
    )
    
    # Step 2: Generate feature tables from feature selection results
    print("Step 2: Generating feature tables from feature selection results...")
    filter_table_dir = os.path.join(base_fs_output_dir, 'filtered_feature_tables')
    generate_feature_tables(
        model_testing_dir=base_fs_output_dir,
        full_feature_table_file=full_feature_table,
        filter_table_dir=filter_table_dir,
        phenotype_column=phenotype_column,
        sample_column=sample_column,
        cut_offs=[3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50]
    )

    # Step 3: Modeling
    print("Step 3: Running modeling experiments...")
    run_experiments(
        input_dir=filter_table_dir,
        base_output_dir=os.path.join(output_dir, 'modeling_results'),
        threads=threads,
        num_runs=num_runs_modeling,
        set_filter=filter_type,
        sample_column=sample_column,
        phenotype_column=phenotype_column
    )

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Workflow for feature selection and modeling with a pre-existing feature table.')

    # Input arguments
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-i', '--full_feature_table', type=str, required=True, help='Path to the full feature table.')

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save results.')

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--num_features', type=int, default=100, help='Number of features to select (default: 100).')
    fs_modeling_group.add_argument('--filter_type', type=str, default='none', help="Filter type for the input data ('none', 'strain', 'phage').")
    fs_modeling_group.add_argument('--method', type=str, default='rfe', choices=['rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                                   help="Feature selection method ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'; default: rfe).")
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection iterations to run (default: 10).')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=10, help='Number of runs per feature table for modeling (default: 10).')

    # Optional column parameters
    optional_columns_group = parser.add_argument_group('Optional columns')
    optional_columns_group.add_argument('--sample_column', type=str, default='strain', help='Column name for the sample identifier (default: strain).')
    optional_columns_group.add_argument('--phenotype_column', type=str, help='Column name for the phenotype (optional).')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')

    args = parser.parse_args()

    # Run the feature selection and modeling workflow
    run_modeling_workflow_from_feature_table(
        full_feature_table=args.full_feature_table,
        output_dir=args.output_dir,
        threads=args.threads,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        method=args.method
    )

if __name__ == "__main__":
    main()
