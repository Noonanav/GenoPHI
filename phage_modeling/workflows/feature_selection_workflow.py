import os
import argparse
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables

def run_feature_selection_workflow(input_path, base_output_dir, threads=4, num_features=500, filter_type='none', num_runs=50, method='rfe'):
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
    """
    # Run multiple iterations of feature selection
    print("Running feature selection iterations...")
    run_feature_selection_iterations(
        input_path=input_path,
        base_output_dir=base_output_dir,
        threads=threads,
        num_features=num_features,
        filter_type=filter_type,
        num_runs=num_runs,
        method=method
    )
    
    # Generate feature tables based on the results
    print("Generating feature tables...")
    filter_table_dir = os.path.join(base_output_dir, 'filtered_feature_tables')
    generate_feature_tables(
        model_testing_dir=base_output_dir,
        full_feature_table_file=input_path,
        filter_table_dir=filter_table_dir,
        cut_offs=[3, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    )

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run feature selection workflow.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input path for the full feature table.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Base output directory for the results.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use.')
    parser.add_argument('--num_features', type=int, default=500, help='Number of features to select during feature selection.')
    parser.add_argument('--filter_type', type=str, default='none', help="Type of filtering to use ('none', 'strain', 'phage').")
    parser.add_argument('--num_runs', type=int, default=50, help='Number of feature selection iterations to run.')
    parser.add_argument('--method', type=str, default='rfe', choices=['rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                        help="Feature selection method to use ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').")
    
    args = parser.parse_args()

    # Run the full feature selection workflow
    run_feature_selection_workflow(
        input_path=args.input,
        base_output_dir=args.output,
        threads=args.threads,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs=args.num_runs,
        method=args.method
    )

if __name__ == "__main__":
    main()
