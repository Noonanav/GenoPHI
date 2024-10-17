import os
import argparse
from phage_modeling.select_feature_modeling import run_experiments

def run_modeling_workflow(input_dir, base_output_dir, threads=4, num_runs=100, set_filter='none', sample_column=None, phenotype_column=None):
    """
    Workflow to run experiments on selected feature tables using grid search and MCC optimization.
    
    Args:
        input_dir (str): Directory containing the feature tables to model.
        base_output_dir (str): Directory to save the results of each experiment.
        threads (int): Number of threads to use.
        num_runs (int): Number of runs to perform per feature table.
        set_filter (str): Filter type for dataset ('none', 'strain', 'phage', 'dataset').
        sample_column (str): Column for sample identification.
        phenotype_column (str): Column for phenotype data.
    """
    # Run the experiments on each feature table in the input directory
    print("Running modeling experiments on feature tables...")
    run_experiments(
        input_dir=input_dir,
        base_output_dir=base_output_dir,
        threads=threads,
        num_runs=num_runs,
        set_filter=set_filter,
        sample_column=sample_column,
        phenotype_column=phenotype_column
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

    args = parser.parse_args()

    # Run the full modeling workflow
    run_modeling_workflow(
        input_dir=args.input_dir,
        base_output_dir=args.output_dir,
        threads=args.threads,
        num_runs=args.num_runs,
        set_filter=args.set_filter,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column
    )

if __name__ == "__main__":
    main()
