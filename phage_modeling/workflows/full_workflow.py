import os
import argparse
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments

def run_full_workflow(input_path_strain, input_path_phage, interaction_matrix, output_dir, tmp_dir="tmp", min_seq_id=0.6, coverage=0.8, sensitivity=7.5, suffix='faa', threads=4, strain_list='none', phage_list='none', strain_column='strain', phage_column='phage', compare=False, source_strain='strain', source_phage='phage', num_features=100, filter_type='none', num_runs_fs=10, num_runs_modeling=10, sample_column='strain', phenotype_column=None, method='rfe'):
    """
    Complete workflow: Feature table generation, feature selection, and modeling.

    Args:
    Input data:
        input_path_strain (str): Path to the input directory or file for strain clustering.
        input_path_phage (str): Path to the input directory or file for phage clustering.
        interaction_matrix (str): Path to the interaction matrix.

    Optional input arguments:
        suffix (str): Suffix for input FASTA files.
        strain_list (str or None): Path to a strain list file, or None for no filtering.
        phage_list (str or None): Path to a phage list file, or None for no filtering.
        strain_column (str): Column in the strain list file containing strain names.
        phage_column (str): Column in the phage list file containing phage names.
        source_strain (str): Prefix for naming selected features for strain in the assignment step.
        source_phage (str): Prefix for naming selected features for phage in the assignment step.
        sample_column (str): Column name for the sample identifier.
        phenotype_column (str): Column name for the phenotype.

    Output arguments:
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.

    Clustering:
        min_seq_id (float): Minimum sequence identity for clustering.
        coverage (float): Minimum coverage for clustering.
        sensitivity (float): Sensitivity for clustering.
        compare (bool): Whether to compare original clusters with assigned clusters.

    Feature selection and modeling:
        filter_type (str): Filter type for the input data ('strain', 'phage', 'none').
        method (str): Feature selection method ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').
        num_features (int): Number of features to select.
        num_runs_fs (int): Number of feature selection iterations.
        num_runs_modeling (int): Number of runs per feature table for modeling.
    
    General:
        threads (int): Number of threads to use.
    """

    # Step 1: Feature table generation for strain and phage
    print("Step 1: Running feature table generation...")
    strain_output_dir = os.path.join(output_dir, "strain")
    phage_output_dir = os.path.join(output_dir, "phage")
    merged_output_dir = os.path.join(output_dir, "merged")

    strain_tmp_dir = os.path.join(tmp_dir, "strain")
    phage_tmp_dir = os.path.join(tmp_dir, "phage")

    # Ensure output directories exist
    if not os.path.exists(merged_output_dir):
        os.makedirs(merged_output_dir)

    # Run clustering and feature assignment for strain
    run_clustering_workflow(input_path_strain, strain_output_dir, strain_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list, strain_column, compare)
    run_feature_assignment(os.path.join(strain_output_dir, "presence_absence_matrix.csv"), os.path.join(strain_output_dir, "features"), source=source_strain, select=strain_list, select_column=strain_column)

    # Run clustering and feature assignment for phage
    run_clustering_workflow(input_path_phage, phage_output_dir, phage_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, phage_list, phage_column, compare)
    run_feature_assignment(os.path.join(phage_output_dir, "presence_absence_matrix.csv"), os.path.join(phage_output_dir, "features"), source=source_phage, select=phage_list, select_column=phage_column)

    # Merge strain and phage feature tables
    merged_feature_table = merge_feature_tables(
        os.path.join(strain_output_dir, "features", "feature_table.csv"),
        os.path.join(phage_output_dir, "features", "feature_table.csv"),
        interaction_matrix,
        merged_output_dir,
        remove_suffix=False
    )
    
    print(f"Merged feature table saved in: {merged_output_dir}")

    # Step 2: Feature selection
    print("Step 2: Running feature selection iterations...")
    base_fs_output_dir = os.path.join(output_dir, 'feature_selection')
    run_feature_selection_iterations(
        input_path=merged_feature_table,
        base_output_dir=base_fs_output_dir,
        threads=threads,
        num_features=num_features,
        filter_type=filter_type,
        num_runs=num_runs_fs,
        method=method,
        sample_column=sample_column,
        phenotype_column=phenotype_column
    )
    
    # Generate feature tables from feature selection
    print("Generating feature tables from feature selection results...")
    filter_table_dir = os.path.join(base_fs_output_dir, 'filtered_feature_tables')
    generate_feature_tables(
        model_testing_dir=base_fs_output_dir,
        full_feature_table_file=merged_feature_table,
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
    parser = argparse.ArgumentParser(description='Complete workflow: Feature table generation, feature selection, and modeling.')
    
    # Input data
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-ih', '--input_strain', type=str, required=True, help='Path to the input directory or file for strain clustering.')
    input_group.add_argument('-ip', '--input_phage', type=str, required=True, help='Path to the input directory or file for phage clustering.')
    input_group.add_argument('-im', '--interaction_matrix', type=str, required=True, help='Path to the interaction matrix.')

    # Optional input arguments
    optional_input_group = parser.add_argument_group('Optional input arguments')
    optional_input_group.add_argument('--suffix', type=str, default='faa', help='Suffix for input FASTA files (default: faa).')
    optional_input_group.add_argument('--strain_list', type=str, default='none', help='Path to a strain list file for filtering (default: none).')
    optional_input_group.add_argument('--phage_list', type=str, default='none', help='Path to a phage list file for filtering (default: none).')
    optional_input_group.add_argument('--strain_column', type=str, default='strain', help='Column in the strain list containing strain names (default: strain).')
    optional_input_group.add_argument('--phage_column', type=str, default='phage', help='Column in the phage list containing phage names (default: phage).')
    optional_input_group.add_argument('--source_strain', type=str, default='strain', help='Prefix for naming selected features for strain in the assignment step (default: strain).')
    optional_input_group.add_argument('--source_phage', type=str, default='phage', help='Prefix for naming selected features for phage in the assignment step (default: phage).')
    optional_input_group.add_argument('--sample_column', type=str, default='strain', help='Column name for the sample identifier (default: strain).')
    optional_input_group.add_argument('--phenotype_column', type=str, help='Column name for the phenotype (optional).')

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('-o', '--output', type=str, required=True, help='Output directory to save results.')
    output_group.add_argument('--tmp', type=str, default="tmp", help='Temporary directory for intermediate files (default: tmp).')

    # Clustering parameters
    clustering_group = parser.add_argument_group('Clustering')
    clustering_group.add_argument('--min_seq_id', type=float, default=0.6, help='Minimum sequence identity for clustering (default: 0.6).')
    clustering_group.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering (default: 0.8).')
    clustering_group.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering (default: 7.5).')
    clustering_group.add_argument('--compare', action='store_true', help='Compare original clusters with assigned clusters.')

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--filter_type', type=str, default='none', help="Filter type for the input data ('none', 'strain', 'phage').")
    fs_modeling_group.add_argument('--method', type=str, default='rfe', choices=['rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                                   help="Feature selection method ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'; default: rfe).")
    fs_modeling_group.add_argument('--num_features', type=int, default=100, help='Number of features to select (default: 100).')
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection iterations to run (default: 10).')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=10, help='Number of runs per feature table for modeling (default: 10).')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')
    
    args = parser.parse_args()

    # Run the full workflow
    run_full_workflow(
        input_path_strain=args.input_strain,
        input_path_phage=args.input_phage,
        interaction_matrix=args.interaction_matrix,
        output_dir=args.output,
        tmp_dir=args.tmp,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        compare=args.compare,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
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
