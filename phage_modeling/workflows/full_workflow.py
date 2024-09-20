import os
import argparse
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments

def run_full_workflow(input_path_strain, input_path_phage, interaction_matrix, output_dir, tmp_dir="tmp", min_seq_id=0.6, coverage=0.8, sensitivity=7.5, suffix='faa', threads=4, strain_list='none', strain_column='strain', compare=False, source_strain='strain', source_phage='phage', num_features=100, filter_type='none', num_runs_fs=10, num_runs_modeling=10, sample_column=None, phenotype_column=None):
    """
    Complete workflow: Feature table generation, feature selection, and modeling.

    Args:
        input_path_strain (str): Path to the input directory or file for strain clustering.
        input_path_phage (str): Path to the input directory or file for phage clustering.
        interaction_matrix (str): Path to the interaction matrix.
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.
        min_seq_id (float): Minimum sequence identity for clustering.
        coverage (float): Minimum coverage for clustering.
        sensitivity (float): Sensitivity for clustering.
        suffix (str): Suffix for input FASTA files.
        threads (int): Number of threads to use.
        strain_list (str or None): Path to a strain list file, or None for no filtering.
        strain_column (str): Column in the strain list file containing strain names.
        compare (bool): Whether to compare original clusters with assigned clusters.
        source_strain (str): Prefix for naming selected features for strain in the assignment step.
        source_phage (str): Prefix for naming selected features for phage in the assignment step.
        num_features (int): Number of features to select during RFE.
        filter_type (str): Filter type for the input data ('strain', 'phage', 'none').
        num_runs_fs (int): Number of feature selection iterations.
        num_runs_modeling (int): Number of runs per feature table for modeling.
        sample_column (str): Column name for the sample identifier.
        phenotype_column (str): Column name for the phenotype.
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
    run_feature_assignment(os.path.join(strain_output_dir, "presence_absence_matrix.csv"), os.path.join(strain_output_dir, "features"), source=source_strain, select=strain_list, select_column='strain')

    # Run clustering and feature assignment for phage
    run_clustering_workflow(input_path_phage, phage_output_dir, phage_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list, strain_column, compare)
    run_feature_assignment(os.path.join(phage_output_dir, "presence_absence_matrix.csv"), os.path.join(phage_output_dir, "features"), source=source_phage, select=strain_list, select_column='phage')

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
        num_runs=num_runs_fs
    )
    
    # Generate feature tables from feature selection
    filter_table_dir = os.path.join(base_fs_output_dir, 'filtered_feature_tables')
    generate_feature_tables(
        model_testing_dir=base_fs_output_dir,
        full_feature_table_file=merged_feature_table,
        filter_table_dir=filter_table_dir,
        cut_offs=[3, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    )

    # Step 3: Modeling
    print("Step 3: Running modeling experiments...")
    run_experiments(
        input_dir=filter_table_dir,
        base_output_dir=os.path.join(output_dir, 'modeling_results'),
        threads=threads,
        num_runs=num_runs_modeling,
        set_filter=filter_type,  # Replaced set_filter with filter_type
        sample_column=sample_column,
        phenotype_column=phenotype_column
    )

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run the full workflow: feature table generation, feature selection, and modeling.')
    parser.add_argument('-ih', '--input_strain', type=str, required=True, help='Input path for strain clustering (directory or file).')
    parser.add_argument('-ip', '--input_phage', type=str, required=True, help='Input path for phage clustering (directory or file).')
    parser.add_argument('-im', '--interaction_matrix', type=str, required=True, help='Path to the interaction matrix.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory to save results.')
    parser.add_argument('--tmp', type=str, default="tmp", help='Temporary directory for intermediate files (default: tmp).')
    parser.add_argument('--min_seq_id', type=float, default=0.6, help='Minimum sequence identity for clustering (default: 0.6).')
    parser.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering (default: 0.8).')
    parser.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering (default: 7.5).')
    parser.add_argument('--suffix', type=str, default='faa', help='Suffix for input FASTA files (default: faa).')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')
    parser.add_argument('--strain_list', type=str, default='none', help='Path to a strain list file for filtering (default: none).')
    parser.add_argument('--strain_column', type=str, default='strain', help='Column in the strain list containing strain names (default: strain).')
    parser.add_argument('--compare', action='store_true', help='Compare original clusters with assigned clusters.')
    parser.add_argument('--source_strain', type=str, default='strain', help='Prefix for naming selected features for strain in the assignment step (default: strain).')
    parser.add_argument('--source_phage', type=str, default='phage', help='Prefix for naming selected features for phage in the assignment step (default: phage).')
    
    parser.add_argument('--num_features', type=int, default=100, help='Number of features to select during feature selection (default: 100).')
    parser.add_argument('--filter_type', type=str, default='none', help="Type of filtering to use during feature selection and modeling ('none', 'strain', 'phage', 'dataset'; default: none).")
    parser.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection iterations to run (default: 10).')
    
    parser.add_argument('--num_runs_modeling', type=int, default=10, help='Number of runs per feature table for modeling (default: 10).')
    parser.add_argument('--sample_column', type=str, help='Column name for the sample identifier (optional).')
    parser.add_argument('--phenotype_column', type=str, help='Column name for the phenotype (optional).')

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
        strain_column=args.strain_column,
        compare=args.compare,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        num_features=args.num_features,
        filter_type=args.filter_type,  # Only one filter type now
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column
    )

if __name__ == "__main__":
    main()
