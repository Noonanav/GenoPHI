import os
import argparse
import pandas as pd
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments
from phage_modeling.workflows.feature_annotations_workflow import run_predictive_proteins_workflow

def run_modeling_workflow_from_feature_table(
    full_feature_table, output_dir, threads=4, num_features=100, filter_type='none',
    num_runs_fs=10, num_runs_modeling=10, sample_column='strain', phenotype_column='interaction',
    method='rfe', annotation_table_path=None, protein_id_col="protein_ID",
    feature2cluster_path=None, cluster2protein_path=None, fasta_dir_or_file=None,
    run_predictive_proteins=False, phage_feature2cluster_path=None, phage_cluster2protein_path=None,
    phage_fasta_dir_or_file=None, task_type='classification', binary_data=False, max_features='none'
):
    """
    Workflow for feature selection, modeling, and predictive protein extraction starting from a pre-generated full feature table.

    Args:
        full_feature_table (str): Path to the full feature table.
        output_dir (str): Directory to save results.
        threads (int): Number of threads to use.
        num_features (int): Number of features to select.
        filter_type (str): Filter type for the input data (default: 'none').
        num_runs_fs (int): Number of feature selection iterations.
        num_runs_modeling (int): Number of runs per feature table for modeling.
        sample_column (str): Column name for the sample identifier.
        phenotype_column (str): Column name for the phenotype.
        method (str): Feature selection method.
        annotation_table_path (str, optional): Path to an optional annotation table for merging predictive protein annotations.
        protein_id_col (str): Column name for protein IDs in the predictive_proteins DataFrame.
        feature2cluster_path (str, optional): Path to the feature-to-cluster mapping file for strains.
        cluster2protein_path (str, optional): Path to the cluster-to-protein mapping file for strains.
        fasta_dir_or_file (str, optional): Path to a FASTA file or directory of FASTA files for strains.
        run_predictive_proteins (bool): Whether to run the predictive proteins extraction step.
        phage_feature2cluster_path (str, optional): Path to the feature-to-cluster mapping file for phages.
        phage_cluster2protein_path (str, optional): Path to the cluster-to-protein mapping file for phages.
        phage_fasta_dir_or_file (str, optional): Path to a FASTA file or directory for phages.
        task_type (str): Either 'classification' or 'regression' (default: classification).
        binary_data (bool): If True, converts feature values to binary (0/1).
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
        phenotype_column=phenotype_column,
        task_type=task_type
    )

    # Step 2: Generate feature tables from feature selection results
    print("Step 2: Generating feature tables from feature selection results...")
    max_features = None if max_features == 'none' else int(max_features)
    filter_table_dir = os.path.join(base_fs_output_dir, 'filtered_feature_tables')
    generate_feature_tables(
        model_testing_dir=base_fs_output_dir,
        full_feature_table_file=full_feature_table,
        filter_table_dir=filter_table_dir,
        phenotype_column=phenotype_column,
        sample_column=sample_column,
        cut_offs=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50],
        binary_data=binary_data,
        max_features=max_features,
        filter_type=filter_type
    )

    # Step 3: Modeling
    print("Step 3: Running modeling experiments...")
    modeling_output_dir = os.path.join(output_dir, 'modeling_results')
    run_experiments(
        input_dir=filter_table_dir,
        base_output_dir=modeling_output_dir,
        threads=threads,
        num_runs=num_runs_modeling,
        set_filter=filter_type,
        sample_column=sample_column,
        phenotype_column=phenotype_column,
        task_type=task_type,
        binary_data=binary_data
    )

    # Conditional Step 4: Predictive Proteins Workflow
    if run_predictive_proteins:
        print("Step 4: Selecting top-performing cutoff and running predictive proteins workflow...")
        
        metrics_file = os.path.join(modeling_output_dir, 'model_performance', 'model_performance_metrics.csv')
        performance_df = pd.read_csv(metrics_file)
        
        # Choose the top cutoff based on metric (MCC for classification, R2 for regression)
        if task_type == 'classification':
            top_cutoff = performance_df.loc[performance_df['MCC'].idxmax(), 'cut_off'].split('_')[-1]
        elif task_type == 'regression':
            top_cutoff = performance_df.loc[performance_df['r2'].idxmax(), 'cut_off'].split('_')[-1]

        # Define paths based on selected top cutoff
        feature_file_path = os.path.join(base_fs_output_dir, 'filtered_feature_tables', f'select_feature_table_cutoff_{top_cutoff}.csv')
        predictive_proteins_output_dir = os.path.join(modeling_output_dir, 'model_performance', 'predictive_proteins')

        # Run predictive proteins workflow for strain
        if feature2cluster_path and cluster2protein_path:
            run_predictive_proteins_workflow(
                feature_file_path=feature_file_path,
                feature2cluster_path=feature2cluster_path,
                cluster2protein_path=cluster2protein_path,
                fasta_dir_or_file=fasta_dir_or_file,
                modeling_dir=os.path.join(modeling_output_dir, f'cutoff_{top_cutoff}'),
                output_dir=predictive_proteins_output_dir,
                output_fasta='predictive_AA_seqs_strain.faa',
                protein_id_col=protein_id_col,
                annotation_table_path=annotation_table_path,
                feature_assignments_path=os.path.join(output_dir, 'strain', 'features', 'feature_assignments.csv'),
                strain_column=sample_column
            )

        # Run predictive proteins workflow for phage, if phage parameters are provided
        if phage_feature2cluster_path and phage_cluster2protein_path:
            run_predictive_proteins_workflow(
                feature_file_path=feature_file_path,
                feature2cluster_path=phage_feature2cluster_path,
                cluster2protein_path=phage_cluster2protein_path,
                fasta_dir_or_file=phage_fasta_dir_or_file,
                modeling_dir=os.path.join(modeling_output_dir, f'cutoff_{top_cutoff}'),
                output_dir=predictive_proteins_output_dir,
                output_fasta='predictive_AA_seqs_phage.faa',
                protein_id_col=protein_id_col,
                annotation_table_path=annotation_table_path,
                feature_assignments_path=os.path.join(output_dir, 'phage', 'features', 'feature_assignments.csv'),
                strain_column='phage'
            )
    else:
        print("Step 4: Predictive proteins workflow skipped.")

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Workflow for feature selection, modeling, and optional predictive protein extraction from a pre-existing feature table.')

    # Input arguments
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-i', '--full_feature_table', type=str, required=True, help='Path to the full feature table.')

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save results.')

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--num_features', type=int, default=100, help='Number of features to select (default: 100).')
    fs_modeling_group.add_argument('--filter_type', type=str, default='none', help="Filter column for the input data.")
    fs_modeling_group.add_argument('--method', type=str, default='rfe', choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                                   help="Feature selection method ('rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'; default: rfe).")
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection iterations to run (default: 10).')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=10, help='Number of runs per feature table for modeling (default: 10).')
    fs_modeling_group.add_argument('--task_type', type=str, default='classification', choices=['classification', 'regression'], help="Specify 'classification' or 'regression' task (default: classification).")
    fs_modeling_group.add_argument('--binary_data', action='store_true', help='If set, converts feature values to binary (1/0); otherwise, continuous values are kept.')
    fs_modeling_group.add_argument('--max_features', default='none', help='Maximum number of features to include in the feature tables.')

    # Predictive proteins and annotations
    predictive_proteins_group = parser.add_argument_group('Predictive Proteins and Annotations')
    predictive_proteins_group.add_argument('--annotation_table_path', type=str, help="Path to an optional annotation table for merging predictive protein annotations (CSV/TSV).")
    predictive_proteins_group.add_argument('--protein_id_col', type=str, default="protein_ID", help="Column name for protein IDs in the predictive_proteins DataFrame.")
    predictive_proteins_group.add_argument('--feature2cluster_path', type=str, help="Path to the strain feature-to-cluster mapping file.")
    predictive_proteins_group.add_argument('--cluster2protein_path', type=str, help="Path to the strain cluster-to-protein mapping file.")
    predictive_proteins_group.add_argument('--fasta_dir_or_file', type=str, help="Path to the strain FASTA file or directory.")
    predictive_proteins_group.add_argument('--phage_feature2cluster_path', type=str, help="Path to the phage feature-to-cluster mapping file.")
    predictive_proteins_group.add_argument('--phage_cluster2protein_path', type=str, help="Path to the phage cluster-to-protein mapping file.")
    predictive_proteins_group.add_argument('--phage_fasta_dir_or_file', type=str, help="Path to the phage FASTA file or directory.")
    predictive_proteins_group.add_argument('--run_predictive_proteins', action='store_true', help="Include to run predictive proteins extraction workflow.")

    # Optional column parameters
    optional_columns_group = parser.add_argument_group('Optional columns')
    optional_columns_group.add_argument('--sample_column', type=str, default='strain', help='Column name for the sample identifier (default: strain).')
    optional_columns_group.add_argument('--phenotype_column', type=str, default='interaction', help='Column name for the phenotype (optional).')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')

    args = parser.parse_args()

    # Run the feature selection, modeling, and optional predictive proteins extraction workflow
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
        method=args.method,
        annotation_table_path=args.annotation_table_path,
        protein_id_col=args.protein_id_col,
        feature2cluster_path=args.feature2cluster_path,
        cluster2protein_path=args.cluster2protein_path,
        fasta_dir_or_file=args.fasta_dir_or_file,
        run_predictive_proteins=args.run_predictive_proteins,
        phage_feature2cluster_path=args.phage_feature2cluster_path,
        phage_cluster2protein_path=args.phage_cluster2protein_path,
        phage_fasta_dir_or_file=args.phage_fasta_dir_or_file,
        task_type=args.task_type,
        binary_data=args.binary_data,
        max_features=args.max_features
    )

if __name__ == "__main__":
    main()
