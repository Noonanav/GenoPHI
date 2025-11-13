import os
import argparse
import logging
import pandas as pd
from genophi.feature_annotations import (
    get_predictive_features,
    get_predictive_proteins, 
    parse_and_filter_aa_sequences,
    parse_feature_information, 
    merge_annotation_table,  # Optional function call
    output_predictive_feature_overview,
    merge_importance_table
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_predictive_proteins_workflow(
    feature_file_path, 
    feature2cluster_path, 
    cluster2protein_path, 
    fasta_dir_or_file, 
    modeling_dir,  
    output_dir=".",
    output_fasta="predictive_AA_seqs.faa", 
    protein_id_col="protein_ID",
    annotation_table_path=None,  # Annotation table is optional
    feature_assignments_path=None, # Feature assignments are optional
    strain_column='strain',  # Column to use for strain
    feature_type='strain',  # Default feature type
    phenotype_column='interaction'  # Default phenotype column
):
    """
    Runs the full workflow to retrieve predictive proteins, optionally merge with annotation table, filter AA sequences, 
    and output a predictive feature overview CSV.

    Args:
        feature_file_path (str): Path to the file containing predictive features.
        feature2cluster_path (str): Path to the file containing the mapping of features to clusters.
        cluster2protein_path (str): Path to the file containing the mapping of clusters to proteins.
        fasta_dir_or_file (str): Path to a FASTA file or directory containing FASTA files.
        modeling_dir (str): Path to the directory containing modeling runs for parsing feature information.
        output_dir (str): Directory to save outputs (default: '.').
        output_fasta (str): Name of the output FASTA file for filtered sequences.
        protein_id_col (str): Column name for protein IDs (default: 'protein_ID').
        annotation_table_path (str, optional): Path to an annotation table (CSV/TSV format). If provided, proteins will be merged with annotation data.
        feature_assignments_path (str, optional): Path to feature assignments CSV file for merging strain information.
        strain_column (str): Column to use for strain information (default: 'strain').
        feature_type (str): Type of features to extract ('strain' or 'phage').
        phenotype_column (str): Column name for the phenotype or target variable (default: 'interaction').

    Returns:
        None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Get predictive features
    logging.info("Step 1: Extracting predictive features.")
    select_features = get_predictive_features(feature_file_path, feature_type=feature_type, sample_column=strain_column, phenotype_column=phenotype_column)

    # Step 2: Get predictive proteins based on selected features
    logging.info("Step 2: Retrieving predictive proteins.")
    predictive_proteins = get_predictive_proteins(
        select_features=select_features,
        feature2cluster_path=feature2cluster_path,
        cluster2protein_path=cluster2protein_path
    )

    # Step 3: Parse feature importance information
    logging.info("Step 3: Parsing feature importance information.")
    feature_importance_df = parse_feature_information(modeling_dir, output_dir)

    # Step 4: Merge predictive proteins with feature importance information
    filtered_proteins = merge_importance_table(
        predictive_proteins=predictive_proteins,
        feature_importance_df=feature_importance_df,
        output_dir=output_dir,
        protein_id_col=protein_id_col,
        prefix=strain_column
    )

    # Optional Step 5: Merge predictive proteins with annotation table if provided
    if annotation_table_path:
        logging.info("Step 4: Merging predictive proteins with annotation table.")
        merge_annotation_table(
            annotation_table_path=annotation_table_path,
            merged_df=filtered_proteins,
            output_dir=output_dir,
            protein_id_col=protein_id_col,
            prefix=strain_column
        )
    else:
        logging.info("Skipping annotation merge as no annotation table was provided.")

    # Step 6: Parse AA sequences from the FASTA file(s)
    logging.info("Step 5: Parsing and filtering AA sequences.")
    genome_protein_df = parse_and_filter_aa_sequences(
        fasta_dir_or_file=fasta_dir_or_file,
        filtered_proteins=filtered_proteins,
        protein_id_col=protein_id_col,
        output_dir=output_dir,
        output_fasta=output_fasta,
        prefix=strain_column
    )

    logging.info("Predictive protein workflow completed.")

    # Step 6: Generate predictive feature overview (CSV)
    if feature_assignments_path:
        logging.info("Step 6: Generating predictive feature overview CSV.")
        feature_assignments_df = pd.read_csv(feature_assignments_path)
        output_predictive_feature_overview(
            predictive_proteins=predictive_proteins,
            feature_assignments_df=feature_assignments_df,
            genome_protein_df=genome_protein_df,
            strain_column=strain_column,
            output_dir=output_dir
        )

# Command-line interface for the workflow
def main():
    parser = argparse.ArgumentParser(description="Run predictive protein workflow.")
    parser.add_argument('--feature_file_path', required=True, help="Path to the file containing predictive features.")
    parser.add_argument('--feature2cluster_path', required=True, help="Path to the file containing feature-to-cluster mappings.")
    parser.add_argument('--cluster2protein_path', required=True, help="Path to the file containing cluster-to-protein mappings.")
    parser.add_argument('--fasta_dir_or_file', required=True, help="Path to the FASTA file or directory containing FASTA files.")
    parser.add_argument('--modeling_dir', required=True, help="Path to the directory containing modeling runs for parsing feature importance.")  
    parser.add_argument('--output_dir', default=".", help="Output directory for output tables and filtered sequences.")
    parser.add_argument('--output_fasta', default="predictive_AA_seqs.faa", help="Name of the output FASTA file for filtered sequences.")
    parser.add_argument('--protein_id_col', default="protein_ID", help="Column name for protein IDs in the predictive_proteins DataFrame.")
    parser.add_argument('--annotation_table_path', help="Path to an optional annotation table (CSV/TSV).")
    parser.add_argument('--feature_assignments_path', help="Path to feature assignments CSV file for merging strain information.")
    parser.add_argument('--strain_column', default='strain', help="Column to use for strain information (default: 'strain').")
    parser.add_argument('--feature_type', default='strain', help="Type of features to extract ('strain' or 'phage').")
    parser.add_argument('--phenotype_column', default='interaction', help="Column name for the phenotype or target variable.")

    args = parser.parse_args()

    # Run the workflow with parsed arguments
    run_predictive_proteins_workflow(
        feature_file_path=args.feature_file_path,
        feature2cluster_path=args.feature2cluster_path,
        cluster2protein_path=args.cluster2protein_path,
        fasta_dir_or_file=args.fasta_dir_or_file,
        modeling_dir=args.modeling_dir,
        output_dir=args.output_dir,
        output_fasta=args.output_fasta,
        protein_id_col=args.protein_id_col,
        annotation_table_path=args.annotation_table_path,
        feature_assignments_path=args.feature_assignments_path,  # Optional
        strain_column=args.strain_column,
        feature_type=args.feature_type,
        phenotype_column=args.phenotype_column
    )

if __name__ == "__main__":
    main()
