import os
import argparse
import logging
from phage_modeling.feature_annotations import (
    get_predictive_features,
    get_predictive_proteins, 
    merge_annotation_table, 
    parse_and_filter_aa_sequences,
    parse_feature_information,
    load_annotation_table  # Added import for load_annotation_table
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_predictive_proteins_workflow(
    feature_file_path, 
    feature2cluster_path, 
    cluster2protein_path, 
    annotation_table_path, 
    fasta_dir_or_file, 
    modeling_dir,  
    output_dir=".",
    output_fasta="predictive_AA_seqs.faa", 
    protein_id_col="protein_ID"
):
    """
    Runs the full workflow to retrieve predictive proteins, merge with annotation table, and filter AA sequences.

    Args:
        feature_file_path (str): Path to the file containing predictive features.
        feature2cluster_path (str): Path to the file containing the mapping of features to clusters.
        cluster2protein_path (str): Path to the file containing the mapping of clusters to proteins.
        annotation_table_path (str): Path to the annotation table (CSV/TSV format).
        fasta_dir_or_file (str): Path to a FASTA file or directory containing FASTA files.
        modeling_dir (str): Path to the directory containing modeling runs for parsing feature information.
        output_fasta (str): Name of the output FASTA file for filtered sequences.
        protein_id_col (str): Column name for protein IDs (default: 'protein_ID').

    Returns:
        None
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Step 1: Get predictive features
    logging.info("Step 1: Extracting predictive features.")
    select_features = get_predictive_features(feature_file_path)

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

    # Step 5: Merge predictive proteins with annotation table
    logging.info("Step 4: Merging predictive proteins with annotation table.")
    merged_proteins = merge_annotation_table(
        annotation_table_path=annotation_table_path,  # Pass the loaded annotation DataFrame
        predictive_proteins=predictive_proteins,
        feature_importance_df=feature_importance_df,
        output_dir=output_dir,
        protein_id_col=protein_id_col
    )

    # Step 6: Parse AA sequences from the FASTA file(s)
    logging.info("Step 5: Parsing and filtering AA sequences.")
    parse_and_filter_aa_sequences(
        fasta_dir_or_file=fasta_dir_or_file,
        filtered_proteins=predictive_proteins,
        protein_id_col=protein_id_col,
        output_dir=output_dir,
        output_fasta=output_fasta
    )

    logging.info("Predictive protein workflow completed.")


# Command-line interface for the workflow
def main():
    parser = argparse.ArgumentParser(description="Run predictive protein workflow.")
    parser.add_argument('--feature_file_path', required=True, help="Path to the file containing predictive features.")
    parser.add_argument('--feature2cluster_path', required=True, help="Path to the file containing feature-to-cluster mappings.")
    parser.add_argument('--cluster2protein_path', required=True, help="Path to the file containing cluster-to-protein mappings.")
    parser.add_argument('--annotation_table_path', required=True, help="Path to the annotation table (CSV or TSV).")
    parser.add_argument('--fasta_dir_or_file', required=True, help="Path to the FASTA file or directory containing FASTA files.")
    parser.add_argument('--modeling_dir', required=True, help="Path to the directory containing modeling runs for parsing feature importance.")  
    parser.add_argument('--output_dir', default=".", help="Output directory for output tables and filtered sequences.")
    parser.add_argument('--output_fasta', default="predictive_AA_seqs.faa", help="Name of the output FASTA file for filtered sequences.")
    parser.add_argument('--protein_id_col', default="protein_ID", help="Column name for protein IDs in the predictive_proteins DataFrame.")

    args = parser.parse_args()

    # Run the workflow with parsed arguments
    run_predictive_proteins_workflow(
        feature_file_path=args.feature_file_path,
        feature2cluster_path=args.feature2cluster_path,
        cluster2protein_path=args.cluster2protein_path,
        annotation_table_path=args.annotation_table_path,
        fasta_dir_or_file=args.fasta_dir_or_file,
        modeling_dir=args.modeling_dir,  
        output_dir=args.output_dir,
        output_fasta=args.output_fasta,
        protein_id_col=args.protein_id_col
    )

if __name__ == "__main__":
    main()
