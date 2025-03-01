import os
import logging
from argparse import ArgumentParser
from phage_modeling.workflows.kmer_assign_features_workflow import run_kmer_assign_features_workflow
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

def kmer_assign_predict_workflow(input_dir, mmseqs_db, clusters_tsv, feature_map, filtered_kmers, aa_sequence_file,
                            tmp_dir, output_dir, model_dir, feature_table=None, phage_feature_table_path=None, genome_type='strain', 
                            genome_list=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4, suffix='faa', 
                            threshold=0.5, reuse_existing=True):
    """
    Runs the full assignment and prediction workflow in sequence by calling functions from specific modules.
    Can reuse existing output files to avoid unnecessary computation.

    Args:
        input_dir (str): Directory containing genome FASTA files for assignment.
        mmseqs_db (str): Path to the existing MMseqs2 database.
        clusters_tsv (str): Path to the clusters TSV file.
        feature_map (str): Path to the feature mapping CSV file.
        filtered_kmers (str): Path to the filtered kmers CSV file.
        aa_sequence_file (str): Path to the FASTA file containing amino acid sequences.
        tmp_dir (str): Temporary directory for intermediate files.
        output_dir (str): Directory to save results.
        feature_table (str, optional): Path to the combined feature table for prediction.
        model_dir (str): Directory containing trained models for prediction.
        phage_feature_table_path (str, optional): Path to the phage feature table for prediction.
        genome_type (str): Type of genomes being processed ('strain' or 'phage').
        genome_list (str, optional): Path to a file with list of genomes to process.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
        suffix (str): Suffix for FASTA files.
        threshold (float): Threshold for kmer matching percentage.
        reuse_existing (bool): Whether to reuse existing output files.
    """
    logging.info("Starting assignment workflow with kmers...")

    # Directory for assignment results
    assign_output_dir = os.path.join(output_dir, "assign_results")
    os.makedirs(assign_output_dir, exist_ok=True)

    # Path to the feature table that will be generated
    strain_feature_table_path = os.path.join(assign_output_dir, f'{genome_type}_combined_feature_table.csv')
    
    # Check if we can skip assignment
    assignment_needed = not (reuse_existing and os.path.exists(strain_feature_table_path))
    
    if assignment_needed:
        # Run assignment workflow with kmers
        run_kmer_assign_features_workflow(
            input_dir=input_dir,
            mmseqs_db=mmseqs_db,
            tmp_dir=tmp_dir,
            output_dir=assign_output_dir,
            feature_map=feature_map,
            filtered_kmers=filtered_kmers,
            aa_sequence_file=aa_sequence_file,
            clusters_tsv=clusters_tsv,
            genome_type=genome_type,
            genome_list=genome_list,
            sensitivity=sensitivity,
            coverage=coverage,
            min_seq_id=min_seq_id,
            threads=threads,
            suffix=suffix,
            threshold=threshold,
            reuse_existing=reuse_existing
        )
    else:
        logging.info(f"Found existing feature table: {strain_feature_table_path}. Skipping assignment.")

    # Skip prediction if the feature table still doesn't exist
    if not os.path.exists(strain_feature_table_path):
        logging.error(f"Feature table not found: {strain_feature_table_path}. Cannot proceed to prediction workflow.")
        return

    # Directory for prediction results
    predict_output_dir = os.path.join(output_dir, "predict_results")
    os.makedirs(predict_output_dir, exist_ok=True)
    
    # Path to the prediction results
    prediction_results_path = os.path.join(predict_output_dir, f"{genome_type}_median_predictions.csv")
    
    # Check if we can skip prediction
    if reuse_existing and os.path.exists(prediction_results_path):
        logging.info(f"Found existing prediction results: {prediction_results_path}. Skipping prediction.")
    else:
        logging.info("Starting prediction workflow...")
        # Run prediction workflow
        run_prediction_workflow(
            input_dir=assign_output_dir,  # The output of assign is now the input for predict
            phage_feature_table_path=phage_feature_table_path,
            model_dir=model_dir,
            feature_table=feature_table,
            output_dir=predict_output_dir
        )

    logging.info("Combined assignment and prediction workflow completed successfully.")

def main():
    parser = ArgumentParser(description="Run both assignment and prediction workflows sequentially with kmer-based assignment.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the existing MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map (selected_features.csv).")
    parser.add_argument('--filtered_kmers', type=str, required=True, help="Path to the filtered kmers CSV file.")
    parser.add_argument('--aa_sequence_file', type=str, required=True, help="Path to the FASTA file containing amino acid sequences.")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory with models for prediction.")
    parser.add_argument('--feature_table', type=str, help="Path to the combined feature table. Optional for single-strain mode.")
    parser.add_argument('--phage_feature_table', type=str, help="Path to the phage feature table. Optional for single-strain mode.")
    parser.add_argument('--genome_type', type=str, choices=['strain', 'phage'], default='strain', help="Type of genome to process.")
    parser.add_argument('--genome_list', type=str, help="Path to file with list of genomes to process.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.6, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")
    parser.add_argument('--suffix', type=str, default='faa', help="Suffix for FASTA files.")
    parser.add_argument('--threshold', type=float, default=0.5, help="Threshold for kmer matching percentage.")
    parser.add_argument('--reuse_existing', action='store_true', help="Reuse existing output files if available.")

    args = parser.parse_args()

    kmer_assign_predict_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        feature_map=args.feature_map,
        filtered_kmers=args.filtered_kmers,
        aa_sequence_file=args.aa_sequence_file,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        model_dir=args.model_dir,
        feature_table=args.feature_table,
        phage_feature_table_path=args.phage_feature_table,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        threshold=args.threshold,
        reuse_existing=args.reuse_existing
    )

if __name__ == "__main__":
    main()
