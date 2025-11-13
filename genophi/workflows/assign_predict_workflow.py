import os
import shutil
import logging
from argparse import ArgumentParser
from genophi.workflows.assign_features_workflow import run_assign_features_workflow
from genophi.workflows.prediction_workflow import run_prediction_workflow

def assign_predict_workflow(input_dir, mmseqs_db, clusters_tsv, feature_map, tmp_dir, output_dir, model_dir, 
                            feature_table=None, strain_feature_table_path=None, phage_feature_table_path=None, 
                            genome_type='strain', genome_list=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, 
                            threads=4, suffix='faa', duplicate_all=False):
    """
    Runs the full assignment and prediction workflow in sequence.
    
    For STRAIN prediction:
        - Assigns features to new strains
        - Uses assigned strain features + existing phage features for prediction
    
    For PHAGE prediction:
        - Assigns features to new phages
        - Uses existing strain features + assigned phage features for prediction

    Args:
        input_dir (str): Directory containing genome FASTA files for assignment.
        mmseqs_db (str): Path to the existing MMseqs2 database.
        clusters_tsv (str): Path to the clusters TSV file.
        feature_map (str): Path to the feature mapping CSV file.
        tmp_dir (str): Temporary directory for intermediate files.
        output_dir (str): Directory to save results.
        model_dir (str): Directory containing trained models for prediction.
        feature_table (str, optional): Path to the combined feature table for filtering prediction features.
        strain_feature_table_path (str, optional): Path to existing strain features (required for phage prediction).
        phage_feature_table_path (str, optional): Path to existing phage features (required for strain prediction).
        genome_type (str): Type of genomes being processed ('strain' or 'phage').
        genome_list (str, optional): Path to a file with list of genomes to process.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
        suffix (str): Suffix for FASTA files.
        duplicate_all (bool): Duplicate all genomes in the feature table for prediction.
    """
    logging.info(f"Starting assignment workflow for {genome_type} prediction...")

    # Validate inputs based on genome type
    if genome_type == 'phage' and strain_feature_table_path is None:
        raise ValueError("strain_feature_table_path is required for phage prediction")
    elif genome_type == 'strain' and phage_feature_table_path is None:
        raise ValueError("phage_feature_table_path is required for strain prediction")

    # Directory for assignment results
    assign_output_dir = os.path.join(output_dir, "assign_results")
    os.makedirs(assign_output_dir, exist_ok=True)

    # The output of assignment - newly assigned features
    assigned_feature_table_path = os.path.join(assign_output_dir, f'{genome_type}_combined_feature_table.csv')

    if not os.path.exists(assigned_feature_table_path):
        logging.error(f"Assigned feature table not found: {assigned_feature_table_path}. Assigning features...")

        # Run assignment workflow
        logging.info(f"Assigning features to new {genome_type}s...")
        run_assign_features_workflow(
            input_dir=input_dir,
            mmseqs_db=mmseqs_db,
            tmp_dir=tmp_dir,
            output_dir=assign_output_dir,
            feature_map=feature_map,
            clusters_tsv=clusters_tsv,
            genome_type=genome_type,
            genome_list=genome_list,
            sensitivity=sensitivity,
            coverage=coverage,
            min_seq_id=min_seq_id,
            threads=threads,
            suffix=suffix,
            duplicate_all=duplicate_all
        )

        logging.info("Assignment completed. Starting prediction workflow...")

    # Directory for prediction results
    predict_output_dir = os.path.join(output_dir, "predict_results")
    os.makedirs(predict_output_dir, exist_ok=True)

    # Set up prediction input directory with the correct feature tables
    predict_input_dir = os.path.join(output_dir, "prediction_input")
    os.makedirs(predict_input_dir, exist_ok=True)

    if genome_type == 'phage':
        # For phage prediction: use existing strain features + newly assigned phage features
        prediction_strain_features = os.path.join(predict_input_dir, "strain_feature_table.csv")
        prediction_phage_features = assigned_feature_table_path
        
        # Copy existing strain features to prediction input
        if not os.path.exists(prediction_strain_features):
            shutil.copy2(strain_feature_table_path, prediction_strain_features)
            logging.info(f"Copied existing strain features for prediction")
        
        # Run prediction workflow
        run_prediction_workflow(
            input_dir=predict_input_dir,  # Contains strain features
            phage_feature_table_path=prediction_phage_features,  # Newly assigned phage features
            model_dir=model_dir,
            output_dir=predict_output_dir,
            feature_table=feature_table,
            strain_source='strain',
            phage_source='phage',
            threads=threads
        )
        
    elif genome_type == 'strain':
        # For strain prediction: use newly assigned strain features + existing phage features
        prediction_strain_features = os.path.join(predict_input_dir, "strain_feature_table.csv")
        prediction_phage_features = phage_feature_table_path
        
        # Copy newly assigned strain features to prediction input
        if not os.path.exists(prediction_strain_features):
            shutil.copy2(assigned_feature_table_path, prediction_strain_features)
            logging.info(f"Copied newly assigned strain features for prediction")
        
        # Run prediction workflow
        run_prediction_workflow(
            input_dir=predict_input_dir,  # Contains newly assigned strain features
            phage_feature_table_path=prediction_phage_features,  # Existing phage features
            model_dir=model_dir,
            output_dir=predict_output_dir,
            feature_table=feature_table,
            strain_source='strain',
            phage_source='phage',
            threads=threads
        )

    logging.info(f"Combined assignment and prediction workflow for {genome_type}s completed successfully.")

def main():
    parser = ArgumentParser(description="Run both assignment and prediction workflows sequentially.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the existing MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map (selected_features.csv).")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory with models for prediction.")
    parser.add_argument('--feature_table', type=str, default=None, help="Path to the feature selection table for filtering.")
    
    # Updated arguments for clearer separation
    parser.add_argument('--strain_feature_table', type=str, help="Path to existing strain features (required for phage prediction).")
    parser.add_argument('--phage_feature_table', type=str, help="Path to existing phage features (required for strain prediction).")
    
    parser.add_argument('--genome_type', type=str, choices=['strain', 'phage'], default='strain', 
                       help="Type of genome to assign features to ('strain' for new strains, 'phage' for new phages).")
    parser.add_argument('--genome_list', type=str, help="Path to file with list of genomes to process.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.4, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")
    parser.add_argument('--suffix', type=str, default='faa', help="Suffix for FASTA files.")
    parser.add_argument('--duplicate_all', action='store_true', help="Duplicate all genomes in the feature table for prediction.")

    args = parser.parse_args()

    assign_predict_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        feature_map=args.feature_map,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        model_dir=args.model_dir,
        feature_table=args.feature_table,
        strain_feature_table_path=args.strain_feature_table,
        phage_feature_table_path=args.phage_feature_table,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        duplicate_all=args.duplicate_all
    )

if __name__ == "__main__":
    main()