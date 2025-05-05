import os
import pandas as pd
import argparse
import logging
import pickle
import numpy as np
from catboost import CatBoostClassifier
from phage_modeling.workflows.assign_features_workflow import run_assign_features_workflow
from collections import defaultdict

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_model(model_file):
    """
    Loads a predictive model from a file (CatBoostClassifier).
    Supports both .pkl and .cbm formats.
    """
    if model_file.endswith('.pkl'):
        with open(model_file, 'rb') as f:
            model = pickle.load(f)
    elif model_file.endswith('.cbm'):
        model = CatBoostClassifier()
        model.load_model(model_file)
    else:
        raise ValueError(f"Unsupported model file format: {model_file}")
    return model

def align_feature_names(model, target_features_testing):
    """
    Align features in the target_features_testing DataFrame with the features expected by the model.
    Raises an error if any required features are missing.
    """
    # Get expected feature names from the model
    model_feature_names = model.feature_names_
    
    if model_feature_names is None:
        logging.warning("Model doesn't have feature names. Cannot align features.")
        return target_features_testing
    
    # Check if all expected features are present
    missing_features = [f for f in model_feature_names if f not in target_features_testing.columns]
    
    if missing_features:
        error_msg = f"Missing required features: {missing_features}"
        logging.error(error_msg)
        raise ValueError(error_msg)
    
    # If we got here, all required features exist, so filter to only the needed ones
    logging.info(f"Filtering input data from {len(target_features_testing.columns)} features to {len(model_feature_names)} features required by the model")
    return target_features_testing[model_feature_names]

def generate_full_feature_table(input_table, phage_feature_table=None, strain_source='strain', phage_source='phage', output_dir=None):
    """
    Generate a full feature table for prediction, combining strain and phage features if phage_feature_table is provided.
    """
    strain_features = ['strain'] + [col for col in input_table.columns if col not in ['strain', 'phage']]

    if phage_feature_table is None:
        logging.info("No phage feature table provided. Using strain features only.")
        return input_table[strain_features]

    # Include all columns from phage_feature_table
    phage_features_table = phage_feature_table

    # Repeat the phage table for each strain
    repeated_phage = pd.DataFrame(
        np.repeat(phage_features_table.values, len(input_table), axis=0),
        columns=phage_features_table.columns
    )

    # Repeat the strain table for each phage
    repeated_strain = pd.DataFrame(
        np.tile(input_table[strain_features].values, (len(phage_features_table), 1)),
        columns=strain_features
    )

    # Combine the repeated data
    prediction_feature_table = pd.concat([repeated_phage.reset_index(drop=True), repeated_strain.reset_index(drop=True)], axis=1)

    # Save the merged table if an output_dir is provided
    if output_dir:
        merged_table_path = os.path.join(output_dir, 'merged_feature_table.csv')
        prediction_feature_table.to_csv(merged_table_path, index=False)
        logging.info(f'Merged feature table saved to: {merged_table_path}')

    return prediction_feature_table

def predict_with_feature_selection_models(feature_selection_dir, prediction_feature_table, single_strain_mode=False, threads=4):
    """
    Runs predictions using models from all feature selection runs and returns predictions with confidence scores.
    
    Args:
        feature_selection_dir (str): Directory containing feature selection run subdirectories.
        prediction_feature_table (pd.DataFrame): Table containing features for prediction.
        single_strain_mode (bool): If True, only strain features are used.
        threads (int): Number of threads to use for prediction.
        
    Returns:
        pd.DataFrame: DataFrame containing all predictions with confidence scores.
    """
    all_predictions_df = pd.DataFrame(columns=['Prediction', 'strain', 'phage', 'run', 'Confidence'] if not single_strain_mode else ['Prediction', 'strain', 'run', 'Confidence'])
    
    # Get all run directories
    run_dirs = [d for d in os.listdir(feature_selection_dir) if d.startswith('run_') and os.path.isdir(os.path.join(feature_selection_dir, d))]
    
    if not run_dirs:
        logging.error(f"No run directories found in {feature_selection_dir}")
        return all_predictions_df
    
    logging.info(f"Found {len(run_dirs)} feature selection runs")
    
    for run_dir in run_dirs:
        run_path = os.path.join(feature_selection_dir, run_dir)
        model_file = os.path.join(run_path, "best_model.pkl")
        
        if not os.path.exists(model_file):
            model_file = os.path.join(run_path, "best_model.cbm")
            if not os.path.exists(model_file):
                logging.warning(f"Model file not found for {run_dir}: Skipping.")
                continue
        
        try:
            model = load_model(model_file)
        except Exception as e:
            logging.error(f"Error loading model for {run_dir}: {e}")
            continue
        
        # Extract identifiers before dropping them
        target_features = prediction_feature_table[['strain', 'phage']] if not single_strain_mode else prediction_feature_table[['strain']]
        
        # Drop identifiers for prediction
        target_features_testing = prediction_feature_table.drop(columns=['strain', 'phage'] if not single_strain_mode else ['strain'])
        
        try:
            # Align features with what the model expects
            aligned_features = align_feature_names(model, target_features_testing)
            
            # Make predictions
            predictions = model.predict(aligned_features, thread_count=threads)
            y_proba = model.predict_proba(aligned_features, thread_count=threads)[:, 1]
            
            predictions_df_temp = pd.DataFrame({
                'Prediction': predictions,
                'Confidence': y_proba,
                'strain': target_features['strain'],
                'run': run_dir
            })
            
            if not single_strain_mode:
                predictions_df_temp['phage'] = target_features['phage']
            
            all_predictions_df = pd.concat([all_predictions_df, predictions_df_temp], ignore_index=True)
            logging.info(f"Made predictions using model from {run_dir}")
        except Exception as e:
            logging.error(f"Error making predictions with model from {run_dir}: {e}")
            continue
    
    return all_predictions_df

def calculate_median_predictions(all_predictions_df, single_strain_mode=False):
    """
    Calculate median confidence score and generate final predictions based on the median confidence.
    """
    group_cols = ['strain', 'phage'] if not single_strain_mode else ['strain']
    median_conf_df = all_predictions_df.groupby(group_cols).agg({
        'Confidence': 'median'
    }).reset_index()
    
    median_conf_df['Final_Prediction'] = (median_conf_df['Confidence'] > 0.5).astype(int)
    
    return median_conf_df

def find_validation_strains(bootstrap_dir, iteration=None):
    """
    Find validation strain lists from bootstrap iterations.
    
    Args:
        bootstrap_dir (str): Directory containing bootstrap iterations.
        iteration (int or None): Specific iteration to use, or None to check all.
        
    Returns:
        str: Path to validation strains file, or None if not found.
    """
    if iteration:
        # Check specific iteration
        iteration_dir = os.path.join(bootstrap_dir, f"iteration_{iteration}")
        validation_path = os.path.join(iteration_dir, "validation_strains.csv")
        if os.path.exists(validation_path):
            return validation_path
    else:
        # Check all iterations
        iteration_dirs = [d for d in os.listdir(bootstrap_dir) if d.startswith("iteration_")]
        for iter_dir in sorted(iteration_dirs):
            validation_path = os.path.join(bootstrap_dir, iter_dir, "validation_strains.csv")
            if os.path.exists(validation_path):
                return validation_path
    
    return None

def run_feature_selection_predict_workflow(
    input_dir, 
    output_dir, 
    bootstrap_dir=None,
    iteration=None,
    use_existing_features=False,
    threads=4, 
    suffix='faa', 
    sensitivity=7.5, 
    coverage=0.8, 
    min_seq_id=0.6, 
    duplicate_all=False
):
    """
    Full workflow for predicting interactions using feature selection models.
    
    Args:
        input_dir (str): Directory containing genome FASTA files.
        output_dir (str): Directory to save results (optional when bootstrap_dir/iteration provided).
        bootstrap_dir (str, optional): Directory containing bootstrap iterations.
        iteration (int, optional): Specific bootstrap iteration to use.
        use_existing_features (bool): Whether to use existing feature assignments.
        threads (int): Number of threads to use.
        suffix (str): Suffix for FASTA files.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        duplicate_all (bool): Duplicate all genomes in the feature table for prediction.
    """
    # Determine the iteration directory and model paths
    if bootstrap_dir and iteration:
        iteration_dir = os.path.join(bootstrap_dir, f"iteration_{iteration}")
        if not os.path.exists(iteration_dir):
            logging.error(f"Iteration directory not found: {iteration_dir}")
            return
            
        # Output directory is within the iteration's model_validation directory
        base_validation_dir = os.path.join(iteration_dir, "model_validation")
        fs_predictions_dir = os.path.join(base_validation_dir, "fs_predictions")
        
        # Model component paths are within the iteration directory
        mmseqs_db = os.path.join(iteration_dir, "tmp", "strain", "mmseqs_db")
        clusters_tsv = os.path.join(iteration_dir, "strain", "clusters.tsv")
        feature_map = os.path.join(iteration_dir, "strain", "features", "selected_features.csv")
        feature_selection_dir = os.path.join(iteration_dir, "feature_selection")
        
        # Check for phage feature table within the iteration
        auto_phage_feature_path = os.path.join(iteration_dir, "phage", "features", "feature_table.csv")
        if os.path.exists(auto_phage_feature_path):
            logging.info(f"Found phage feature table: {auto_phage_feature_path}")
            phage_feature_table_path = auto_phage_feature_path
        else:
            logging.info("No phage feature table found in the iteration directory.")
            phage_feature_table_path = None
            
    else:
        # Use provided output directory or create a default one
        if not output_dir:
            from datetime import datetime
            date_str = datetime.now().strftime("%Y-%m-%d")
            output_dir = f"./fs_predictions_{date_str}"
            logging.info(f"No output directory specified. Using default: {output_dir}")
        fs_predictions_dir = output_dir
        
        # Error out since we need the model components
        logging.error("When not using bootstrap_dir and iteration, you must specify model component paths.")
        return
    
    # Create subdirectories within the fs_predictions_dir
    assign_output_dir = os.path.join(fs_predictions_dir, "assign_results")
    predict_output_dir = os.path.join(fs_predictions_dir, "predict_results")
    tmp_dir = os.path.join(fs_predictions_dir, "tmp")
    
    os.makedirs(fs_predictions_dir, exist_ok=True)
    os.makedirs(assign_output_dir, exist_ok=True)
    os.makedirs(predict_output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Set up a log file for this run
    log_file = os.path.join(fs_predictions_dir, "fs_predict_workflow.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)
    
    logging.info(f"Output will be saved to: {fs_predictions_dir}")
    
    # Check if the necessary files exist
    for file_path in [mmseqs_db, clusters_tsv, feature_map]:
        if not os.path.exists(file_path):
            logging.error(f"Required file not found: {file_path}")
            return
    
    if not os.path.exists(feature_selection_dir) or not os.path.isdir(feature_selection_dir):
        logging.error(f"Feature selection directory not found: {feature_selection_dir}")
        return
    
    # Check for validation strains if not explicitly provided
    genome_list_path = None
    if bootstrap_dir and iteration and not genome_list_path:
        validation_file = os.path.join(iteration_dir, "validation_strains.csv")
        if os.path.exists(validation_file):
            genome_list_path = validation_file
            logging.info(f"Using validation strains from: {validation_file}")
        else:
            logging.warning(f"No validation strains found in: {validation_file}")
    
    # Check for modified_AAs directory
    modified_aa_dir = os.path.join(iteration_dir, 'strain', 'modified_AAs', 'strain') if bootstrap_dir and iteration else None
    if modified_aa_dir and os.path.exists(modified_aa_dir):
        logging.info(f"Found modified_AAs directory: {modified_aa_dir}")
        logging.info("Using modified_AAs directory for feature assignment")
        input_dir_to_use = modified_aa_dir
    else:
        input_dir_to_use = input_dir
        logging.info("No modified_AAs directory found. Using original input directory.")
    
    # Decide whether to use existing feature assignments or create new ones
    existing_feature_table = None
    if use_existing_features and bootstrap_dir and iteration:
        # Path to existing feature table from validation
        validation_dir = os.path.join(iteration_dir, "model_validation")
        existing_feature_path = os.path.join(validation_dir, "assign_results", "strain_combined_feature_table.csv")
        
        if os.path.exists(existing_feature_path):
            logging.info(f"Using existing feature assignments from: {existing_feature_path}")
            existing_feature_table = existing_feature_path
    
    # Run feature assignment or use existing
    strain_feature_path = os.path.join(assign_output_dir, "strain_combined_feature_table.csv")
    
    if existing_feature_table and use_existing_features:
        # Copy existing feature table to our output directory
        import shutil
        shutil.copy2(existing_feature_table, strain_feature_path)
        logging.info(f"Copied existing feature table to: {strain_feature_path}")
    else:
        # Run new feature assignment
        logging.info("Assigning features to input genomes...")
        run_assign_features_workflow(
            input_dir=input_dir_to_use,
            mmseqs_db=mmseqs_db,
            tmp_dir=tmp_dir,
            output_dir=assign_output_dir,
            feature_map=feature_map,
            clusters_tsv=clusters_tsv,
            genome_type='strain',
            genome_list=genome_list_path,
            sensitivity=sensitivity,
            coverage=coverage,
            min_seq_id=min_seq_id,
            threads=threads,
            suffix=suffix,
            duplicate_all=duplicate_all
        )
    
    # Check if feature assignment was successful
    if not os.path.exists(strain_feature_path):
        logging.error(f"Feature table not found: {strain_feature_path}")
        return
    
    # Load assigned features
    input_table = pd.read_csv(strain_feature_path)
    
    # Load phage feature table if available
    single_strain_mode = phage_feature_table_path is None
    phage_feature_table = None
    
    if not single_strain_mode:
        if not os.path.exists(phage_feature_table_path):
            logging.error(f"Phage feature table not found: {phage_feature_table_path}")
            return
        phage_feature_table = pd.read_csv(phage_feature_table_path)
        logging.info(f"Using phage feature table with {len(phage_feature_table)} entries")
    
    # Generate full feature table
    logging.info("Generating full feature table...")
    prediction_feature_table = generate_full_feature_table(
        input_table, 
        phage_feature_table, 
        output_dir=predict_output_dir
    )
    
    # Run predictions using feature selection models
    logging.info("Running predictions using feature selection models...")
    all_predictions_df = predict_with_feature_selection_models(
        feature_selection_dir, 
        prediction_feature_table, 
        single_strain_mode, 
        threads
    )
    
    if all_predictions_df.empty:
        logging.error("No predictions were made. Check the feature selection models.")
        return
    
    # Save all predictions
    all_predictions_path = os.path.join(predict_output_dir, "fs_models_all_predictions.csv")
    all_predictions_df.to_csv(all_predictions_path, index=False)
    logging.info(f"All predictions saved to: {all_predictions_path}")
    
    # Calculate median predictions
    logging.info("Calculating median predictions...")
    median_predictions_df = calculate_median_predictions(all_predictions_df, single_strain_mode)
    
    # Save median predictions
    median_predictions_path = os.path.join(predict_output_dir, "fs_models_median_predictions.csv")
    median_predictions_df.to_csv(median_predictions_path, index=False)
    logging.info(f"Median predictions saved to: {median_predictions_path}")
    
    # Compare predictions by threshold
    run_counts = defaultdict(int)
    for run in all_predictions_df['run'].unique():
        run_counts[run] += 1
    
    logging.info(f"Used {len(run_counts)} feature selection models for prediction")
    
    # Calculate prediction statistics at different thresholds
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    threshold_stats = []
    
    for threshold in thresholds:
        median_predictions_df['Prediction_Threshold'] = (median_predictions_df['Confidence'] > threshold).astype(int)
        positive_count = median_predictions_df['Prediction_Threshold'].sum()
        total_count = len(median_predictions_df)
        
        threshold_stats.append({
            'Threshold': threshold,
            'Positive_Count': positive_count,
            'Total_Count': total_count,
            'Positive_Percentage': (positive_count / total_count) * 100 if total_count > 0 else 0
        })
    
    threshold_df = pd.DataFrame(threshold_stats)
    threshold_path = os.path.join(fs_predictions_dir, "threshold_summary.csv")
    threshold_df.to_csv(threshold_path, index=False)
    logging.info(f"Threshold summary saved to: {threshold_path}")
    
    logging.info("Feature selection prediction workflow completed successfully")


def run_all_iterations(
    input_dir, 
    bootstrap_dir, 
    use_existing_features=False,
    threads=4, 
    suffix='faa', 
    sensitivity=7.5, 
    coverage=0.8, 
    min_seq_id=0.6, 
    duplicate_all=False
):
    """
    Run the feature selection prediction workflow for all iterations in the bootstrap directory.
    
    Args:
        input_dir (str): Directory containing genome FASTA files.
        bootstrap_dir (str): Directory containing bootstrap iterations.
        use_existing_features (bool): Whether to use existing feature assignments.
        threads (int): Number of threads to use.
        suffix (str): Suffix for FASTA files.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        duplicate_all (bool): Duplicate all genomes in the feature table for prediction.
    """
    # Set up a log file for the entire run
    log_file = os.path.join(bootstrap_dir, "fs_all_iterations.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)
    
    # Find all iteration directories
    iteration_dirs = [d for d in os.listdir(bootstrap_dir) if d.startswith("iteration_") and os.path.isdir(os.path.join(bootstrap_dir, d))]
    
    if not iteration_dirs:
        logging.error(f"No iteration directories found in {bootstrap_dir}")
        return
    
    logging.info(f"Found {len(iteration_dirs)} iterations in {bootstrap_dir}")
    
    # Process each iteration
    for iter_dir in sorted(iteration_dirs):
        try:
            # Extract iteration number
            iteration = int(iter_dir.split('_')[1])
            logging.info(f"Processing {iter_dir} (iteration {iteration})...")
            
            # Run workflow for this iteration
            run_feature_selection_predict_workflow(
                input_dir=input_dir,
                output_dir=None,  # Will be determined based on bootstrap_dir and iteration
                bootstrap_dir=bootstrap_dir,
                iteration=iteration,
                use_existing_features=use_existing_features,
                threads=threads,
                suffix=suffix,
                sensitivity=sensitivity,
                coverage=coverage,
                min_seq_id=min_seq_id,
                duplicate_all=duplicate_all
            )
        except Exception as e:
            logging.error(f"Error processing {iter_dir}: {e}")
            continue
    
    # Compile summary of all iterations
    compile_cross_iteration_summary(bootstrap_dir)
    
    logging.info(f"Completed processing all iterations in {bootstrap_dir}")


def compile_cross_iteration_summary(bootstrap_dir):
    """
    Compile a summary of predictions across all iterations.
    
    Args:
        bootstrap_dir (str): Directory containing bootstrap iterations.
    """
    all_median_predictions = []
    
    # Find all iteration directories
    iteration_dirs = [d for d in os.listdir(bootstrap_dir) if d.startswith("iteration_") and os.path.isdir(os.path.join(bootstrap_dir, d))]
    
    for iter_dir in sorted(iteration_dirs):
        # Path to median predictions file
        median_file = os.path.join(bootstrap_dir, iter_dir, "model_validation", "fs_predictions", 
                                   "predict_results", "fs_models_median_predictions.csv")
        
        if os.path.exists(median_file):
            try:
                median_df = pd.read_csv(median_file)
                median_df['iteration'] = iter_dir
                all_median_predictions.append(median_df)
            except Exception as e:
                logging.warning(f"Error reading predictions from {iter_dir}: {e}")
    
    if not all_median_predictions:
        logging.warning("No prediction results found across iterations")
        return
    
    # Combine results from all iterations
    combined_df = pd.concat(all_median_predictions, ignore_index=True)
    
    # Save combined results
    output_file = os.path.join(bootstrap_dir, "fs_models_all_iterations_summary.csv")
    combined_df.to_csv(output_file, index=False)
    logging.info(f"Compiled summary from {len(all_median_predictions)} iterations saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Predict interactions using feature selection models.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--output_dir', type=str, help="Directory to save results. Optional when bootstrap_dir and iteration are provided.")
    parser.add_argument('--bootstrap_dir', type=str, help="Directory containing bootstrap iterations.")
    parser.add_argument('--iteration', type=int, help="Specific bootstrap iteration to use.")
    parser.add_argument('--all_iterations', action='store_true', help="Process all iterations in bootstrap directory.")
    parser.add_argument('--use_existing_features', action='store_true', help="Use existing feature assignments from validation.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads to use.")
    parser.add_argument('--suffix', type=str, default='faa', help="Suffix for FASTA files.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.6, help="Minimum sequence identity for assignment.")
    parser.add_argument('--duplicate_all', action='store_true', help="Duplicate all genomes in the feature table for prediction.")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.all_iterations and args.iteration is not None:
        parser.error("Cannot use --all_iterations and --iteration together")
    
    if args.all_iterations and not args.bootstrap_dir:
        parser.error("--all_iterations requires --bootstrap_dir to be specified")
    
    if not args.bootstrap_dir and args.iteration is not None:
        parser.error("--iteration requires --bootstrap_dir to be specified")
    
    if not args.output_dir and not args.bootstrap_dir:
        parser.error("--output_dir is required when not using --bootstrap_dir")
    
    if args.all_iterations:
        run_all_iterations(
            input_dir=args.input_dir,
            bootstrap_dir=args.bootstrap_dir,
            use_existing_features=args.use_existing_features,
            threads=args.threads,
            suffix=args.suffix,
            sensitivity=args.sensitivity,
            coverage=args.coverage,
            min_seq_id=args.min_seq_id,
            duplicate_all=args.duplicate_all
        )
    else:
        run_feature_selection_predict_workflow(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            bootstrap_dir=args.bootstrap_dir,
            iteration=args.iteration,
            use_existing_features=args.use_existing_features,
            threads=args.threads,
            suffix=args.suffix,
            sensitivity=args.sensitivity,
            coverage=args.coverage,
            min_seq_id=args.min_seq_id,
            duplicate_all=args.duplicate_all
        )

if __name__ == "__main__":
    main()