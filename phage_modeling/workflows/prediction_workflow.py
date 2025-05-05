import pandas as pd
import os
import pickle
from catboost import CatBoostClassifier
from argparse import ArgumentParser
import numpy as np
import logging

def generate_full_feature_table(input_table, phage_feature_table=None, strain_source='strain', phage_source='phage', output_dir=None):
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

def align_feature_names(model, target_features_testing):
    """
    Align features in the target_features_testing DataFrame with the features expected by the model.
    Filters out extra features and raises an error if any required features are missing.
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

from catboost import CatBoostClassifier
import pandas as pd
import os
import logging

def predict_interactions(model_dir, prediction_feature_table, single_strain_mode=False, threads=4):
    """
    Runs predictions using all models in the model_dir and returns predictions with confidence scores.
    Uses multiple threads for prediction to speed up computations.

    Args:
        model_dir (str): Directory with model subdirectories.
        prediction_feature_table (pd.DataFrame): Table containing features for prediction.
        single_strain_mode (bool): If True, only strain features are used.
        threads (int): Number of threads to use for prediction. Default is 4.

    Returns:
        pd.DataFrame: DataFrame containing all predictions with confidence scores.
    """
    all_predictions_df = pd.DataFrame(columns=['Prediction', 'strain', 'phage', 'run', 'Confidence'] if not single_strain_mode else ['Prediction', 'strain', 'run', 'Confidence'])

    subdirs = [subdir for subdir in os.listdir(model_dir) if subdir.startswith('run')]

    for subdir in subdirs:
        model_subdir_path = os.path.join(model_dir, subdir)
        model_file = os.path.join(model_subdir_path, "best_model.pkl")

        if not os.path.exists(model_file):
            model_file = os.path.join(model_subdir_path, "best_model.cbm")
            if not os.path.exists(model_file):
                logging.warning(f"Model file not found for run {subdir}: Skipping.")
                continue

        try:
            model = load_model(model_file)
        except Exception as e:
            logging.error(f"Error loading model for run {subdir}: {e}")
            continue

        # Extract identifiers before dropping them
        target_features = prediction_feature_table[['strain', 'phage']] if not single_strain_mode else prediction_feature_table[['strain']]
        
        # Drop identifiers for prediction
        target_features_testing = prediction_feature_table.drop(columns=['strain', 'phage'] if not single_strain_mode else ['strain'])
        
        try:
            # Align features with what the model expects
            aligned_features = align_feature_names(model, target_features_testing)
            
            # Using thread_count for prediction
            predictions = model.predict(aligned_features, thread_count=threads)
            y_proba = model.predict_proba(aligned_features, thread_count=threads)[:, 1]

            predictions_df_temp = pd.DataFrame({
                'Prediction': predictions,
                'Confidence': y_proba,
                'strain': target_features['strain'],
                'run': subdir
            })

            if not single_strain_mode:
                predictions_df_temp['phage'] = target_features['phage']

            all_predictions_df = pd.concat([all_predictions_df, predictions_df_temp], ignore_index=True)
        except Exception as e:
            logging.error(f"Error making predictions with model from run {subdir}: {e}")
            continue

    return all_predictions_df

def calculate_median_predictions(all_predictions_df, single_strain_mode=False):
    """
    Calculate median confidence score and generate final predictions based on the average confidence.
    """
    group_cols = ['strain', 'phage'] if not single_strain_mode else ['strain']
    median_conf_df = all_predictions_df.groupby(group_cols).agg({
        'Confidence': 'median'
    }).reset_index()

    median_conf_df['Final_Prediction'] = (median_conf_df['Confidence'] > 0.5).astype(int)

    return median_conf_df

def run_prediction_workflow(input_dir, phage_feature_table_path, model_dir, output_dir, feature_table=None, strain_source='strain', phage_source='phage', threads=4):
    """
    Full workflow for predicting interactions using a combined strain feature table and optionally phage-specific features.

    Args:
        input_dir (str): Directory containing the combined strain feature table.
        phage_feature_table_path (str or None): Path to the phage feature table. If None, single-strain mode is used.
        model_dir (str): Directory with model subdirectories.
        output_dir (str): Output directory for saving predictions.
        strain_source (str): Prefix used for strain features.
        phage_source (str): Prefix used for phage features.
    """
    os.makedirs(output_dir, exist_ok=True)

    phage_feature_table = pd.read_csv(phage_feature_table_path) if phage_feature_table_path else None

    input_files = [f for f in os.listdir(input_dir) if f.endswith('_feature_table.csv')]
    if not input_files:
        logging.error(f"No strain feature table found in {input_dir}")
        return  # or raise an exception
    input_file = input_files[0]

    logging.info(f'Processing file: {input_file}')

    input_table = pd.read_csv(os.path.join(input_dir, input_file))

    single_strain_mode = phage_feature_table is None

    print('Generating full feature table...')
    prediction_feature_table = generate_full_feature_table(input_table, phage_feature_table, strain_source, phage_source, output_dir=output_dir)

    if feature_table:
        logging.info(f'Loading feature table from {feature_table}')
        feature_table = pd.read_csv(feature_table)

        select_columns = feature_table.columns
        select_columns = [col for col in select_columns if col != 'interaction']

        prediction_feature_table = prediction_feature_table[select_columns]

    logging.info('Running predictions...')
    all_predictions_df = predict_interactions(model_dir, prediction_feature_table, single_strain_mode, threads)

    all_predictions_df.to_csv(os.path.join(output_dir, f'{strain_source}_all_predictions.csv'), index=False)

    logging.info('Calculating median predictions...')
    median_predictions_df = calculate_median_predictions(all_predictions_df, single_strain_mode)
    median_predictions_df.to_csv(os.path.join(output_dir, f'{strain_source}_median_predictions.csv'), index=False)

    logging.info('Workflow completed successfully.')

def main():
    parser = ArgumentParser(description="Predict interactions using strain-specific and optionally phage-specific feature tables.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory with strain-specific feature tables.")
    parser.add_argument('--phage_feature_table', type=str, help="Path to the phage feature table. Optional for single-strain mode.")
    parser.add_argument('--feature_table', type=str, default=None, help="Path to the combined feature table. Optional for single-strain mode.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory with models (run_* subdirectories).")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save predictions.")
    parser.add_argument('--strain_source', type=str, default='strain', help="Prefix used for strain features.")
    parser.add_argument('--phage_source', type=str, default='phage', help="Prefix used for phage features.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads to use for prediction.")

    args = parser.parse_args()

    run_prediction_workflow(
        input_dir=args.input_dir,
        phage_feature_table_path=args.phage_feature_table,
        feature_table=args.feature_table,
        model_dir=args.model_dir,
        output_dir=args.output_dir,
        strain_source=args.strain_source,
        phage_source=args.phage_source,
        threads=args.threads
    )

if __name__ == '__main__':
    main()
