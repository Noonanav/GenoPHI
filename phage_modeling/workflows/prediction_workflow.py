import pandas as pd
import os
import pickle
from catboost import CatBoostClassifier
from argparse import ArgumentParser

def generate_full_feature_table(input_table, phage_feature_table, strain_source='strain', phage_source='phage'):
    """
    Concatenates strain feature table with phage feature table for predictions.
    Each row in the final table is a combination of one strain row and one phage row.
    The strain and phage features are identified based on the prefixes defined during feature selection.
    """
    # strain features (those starting with strain_source prefix, e.g., 'hc_' for strain)
    strain_features = ['strain'] + [col for col in input_table.columns if col.startswith(f'{strain_source[0]}c_')]

    # Phage features (those starting with phage_source prefix, e.g., 'pc_' for phage)
    phage_features_table = phage_feature_table[['phage'] + [col for col in phage_feature_table.columns if col.startswith(f'{phage_source[0]}c_')]]

    # Create repeated rows of strain features for each phage entry
    repeated_strain = pd.concat([input_table[strain_features]] * len(phage_features_table), ignore_index=True)

    # Create the prediction feature table by concatenating strain and phage rows
    prediction_feature_table = pd.concat([phage_features_table.reset_index(drop=True), repeated_strain], axis=1)

    return prediction_feature_table

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

def predict_interactions(model_dir, prediction_feature_table):
    """
    Runs predictions using all models in the model_dir and returns individual predictions
    as well as the mean confidence score.
    """
    all_predictions_df = pd.DataFrame(columns=['Prediction', 'strain', 'phage', 'run', 'Confidence'])

    subdirs = [subdir for subdir in os.listdir(model_dir) if subdir.startswith('run')]

    for subdir in subdirs:
        model_subdir_path = os.path.join(model_dir, subdir)
        model_file = os.path.join(model_subdir_path, "best_model.pkl")

        # If pkl model not found, use cbm
        if not os.path.exists(model_file):
            model_file = os.path.join(model_subdir_path, "best_model.cbm")

        model = load_model(model_file)

        # Get strain and phage columns
        target_features_testing = prediction_feature_table.drop(columns=['strain', 'phage'])
        target_features = prediction_feature_table[['strain', 'phage']]

        # Run predictions
        predictions = model.predict(target_features_testing)
        y_proba = model.predict_proba(target_features_testing)[:, 1]  # Confidence score for the positive class

        # Store results
        predictions_df_temp = pd.DataFrame({
            'Prediction': predictions,
            'Confidence': y_proba,
            'strain': target_features['strain'],
            'phage': target_features['phage'],
            'run': subdir
        })

        all_predictions_df = pd.concat([all_predictions_df, predictions_df_temp], ignore_index=True)

    return all_predictions_df

def calculate_mean_predictions(all_predictions_df):
    """
    Calculate mean confidence score and generate final predictions based on the average confidence.
    """
    mean_conf_df = all_predictions_df.groupby(['strain', 'phage']).agg({
        'Confidence': 'mean'
    }).reset_index()

    # Final prediction based on mean confidence > 0.5
    mean_conf_df['Final_Prediction'] = (mean_conf_df['Confidence'] > 0.5).astype(int)

    return mean_conf_df

def run_full_prediction_workflow(input_dir, phage_feature_table_path, model_dir, output_dir, strain_source='strain', phage_source='phage'):
    """
    Full workflow for predicting phage-strain interactions using multiple models and strain-specific feature tables.

    Args:
        input_dir (str): Directory with strain-specific feature tables.
        phage_feature_table_path (str): Path to the phage feature table.
        model_dir (str): Directory with model subdirectories.
        output_dir (str): Output directory for saving predictions.
        strain_source (str): Source used for strain feature naming (prefix).
        phage_source (str): Source used for phage feature naming (prefix).
    """
    # Load phage feature table
    print('Loading phage feature table...')
    phage_feature_table = pd.read_csv(phage_feature_table_path)

    # Prepare output directory
    os.makedirs(output_dir, exist_ok=True)

    # Load strain feature tables and generate full feature table
    print('Generating full feature table...')
    input_files = [f for f in os.listdir(input_dir) if f.endswith('_feature_table.csv')]
    full_predictions_df = pd.DataFrame()

    for input_file in input_files:
        strain_name = input_file.split('_feature_table')[0]
        print(f'Processing strain: {strain_name}')

        # Load the strain-specific feature table
        input_table = pd.read_csv(os.path.join(input_dir, input_file))

        # Generate the full prediction table (strain + phage)
        prediction_feature_table = generate_full_feature_table(input_table, phage_feature_table, strain_source, phage_source)

        # Run predictions using all models
        print('Running predictions using all models...')
        all_predictions_df = predict_interactions(model_dir, prediction_feature_table)

        # Save all predictions for this strain
        all_predictions_df.to_csv(os.path.join(output_dir, f'{strain_name}_all_predictions.csv'), index=False)

        # Concatenate predictions for the final output
        full_predictions_df = pd.concat([full_predictions_df, all_predictions_df], ignore_index=True)

    # Save all predictions in one CSV
    print('Saving all predictions...')
    full_predictions_df.to_csv(os.path.join(output_dir, 'all_predictions.csv'), index=False)

    # Calculate mean confidence and final predictions
    print('Calculating mean predictions...')
    mean_predictions_df = calculate_mean_predictions(full_predictions_df)

    # Save the mean predictions
    print('Saving mean predictions...')
    mean_predictions_df.to_csv(os.path.join(output_dir, 'mean_predictions.csv'), index=False)

    print('Workflow completed successfully.')

def main():
    """
    Command-line interface (CLI) for running the full prediction workflow.
    """
    parser = ArgumentParser(description="Predict phage-strain interactions using strain-specific feature tables and multiple models.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory with strain-specific feature tables.")
    parser.add_argument('--phage_feature_table', type=str, required=True, help="Path to the phage feature table.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory with models (run_* subdirectories).")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save predictions.")
    parser.add_argument('--strain_source', type=str, default='strain', help="Prefix used for strain features.")
    parser.add_argument('--phage_source', type=str, default='phage', help="Prefix used for phage features.")

    args = parser.parse_args()

    # Call the workflow function
    run_full_prediction_workflow(
        input_dir=args.input_dir,
        phage_feature_table_path=args.phage_feature_table,
        model_dir=args.model_dir,
        output_dir=args.output_dir,
        strain_source=args.strain_source,
        phage_source=args.phage_source
    )

if __name__ == '__main__':
    main()
