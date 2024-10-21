import os
import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef, roc_curve, precision_recall_curve
from plotnine import ggplot, aes, geom_line, geom_abline, labs, theme, guides, guide_legend, element_rect, element_blank, element_line, geom_vline, geom_jitter
from tqdm import tqdm
import time
import joblib
import re
from phage_modeling.feature_selection import load_and_prepare_data, filter_data, train_and_evaluate, grid_search, save_feature_importances
import shap
import matplotlib.pyplot as plt

# Set environment variables to control threading
os.environ['OMP_NUM_THREADS'] = '12'
os.environ['OPENBLAS_NUM_THREADS'] = '12'
os.environ['MKL_NUM_THREADS'] = '12'
os.environ['VECLIB_MAXIMUM_THREADS'] = '12'
os.environ['NUMEXPR_NUM_THREADS'] = '12'

def plot_custom_shap_beeswarm(shap_values_df, output_dir, prefix=None):
    """
    Plots and saves a custom SHAP beeswarm plot using plotnine.

    Args:
        shap_values_df (DataFrame): DataFrame containing SHAP values, feature values, and feature importance.
        output_dir (str): Directory to save the plot.
    """
    shap_values_df['value'] = shap_values_df['value'].astype(str)
    shap_values_df['shap_importance'] = shap_values_df.groupby('feature')['shap_value'].transform(lambda x: abs(x).mean())

    # Select the top 20 features by SHAP importance
    top_20_shap_df = shap_values_df.sort_values('shap_importance', ascending=False).drop_duplicates('feature').head(20)
    full_shap_values_df_top20 = shap_values_df[shap_values_df['feature'].isin(top_20_shap_df['feature'])]
    full_shap_values_df_top20['value'] = full_shap_values_df_top20['value'].astype(str)

    # Custom beeswarm plot using plotnine
    shap_plot = (
        ggplot(full_shap_values_df_top20, aes(y='reorder(feature, shap_importance)', x='shap_value', fill='value', group='feature')) +
        geom_vline(xintercept=0, color='black', size=0.5) +
        geom_jitter(height=0.2, alpha=0.3, stroke=0) +
        labs(x='SHAP Value', y='Feature', fill='Feature\nValue') +
        guides(fill=guide_legend(override_aes={'size': 6})) +
        theme(figure_size=(6, 8),
              panel_background=element_rect(fill='white'),
              panel_grid_major_x=element_blank(),
              panel_grid_minor_x=element_blank(),
              panel_grid_major_y=element_line(color='lightgrey', size=0.5),
              panel_grid_minor_y=element_blank(),
              axis_line_x=element_line(color='black', size=0.8),
              axis_line_y=element_blank())
    )

    # Save the custom beeswarm plot
    if prefix:
        shap_plot_path = os.path.join(output_dir, f"{prefix}_shap_summary_jitter.png")
    else:
        shap_plot_path = os.path.join(output_dir, "shap_summary_jitter.png")
    shap_plot.save(shap_plot_path)
    print(f"Custom SHAP beeswarm plot saved to {shap_plot_path}")

def model_testing_select_MCC(input, output_dir, threads, random_state, set_filter='none', sample_column=None, phenotype_column=None):
    """
    Runs a single experiment for feature table, training a CatBoost model with grid search and saving results.
    
    Args:
        input (str): Path to input feature table.
        output_dir (str): Directory to store results.
        threads (int): Number of threads for training.
        random_state (int): Seed for reproducibility.
        set_filter (str): Filter type for the dataset ('none', 'strain', 'phage', 'dataset').
        sample_column (str): Name of the sample column (optional).
        phenotype_column (str): Name of the phenotype column (optional).
    """
    start_time = time.time()
    
    # Load and prepare data
    X, y, full_feature_table = load_and_prepare_data(input, sample_column, phenotype_column)
    X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids = filter_data(
        X, y, full_feature_table, set_filter, sample_column=sample_column if sample_column else 'strain', random_state=random_state
    )
    print(f"Training data shape: {X_train.shape}, Testing data shape: {X_test.shape}")

    # Define grid search parameters
    param_grid = {
        'iterations': [500, 1000],
        'learning_rate': [0.05, 0.1],
        'depth': [4, 6],
        'loss_function': ['Logloss'],
        'thread_count': [threads]
    }

    # Run grid search to get the best model
    best_model, best_params, best_mcc = grid_search(X_train, y_train, X_test, y_test, X_test_sample_ids, param_grid, output_dir, phenotype_column)

    if best_model is None:
        print(f"Skipping iteration: Best MCC is {best_mcc}, no model found.")
        return

    print(f"Best Model Parameters: {best_params}, MCC: {best_mcc}")

    # Save feature importances
    feature_importances_path = os.path.join(output_dir, "feature_importances.csv")
    save_feature_importances(best_model, X_train, feature_importances_path)
    
    # Save model
    best_model_path = os.path.join(output_dir, "best_model.pkl")
    with open(best_model_path, 'wb') as f:
        joblib.dump(best_model, f)
    print(f"Best model saved to {best_model_path}")

    # Calculate and save SHAP values
    explainer = shap.TreeExplainer(best_model)
    shap_values = explainer.shap_values(X_train)

    # Save SHAP values as .npy
    shap_values_npy_path = os.path.join(output_dir, "shap_values.npy")
    np.save(shap_values_npy_path, shap_values)
    print(f"SHAP values saved as .npy to {shap_values_npy_path}")
    
    # Save SHAP values as CSV for each feature
    X_train_sample_ids = X_train_sample_ids.reset_index(drop=True)
    X_train_df = X_train.copy().reset_index(drop=True)
    X_train_id_columns = list(X_train_sample_ids.columns)

    shap_values_df = pd.DataFrame(shap_values, columns=X_train.columns)
    shap_values_df[X_train_id_columns] = X_train_sample_ids.reset_index(drop=True)
    shap_values_df = shap_values_df.melt(id_vars=X_train_id_columns, var_name='feature', value_name='shap_value')

    X_train_df[X_train_id_columns] = X_train_sample_ids.reset_index(drop=True)
    X_train_df = X_train_df.melt(id_vars=X_train_id_columns, var_name='feature', value_name='value')

    shap_values_df = shap_values_df.merge(X_train_df, on=X_train_id_columns + ['feature'], how='left')

    shap_values_csv_path = os.path.join(output_dir, "shap_importances.csv")
    shap_values_df.to_csv(shap_values_csv_path, index=False)
    print(f"SHAP values saved to {shap_values_csv_path}")

    # Plot and save SHAP summary bar plot
    shap_summary_bar_path = os.path.join(output_dir, "shap_summary_bar.png")
    shap.summary_plot(shap_values, X_train, plot_type="bar", show=False)
    plt.savefig(shap_summary_bar_path)
    plt.close()
    print(f"SHAP summary bar plot saved to {shap_summary_bar_path}")

    # Custom SHAP beeswarm plot using plotnine
    plot_custom_shap_beeswarm(shap_values_df, output_dir)

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")

def extract_cutoff_from_filename(filename):
    """
    Extracts the numeric cutoff value from the filename (e.g., 7 from select_feature_table_cutoff_7.csv).

    Args:
        filename (str): The filename of the feature table.

    Returns:
        str: The extracted cutoff value.
    """
    match = re.search(r'cutoff_(\d+)', filename)
    if match:
        return match.group(1)  # Return the cutoff value as a string
    else:
        raise ValueError(f"Could not extract cutoff from filename: {filename}")

def parse_model_predictions_and_performance(model_dir):
    """
    Parses the prediction and performance data from the run directories within the specified model directory.

    Args:
        model_dir (str): Base directory for the select feature modeling containing cutoff subdirectories.

    Returns:
        model_predictions_df (DataFrame): Combined DataFrame of all model predictions.
        model_performance_df (DataFrame): Combined DataFrame of model performance (top models summary).
    """
    model_predictions_df = pd.DataFrame()
    model_performance_df = pd.DataFrame()

    # Get all cutoff subdirectories (excluding CSV files)
    cut_offs = [x for x in os.listdir(model_dir) if '.csv' not in x]

    # Loop through each cutoff directory and parse run directories
    for cut_off in cut_offs:
        print(f"Parsing cut-off: {cut_off}")
        cut_off_dir = os.path.join(model_dir, cut_off)
        run_dirs = [x for x in os.listdir(cut_off_dir) if 'run' in x]

        # Parse predictions from each run
        for run in run_dirs:
            predictions_file = os.path.join(cut_off_dir, run, 'best_model_predictions.csv')
            if os.path.exists(predictions_file):
                predictions_temp = pd.read_csv(predictions_file)
                predictions_temp['run'] = run
                predictions_temp['cut_off'] = cut_off
                model_predictions_df = pd.concat([model_predictions_df, predictions_temp])

        # Parse top models summary for this cutoff
        top_models_file = os.path.join(cut_off_dir, 'top_models_summary.csv')
        if os.path.exists(top_models_file):
            top_models_temp = pd.read_csv(top_models_file)
            top_models_temp = top_models_temp[['mcc']].reset_index()
            top_models_temp['index'] = ['_'.join(['run', str(x)]) for x in top_models_temp['index']]
            top_models_temp = top_models_temp.rename(columns={'index': 'run'})
            top_models_temp['cut_off'] = cut_off
            model_performance_df = pd.concat([model_performance_df, top_models_temp])

    return model_predictions_df, model_performance_df

def evaluate_model_performance(predictions_file, output_dir, sample_column='strain', phenotype_column='interaction'):
    """
    Evaluates model performances from the prediction data and generates performance plots grouped by cutoff.

    Args:
        predictions_file (str): Path to the CSV file containing model predictions.
        output_dir (str): Directory to save performance plots and evaluation metrics.
        sample_column (str): Column name for the sample identifier.
        phenotype_column (str): Column name for the phenotype.
    """
    # Create the output directory if it doesn't exist
    model_performance_dir = os.path.join(output_dir, 'model_performance')
    os.makedirs(model_performance_dir, exist_ok=True)

    # Load model predictions
    model_predictions_df_full = pd.read_csv(predictions_file)
    model_predictions_df_full['cut_off'] = model_predictions_df_full['cut_off'].astype(str)

    # Determine grouping columns
    grouping_columns = ['cut_off', sample_column, phenotype_column]
    if 'phage' in model_predictions_df_full.columns:
        grouping_columns.insert(2, 'phage')  # Include phage in grouping if present

    # Calculate average prediction and confidence per group
    model_predictions_df_calcs = model_predictions_df_full.groupby(grouping_columns).agg({
        'Prediction': 'mean',
        'Confidence': 'mean'
    }).reset_index()

    # Convert confidence to binary predictions
    model_predictions_df_calcs['Prediction'] = (model_predictions_df_calcs['Confidence'] > 0.5).astype(int)

    # Function to calculate metrics
    def calculate_metrics(df):
        y_true = df[phenotype_column]
        y_pred_prob = df['Confidence']
        y_pred = df['Prediction']
        
        metrics = {
            'AUC': roc_auc_score(y_true, y_pred_prob),
            'Accuracy': accuracy_score(y_true, y_pred),
            'Precision': precision_score(y_true, y_pred),
            'Recall': recall_score(y_true, y_pred),
            'F1': f1_score(y_true, y_pred),
            'MCC': matthews_corrcoef(y_true, y_pred)
        }
        
        fpr, tpr, roc_thresholds = roc_curve(y_true, y_pred_prob)
        roc_df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'thresholds': roc_thresholds, 'cut_off': df['cut_off'].iloc[0]})
        
        precision, recall, pr_thresholds = precision_recall_curve(y_true, y_pred_prob)
        pr_df = pd.DataFrame({'precision': precision, 'recall': recall, 'thresholds': list(pr_thresholds) + [1.0], 'cut_off': df['cut_off'].iloc[0]})
        
        return metrics, roc_df, pr_df

    # Group by model and cut_off and calculate metrics
    metrics_list = []
    roc_list = []
    pr_list = []

    for cut_off, group in model_predictions_df_calcs.groupby('cut_off'):
        metrics, roc_df, pr_df = calculate_metrics(group)
        metrics['cut_off'] = cut_off
        metrics_list.append(metrics)
        roc_list.append(roc_df)
        pr_list.append(pr_df)

    # Convert results to dataframes
    metrics_df = pd.DataFrame(metrics_list)
    roc_df = pd.concat(roc_list).reset_index(drop=True)
    pr_df = pd.concat(pr_list).reset_index(drop=True)

    # Plot and save ROC curve
    roc_curve_plot = (
        ggplot(roc_df, aes(x='fpr', y='tpr', color='cut_off')) +
        geom_line() +
        geom_abline(linetype='dashed') +
        labs(x='False Positive Rate', y='True Positive Rate', color='Cut Off') +
        theme(figure_size=(5, 4))
    )
    roc_curve_plot.save(os.path.join(model_performance_dir, 'roc_curve.png'))

    # Plot and save Precision-Recall curve
    pr_curve_plot = (
        ggplot(pr_df, aes(x='recall', y='precision', color='cut_off')) +
        geom_line() +
        labs(x='Recall', y='Precision', color='Cut Off') +
        theme(figure_size=(5, 4))
    )
    pr_curve_plot.save(os.path.join(model_performance_dir, 'pr_curve.png'))

    # Calculate and plot hit rate and hit ratio
    def calculate_hit_rate(df):
        df = df.sort_values(by='Confidence', ascending=False).reset_index(drop=True)
        df['cumulative_hits'] = df[phenotype_column].cumsum()
        df['cumulative_total'] = np.arange(1, len(df) + 1)
        df['hit_rate'] = df['cumulative_hits'] / df['cumulative_total']
        df['fraction_of_samples'] = df['cumulative_total'] / len(df)
        return df

    def calculate_hit_ratio(df):
        df = df.sort_values(by='Confidence', ascending=False).reset_index(drop=True)
        df['cumulative_true_positives'] = df[phenotype_column].cumsum()
        df['total_true_positives'] = df[phenotype_column].sum()
        df['hit_ratio'] = df['cumulative_true_positives'] / df['total_true_positives']
        df['fraction_of_samples'] = np.arange(1, len(df) + 1) / len(df)
        return df

    hit_rate_list = []
    hit_ratio_list = []
    for cut_off, group in model_predictions_df_calcs.groupby('cut_off'):
        hit_rate_df = calculate_hit_rate(group)
        hit_rate_df['cut_off'] = cut_off
        hit_rate_list.append(hit_rate_df)

        hit_ratio_df = calculate_hit_ratio(group)
        hit_ratio_df['cut_off'] = cut_off
        hit_ratio_list.append(hit_ratio_df)

    hit_rate_df = pd.concat(hit_rate_list).reset_index(drop=True)
    hit_ratio_df = pd.concat(hit_ratio_list).reset_index(drop=True)

    # Plot and save Hit Rate curve
    hit_rate_curve_plot = (
        ggplot(hit_rate_df, aes(x='fraction_of_samples', y='hit_rate', color='cut_off')) +
        geom_line() +
        labs(x='Fraction of Samples', y='Hit Rate', color='Cut Off') +
        theme(figure_size=(5, 4))
    )
    hit_rate_curve_plot.save(os.path.join(model_performance_dir, 'hit_rate_curve.png'))

    # Plot and save Hit Ratio curve
    hit_ratio_curve_plot = (
        ggplot(hit_ratio_df, aes(x='fraction_of_samples', y='hit_ratio', color='cut_off')) +
        geom_line() +
        labs(x='Fraction of Samples', y='Hit Ratio', color='Cut Off') +
        theme(figure_size=(5, 4))
    )
    hit_ratio_curve_plot.save(os.path.join(model_performance_dir, 'hit_ratio_curve.png'))

    # Save the metrics DataFrame
    metrics_df.to_csv(os.path.join(model_performance_dir, 'model_performance_metrics.csv'), index=False)
    print(f"Metrics saved to {os.path.join(model_performance_dir, 'model_performance_metrics.csv')}")

def run_experiments(input_dir, base_output_dir, threads, num_runs, set_filter='none', sample_column='strain', phenotype_column='interaction'):
    """
    Iterates through feature tables in a directory, running the model testing process for each.
    
    Args:
        input_dir (str): Directory containing feature tables.
        base_output_dir (str): Base directory to store results.
        threads (int): Number of threads for training.
        num_runs (int): Number of runs to perform per table.
        set_filter (str): Filter type for the dataset ('none', 'strain', 'phage', 'dataset').
        sample_column (str): Name of the sample column (optional).
        phenotype_column (str): Name of the phenotype column (optional).
    """
    start_total_time = time.time()
    feature_tables = os.listdir(input_dir)

    for feature_table in feature_tables:
        feature_table_path = os.path.join(input_dir, feature_table)
        
        # Use the extract_cutoff_from_filename to generate the directory name
        cutoff_value = extract_cutoff_from_filename(feature_table)
        model_output_dir = os.path.join(base_output_dir, f'cutoff_{cutoff_value}')  # Naming with cutoff

        if not os.path.exists(model_output_dir):
            os.makedirs(model_output_dir)

        top_models_df = pd.DataFrame()
        top_models_shap_df = pd.DataFrame()

        for i in tqdm(range(num_runs), desc=f"Running Experiments for cutoff {cutoff_value}"):
            output_dir = os.path.join(model_output_dir, f'run_{i}')
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            random_state = i

            model_performance_path = os.path.join(output_dir, 'model_performance.csv')
            if not os.path.exists(model_performance_path):
                model_testing_select_MCC(feature_table_path, output_dir, threads, random_state, set_filter, sample_column, phenotype_column)

            if os.path.exists(model_performance_path):
                run_results = pd.read_csv(model_performance_path)
                top_model = run_results.nlargest(1, 'mcc')
                top_models_df = pd.concat([top_models_df, top_model])

            shap_values_csv_path = os.path.join(output_dir, "shap_importances.csv")
            if os.path.exists(shap_values_csv_path):
                # Load SHAP values and append to top_models_shap_df
                shap_values_temp = pd.read_csv(shap_values_csv_path)
                top_models_shap_df = pd.concat([top_models_shap_df, shap_values_temp])

        # Generate SHAP summary plot for the top models
        model_performance_dir = os.path.join(base_output_dir, 'model_performance')
        os.makedirs(model_performance_dir, exist_ok=True)

        plot_custom_shap_beeswarm(top_models_shap_df, model_performance_dir, prefix=f'cutoff_{cutoff_value}')

        top_models_summary_path = os.path.join(model_output_dir, 'top_models_summary.csv')
        top_models_df.to_csv(top_models_summary_path, index=False)
        print(f"Top models saved to {top_models_summary_path}")

    # After running all experiments, parse the predictions and performance
    model_predictions_output = os.path.join(base_output_dir, 'select_features_model_predictions.csv')
    
    # Parse all predictions from the run directories
    model_predictions_df, model_performance_df = parse_model_predictions_and_performance(base_output_dir)

    # Save parsed predictions and performance data
    model_predictions_df.to_csv(model_predictions_output, index=False)
    model_performance_output = os.path.join(base_output_dir, 'select_features_model_performance.csv')
    model_performance_df.to_csv(model_performance_output, index=False)

    print(f"Model predictions saved to {model_predictions_output}")
    print(f"Model performance saved to {model_performance_output}")

    # Now evaluate model performance and generate performance plots
    evaluate_model_performance(
        predictions_file=model_predictions_output,
        output_dir=base_output_dir,
        sample_column=sample_column,
        phenotype_column=phenotype_column
    )

    end_total_time = time.time()
    print(f"All experiments completed in {end_total_time - start_total_time:.2f} seconds.")

# To run the experiments
# Example usage:
# run_experiments(input_dir='path/to/feature/tables', base_output_dir='path/to/save/results', threads=4, num_runs=50)
