import os
import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import time
import joblib
from feature_selection import load_and_prepare_data, filter_data, train_and_evaluate, grid_search, save_feature_importances

# Set environment variables to control threading
os.environ['OMP_NUM_THREADS'] = '6'
os.environ['OPENBLAS_NUM_THREADS'] = '6'
os.environ['MKL_NUM_THREADS'] = '6'
os.environ['VECLIB_MAXIMUM_THREADS'] = '6'
os.environ['NUMEXPR_NUM_THREADS'] = '6'

def model_testing_select_MCC(input, output_dir, threads, random_state, set_filter='none', sample_column=None, phenotype_column=None):
    """
    Runs a single experiment for feature table, training a CatBoost model with grid search and saving results.
    
    Args:
        input (str): Path to input feature table.
        output_dir (str): Directory to store results.
        threads (int): Number of threads for training.
        random_state (int): Seed for reproducibility.
        set_filter (str): Filter type for the dataset ('none', 'host', 'phage', 'dataset').
        sample_column (str): Name of the sample column (optional).
        phenotype_column (str): Name of the phenotype column (optional).
    """
    start_time = time.time()
    
    X, y, full_feature_table = load_and_prepare_data(input, sample_column, phenotype_column)
    
    X_train, X_test, y_train, y_test, X_test_sample_ids = filter_data(X, y, full_feature_table, set_filter, sample_column=sample_column if sample_column else 'strain', random_state=random_state)
    print(f"Training data shape: {X_train.shape}, Testing data shape: {X_test.shape}")

    param_grid = {
        'iterations': [500, 1000],
        'learning_rate': [0.05, 0.1],
        'depth': [4, 6],
        'loss_function': ['Logloss'],
        'thread_count': [threads]
    }

    best_model, best_params, best_mcc = grid_search(X_train, y_train, X_test, y_test, X_test_sample_ids, param_grid, output_dir)

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

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")

def run_experiments(input_dir, base_output_dir, threads, num_runs, set_filter='none', sample_column=None, phenotype_column=None):
    """
    Iterates through feature tables in a directory, running the model testing process for each.

    Args:
        input_dir (str): Directory containing feature tables.
        base_output_dir (str): Base directory to store results.
        threads (int): Number of threads for training.
        num_runs (int): Number of runs to perform per table.
        set_filter (str): Filter type for the dataset ('none', 'host', 'phage', 'dataset').
        sample_column (str): Name of the sample column (optional).
        phenotype_column (str): Name of the phenotype column (optional).
    """
    start_total_time = time.time()
    
    feature_tables = os.listdir(input_dir)
    for feature_table in feature_tables:
        feature_table_path = os.path.join(input_dir, feature_table)
        model_output_dir = os.path.join(base_output_dir, os.path.splitext(feature_table)[0])

        if not os.path.exists(model_output_dir):
            os.makedirs(model_output_dir)

        top_models_df = pd.DataFrame()

        for i in tqdm(range(num_runs), desc=f"Running Experiments for {feature_table}"):
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

        top_models_summary_path = os.path.join(model_output_dir, 'top_models_summary.csv')
        top_models_df.to_csv(top_models_summary_path, index=False)
        print(f"Top models saved to {top_models_summary_path}")
    
    end_total_time = time.time()
    print(f"All experiments completed in {end_total_time - start_total_time:.2f} seconds.")

# To run the experiments
# Example usage:
# run_experiments(input_dir='path/to/feature/tables', base_output_dir='path/to/save/results', threads=4, num_runs=50)
