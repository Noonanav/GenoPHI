import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from catboost import CatBoostClassifier
from sklearn.feature_selection import RFE
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, confusion_matrix, roc_curve, auc, precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from tqdm import tqdm
import time

# Function to load and prepare data
def load_and_prepare_data(input_path):
    """
    Loads the input feature table, drops unnecessary columns, and splits into features and target.

    Args:
        input_path (str): Path to the input CSV file containing the full feature table.

    Returns:
        X (DataFrame): Features for modeling.
        y (Series): Target variable (interaction).
        full_feature_table (DataFrame): The complete feature table after cleaning.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"The input file {input_path} does not exist.")
    
    full_feature_table = pd.read_csv(input_path)
    if full_feature_table.empty:
        raise ValueError("Input data is empty.")
    
    full_feature_table = full_feature_table.dropna()
    X = full_feature_table.drop(['strain', 'phage', 'interaction', 'header', 'contig_id', 'orf_ko'], axis=1, errors='ignore')
    y = full_feature_table['interaction']
    
    print(f"Number of positive samples: {y.sum()}")
    print(f"Number of negative samples: {len(y) - y.sum()}")
    print("Data loaded and prepared, split into features and target.")
    
    return X, y, full_feature_table

# Function to filter the data based on strain or phage
def filter_data(X, y, full_feature_table, filter_type, random_state=42):
    """
    Filters the data by strain or phage and splits into training and testing sets.

    Args:
        X (DataFrame): Features.
        y (Series): Target variable.
        full_feature_table (DataFrame): The full feature table with metadata.
        filter_type (str): 'none', 'host', or 'phage' to determine how the data should be filtered.
        random_state (int): Seed for reproducibility.

    Returns:
        X_train (DataFrame): Training features.
        X_test (DataFrame): Testing features.
        y_train (Series): Training target.
        y_test (Series): Testing target.
        X_test_sample_ids (DataFrame): Metadata of the test set samples.
    """
    if filter_type == 'none':
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
        test_idx = X_test.index
        X_test_sample_ids = full_feature_table.loc[test_idx, ['strain', 'phage']]
    else:
        if filter_type == 'host':
            group = 'strain'
        elif filter_type == 'phage':
            group = 'phage'
        else:
            raise ValueError("Invalid filter type.")
        
        groups = full_feature_table[group].unique()
        np.random.seed(random_state)
        train_groups = np.random.choice(groups, size=int(0.8 * len(groups)), replace=False)
        test_groups = np.setdiff1d(groups, train_groups)
        
        train_idx = full_feature_table[group].isin(train_groups)
        test_idx = full_feature_table[group].isin(test_groups)
        
        X_train = X[train_idx]
        X_test = X[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]
        X_test_sample_ids = full_feature_table.loc[test_idx, ['strain', 'phage']]

    return X_train, X_test, y_train, y_test, X_test_sample_ids

# Function to perform Recursive Feature Elimination (RFE)
def perform_rfe(X_train, y_train, num_features, threads, output_dir):
    """
    Performs Recursive Feature Elimination (RFE) to select the top features.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Number of features to select.
        threads (int): Number of threads to use for CatBoost.
        output_dir (str): Directory to store intermediate CatBoost information.

    Returns:
        rfe (RFE object): Fitted RFE model.
        selected_features (Index): List of selected features.
    """
    total_features = X_train.shape[1]
    step_size = max(1, int((total_features - num_features) / 10))  # Ensure step_size is at least 1

    model = CatBoostClassifier(
        iterations=500, 
        learning_rate=0.1, 
        depth=4, 
        verbose=10, 
        thread_count=threads, 
        train_dir=os.path.join(output_dir, '..', 'catboost_info')
    )
    
    print(f"Performing Recursive Feature Elimination (RFE) with step_size: {step_size}...")

    rfe = RFE(estimator=model, n_features_to_select=num_features, step=step_size, verbose=10)
    rfe.fit(X_train, y_train)
    
    selected_features = X_train.columns[rfe.support_]
    print(f"RFE selected {len(selected_features)} features.")
    
    return rfe, selected_features

# Function to perform grid search
def grid_search(X_train, y_train, X_test, y_test, X_test_sample_ids, param_grid, output_dir):
    """
    Performs grid search to find the best hyperparameters for CatBoost.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        X_test (DataFrame): Testing features.
        y_test (Series): Testing target.
        X_test_sample_ids (DataFrame): Metadata for the test set samples.
        param_grid (dict): Dictionary of hyperparameters for grid search.
        output_dir (str): Directory to save results.

    Returns:
        best_model (CatBoostClassifier): The model with the best performance.
        best_params (dict): The hyperparameters of the best model.
        best_mcc (float): The highest MCC score achieved during grid search.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    best_mcc = 0
    best_model = None
    best_params = None
    results = []
    
    print("Starting grid search...")
    for idx, params in enumerate(itertools.product(*param_grid.values()), start=1):
        params = dict(zip(param_grid.keys(), params))
        model, accuracy, f1, mcc, y_pred = train_and_evaluate(X_train, y_train, X_test, y_test, params, output_dir)
        
        results.append({**params, 'accuracy': accuracy, 'f1_score': f1, 'mcc': mcc})
        
        # Plot confusion matrix
        plot_confusion_matrix(y_test, y_pred, f"{output_dir}/conf_matrix_{idx}.png")
        # Plot ROC curve
        plot_roc_curve(y_test, model.predict_proba(X_test), f"{output_dir}/roc_curve_{idx}.png")
        # Plot precision-recall curve
        plot_precision_recall_curve(y_test, model.predict_proba(X_test), f"{output_dir}/precision_recall_curve_{idx}.png")
        
        if mcc >= best_mcc:
            best_mcc = mcc
            best_model = model
            best_params = params

            best_predictions_df = X_test_sample_ids.copy()
            best_predictions_df['Prediction'] = y_pred
            best_predictions_df['Confidence'] = model.predict_proba(X_test)[:, 1]
            best_predictions_df['interaction'] = y_test
            best_predictions_df.to_csv(f"{output_dir}/best_model_predictions.csv", index=False)

    pd.DataFrame(results).to_csv(f"{output_dir}/model_performance.csv", index=False)
    
    return best_model, best_params, best_mcc

# Utility functions to plot graphs and save feature importances
def plot_confusion_matrix(y_test, y_pred, output_path):
    """
    Plots and saves a confusion matrix.

    Args:
        y_test (Series): True labels for the test set.
        y_pred (Series): Predicted labels for the test set.
        output_path (str): Path to save the confusion matrix plot.
    """
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.title('Confusion Matrix')
    plt.savefig(output_path)
    plt.close()

def plot_roc_curve(y_test, y_scores, output_path):
    """
    Plots and saves a ROC curve.

    Args:
        y_test (Series): True labels for the test set.
        y_scores (ndarray): Predicted probabilities for the test set.
        output_path (str): Path to save the ROC curve plot.
    """
    fpr, tpr, _ = roc_curve(y_test, y_scores[:, 1])
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (area = {auc(fpr, tpr):.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.savefig(output_path)
    plt.close()

def plot_precision_recall_curve(y_test, y_scores, output_path):
    """
    Plots and saves a precision-recall curve.

    Args:
        y_test (Series): True labels for the test set.
        y_scores (ndarray): Predicted probabilities for the test set.
        output_path (str): Path to save the precision-recall curve plot.
    """
    precision, recall, _ = precision_recall_curve(y_test, y_scores[:, 1])
    plt.figure()
    plt.step(recall, precision, where='post')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.savefig(output_path)
    plt.close()

def save_feature_importances(model, X_train_selected, output_path):
    """
    Saves the feature importances from the model.

    Args:
        model (CatBoostClassifier): The trained model.
        X_train_selected (DataFrame): The training data with selected features.
        output_path (str): Path to save the feature importances.
    """
    feature_importances = model.get_feature_importance()
    selected_features = X_train_selected.columns
    importance_df = pd.DataFrame({
        'Feature': selected_features,
        'Importance': feature_importances
    }).sort_values(by='Importance', ascending=False)
    
    importance_df.to_csv(output_path, index=False)
    print(f"Feature importances saved to {output_path}")

def run_feature_selection_iterations(
    input_path, base_output_dir, threads, num_features, filter_type, num_runs, select_cols=False, sample_column=None, phenotype_column=None
):
    """
    Runs multiple iterations of feature selection, saves the results in `run_*` directories, and tracks feature occurrences.
    
    Args:
        input_path (str): Path to the input feature table.
        base_output_dir (str): Base output directory where results for each run will be stored.
        threads (int): Number of threads to use for feature selection.
        num_features (int): Number of features to select.
        filter_type (str): Filter type for the input data ('host', 'phage', 'none').
        num_runs (int): Number of runs to perform.
        select_cols (bool): Whether to run with selected columns.
        sample_column (str): Column name for the sample/strain (if using selected columns).
        phenotype_column (str): Column name for the phenotype (if using selected columns).
    """
    
    # Ensure base output directory exists
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    # Dictionary to keep track of feature occurrences across runs
    features_occurrence = {}

    start_total_time = time.time()

    # Run feature selection for each iteration
    for i in tqdm(range(num_runs), desc="Running Feature Selection Iterations"):
        output_dir = os.path.join(base_output_dir, f'run_{i}')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        random_state = i

        # Load and prepare data
        X, y, full_feature_table = load_and_prepare_data(input_path)
        X_train, X_test, y_train, y_test, X_test_sample_ids = filter_data(X, y, full_feature_table, filter_type)

        # Perform Recursive Feature Elimination (RFE)
        rfe_model, selected_features = perform_rfe(X_train, y_train, num_features, threads, output_dir)

        # Define parameter grid for grid search
        param_grid = {
            'iterations': [500, 1000],
            'learning_rate': [0.05, 0.1],
            'depth': [4, 6],
            'loss_function': ['Logloss'],
            'thread_count': [threads]
        }

        # Perform grid search for hyperparameter tuning
        best_model, best_params, best_mcc = grid_search(X_train, y_train, X_test, y_test, X_test_sample_ids, param_grid, output_dir)

        # Save feature importances and track feature occurrences
        feature_importances_path = os.path.join(output_dir, 'feature_importances.csv')
        save_feature_importances(best_model, pd.DataFrame(X_train, columns=selected_features), feature_importances_path)

        # Load feature importances and update occurrences
        features_df = pd.read_csv(feature_importances_path)
        for feature in features_df['Feature'].values:
            features_occurrence[feature] = features_occurrence.get(feature, 0) + 1

    # Save feature occurrences summary
    features_occurrence_df = pd.DataFrame(list(features_occurrence.items()), columns=['Feature', 'Occurrence'])
    features_occurrence_df.sort_values(by='Occurrence', ascending=False, inplace=True)
    features_occurrence_path = os.path.join(base_output_dir, 'features_occurrence.csv')
    features_occurrence_df.to_csv(features_occurrence_path, index=False)

    end_total_time = time.time()
    print(f"Feature selection iterations completed in {end_total_time - start_total_time:.2f} seconds.")

def generate_feature_tables(model_testing_dir, full_feature_table_file, filter_table_dir, cut_offs=[3, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50]):
    """
    Generate and save feature tables based on feature selection results from multiple runs in the main directory.
    
    Args:
        model_testing_dir (str): Directory containing feature selection runs.
        full_feature_table_file (str): Path to the full feature table CSV.
        filter_table_dir (str): Directory where filtered feature tables will be saved.
        cut_offs (list): List of thresholds for feature occurrences to be used for filtering.
    """
    # Load the full feature table
    full_feature_table = pd.read_csv(full_feature_table_file)
    interaction_count = full_feature_table.shape[0]
    print('Interaction count:', interaction_count)

    # Get all run directories in the main model testing directory
    run_dirs = [x for x in os.listdir(model_testing_dir) if 'run' in x]

    # Dictionary to store feature occurrences across all runs
    features_occurrence = {}

    # Iterate over each run directory
    for run in run_dirs:
        feature_importances_path = os.path.join(model_testing_dir, run, 'feature_importances.csv')
        if os.path.exists(feature_importances_path):
            features_df = pd.read_csv(feature_importances_path)
            
            # Update feature occurrence counts
            for feature in features_df['Feature'].values:
                features_occurrence[feature] = features_occurrence.get(feature, 0) + 1

    # Convert the feature occurrences into a DataFrame
    features_occurrence_df = pd.DataFrame(list(features_occurrence.items()), columns=['Feature', 'Occurrence'])
    features_occurrence_df.sort_values(by='Occurrence', ascending=False, inplace=True)

    # Create occurrence counts table for analysis (optional)
    occurrence_counts = features_occurrence_df['Occurrence'].value_counts().reset_index()
    occurrence_counts = occurrence_counts.rename(columns={'index': 'Occurrence', 'Occurrence': 'Count'})

    # Process each cut-off value
    for cut_off in cut_offs:
        features_occurrence_filter = features_occurrence_df[features_occurrence_df['Occurrence'] >= cut_off]

        # Ensure the filtered feature set is within reasonable bounds
        if 20 < len(features_occurrence_filter) < interaction_count / 20:
            print(f'Cut-off: {cut_off} - Features: {len(features_occurrence_filter)}')
            select_features = features_occurrence_filter['Feature'].tolist()

            # Select the relevant features from the full feature table
            select_feature_table = full_feature_table[['strain', 'phage', 'interaction'] + select_features]
            select_feature_table = select_feature_table.melt(id_vars=['strain', 'phage', 'interaction'], var_name='Feature', value_name='Value')
            select_feature_table['Value'] = select_feature_table['Value'].apply(lambda x: 1 if x > 0 else 0)
            select_feature_table = select_feature_table.pivot_table(index=['strain', 'phage', 'interaction'], columns='Feature', values='Value').reset_index()

            # Create directory for the filtered feature tables if it doesn't exist
            if not os.path.exists(filter_table_dir):
                os.makedirs(filter_table_dir)

            # Save the filtered feature table
            select_feature_table_path = os.path.join(filter_table_dir, f'select_feature_table_cutoff_{cut_off}.csv')
            print(f"Saving feature table to {select_feature_table_path}")
            select_feature_table.to_csv(select_feature_table_path, index=False)
