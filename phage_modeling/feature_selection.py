import os
import logging
import pandas as pd
import numpy as np
from sklearn.feature_selection import RFE, SelectKBest, f_classif, f_regression, chi2, SelectFromModel
from sklearn.linear_model import LogisticRegression, Lasso, LinearRegression
from sklearn.cluster import AgglomerativeClustering
import shap
from sklearn.model_selection import train_test_split
from catboost import CatBoostClassifier, CatBoostRegressor
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, confusion_matrix, roc_curve, auc, precision_recall_curve, average_precision_score, mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from tqdm import tqdm
from hdbscan import HDBSCAN
import time

# Function to load and prepare data
def load_and_prepare_data(input_path, sample_column=None, phenotype_column=None, filter_type='none'):
    """
    Loads the input feature table, drops unnecessary columns, and splits into features and target.

    Args:
        input_path (str): Path to the input CSV file containing the full feature table.
        sample_column (str): Optional name of the column to retain for sample identifiers.
        phenotype_column (str): Optional name of the column to retain for phenotype information.

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

    # Prepare the feature set and drop unnecessary columns
    drop_columns = ['strain', 'phage', 'interaction', 'header', 'contig_id', 'orf_ko', filter_type]
    
    # Ensure the sample and phenotype columns are retained if specified
    if sample_column:
        drop_columns.remove('strain')  # Keep 'strain' or replace with sample_column
        drop_columns.append(sample_column)  # Add custom sample column if provided
    
    if phenotype_column:
        drop_columns.remove('interaction')  # Keep 'interaction' or replace with phenotype_column
        drop_columns.append(phenotype_column)  # Add custom phenotype column if provided

    X = full_feature_table.drop(drop_columns, axis=1, errors='ignore')
    
    # Determine the target variable (default 'interaction' or custom phenotype_column)
    target_column = phenotype_column if phenotype_column else 'interaction'
    y = full_feature_table[target_column]

    print(f"Number of positive samples: {y.sum()}")
    print(f"Number of negative samples: {len(y) - y.sum()}")
    print("Data loaded and prepared, split into features and target.")
    
    return X, y, full_feature_table

# Function to filter the data based on strain or phage
def filter_data(
    X, y, 
    full_feature_table, 
    filter_type, 
    random_state=42, 
    sample_column='strain', 
    output_dir=None,
    use_clustering=False, 
    cluster_method='hdbscan',
    n_clusters=20,
    check_feature_presence=False,
    ensure_balanced_split=True,  # New parameter to enable/disable this check
    max_attempts=10,             # Maximum attempts to find a valid split
    **kwargs
):
    """
    Filters the data by strain or phage and splits into training and testing sets. 
    If use_clustering=True, it first clusters based on feature content before splitting.

    Args:
        X (DataFrame): Features.
        y (Series): Target variable.
        full_feature_table (DataFrame): The full feature table with metadata.
        filter_type (str): 'none', 'strain', 'phage' to determine how the data should be filtered.
        random_state (int): Seed for reproducibility.
        sample_column (str): Column to use as the sample identifier (default: 'strain').
        output_dir (str): Directory to store intermediate files.
        use_clustering (bool): Whether to apply clustering before splitting. Default is False.
        cluster_method (str): Clustering method to use when use_clustering=True. Options: 'hdbscan' or 'hierarchical'.
        n_clusters (int): Number of clusters for hierarchical clustering (default: 20).
        check_feature_presence (bool): If True, only include features present in both train and test sets.
        ensure_balanced_split (bool): If True, ensures both train and test sets contain at least one positive sample.
                                     Only applies when y contains only 0s and 1s.
        max_attempts (int): Maximum number of attempts to find a valid split when ensure_balanced_split is True.
        **kwargs: Additional parameters for the clustering method.

    Returns:
        X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids
    """
    # Check if y is binary (contains only 0s and 1s)
    is_binary = ((y == 0) | (y == 1)).all()
    
    # Only apply balanced split logic if y is binary and ensure_balanced_split is True
    apply_balanced_split = ensure_balanced_split and is_binary
    
    if ensure_balanced_split and not is_binary:
        logging.warning("ensure_balanced_split is set to True but y is not binary (0s and 1s only). Balanced split logic will be skipped.")
    
    # ---- 1️⃣ Standard Random Split if No Clustering ----
    if not use_clustering or filter_type == 'none':
        # If we need to ensure a valid split with positives in both sets
        if apply_balanced_split:
            for attempt in range(max_attempts):
                # Use a more diverse seed generation strategy
                current_seed = random_state * 1000 + attempt * 17
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.2, random_state=current_seed
                )
                
                # Check if both train and test sets have at least one positive
                if y_train.sum() > 0 and y_test.sum() > 0:
                    logging.info(f"Found valid split on attempt {attempt+1} with seed {current_seed}")
                    break
                
                if attempt == max_attempts - 1:
                    logging.warning(f"Failed to find split with positives in both sets after {max_attempts} attempts. Using last attempt.")
        else:
            # Regular splitting without the check
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.2, random_state=random_state
            )
            
        test_idx = X_test.index
        train_idx = X_train.index

        meta_columns = [sample_column]
        if 'phage' in full_feature_table.columns:
            meta_columns.append('phage')
        X_test_sample_ids = full_feature_table.loc[test_idx, meta_columns]
        X_train_sample_ids = full_feature_table.loc[train_idx, meta_columns]

        if check_feature_presence:
            # Find features present in both train and test sets
            train_features = (X_train > 0).any(axis=0)
            test_features = (X_test > 0).any(axis=0)
            common_features = train_features & test_features
            
            # Filter features
            X_train = X_train.loc[:, common_features]
            X_test = X_test.loc[:, common_features]
            
            logging.info(f"Original features: {len(train_features)}")
            logging.info(f"Features present in both sets: {common_features.sum()}")
            logging.info(f"Features removed: {len(train_features) - common_features.sum()}")

        return X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids

    # ---- 2️⃣ Validate filter_type Column ----
    if filter_type not in full_feature_table.columns:
        raise ValueError(f"Filter type '{filter_type}' must be a column in full_feature_table.")

    # ---- 3️⃣ Prepare Feature Table for Clustering ----
    feature_columns = [col for col in full_feature_table.columns if col.startswith(f"{filter_type[0]}c_")]
    filter_type_feature_table = full_feature_table[[filter_type] + feature_columns].drop_duplicates()

    if feature_columns == []:
        logging.warning(f"No feature columns found for clustering. Falling back to group-based split.")
        use_clustering = False
        return filter_data(X, y, full_feature_table, filter_type, random_state, sample_column, output_dir, 
                          use_clustering=False, ensure_balanced_split=ensure_balanced_split)

    # ---- 4️⃣ Apply Clustering Based on Selected Method ----
    if use_clustering:
        if cluster_method == 'hdbscan':
            logging.info(f"Clustering using HDBSCAN with parameters: {kwargs}")
            clusterer = HDBSCAN(**kwargs)
            cluster_labels = clusterer.fit_predict(filter_type_feature_table[feature_columns])

            print(f"Number of clusters: {len(np.unique(cluster_labels[cluster_labels != -1]))}")
            print(f"Number of noise points: {np.sum(cluster_labels == -1)}")
            
            # Count the number of clusters INCLUDING noise points labeled as -1
            num_clusters = len(np.unique(cluster_labels))
            
            # If fewer than 5 clusters, switch to hierarchical clustering
            if num_clusters < 5:
                logging.info(f"HDBSCAN found only {num_clusters} clusters (including noise). Switching to hierarchical clustering with 5 clusters.")
                cluster_method = 'hierarchical'
                n_clusters = 5
                clusterer = AgglomerativeClustering(n_clusters=n_clusters)
                cluster_labels = clusterer.fit_predict(filter_type_feature_table[feature_columns])
                print(f"Number of clusters after switching to hierarchical: {len(np.unique(cluster_labels))}")
            else:
                # Handle noise points (label -1) by giving them unique cluster IDs
                max_cluster_label = cluster_labels.max()
                for i, label in enumerate(cluster_labels):
                    if label == -1:  # If it's noise
                        max_cluster_label += 1
                        cluster_labels[i] = max_cluster_label

        elif cluster_method == 'hierarchical':
            logging.info(f"Clustering using Hierarchical Clustering with {n_clusters} clusters")
            clusterer = AgglomerativeClustering(n_clusters=n_clusters)
            cluster_labels = clusterer.fit_predict(filter_type_feature_table[feature_columns])
            
            print(f"Number of clusters: {len(np.unique(cluster_labels))}")

        else:
            raise ValueError(f"Unsupported clustering method: {cluster_method}. Choose 'hdbscan' or 'hierarchical'.")

        # Assign clusters to samples
        filter_type_feature_table["cluster"] = cluster_labels
        print(f"Number of clusters for splitting: {len(filter_type_feature_table['cluster'].unique())}")
        full_feature_table = full_feature_table.merge(
            filter_type_feature_table[[filter_type, "cluster"]], on=filter_type, how="left"
        )

        if output_dir:
            cluster_file = os.path.join(output_dir, f"{filter_type}_clusters.csv")
            cluster_df = filter_type_feature_table[[filter_type, "cluster"]].drop_duplicates()
            cluster_df.to_csv(cluster_file, index=False)
            logging.info(f"Saved cluster labels to: {cluster_file}")

        group_col = "cluster"  # Use cluster labels for splitting
    else:
        group_col = filter_type  # Use original filter type column

    # ---- 5️⃣ Perform Group-Based Splitting (Clustering or Normal) ----
    groups = full_feature_table[group_col].unique()

    if apply_balanced_split:
        # Track which groups contain positive samples
        group_positive = {}
        for group in groups:
            group_mask = full_feature_table[group_col] == group
            group_indices = full_feature_table.index[group_mask]
            if len(group_indices) > 0:  # Check if the group has any samples
                group_y = y.loc[group_indices]
                group_positive[group] = group_y.sum() > 0
            else:
                group_positive[group] = False
        
        # Find groups with positives
        positive_groups = [g for g, has_pos in group_positive.items() if has_pos]
        
        if len(positive_groups) < 2:
            logging.warning(f"Not enough groups with positive samples (found {len(positive_groups)}). Cannot ensure balanced split.")
            # Fall back to random splitting
            np.random.seed(random_state)
            train_groups = np.random.choice(groups, size=int(0.8 * len(groups)), replace=False)
            test_groups = np.setdiff1d(groups, train_groups)
        else:
            # Try multiple random seeds to find a valid split
            found_valid_split = False
            for attempt in range(max_attempts):
                # Use a more diverse seed generation strategy
                current_seed = random_state * 1000 + attempt * 17
                np.random.seed(current_seed)
                
                # Ensure at least one positive group in train and test
                np.random.shuffle(positive_groups)
                pos_train_groups = positive_groups[:max(1, int(0.8 * len(positive_groups)))]
                pos_test_groups = np.setdiff1d(positive_groups, pos_train_groups)
                
                # Fill remaining train groups with other groups
                remaining_groups = np.setdiff1d(groups, positive_groups)
                if len(remaining_groups) > 0:
                    other_train_count = int(0.8 * len(groups)) - len(pos_train_groups)
                    if other_train_count > 0:
                        other_train_groups = np.random.choice(
                            remaining_groups, 
                            size=min(other_train_count, len(remaining_groups)), 
                            replace=False
                        )
                        train_groups = np.concatenate([pos_train_groups, other_train_groups])
                    else:
                        train_groups = pos_train_groups
                else:
                    train_groups = pos_train_groups
                
                test_groups = np.setdiff1d(groups, train_groups)
                
                # Check that both sets will have positive samples
                train_idx = full_feature_table[group_col].isin(train_groups)
                test_idx = full_feature_table[group_col].isin(test_groups)
                
                train_pos = y[train_idx].sum() > 0
                test_pos = y[test_idx].sum() > 0
                
                if train_pos and test_pos:
                    found_valid_split = True
                    logging.info(f"Found valid group-based split on attempt {attempt+1} with seed {current_seed}")
                    break
            
            if not found_valid_split:
                logging.warning(f"Failed to find a valid group-based split after {max_attempts} attempts. Using last attempt.")
    else:
        # Standard group-based split without ensuring positive samples
        np.random.seed(random_state)
        train_groups = np.random.choice(groups, size=int(0.8 * len(groups)), replace=False)
        test_groups = np.setdiff1d(groups, train_groups)

    train_idx = full_feature_table[group_col].isin(train_groups)
    test_idx = full_feature_table[group_col].isin(test_groups)

    X_train = X[train_idx]
    X_test = X[test_idx]
    y_train = y[train_idx]
    y_test = y[test_idx]

    # Print distribution information for debugging
    if is_binary:
        logging.info(f"Train set positive samples: {y_train.sum()} out of {len(y_train)}")
        logging.info(f"Test set positive samples: {y_test.sum()} out of {len(y_test)}")
    else:
        logging.info(f"Train set size: {len(y_train)}, Test set size: {len(y_test)}")
        logging.info(f"Train set distribution: min={y_train.min()}, mean={y_train.mean():.2f}, max={y_train.max()}")
        logging.info(f"Test set distribution: min={y_test.min()}, mean={y_test.mean():.2f}, max={y_test.max()}")

    # ---- 6️⃣ Check Feature Presence if Requested ----
    if check_feature_presence:
        # Find features present in both train and test sets
        train_features = (X_train > 0).any(axis=0)
        test_features = (X_test > 0).any(axis=0)
        common_features = train_features & test_features
        
        # Filter features
        X_train = X_train.loc[:, common_features]
        X_test = X_test.loc[:, common_features]
        
        logging.info(f"Original features: {len(train_features)}")
        logging.info(f"Features present in both sets: {common_features.sum()}")
        logging.info(f"Features removed: {len(train_features) - common_features.sum()}")

    # ---- 7️⃣ Prepare Metadata for Train/Test Sets ----
    meta_columns = [sample_column]
    if 'phage' in full_feature_table.columns:
        meta_columns.append('phage')
    X_test_sample_ids = full_feature_table.loc[test_idx, meta_columns]
    X_train_sample_ids = full_feature_table.loc[train_idx, meta_columns]

    # ---- 8️⃣ Check for Valid Training Set ----
    unique_values = y_train.nunique()
    if unique_values < 2:
        logging.warning(
            f"Training set contains only one unique target value ({y_train.unique()[0]}). "
            f"Skipping this split."
        )
        return None, None, None, None, None, None

    return X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids

def compute_phage_weights(X_train_sample_ids, y_train, phage_column, method='log10', smoothing=1.0):
    """
    Compute class weights for each phage using only X_train_sample_ids and y_train based on the specified method.
    Args:
      - X_train_sample_ids (DataFrame): Contains sample identifiers, including `phage_column`.
      - y_train (Series): Training target labels (0 or 1).
      - phage_column (str): Column name representing the phage.
      - method (str): Method to calculate weights. Options are 'log10' (default), 'inverse_frequency', or 'balanced'.
      - smoothing (float): Small constant to prevent extreme weight ratios.
    Returns:
      Dict[phage] = {0: weight_for_neg, 1: weight_for_pos}
    """
    if phage_column not in X_train_sample_ids.columns:
        raise ValueError(f"Column '{phage_column}' not found in X_train_sample_ids.")

    # Merge y_train with X_train_sample_ids to ensure alignment
    phage_data = X_train_sample_ids.copy()
    phage_data["interaction"] = y_train.values  # Ensure y_train aligns with indices

    # Compute positive and negative counts per phage
    phage_counts = phage_data.groupby(phage_column)["interaction"].agg(["sum", "count"])
    phage_counts.rename(columns={"sum": "pos_count", "count": "total_count"}, inplace=True)
    phage_counts["neg_count"] = phage_counts["total_count"] - phage_counts["pos_count"]

    # Compute weights based on the specified method
    phage_weights = {}
    for phage, row in phage_counts.iterrows():
        pos_count = row["pos_count"]
        neg_count = row["neg_count"]
        if pos_count == 0:
            phage_weights[phage] = {0: 1.0, 1: smoothing}  # Avoid infinite scaling
        else:
            if method == 'log10':
                # Log10 method (default, similar to log1p but using base 10)
                weight_1 = max(1.0, np.log10(neg_count / (pos_count + smoothing)) + 1)
            elif method == 'inverse_frequency':
                # Inverse frequency method
                weight_1 = (row["total_count"] / (pos_count + smoothing))
            elif method == 'balanced':
                # Balanced method
                weight_1 = (row["total_count"] - pos_count) / (pos_count + smoothing)
            else:
                raise ValueError(f"Unsupported method: {method}. Choose 'log10', 'inverse_frequency', or 'balanced'.")
            phage_weights[phage] = {0: 1.0, 1: weight_1}

    return phage_weights

def build_row_weights(y_train, X_train_sample_ids, phage_weights, phage_column):
    """
    Assign sample weights based on phage weights using X_train_sample_ids.

    Args:
      - y_train (Series): Training target labels (0 or 1).
      - X_train_sample_ids (DataFrame): Contains sample identifiers, including `phage_column`.
      - phage_weights (dict): Precomputed weights per phage.
      - phage_column (str): Column name representing phage ID.

    Returns:
      - sample_weights (numpy array): Array of per-row sample weights.
    """
    if phage_column not in X_train_sample_ids.columns:
        raise ValueError(f"Column '{phage_column}' not found in X_train_sample_ids.")

    # Map phage_column from X_train_sample_ids
    phage_series = X_train_sample_ids[phage_column]

    # Initialize weight array
    sample_weights = np.zeros(len(y_train), dtype=float)

    for i, idx in enumerate(y_train.index):
        row_phage = phage_series.loc[idx]  # Get phage ID
        row_label = y_train.loc[idx]  # Get label
        sample_weights[i] = phage_weights.get(row_phage, {0: 1.0, 1: 1.0}).get(row_label, 1.0)

    return sample_weights

# Function to perform Recursive Feature Elimination (RFE)
def perform_rfe(
    X_train, 
    y_train, 
    X_train_sample_ids,
    num_features, 
    threads, 
    output_dir, 
    task_type='classification', 
    phage_column='phage',
    use_dynamic_weights=False,
    weights_method='log10',
    max_ram=8
):
    """
    Performs Recursive Feature Elimination (RFE) to select the top features.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Number of features to select.
        threads (int): Number of threads to use for CatBoost.
        output_dir (str): Directory to store intermediate CatBoost information.
        task_type (str): Task type for model ('classification' or 'regression').

    Returns:
        rfe (RFE object): Fitted RFE model.
        selected_features (Index): List of selected features.
    """
    total_features = X_train.shape[1]
    step_size = max(1, int((total_features - num_features) / 10))  # Ensure step_size is at least 1

    # Initialize model based on task type
    if task_type == 'classification':
        model = CatBoostClassifier(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=10,
            thread_count=threads,
            train_dir=os.path.join(output_dir, '..', 'catboost_info'),
            used_ram_limit=f"{max_ram}gb"  # Set the RAM limit
        )
    elif task_type == 'regression':
        model = CatBoostRegressor(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=10,
            thread_count=threads,
            train_dir=os.path.join(output_dir, '..', 'catboost_info'),
            used_ram_limit=f"{max_ram}gb"  # Set the RAM limit
        )
    else:
        raise ValueError("task_type must be 'classification' or 'regression'")

    print(f"Performing Recursive Feature Elimination (RFE) with step_size: {step_size} for {task_type}...")

    # If requested, build sample weights for each row
    sample_weight = None
    if use_dynamic_weights and phage_column and (phage_column in X_train_sample_ids.columns):
        print('Length X_train_sample_ids: ', len(X_train_sample_ids))
        print('Length y_train: ', len(y_train))
        phage_weights = compute_phage_weights(X_train_sample_ids, y_train, phage_column, method=weights_method)
        sample_weight = build_row_weights(y_train, X_train_sample_ids, phage_weights, phage_column)
        print(phage_weights)
        phage_weights_df = pd.DataFrame(phage_weights).T
        phage_weights_df.to_csv(f"{output_dir}/phage_weights.csv")
        logging.info("Using dynamic phage-based sample weights during RFE.")
    else:
        print("RFE without dynamic sample weights.")

    # Set up and fit RFE
    rfe = RFE(estimator=model, n_features_to_select=num_features, step=step_size, verbose=10)
    rfe.fit(X_train, y_train, sample_weight=sample_weight)
    
    selected_features = X_train.columns[rfe.support_]
    print(f"RFE selected {len(selected_features)} features.")
    
    return rfe, selected_features

def shap_rfe(X_train, y_train, num_features, threads, task_type='classification', max_ram=8):
    """
    Performs Recursive Feature Elimination (RFE) based on SHAP feature importances.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Desired number of features to select.
        threads (int): Number of threads to use for CatBoost training.
        task_type (str): Task type for model ('classification' or 'regression').

    Returns:
        X_train_selected (DataFrame): Training features with the selected top features.
        selected_features (Index): List of selected feature names.
    """
    total_features = X_train.shape[1]
    step_size = max(1, int((total_features - num_features) / 10))  # Ensure step_size is at least 1
    current_features = X_train.columns.tolist()

    # Select the model based on task type
    if task_type == 'classification':
        model = CatBoostClassifier(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=0,
            thread_count=threads,
            used_ram_limit=f"{max_ram}gb"  # Set the RAM limi
        )
    elif task_type == 'regression':
        model = CatBoostRegressor(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=0,
            thread_count=threads,
            used_ram_limit=f"{max_ram}gb"  # Set the RAM limi
        )
    else:
        raise ValueError("task_type must be 'classification' or 'regression'")
    
    # Perform SHAP-based RFE
    while len(current_features) > num_features:
        # Train the model with current features
        model.fit(X_train[current_features], y_train)

        # Calculate SHAP values
        explainer = shap.TreeExplainer(model, approximate=True)
        shap_values = explainer.shap_values(X_train[current_features])
        
        # Calculate mean absolute SHAP values for each feature
        shap_importances = np.abs(shap_values).mean(axis=0)
        feature_importances_df = pd.DataFrame({
            'Feature': current_features,
            'Importance': shap_importances
        })
        
        # Sort features by SHAP importance
        feature_importances_df = feature_importances_df.sort_values(by='Importance', ascending=False)
        
        # Dynamically adjust step size if close to target
        remaining_features = len(current_features)
        step_size = min(step_size, remaining_features - num_features)
        
        # Remove the bottom features according to adjusted step size
        to_remove = feature_importances_df.tail(step_size)['Feature'].tolist()
        current_features = [f for f in current_features if f not in to_remove]
        
        print(f"Removed {len(to_remove)} features. Remaining features: {len(current_features)}")
    
    print(f"SHAP-RFE selected {len(current_features)} features.")
    
    # Return selected features and transformed X_train
    X_train_selected = X_train[current_features]
    
    return X_train_selected, current_features

def select_k_best_feature_selection(X_train, y_train, num_features, task_type='classification'):
    """
    Selects the top features using the SelectKBest method with ANOVA F-test for classification 
    or f_regression for regression.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Number of top features to select.
        task_type (str): Task type ('classification' or 'regression').

    Returns:
        X_train_selected (DataFrame): Training features with the selected top features.
        selected_features (Index): List of selected feature names.
    """
    print("Selecting features using SelectKBest...")

    if task_type == 'classification':
        skb = SelectKBest(score_func=f_classif, k=num_features)
    elif task_type == 'regression':
        skb = SelectKBest(score_func=f_regression, k=num_features)
    else:
        raise ValueError("task_type must be 'classification' or 'regression'")

    skb.fit(X_train, y_train)

    # Get the boolean mask of selected features
    support = skb.get_support()
    selected_features = X_train.columns[support]

    # Transform X_train to a DataFrame with selected features' column names
    X_train_selected = pd.DataFrame(skb.transform(X_train), columns=selected_features, index=X_train.index)
    
    print(f"SelectKBest selected {len(selected_features)} features.")
    return X_train_selected, selected_features


def chi_squared_feature_selection(X_train, y_train, num_features):
    """
    Selects top features using the Chi-Squared Test.

    Args:
        X_train (DataFrame): Training features (must be non-negative).
        y_train (Series): Training target.
        num_features (int): Number of top features to select.

    Returns:
        X_train_selected (DataFrame): Transformed training features with selected features.
        selected_features (Index): List of selected feature names.
    """
    print("Selecting features using Chi-Squared Test...")
    chi2_selector = SelectKBest(score_func=chi2, k=num_features)
    chi2_selector.fit(X_train, y_train)
    
    # Get the boolean mask of selected features
    support = chi2_selector.get_support()
    selected_features = X_train.columns[support]
    
    # Transform X_train to a DataFrame with selected features' column names
    X_train_selected = pd.DataFrame(chi2_selector.transform(X_train), columns=selected_features, index=X_train.index)
    
    print(f"Chi-Squared Test selected {len(selected_features)} features.")
    return X_train_selected, selected_features

def lasso_feature_selection(X_train, y_train, num_features, task_type='classification'):
    """
    Selects top features using Lasso regularization for classification or regression.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Number of top features to select.
        task_type (str): Task type ('classification' or 'regression').

    Returns:
        X_train_selected (DataFrame): Transformed training features with selected features.
        selected_features (Index): List of selected feature names.
    """
    print("Selecting features using Lasso regularization...")

    # Choose model based on task type
    if task_type == 'classification':
        model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000)
    elif task_type == 'regression':
        model = Lasso(max_iter=1000)
    else:
        raise ValueError("task_type must be 'classification' or 'regression'")
    
    # Fit the model and select features
    model.fit(X_train, y_train)
    selector = SelectFromModel(model, max_features=num_features, prefit=True)
    support = selector.get_support()
    selected_features = X_train.columns[support]
    
    # Transform X_train to a DataFrame with selected features' column names
    X_train_selected = pd.DataFrame(selector.transform(X_train), columns=selected_features, index=X_train.index)
    
    print(f"Lasso selected {len(selected_features)} features.")
    return X_train_selected, selected_features

def shap_feature_selection(X_train, y_train, num_features, threads, task_type='classification', max_ram=8):
    """
    Selects top features based on SHAP values for classification or regression.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Number of top features to select.
        threads (int): Number of threads to use for CatBoost training.
        task_type (str): Task type ('classification' or 'regression').

    Returns:
        X_train_selected (DataFrame): Transformed training features with selected features.
        selected_features (Index): List of selected feature names.
    """
    # Choose the model based on the task type
    if task_type == 'classification':
        model = CatBoostClassifier(iterations=500, learning_rate=0.1, depth=4, verbose=0, thread_count=threads, used_ram_limit=f"{max_ram}gb")
    elif task_type == 'regression':
        model = CatBoostRegressor(iterations=500, learning_rate=0.1, depth=4, verbose=0, thread_count=threads, used_ram_limit=f"{max_ram}gb")
    else:
        raise ValueError("task_type must be 'classification' or 'regression'")
    
    # Train the model
    model.fit(X_train, y_train)
    
    # Calculate SHAP values
    explainer = shap.TreeExplainer(model, approximate=True)
    shap_values = explainer.shap_values(X_train)
    
    # Calculate mean absolute SHAP values for each feature
    shap_importances = np.abs(shap_values).mean(axis=0)
    top_indices = np.argsort(shap_importances)[-num_features:]
    selected_features = X_train.columns[top_indices]
    
    print(f"SHAP selected {len(selected_features)} features.")
    return X_train[selected_features], selected_features

def train_and_evaluate(X_train, 
                       y_train, 
                       X_test, 
                       y_test,
                       X_train_sample_ids,
                       params, 
                       output_dir, 
                       max_ram=8,
                       phage_column=None,
                       use_dynamic_weights=False,
                       weights_method='log10'):
    """
    Train a CatBoost model and evaluate it on the test set.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        X_test (DataFrame): Test features.
        y_test (Series): Test target.
        params (dict): Model hyperparameters.
        output_dir (str): Directory to save evaluation results.

    Returns:
        model: Trained CatBoost model.
        accuracy (float): Accuracy on the test set.
        f1 (float): F1 score on the test set.
        mcc (float): Matthews Correlation Coefficient on the test set.
        y_pred (array): Predictions on the test set.
    """
    # Setting up CatBoost's training directory
    train_dir = os.path.join(output_dir, '..', 'catboost_info')
    model = CatBoostClassifier(**params, train_dir=train_dir, used_ram_limit=f"{max_ram}gb")

    print(f"Training with parameters: {params}")

    # Build sample weights if requested
    sample_weight = None
    if use_dynamic_weights and phage_column and (phage_column in X_train_sample_ids.columns):
        print('Length X_train_sample_ids: ', len(X_train_sample_ids))
        print('Length y_train: ', len(y_train))
        phage_weights = compute_phage_weights(X_train_sample_ids, y_train, phage_column, method=weights_method)
        sample_weight = build_row_weights(y_train, X_train_sample_ids, phage_weights, phage_column)
        print(phage_weights)
        phage_weights_df = pd.DataFrame(phage_weights).T
        phage_weights_df.to_csv(f"{output_dir}/phage_weights.csv")
        print("Using dynamic phage-based class weights.")
    else:
        print("Training without dynamic class weights.")
    
    # Training the model with early stopping
    model.fit(
        X_train, 
        y_train,
        sample_weight=sample_weight,
        eval_set=(X_test, y_test), 
        plot=False, 
        verbose=10, 
        early_stopping_rounds=100
    )

    # Make predictions
    y_pred = model.predict(X_test)

    # Calculate evaluation metrics
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)

    print(f"Training completed. Params: {params}, Accuracy: {accuracy}, F1 Score: {f1}, MCC: {mcc}")

    return model, accuracy, f1, mcc, y_pred

def train_and_evaluate_regressor(X_train, y_train, X_test, y_test, params, output_dir, max_ram=8):
    """
    Train a CatBoost regressor and evaluate it on the test set.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        X_test (DataFrame): Test features.
        y_test (Series): Test target.
        params (dict): Model hyperparameters.
        output_dir (str): Directory to save evaluation results.

    Returns:
        model: Trained CatBoost regressor model.
        mse (float): Mean Squared Error on the test set.
        r2 (float): R2 score on the test set.
        y_pred (array): Predictions on the test set.
    """
    train_dir = os.path.join(output_dir, '..', 'catboost_info')
    model = CatBoostRegressor(**params, train_dir=train_dir, used_ram_limit=f"{max_ram}gb")

    print(f"Training regressor with parameters: {params}")
    
    # Train with early stopping
    model.fit(X_train, y_train, eval_set=(X_test, y_test), plot=False, verbose=10, early_stopping_rounds=100)

    # Make predictions
    y_pred = model.predict(X_test)

    # Calculate evaluation metrics
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    print(f"Training completed. Params: {params}, MSE: {mse}, R2: {r2}")

    return model, mse, r2, y_pred

# Function to perform grid search
def grid_search(
    X_train, 
    y_train, 
    X_test, 
    y_test, 
    X_train_sample_ids, 
    X_test_sample_ids,
    param_grid, 
    output_dir, 
    phenotype_column='interaction', 
    phage_column='phage', 
    use_dynamic_weights=False,
    weights_method='log10',
    max_ram=8
):
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
        phenotype_column (str): Column name for the interaction or target variable.

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
        model, accuracy, f1, mcc, y_pred = train_and_evaluate(
            X_train, 
            y_train, 
            X_test, 
            y_test, 
            X_train_sample_ids,
            params, 
            output_dir, 
            max_ram=max_ram,
            phage_column=phage_column,
            use_dynamic_weights=use_dynamic_weights,
            weights_method=weights_method
        )
        
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
            best_predictions_df[phenotype_column] = y_test
            best_predictions_df.to_csv(f"{output_dir}/best_model_predictions.csv", index=False)

    pd.DataFrame(results).to_csv(f"{output_dir}/model_performance.csv", index=False)

    # Return None if no model is found
    if best_model is None:
        logging.warning("No valid model was found in grid search.")
    
    return best_model, best_params, best_mcc

def grid_search_regressor(X_train, y_train, X_test, y_test, X_test_sample_ids, param_grid, output_dir, phenotype_column='interaction', max_ram=8):
    """
    Performs grid search to find the best hyperparameters for CatBoost regression.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        X_test (DataFrame): Testing features.
        y_test (Series): Testing target.
        X_test_sample_ids (DataFrame): Metadata for the test set samples.
        param_grid (dict): Dictionary of hyperparameters for grid search.
        output_dir (str): Directory to save results.
        phenotype_column (str): Column name for the interaction or target variable.

    Returns:
        best_model (CatBoostRegressor): The model with the best performance.
        best_params (dict): The hyperparameters of the best model.
        best_r2 (float): The highest R2 score achieved during grid search.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    best_r2 = -float('inf')
    best_model = None
    best_params = None
    results = []
    
    print("Starting grid search for regression...")
    for idx, params in enumerate(itertools.product(*param_grid.values()), start=1):
        params = dict(zip(param_grid.keys(), params))
        model, mse, r2, y_pred = train_and_evaluate_regressor(X_train, y_train, X_test, y_test, params, output_dir, max_ram=max_ram)
        
        # Save performance results for this iteration
        results.append({**params, 'mse': mse, 'r2': r2})
        
        # Generate plots for each parameter set
        plot_dir = os.path.join(output_dir, f"plots_{idx}")
        os.makedirs(plot_dir, exist_ok=True)
        
        plot_predicted_vs_actual(y_test, y_pred, os.path.join(plot_dir, f"predicted_vs_actual_{idx}.png"))
        plot_residuals(y_test, y_pred, os.path.join(plot_dir, f"residuals_{idx}.png"))
        plot_residual_distribution(y_test, y_pred, os.path.join(plot_dir, f"residual_distribution_{idx}.png"))
        
        # Check if current model has the best R2 score
        if r2 >= best_r2:
            best_r2 = r2
            best_model = model
            best_params = params

            best_predictions_df = X_test_sample_ids.copy()
            best_predictions_df['Prediction'] = y_pred
            best_predictions_df[phenotype_column] = y_test
            best_predictions_df.to_csv(f"{output_dir}/best_model_predictions.csv", index=False)

    # Save the results DataFrame
    pd.DataFrame(results).to_csv(f"{output_dir}/model_performance.csv", index=False)

    # Return None if no model is found
    if best_model is None:
        logging.warning("No valid model was found in grid search.")
    
    return best_model, best_params, best_r2

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

def plot_predicted_vs_actual(y_test, y_pred, output_path):
    """
    Plots and saves a Predicted vs Actual values plot for regression, with a dashed 1-1 line and a linear regression line.

    Args:
        y_test (Series): True values for the test set.
        y_pred (Series): Predicted values for the test set.
        output_path (str): Path to save the Predicted vs Actual plot.
    """
    plt.figure(figsize=(8, 6))
    
    # Scatter plot of actual vs predicted values
    plt.scatter(y_test, y_pred, alpha=0.6, label='Predicted vs Actual')
    
    # Plot the 1-1 line
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2, color='red', label='1-1 Line')
    
    # Fit and plot the linear regression line
    model = LinearRegression().fit(np.array(y_test).reshape(-1, 1), y_pred)
    y_pred_line = model.predict(np.array(y_test).reshape(-1, 1))
    plt.plot(y_test, y_pred_line, color='blue', lw=2, label='Regression Line')
    
    # Labels and title
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title('Predicted vs Actual Values')
    plt.legend()
    
    # Save the plot
    plt.savefig(output_path)
    plt.close()
    print(f"Predicted vs Actual plot saved to {output_path}")

def plot_residuals(y_test, y_pred, output_path):
    """
    Plots and saves a residuals plot for regression.

    Args:
        y_test (Series): True values for the test set.
        y_pred (Series): Predicted values for the test set.
        output_path (str): Path to save the residuals plot.
    """
    residuals = y_test - y_pred
    plt.figure(figsize=(8, 6))
    plt.scatter(y_pred, residuals, alpha=0.6)
    plt.axhline(0, color='red', linestyle='--', linewidth=2)
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Residuals Plot')
    plt.savefig(output_path)
    plt.close()
    print(f"Residual plot saved to {output_path}")

def plot_residual_distribution(y_test, y_pred, output_path):
    """
    Plots and saves a distribution plot for residuals.

    Args:
        y_test (Series): True values for the test set.
        y_pred (Series): Predicted values for the test set.
        output_path (str): Path to save the residual distribution plot.
    """
    residuals = y_test - y_pred
    plt.figure(figsize=(8, 6))
    sns.histplot(residuals, kde=True, color='blue', bins=30)
    plt.axvline(0, color='red', linestyle='--', linewidth=2)
    plt.xlabel('Residuals')
    plt.ylabel('Frequency')
    plt.title('Residuals Distribution')
    plt.savefig(output_path)
    plt.close()
    print(f"Residual distribution plot saved to {output_path}")

def save_feature_importances(best_model, selected_features, feature_importances_path):
    """
    Saves feature importances from the model into a CSV file.
    
    Args:
        best_model: Trained model with feature importances.
        selected_features (DataFrame): DataFrame containing selected features used for training.
        feature_importances_path (str): Path to save the feature importances CSV.
    """
    if not hasattr(best_model, "feature_importances_"):
        logging.warning("No feature importances found on the model. Skipping save.")
        return
        
    feature_importances = best_model.feature_importances_

    # Ensure the lengths match
    if len(selected_features.columns) != len(feature_importances):
        logging.error("Mismatch between the number of selected features and the number of feature importances.")
        logging.info(f"Number of selected features: {len(selected_features.columns)}")
        logging.info(f"Number of feature importances: {len(feature_importances)}")
        return
    
    importance_df = pd.DataFrame({
        'Feature': selected_features.columns,
        'Importance': feature_importances
    })

    importance_df = importance_df.sort_values(by='Importance', ascending=False)

    importance_df.to_csv(feature_importances_path, index=False)
    logging.info(f"Feature importances saved to {feature_importances_path}")

def run_feature_selection_iterations(
    input_path, 
    base_output_dir, 
    threads, 
    num_features, 
    filter_type, 
    num_runs, 
    select_cols=False, 
    sample_column='strain', 
    phage_column='phage', 
    phenotype_column=None,
    use_dynamic_weights=False,
    weights_method='log10',
    method='rfe', 
    task_type='classification', 
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    max_ram=8
):
    """
    Runs multiple iterations of feature selection, saves the results in `run_*` directories, and tracks feature occurrences.
    
    Args:
        input_path (str): Path to the input feature table.
        base_output_dir (str): Base output directory where results for each run will be stored.
        threads (int): Number of threads to use for feature selection.
        num_features (int): Number of features to select.
        filter_type (str): Filter type for the input data ('strain', 'phage', 'none').
        num_runs (int): Number of runs to perform.
        select_cols (bool): Whether to run with selected columns.
        sample_column (str): Column name for the sample/strain (if using selected columns).
        use_dynamic_weights: Whether to use dynamic weights for phage-based samples.
        phenotype_column (str): Column name for the phenotype (if using selected columns).
        method (str): Feature selection method ('rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').
        task_type (str): Task type ('classification' or 'regression').
        use_clustering (bool): Whether to use clustering for filtering.
        cluster_method (str): Clustering method to use.
        n_clusters (int): Number of clusters for clustering.
        min_cluster_size (int): Minimum cluster size for filtering.
        min_samples (int): Minimum number of samples for filtering (default: None for same as min_cluster_size).
        cluster_selection_epsilon (float): Epsilon value for clustering.
        max_ram (int): Maximum RAM to use for CatBoost training.
    """
    
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    features_occurrence = {}
    start_total_time = time.time()

    for i in tqdm(range(num_runs), desc="Running Feature Selection Iterations"):
        output_dir = os.path.join(base_output_dir, f'run_{i}')
        feature_importances_path = os.path.join(output_dir, 'feature_importances.csv')

        if not os.path.exists(feature_importances_path):
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            random_state = i

            X, y, full_feature_table = load_and_prepare_data(input_path, sample_column=sample_column, phenotype_column=phenotype_column, filter_type=filter_type)
            X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids = filter_data(
                X, y, 
                full_feature_table, 
                filter_type, 
                random_state=random_state, 
                sample_column=sample_column, 
                output_dir=output_dir,
                use_clustering=use_clustering, 
                cluster_method=cluster_method,
                n_clusters=n_clusters,
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                cluster_selection_epsilon=cluster_selection_epsilon,
                check_feature_presence=check_feature_presence
            )

            if num_features == 'none':
                num_features = int(len(X_train) / 10) if len(X_train) < 500 else int(len(X_train) / 20)
            else:
                num_features = int(num_features)

            if X_train is None:
                logging.info("Skipping this run due to insufficient training data.")
                continue  # Skip this iteration and proceed to the next

            # Apply selected feature selection method
            if method == 'rfe':
                _, selected_features = perform_rfe(X_train, y_train, X_train_sample_ids, num_features, threads, output_dir, task_type=task_type, phage_column=phage_column, use_dynamic_weights=use_dynamic_weights, max_ram=max_ram)
            elif method == 'shap_rfe':
                X_train, selected_features = shap_rfe(X_train, y_train, num_features, threads, task_type=task_type, max_ram=max_ram)
            elif method == 'select_k_best':
                X_train, selected_features = select_k_best_feature_selection(X_train, y_train, num_features, task_type=task_type)
            elif method == 'chi_squared' and task_type == 'classification':
                X_train, selected_features = chi_squared_feature_selection(X_train, y_train, num_features)
            elif method == 'lasso':
                X_train, selected_features = lasso_feature_selection(X_train, y_train, num_features, task_type=task_type)
            elif method == 'shap':
                X_train, selected_features = shap_feature_selection(X_train, y_train, num_features, threads, task_type=task_type, max_ram=max_ram)
            else:
                raise ValueError(f"Unsupported feature selection method: {method} or incompatible task_type.")

            X_train_selected = X_train[selected_features]
            X_test_selected = X_test[selected_features]

            param_grid = {
                'iterations': [500, 1000],
                'learning_rate': [0.05, 0.1],
                'depth': [4, 6],
                'thread_count': [threads]
            }

            if task_type == 'classification':
                param_grid['loss_function'] = ['Logloss']
                best_model, best_params, best_mcc = grid_search(
                    X_train_selected, 
                    y_train, 
                    X_test_selected, 
                    y_test, 
                    X_train_sample_ids,
                    X_test_sample_ids, 
                    param_grid, 
                    output_dir, 
                    phenotype_column=phenotype_column,
                    phage_column=phage_column,
                    use_dynamic_weights=use_dynamic_weights,
                    weights_method=weights_method,
                    max_ram=max_ram)
                best_metric = best_mcc
            elif task_type == 'regression':
                param_grid['loss_function'] = ['RMSE']
                best_model, best_params, best_r2 = grid_search_regressor(X_train_selected, y_train, X_test_selected, y_test, X_test_sample_ids, param_grid, output_dir, phenotype_column=phenotype_column, max_ram=max_ram)
                best_metric = best_r2
            else:
                raise ValueError("task_type must be 'classification' or 'regression'")

            if best_model is None:
                logging.warning(f"No best model found for iteration {i}. Skipping feature importance saving.")
                continue

            # Save feature importances
            save_feature_importances(best_model, pd.DataFrame(X_train_selected, columns=selected_features), feature_importances_path)

            features_df = pd.read_csv(feature_importances_path)
            for feature in features_df['Feature'].values:
                features_occurrence[feature] = features_occurrence.get(feature, 0) + 1
        else:
            print(f"Feature importances already exist for run {i}. Skipping feature selection.")

    features_occurrence_df = pd.DataFrame(list(features_occurrence.items()), columns=['Feature', 'Occurrence'])
    features_occurrence_df.sort_values(by='Occurrence', ascending=False, inplace=True)
    features_occurrence_path = os.path.join(base_output_dir, 'features_occurrence.csv')
    features_occurrence_df.to_csv(features_occurrence_path, index=False)

    end_total_time = time.time()
    print(f"Feature selection iterations completed in {end_total_time - start_total_time:.2f} seconds.")

def generate_feature_tables(
    model_testing_dir, full_feature_table_file, filter_table_dir, 
    phenotype_column=None, sample_column='strain', cut_offs=[3, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    binary_data=False, max_features=None, filter_type='strain'
):
    """
    Generate and save feature tables based on feature selection results from multiple runs in the main directory.
    
    Args:
        model_testing_dir (str): Directory containing feature selection runs.
        full_feature_table_file (str): Path to the full feature table CSV.
        filter_table_dir (str): Directory where filtered feature tables will be saved.
        phenotype_column (str): Column name for the target variable.
        sample_column (str): Column name for the sample or strain identifier.
        cut_offs (list): List of thresholds for feature occurrences to be used for filtering.
        binary_data (bool): If True, converts feature values to binary (1/0). Default is False for continuous values.
    """
    if phenotype_column is None:
        phenotype_column = 'interaction'

    full_feature_table = pd.read_csv(full_feature_table_file)
    interaction_count = full_feature_table.shape[0]
    print('Interaction count:', interaction_count)

    run_dirs = [x for x in os.listdir(model_testing_dir) if 'run' in x]
    features_occurrence = {}

    for run in run_dirs:
        feature_importances_path = os.path.join(model_testing_dir, run, 'feature_importances.csv')
        if os.path.exists(feature_importances_path):
            features_df = pd.read_csv(feature_importances_path)
            for feature in features_df['Feature'].values:
                features_occurrence[feature] = features_occurrence.get(feature, 0) + 1

    features_occurrence_df = pd.DataFrame(list(features_occurrence.items()), columns=['Feature', 'Occurrence'])
    features_occurrence_df.sort_values(by='Occurrence', ascending=False, inplace=True)


    min_features = 5 if interaction_count < 500 else 20
    if max_features is None:
        max_features = interaction_count / 10 if interaction_count < 500 else interaction_count / 20

    for cut_off in cut_offs:
        features_occurrence_filter = features_occurrence_df[features_occurrence_df['Occurrence'] >= cut_off]
        num_features = len(features_occurrence_filter)
        print(f'Cut-off: {cut_off} - Features: {num_features}')

        if min_features < num_features < max_features:
            select_features = features_occurrence_filter['Feature'].tolist()
            id_vars = [sample_column]
            if 'phage' in full_feature_table.columns:
                id_vars.append('phage')
            if phenotype_column in full_feature_table.columns:
                id_vars.append(phenotype_column)
            if filter_type in full_feature_table.columns:
                id_vars.append(filter_type)

            id_vars = list(set(id_vars))

            select_feature_table = full_feature_table[id_vars + select_features]
            select_feature_table = select_feature_table.melt(
                id_vars=id_vars, var_name='Feature', value_name='Value'
            )
            if binary_data:
                select_feature_table['Value'] = select_feature_table['Value'].apply(lambda x: 1 if x > 0 else 0)
            
            select_feature_table = select_feature_table.pivot_table(
                index=id_vars, columns='Feature', values='Value'
            ).reset_index()

            os.makedirs(filter_table_dir, exist_ok=True)
            select_feature_table_path = os.path.join(filter_table_dir, f'select_feature_table_cutoff_{cut_off}.csv')
            print(f"Saving feature table to {select_feature_table_path}")
            select_feature_table.to_csv(select_feature_table_path, index=False)