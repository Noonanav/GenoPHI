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
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, confusion_matrix, roc_curve, auc, precision_recall_curve, average_precision_score, mean_squared_error, r2_score, hamming_loss, precision_recall_fscore_support
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from tqdm import tqdm
from hdbscan import HDBSCAN
import time


def drop_single_class_columns(y_train, *other_frames):
    """Drop multi-label/-target columns that are single-valued in y_train.

    In a CV fold a rare target can end up all-zeros (or all-ones) in the training
    split. CatBoost's MultiLogloss/MultiRMSE rejects such a column ("All train
    targets are equal"). Drop those columns from y_train (and the same columns
    from any other DataFrames passed, e.g. y_test) so the joint model trains on
    the learnable targets for this fold. A no-op for a Series (single target).

    Returns (y_train, *other_frames, dropped) with the dropped column names.
    """
    if not hasattr(y_train, 'columns'):
        return (y_train, *other_frames, [])
    keep = [c for c in y_train.columns if y_train[c].nunique(dropna=True) > 1]
    dropped = [c for c in y_train.columns if c not in keep]
    if dropped:
        logging.warning(
            f"Dropping single-class target column(s) for this split (unlearnable "
            f"jointly here): {dropped}"
        )
    y_train = y_train[keep]
    trimmed = tuple(f[keep] if hasattr(f, 'columns') else f for f in other_frames)
    return (y_train, *trimmed, dropped)


# Function to load and prepare data
def load_and_prepare_data(input_path, sample_column=None, phenotype_column=None, filter_type='none', task_type='classification'):
    """
    Loads the input feature table, drops unnecessary columns, and splits into features and target.

    Args:
        input_path (str): Path to the input CSV file containing the full feature table.
        sample_column (str): Optional name of the column to retain for sample identifiers.
        phenotype_column (str): Optional name of the column to retain for phenotype information.
        filter_type (str): Filter type for the dataset ('none', 'strain', 'phage', 'dataset').
        task_type (str): Task type ('classification' or 'regression') for validation. Default is 'classification'.

    Returns:
        X (DataFrame): Features for modeling.
        y (Series): Target variable (interaction).
        full_feature_table (DataFrame): The complete feature table after cleaning.
    """
    from genophi.utils import validate_phenotype_task_type

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"The input file {input_path} does not exist.")

    full_feature_table = pd.read_csv(input_path)
    if full_feature_table.empty:
        raise ValueError("Input data is empty.")

    if full_feature_table.isnull().any().any():
        print("Removing rows with missing values...")
        full_feature_table.dropna(inplace=True)
    else:
        print("No missing values found.")
    full_feature_table = full_feature_table.reset_index(drop=True)

    # 1. Identify Target Column(s) First.
    #    phenotype_column may be a single name (str -> Series y, single-target)
    #    or a list of names (-> DataFrame y, multi-target / multi-label).
    if phenotype_column is None:
        target_columns = ['interaction']
        multi_target = False
    elif isinstance(phenotype_column, (list, tuple)):
        target_columns = list(phenotype_column)
        multi_target = True
    else:
        target_columns = [phenotype_column]
        multi_target = False

    missing = [c for c in target_columns if c not in full_feature_table.columns]
    if missing:
        raise ValueError(f"Target column(s) {missing} not found in input data.")

    if multi_target:
        y = full_feature_table[target_columns]            # DataFrame (N x K)
    else:
        y = full_feature_table[target_columns[0]]         # Series (single target)

    # 2. Validate that phenotype data type matches task type (Series or DataFrame).
    validate_phenotype_task_type(y, task_type,
                                 target_columns[0] if not multi_target else None)

    # 3. Define Columns to Drop (Metadata + ALL target columns to prevent leakage).
    drop_columns = ['strain', 'phage', 'header', 'contig_id', 'orf_ko', filter_type] + target_columns

    # Handle sample column logic
    if sample_column and sample_column != 'strain':
        drop_columns.append(sample_column)
    elif 'strain' not in drop_columns: # Ensure strain is dropped if it is the sample column
        drop_columns.append('strain')

    X = full_feature_table.drop(columns=drop_columns, errors='ignore')
    X = X.select_dtypes(include=[np.number])

    # Leakage guard: no target column may survive in the feature matrix.
    leaked = [c for c in target_columns if c in X.columns]
    if leaked:
        raise ValueError(
            f"Target column(s) {leaked} leaked into the feature matrix X. "
            f"This must not happen; aborting to avoid label leakage."
        )

    if multi_target:
        print(f"Targets ({len(target_columns)}): {target_columns}")
        print(f"Per-target positive counts:\n{y.sum().to_string()}")
    else:
        print(f"Number of positive samples: {y.sum()}")
        print(f"Number of negative samples: {len(y) - y.sum()}")
    print(f"Final feature count: {X.shape[1]}")
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
    filter_by_cluster_presence=False,  # NEW: Filter features by cluster/group presence
    min_cluster_presence=2,            # NEW: Minimum clusters/groups a feature must appear in
    ensure_balanced_split=True,
    max_attempts=10,
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
        filter_by_cluster_presence (bool): If True, only include features present in multiple clusters/groups.
        min_cluster_presence (int): Minimum number of clusters/groups a feature must be present in (default: 2).
        ensure_balanced_split (bool): If True, ensures both train and test sets contain at least one positive sample.
                                     Only applies when y contains only 0s and 1s.
        max_attempts (int): Maximum number of attempts to find a valid split when ensure_balanced_split is True.
        **kwargs: Additional parameters for the clustering method.

    Returns:
        X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids
    """
    # Check if y is binary (contains only 0s and 1s). For multi-target y
    # (DataFrame), require EVERY column to be binary for the binary-only
    # balanced-split logic to apply; otherwise skip it. Reduce to a scalar bool
    # so the conditionals below are unambiguous for both Series and DataFrame.
    binary_mask = ((y == 0) | (y == 1)).all()
    is_binary = bool(binary_mask.all()) if hasattr(binary_mask, 'all') else bool(binary_mask)

    # The balanced-split positive-coverage check below is written for a single
    # binary Series; only apply it when y is a single-column (Series) binary target.
    is_multi_target = hasattr(y, 'columns')
    apply_balanced_split = ensure_balanced_split and is_binary and not is_multi_target

    if ensure_balanced_split and not apply_balanced_split:
        if is_multi_target:
            logging.warning("ensure_balanced_split is True but y is multi-target (DataFrame). Per-column balanced split is not applied here; relying on the multi-label split path.")
        else:
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
                          use_clustering=False, ensure_balanced_split=ensure_balanced_split,
                          filter_by_cluster_presence=filter_by_cluster_presence, 
                          min_cluster_presence=min_cluster_presence)

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
            logging.info(f"Clustering using Hierarchical Clustering")
            
            # Ensure n_clusters doesn't exceed the number of samples
            num_samples = len(filter_type_feature_table)
            adjusted_n_clusters = min(n_clusters, num_samples - 1)
            
            if adjusted_n_clusters != n_clusters:
                logging.warning(f"Reduced number of clusters from {n_clusters} to {adjusted_n_clusters} because we only have {num_samples} samples")
            
            clusterer = AgglomerativeClustering(n_clusters=adjusted_n_clusters)
            cluster_labels = clusterer.fit_predict(filter_type_feature_table[feature_columns])
            
            print(f"Number of clusters: {len(np.unique(cluster_labels))}")

        else:
            raise ValueError(f"Unsupported clustering method: {cluster_method}. Choose 'hdbscan' or 'hierarchical'.")

        # Assign clusters to samples
        filter_type_feature_table["cluster"] = cluster_labels
        print(f"Number of clusters for splitting: {len(filter_type_feature_table['cluster'].unique())}")
        
        # Create mapping dictionary with validation
        try:
            cluster_mapping = dict(zip(filter_type_feature_table[filter_type], 
                                    filter_type_feature_table["cluster"]))
            logging.info(f"Created cluster mapping for {len(cluster_mapping)} {filter_type} groups")
            
            # Apply mapping to add cluster column
            full_feature_table["cluster"] = full_feature_table[filter_type].map(cluster_mapping)
            
            # Validate mapping success
            unmapped_rows = full_feature_table["cluster"].isna()
            if unmapped_rows.any():
                unmapped_count = unmapped_rows.sum()
                unmapped_values = full_feature_table.loc[unmapped_rows, filter_type].unique()
                logging.error(f"Failed to map {unmapped_count} rows to clusters. "
                            f"Unmapped {filter_type} values: {unmapped_values}")
                raise ValueError(f"Cluster mapping failed for {unmapped_count} rows")
            
            logging.info("Cluster mapping completed successfully")
            
        except Exception as e:
            logging.error(f"Error during cluster mapping: {e}")
            raise

        if output_dir:
            cluster_file = os.path.join(output_dir, f"{filter_type}_clusters.csv")
            cluster_df = filter_type_feature_table[[filter_type, "cluster"]].drop_duplicates()
            cluster_df.to_csv(cluster_file, index=False)
            logging.info(f"Saved cluster labels to: {cluster_file}")

        group_col = "cluster"  # Use cluster labels for splitting
    else:
        group_col = filter_type  # Use original filter type column

    # ---- NEW: Filter by Cluster/Group Presence if Requested ----
    if filter_by_cluster_presence:
        # Determine what to use as "groups" for filtering
        if use_clustering and 'cluster' in full_feature_table.columns:
            group_col_for_filtering = 'cluster'
            logging.info(f"Filtering features by cluster presence (min_cluster_presence={min_cluster_presence})")
        else:
            group_col_for_filtering = filter_type
            logging.info(f"Filtering features by {filter_type} group presence (min_group_presence={min_cluster_presence})")
        
        # Get feature columns that match those in X
        feature_columns_for_filtering = [col for col in feature_columns if col in X.columns]
        
        if feature_columns_for_filtering:
            # Create a matrix showing feature presence by group/cluster
            group_feature_presence = full_feature_table.groupby(group_col_for_filtering)[feature_columns_for_filtering].apply(
                lambda group: (group > 0).any()
            )
            
            # Count how many groups/clusters each feature appears in
            feature_group_counts = group_feature_presence.sum(axis=0)
            
            # Keep features that appear in at least min_cluster_presence groups/clusters
            valid_features = feature_group_counts[feature_group_counts >= min_cluster_presence].index.tolist()
            
            # Only filter features of the filter_type, keep all other features
            original_feature_count = len(feature_columns_for_filtering)
            
            # Get all features that are NOT of the filter_type (e.g., keep phage features when filtering strain)
            other_features = [col for col in X.columns if col not in feature_columns_for_filtering]
            
            # Combine valid filtered features with all other features
            all_features_to_keep = valid_features + other_features
            X = X[all_features_to_keep]
            
            logging.info(f"Original {filter_type} features: {original_feature_count}")
            logging.info(f"{filter_type} features present in >= {min_cluster_presence} groups/clusters: {len(valid_features)}")
            logging.info(f"{filter_type} features removed: {original_feature_count - len(valid_features)}")
            logging.info(f"Other feature types kept: {len(other_features)}")
        else:
            logging.warning("No matching feature columns found for cluster/group-based filtering")

    # ---- 5️⃣ Perform Group-Based Splitting (Clustering or Normal) ----
    # Sort groups to ensure deterministic order (pandas unique() order is not guaranteed)
    groups = np.sort(full_feature_table[group_col].unique())

    # Calculate group/cluster sizes for sample-based splitting
    group_sizes = {}
    for group in groups:
        group_mask = full_feature_table[group_col] == group
        group_sizes[group] = group_mask.sum()
    
    total_samples = sum(group_sizes.values())
    target_train_size = int(0.8 * total_samples)

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
            # CHANGED: Fall back to sample-based random splitting instead of cluster-count splitting
            np.random.seed(random_state)
            all_groups = list(groups)
            np.random.shuffle(all_groups)
            
            train_groups = []
            current_train_size = 0
            
            for group in all_groups:
                train_groups.append(group)
                current_train_size += group_sizes[group]
                if current_train_size >= target_train_size:
                    break
            
            test_groups = np.setdiff1d(groups, train_groups)
        else:
            # Try multiple random seeds to find a valid split
            found_valid_split = False
            for attempt in range(max_attempts):
                # Use a more diverse seed generation strategy
                current_seed = random_state * 1000 + attempt * 17
                np.random.seed(current_seed)
                
                # Sample-based splitting while ensuring positive groups in both sets
                # Shuffle all groups to randomize selection order
                all_groups = list(groups)
                np.random.shuffle(all_groups)
                
                # Keep adding groups until we reach target, ensuring positive groups are distributed
                train_groups = []
                current_train_size = 0
                used_positive_groups = []
                
                for group in all_groups:
                    # Add the group
                    train_groups.append(group)
                    current_train_size += group_sizes[group]
                    if group in positive_groups:
                        used_positive_groups.append(group)
                    
                    # Stop when we reach target size
                    if current_train_size >= target_train_size:
                        break
                
                test_groups = np.setdiff1d(groups, train_groups)
                remaining_positive_groups = np.setdiff1d(positive_groups, used_positive_groups)
                
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
        # Standard sample-based split without ensuring positive samples
        np.random.seed(random_state)
        all_groups = list(groups)
        np.random.shuffle(all_groups)
        
        train_groups = []
        current_train_size = 0
        
        for group in all_groups:
            train_groups.append(group)
            current_train_size += group_sizes[group]
            if current_train_size >= target_train_size:
                break
        
        test_groups = np.setdiff1d(groups, train_groups)

    train_idx = full_feature_table[group_col].isin(train_groups)
    test_idx = full_feature_table[group_col].isin(test_groups)

    # Post-splitting validation: ensure we have valid train/test sets
    train_sample_count = train_idx.sum() if isinstance(train_idx, pd.Series) else len(train_idx)
    test_sample_count = test_idx.sum() if isinstance(test_idx, pd.Series) else len(test_idx)

    # Fallback to random split if group-based splitting produced invalid results
    if train_sample_count == 0 or test_sample_count == 0 or test_sample_count < 2:
        logging.warning(f"Group-based splitting produced invalid sets (train: {train_sample_count}, test: {test_sample_count}). "
                    f"Falling back to random split to ensure valid train/test sets.")
        
        # Use the same balanced splitting logic as the random split section
        if apply_balanced_split:
            for attempt in range(max_attempts):
                current_seed = random_state * 1000 + attempt * 17
                X_train_temp, X_test_temp, y_train_temp, y_test_temp = train_test_split(
                    X, y, test_size=0.2, random_state=current_seed
                )
                
                # Check if both sets have positives
                if y_train_temp.sum() > 0 and y_test_temp.sum() > 0:
                    # Create boolean masks instead of using indices directly
                    train_idx = pd.Series(False, index=full_feature_table.index)
                    test_idx = pd.Series(False, index=full_feature_table.index)
                    train_idx.loc[X_train_temp.index] = True
                    test_idx.loc[X_test_temp.index] = True
                    logging.info(f"Fallback random split successful on attempt {attempt+1}")
                    break
                
                if attempt == max_attempts - 1:
                    # Use last attempt even if not perfect
                    train_idx = pd.Series(False, index=full_feature_table.index)
                    test_idx = pd.Series(False, index=full_feature_table.index)
                    train_idx.loc[X_train_temp.index] = True
                    test_idx.loc[X_test_temp.index] = True
                    logging.warning(f"Fallback random split used last attempt after {max_attempts} tries")
        else:
            # Simple random split for non-balanced case
            X_train_temp, X_test_temp, y_train_temp, y_test_temp = train_test_split(
                X, y, test_size=0.2, random_state=random_state
            )
            # Create boolean masks instead of using indices directly
            train_idx = pd.Series(False, index=full_feature_table.index)
            test_idx = pd.Series(False, index=full_feature_table.index)
            train_idx.loc[X_train_temp.index] = True
            test_idx.loc[X_test_temp.index] = True
            logging.info("Fallback to simple random split successful")

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
    max_ram=8,
    random_state=42
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

    # Multi-output: drop any target column that is single-class in this split --
    # MultiLogloss/MultiRMSE reject an all-equal column ("All train targets are
    # equal"). Train jointly on the learnable targets for this split.
    y_train, dropped = drop_single_class_columns(y_train)
    if hasattr(y_train, 'columns') and y_train.shape[1] == 0:
        raise ValueError("All target columns are single-class in this split; "
                         "cannot run multi-output RFE.")

    # Detect multi-output y (DataFrame): RFE still applies because the CatBoost
    # MultiLogloss/MultiRMSE estimator exposes a single joint feature-importance
    # vector; only the loss function / y shape differs. Classification ->
    # MultiLogloss; regression -> MultiRMSE.
    is_multi_output = hasattr(y_train, 'columns')
    is_multilabel = is_multi_output and task_type == 'classification'
    is_multitarget = is_multi_output and task_type == 'regression'

    # Initialize model based on task type
    if is_multilabel:
        model = CatBoostClassifier(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            loss_function='MultiLogloss',
            verbose=10,
            thread_count=threads,
            train_dir=os.path.join(output_dir, '..', 'catboost_info'),
            used_ram_limit=f"{max_ram}gb",
            random_seed=random_state
        )
    elif is_multitarget:
        model = CatBoostRegressor(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            loss_function='MultiRMSE',
            verbose=10,
            thread_count=threads,
            train_dir=os.path.join(output_dir, '..', 'catboost_info'),
            used_ram_limit=f"{max_ram}gb",
            random_seed=random_state
        )
    elif task_type == 'classification':
        model = CatBoostClassifier(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=10,
            thread_count=threads,
            train_dir=os.path.join(output_dir, '..', 'catboost_info'),
            used_ram_limit=f"{max_ram}gb",  # Set the RAM limit
            random_seed=random_state
        )
    elif task_type == 'regression':
        model = CatBoostRegressor(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=10,
            thread_count=threads,
            train_dir=os.path.join(output_dir, '..', 'catboost_info'),
            used_ram_limit=f"{max_ram}gb",  # Set the RAM limit
            random_seed=random_state
        )
    else:
        raise ValueError("task_type must be 'classification' or 'regression'")

    print(f"Performing Recursive Feature Elimination (RFE) with step_size: {step_size} for {task_type}...")

    # If requested, build sample weights for each row. Dynamic phage weights are
    # defined for a single binary target; skip them for multi-label y.
    sample_weight = None
    if is_multi_output and use_dynamic_weights:
        logging.warning("Dynamic phage weights are not supported for multi-output RFE; proceeding without them.")
    if (not is_multi_output) and use_dynamic_weights and phage_column and (phage_column in X_train_sample_ids.columns):
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

def shap_rfe(
    X_train, 
    y_train, 
    num_features, 
    threads, 
    task_type='classification', 
    max_ram=8, 
    random_state=42
):
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
            used_ram_limit=f"{max_ram}gb",  # Set the RAM limit
            random_seed=random_state
        )
    elif task_type == 'regression':
        model = CatBoostRegressor(
            iterations=500,
            learning_rate=0.1,
            depth=4,
            verbose=0,
            thread_count=threads,
            used_ram_limit=f"{max_ram}gb",  # Set the RAM limit
            random_seed=random_state
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

def lasso_feature_selection(X_train, y_train, num_features, task_type='classification', random_state=42):
    """
    Selects top features using Lasso regularization for classification or regression.

    Args:
        X_train (DataFrame): Training features.
        y_train (Series): Training target.
        num_features (int): Number of top features to select.
        task_type (str): Task type ('classification' or 'regression').
        random_state (int): Random seed for reproducibility.

    Returns:
        X_train_selected (DataFrame): Transformed training features with selected features.
        selected_features (Index): List of selected feature names.
    """
    print("Selecting features using Lasso regularization...")

    # Choose model based on task type
    if task_type == 'classification':
        model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=1000, random_state=random_state)
    elif task_type == 'regression':
        model = Lasso(max_iter=1000, random_state=random_state)
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

def shap_feature_selection(
    X_train, 
    y_train, 
    num_features, 
    threads, 
    task_type='classification', 
    max_ram=8, 
    random_state=42
):
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
        model = CatBoostClassifier(
            iterations=500, 
            learning_rate=0.1, 
            depth=4, 
            verbose=0, 
            thread_count=threads, 
            used_ram_limit=f"{max_ram}gb",
            random_seed=random_state
        )
    elif task_type == 'regression':
        model = CatBoostRegressor(
            iterations=500, 
            learning_rate=0.1, 
            depth=4, 
            verbose=0, 
            thread_count=threads, 
            used_ram_limit=f"{max_ram}gb",
            random_seed=random_state
        )
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
                       weights_method='log10',
                       random_state=42):
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
    model = CatBoostClassifier(**params, train_dir=train_dir, used_ram_limit=f"{max_ram}gb", random_seed=random_state)

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

def train_and_evaluate_regressor(
    X_train, 
    y_train, 
    X_test, 
    y_test, 
    params, 
    output_dir, 
    max_ram=8, 
    random_state=42
):
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
    model = CatBoostRegressor(**params, train_dir=train_dir, used_ram_limit=f"{max_ram}gb", random_seed=random_state)

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
    max_ram=8,
    random_state=42
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
            weights_method=weights_method,
            random_state=random_state
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


def _tune_per_label_thresholds(y_true, proba, label_names):
    """Pick the decision threshold per label that maximizes F1 on the given set.

    Args:
        y_true (np.ndarray): N x K binary ground-truth matrix.
        proba (np.ndarray): N x K predicted-probability matrix (P(label=1)).
        label_names (list[str]): K label column names.

    Returns:
        dict[str, float]: label name -> chosen threshold (defaults to 0.5).
    """
    thresholds = {}
    candidate_grid = np.linspace(0.05, 0.95, 19)
    for k, name in enumerate(label_names):
        yk = np.asarray(y_true)[:, k]
        pk = np.asarray(proba)[:, k]
        # If a label has no positives in this set, F1 tuning is undefined; use 0.5.
        if yk.sum() == 0:
            thresholds[name] = 0.5
            continue
        best_t, best_f1 = 0.5, -1.0
        for t in candidate_grid:
            preds = (pk >= t).astype(int)
            f1 = f1_score(yk, preds, zero_division=0)
            if f1 > best_f1:
                best_f1, best_t = f1, float(t)
        thresholds[name] = best_t
    return thresholds


def _multilabel_proba(model, X):
    """Extract an N x K P(label=1) matrix from a CatBoost MultiLogloss model.

    CatBoost's predict_proba for MultiLogloss returns a list/array shaped per
    label; normalize it to a dense N x K array of positive-class probabilities.
    """
    proba = model.predict_proba(X)
    proba = np.asarray(proba)
    # Some CatBoost versions return (K, N, 2); others (N, K). Normalize to (N, K).
    if proba.ndim == 3:
        # (K, N, 2) -> take positive class, transpose to (N, K)
        proba = proba[:, :, 1].T
    return proba


def train_and_evaluate_multilabel(
    X_train, y_train, X_test, y_test, params, output_dir,
    max_ram=8, random_state=42,
):
    """Train a single CatBoost MultiLogloss model over all labels jointly.

    Args:
        X_train, X_test (DataFrame): Features.
        y_train, y_test (DataFrame): N x K binary label matrices.
        params (dict): CatBoost hyperparameters (must request MultiLogloss).
        output_dir (str): Directory for CatBoost train_dir.

    Returns:
        model: the trained CatBoostClassifier.
        metrics (dict): aggregate metrics (micro/macro F1, Hamming loss, subset acc).
        thresholds (dict): per-label tuned thresholds (from the TEST set here;
            in production tune on a validation split).
        proba (np.ndarray): N x K probability matrix on X_test.
    """
    label_names = list(y_train.columns)
    train_dir = os.path.join(output_dir, '..', 'catboost_info')
    model = CatBoostClassifier(
        **params, train_dir=train_dir,
        used_ram_limit=f"{max_ram}gb", random_seed=random_state,
    )
    print(f"Training MultiLogloss model ({len(label_names)} labels) with params: {params}")
    model.fit(X_train, y_train, eval_set=(X_test, y_test),
              plot=False, verbose=10, early_stopping_rounds=100)

    proba = _multilabel_proba(model, X_test)
    thresholds = _tune_per_label_thresholds(y_test.values, proba, label_names)

    thresh_vec = np.array([thresholds[n] for n in label_names])
    y_pred = (proba >= thresh_vec).astype(int)
    y_true = y_test.values

    micro_f1 = f1_score(y_true, y_pred, average='micro', zero_division=0)
    macro_f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)
    h_loss = hamming_loss(y_true, y_pred)
    subset_acc = float((y_pred == y_true).all(axis=1).mean())

    metrics = {
        'micro_f1': micro_f1, 'macro_f1': macro_f1,
        'hamming_loss': h_loss, 'subset_accuracy': subset_acc,
    }
    print(f"Multilabel done. micro-F1={micro_f1:.3f} macro-F1={macro_f1:.3f} "
          f"hamming={h_loss:.3f} subset_acc={subset_acc:.3f}")
    return model, metrics, thresholds, proba


def grid_search_multilabel(
    X_train, y_train, X_test, y_test,
    X_train_sample_ids, X_test_sample_ids,
    param_grid, output_dir, phenotype_columns,
    max_ram=8, random_state=42,
):
    """Grid search for a joint multi-label CatBoost (MultiLogloss) model.

    Ranks models on macro-F1 (handles label imbalance better than micro for the
    receptor-profile use case). Returns the single best model plus a metadata
    dict (target_names + per-label thresholds) for the sidecar JSON.

    Returns:
        best_model, best_params, best_macro_f1, metadata
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Drop labels that are single-class in this split (MultiLogloss rejects them);
    # trim y_test to the same labels so the matrices stay aligned.
    y_train, y_test, _ = drop_single_class_columns(y_train, y_test)
    if y_train.shape[1] == 0:
        logging.warning("All labels single-class in this split; no multilabel model.")
        return None, None, -1.0, None

    label_names = list(y_train.columns)
    best_macro_f1 = -1.0
    best_model = None
    best_params = None
    best_thresholds = None
    results = []

    print("Starting MultiLogloss grid search...")
    for idx, values in enumerate(itertools.product(*param_grid.values()), start=1):
        params = dict(zip(param_grid.keys(), values))
        model, metrics, thresholds, proba = train_and_evaluate_multilabel(
            X_train, y_train, X_test, y_test, params, output_dir,
            max_ram=max_ram, random_state=random_state,
        )
        results.append({**params, **metrics})

        if metrics['macro_f1'] >= best_macro_f1:
            best_macro_f1 = metrics['macro_f1']
            best_model = model
            best_params = params
            best_thresholds = thresholds

            # Save per-sample predictions using actual label names.
            thresh_vec = np.array([thresholds[n] for n in label_names])
            y_pred = (proba >= thresh_vec).astype(int)
            preds_df = X_test_sample_ids.copy().reset_index(drop=True)
            for k, name in enumerate(label_names):
                preds_df[f'Prediction_{name}'] = y_pred[:, k]
                preds_df[f'Confidence_{name}'] = proba[:, k]
                preds_df[f'True_{name}'] = y_test[name].values
            preds_df.to_csv(os.path.join(output_dir, 'best_model_predictions.csv'), index=False)

    pd.DataFrame(results).to_csv(os.path.join(output_dir, 'model_performance.csv'), index=False)

    if best_model is None:
        logging.warning("No valid multilabel model was found in grid search.")
        return None, None, best_macro_f1, None

    metadata = {
        'specific_task': 'multilabel_classification',
        'strategy': 'joint',
        'task_type': 'classification',
        'target_names': label_names,
        'per_label_thresholds': best_thresholds,
        'feature_names': list(X_train.columns),
    }
    return best_model, best_params, best_macro_f1, metadata


def train_and_evaluate_multitarget(
    X_train, y_train, X_test, y_test, params, output_dir,
    max_ram=8, random_state=42,
):
    """Train a single CatBoost MultiRMSE model over all targets jointly.

    Args:
        X_train, X_test (DataFrame): Features.
        y_train, y_test (DataFrame): N x T continuous target matrices.
        params (dict): CatBoost hyperparameters (must request MultiRMSE).
        output_dir (str): Directory for CatBoost train_dir.

    Returns:
        model: trained CatBoostRegressor.
        metrics (dict): per-target R2/RMSE/MAE + aggregate mean_r2,
            mean_normalized_rmse.
        y_pred (np.ndarray): N x T prediction matrix on X_test.
    """
    target_names = list(y_train.columns)
    train_dir = os.path.join(output_dir, '..', 'catboost_info')
    model = CatBoostRegressor(
        **params, train_dir=train_dir,
        used_ram_limit=f"{max_ram}gb", random_seed=random_state,
    )
    print(f"Training MultiRMSE model ({len(target_names)} targets) with params: {params}")
    model.fit(X_train, y_train, eval_set=(X_test, y_test),
              plot=False, verbose=10, early_stopping_rounds=100)

    y_pred = np.asarray(model.predict(X_test))
    y_true = y_test.values

    per_target = {}
    norm_rmses = []
    r2s = []
    for t_idx, name in enumerate(target_names):
        yt = y_true[:, t_idx]
        yp = y_pred[:, t_idx]
        rmse = float(np.sqrt(mean_squared_error(yt, yp)))
        mae = float(np.mean(np.abs(yt - yp)))
        r2 = float(r2_score(yt, yp))
        # Normalize RMSE by target spread so targets on different scales are
        # comparable (avoids large-scale targets dominating the aggregate).
        spread = float(yt.max() - yt.min())
        nrmse = rmse / spread if spread > 0 else float('nan')
        per_target[name] = {'r2': r2, 'rmse': rmse, 'mae': mae, 'normalized_rmse': nrmse}
        r2s.append(r2)
        if not np.isnan(nrmse):
            norm_rmses.append(nrmse)

    metrics = {f'r2_{n}': per_target[n]['r2'] for n in target_names}
    metrics.update({f'rmse_{n}': per_target[n]['rmse'] for n in target_names})
    metrics.update({f'mae_{n}': per_target[n]['mae'] for n in target_names})
    metrics.update({f'normalized_rmse_{n}': per_target[n]['normalized_rmse'] for n in target_names})
    metrics['mean_r2'] = float(np.mean(r2s))
    metrics['mean_normalized_rmse'] = float(np.mean(norm_rmses)) if norm_rmses else float('nan')

    print(f"Multitarget done. mean_R2={metrics['mean_r2']:.3f} "
          f"mean_nRMSE={metrics['mean_normalized_rmse']:.3f}")
    return model, metrics, y_pred


def grid_search_multitarget(
    X_train, y_train, X_test, y_test,
    X_train_sample_ids, X_test_sample_ids,
    param_grid, output_dir, phenotype_columns,
    max_ram=8, random_state=42,
):
    """Grid search for a joint multi-target CatBoost (MultiRMSE) model.

    Ranks models on mean per-target R2 (the higher-is-better counterpart of the
    binary path's r2). Returns the single best model plus a metadata dict.

    Returns:
        best_model, best_params, best_mean_r2, metadata
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Drop targets that are single-valued in this split (MultiRMSE rejects them);
    # trim y_test to match.
    y_train, y_test, _ = drop_single_class_columns(y_train, y_test)
    if y_train.shape[1] == 0:
        logging.warning("All targets single-valued in this split; no multitarget model.")
        return None, None, -float('inf'), None

    target_names = list(y_train.columns)
    best_mean_r2 = -float('inf')
    best_model = None
    best_params = None
    results = []

    print("Starting MultiRMSE grid search...")
    for idx, values in enumerate(itertools.product(*param_grid.values()), start=1):
        params = dict(zip(param_grid.keys(), values))
        model, metrics, y_pred = train_and_evaluate_multitarget(
            X_train, y_train, X_test, y_test, params, output_dir,
            max_ram=max_ram, random_state=random_state,
        )
        results.append({**params, **metrics})

        if metrics['mean_r2'] >= best_mean_r2:
            best_mean_r2 = metrics['mean_r2']
            best_model = model
            best_params = params

            preds_df = X_test_sample_ids.copy().reset_index(drop=True)
            for k, name in enumerate(target_names):
                preds_df[f'Prediction_{name}'] = y_pred[:, k]
                preds_df[f'True_{name}'] = y_test[name].values
            preds_df.to_csv(os.path.join(output_dir, 'best_model_predictions.csv'), index=False)

    pd.DataFrame(results).to_csv(os.path.join(output_dir, 'model_performance.csv'), index=False)

    if best_model is None:
        logging.warning("No valid multitarget model was found in grid search.")
        return None, None, best_mean_r2, None

    metadata = {
        'specific_task': 'multitarget_regression',
        'strategy': 'joint',
        'task_type': 'regression',
        'target_names': target_names,
        'feature_names': list(X_train.columns),
    }
    return best_model, best_params, best_mean_r2, metadata


def grid_search_regressor(
    X_train, 
    y_train, 
    X_test, 
    y_test, 
    X_test_sample_ids, 
    param_grid, 
    output_dir, 
    phenotype_column='interaction', 
    max_ram=8, 
    random_state=42
):
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
        model, mse, r2, y_pred = train_and_evaluate_regressor(
            X_train, 
            y_train, 
            X_test, 
            y_test, 
            params, 
            output_dir, 
            max_ram=max_ram,
            random_state=random_state
        )
        
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
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
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
            # Set global random seeds for full reproducibility
            np.random.seed(random_state)

            X, y, full_feature_table = load_and_prepare_data(input_path, sample_column=sample_column, phenotype_column=phenotype_column, filter_type=filter_type, task_type=task_type)
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
                check_feature_presence=check_feature_presence,
                filter_by_cluster_presence=filter_by_cluster_presence,
                min_cluster_presence=min_cluster_presence
            )

            if num_features == 'none':
                num_features = int(len(X_train) / 10) if len(X_train) < 500 else int(len(X_train) / 20)
            else:
                num_features = int(num_features)

            if X_train is None:
                logging.info("Skipping this run due to insufficient training data.")
                continue  # Skip this iteration and proceed to the next

            # Multi-label y (DataFrame): only the default 'rfe' method is
            # supported for joint multi-label feature selection (it uses a
            # CatBoost MultiLogloss estimator whose joint importances drive RFE).
            if hasattr(y_train, 'columns') and method != 'rfe':
                raise NotImplementedError(
                    f"Feature selection method '{method}' is not supported for "
                    f"joint multi-label targets. Use method='rfe' (default)."
                )

            # Apply selected feature selection method
            if method == 'rfe':
                _, selected_features = perform_rfe(X_train, y_train, X_train_sample_ids, num_features, threads, output_dir, task_type=task_type, phage_column=phage_column, use_dynamic_weights=use_dynamic_weights, max_ram=max_ram, random_state=random_state)
            elif method == 'shap_rfe':
                X_train, selected_features = shap_rfe(X_train, y_train, num_features, threads, task_type=task_type, max_ram=max_ram, random_state=random_state)
            elif method == 'select_k_best':
                X_train, selected_features = select_k_best_feature_selection(X_train, y_train, num_features, task_type=task_type)
            elif method == 'chi_squared' and task_type == 'classification':
                X_train, selected_features = chi_squared_feature_selection(X_train, y_train, num_features)
            elif method == 'lasso':
                X_train, selected_features = lasso_feature_selection(X_train, y_train, num_features, task_type=task_type, random_state=random_state)
            elif method == 'shap':
                X_train, selected_features = shap_feature_selection(X_train, y_train, num_features, threads, task_type=task_type, max_ram=max_ram, random_state=random_state)
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

            if hasattr(y_train, 'columns') and task_type == 'classification':
                # Joint multi-label: evaluate selected features with a single
                # MultiLogloss model (ranked on macro-F1).
                param_grid['loss_function'] = ['MultiLogloss']
                best_model, best_params, best_macro_f1, _meta = grid_search_multilabel(
                    X_train_selected,
                    y_train,
                    X_test_selected,
                    y_test,
                    X_train_sample_ids,
                    X_test_sample_ids,
                    param_grid,
                    output_dir,
                    phenotype_columns=list(y_train.columns),
                    max_ram=max_ram,
                    random_state=random_state
                )
                best_metric = best_macro_f1
            elif hasattr(y_train, 'columns') and task_type == 'regression':
                # Joint multi-target: evaluate with a single MultiRMSE model
                # (ranked on mean per-target R2).
                param_grid['loss_function'] = ['MultiRMSE']
                best_model, best_params, best_mean_r2, _meta = grid_search_multitarget(
                    X_train_selected,
                    y_train,
                    X_test_selected,
                    y_test,
                    X_train_sample_ids,
                    X_test_sample_ids,
                    param_grid,
                    output_dir,
                    phenotype_columns=list(y_train.columns),
                    max_ram=max_ram,
                    random_state=random_state
                )
                best_metric = best_mean_r2
            elif task_type == 'classification':
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
                    max_ram=max_ram,
                    random_state=random_state
                )
                best_metric = best_mcc
            elif task_type == 'regression':
                param_grid['loss_function'] = ['RMSE']
                best_model, best_params, best_r2 = grid_search_regressor(
                    X_train_selected, 
                    y_train, 
                    X_test_selected, 
                    y_test, 
                    X_test_sample_ids, 
                    param_grid, 
                    output_dir, 
                    phenotype_column=phenotype_column, 
                    max_ram=max_ram,
                    random_state=random_state
                )
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
    binary_data=False, max_features=None, min_features=None, filter_type='strain'
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
        max_features (int, float, str, or None): Maximum number of features (exclusive upper bound).
            If None, auto-calculated based on interaction count. Can be passed as string 'none'.
        min_features (int, str, or None): Minimum number of features required (inclusive lower bound).
            If None, auto-calculated based on interaction count. Can be passed as string 'none'.
        filter_type (str): Type of filter to apply ('strain', 'phage', etc.).
    """
    if phenotype_column is None:
        phenotype_column = 'interaction'

    full_feature_table = pd.read_csv(full_feature_table_file)
    interaction_count = full_feature_table.shape[0]
    print('Interaction count:', interaction_count)

    run_dirs = sorted([x for x in os.listdir(model_testing_dir) if 'run' in x])
    features_occurrence = {}

    for run in run_dirs:
        feature_importances_path = os.path.join(model_testing_dir, run, 'feature_importances.csv')
        if os.path.exists(feature_importances_path):
            features_df = pd.read_csv(feature_importances_path)
            for feature in features_df['Feature'].values:
                features_occurrence[feature] = features_occurrence.get(feature, 0) + 1

    features_occurrence_df = pd.DataFrame(list(features_occurrence.items()), columns=['Feature', 'Occurrence'])
    features_occurrence_df.sort_values(by='Occurrence', ascending=False, inplace=True)

    # Validate and convert max_features if it is passed as a string (e.g., from CLI args)
    if isinstance(max_features, str):
        if max_features.lower() == 'none':
            max_features = None
        else:
            try:
                max_features = float(max_features)
            except ValueError:
                logging.warning(f"Invalid max_features value: {max_features}. Using default calculation.")
                max_features = None

    # Validate and convert min_features if it is passed as a string (e.g., from CLI args)
    if isinstance(min_features, str):
        if min_features.lower() == 'none':
            min_features = None
        else:
            try:
                min_features = int(min_features)
            except ValueError:
                logging.warning(f"Invalid min_features value: {min_features}. Using default calculation.")
                min_features = None

    if min_features is None:
        min_features = 5 if interaction_count < 500 else 20
    if max_features is None:
        max_features = interaction_count / 10 if interaction_count < 500 else interaction_count / 20

    for cut_off in cut_offs:
        features_occurrence_filter = features_occurrence_df[features_occurrence_df['Occurrence'] >= cut_off]
        num_features = len(features_occurrence_filter)
        print(f'Cut-off: {cut_off} - Features: {num_features}')

        if min_features <= num_features < max_features:
            select_features = features_occurrence_filter['Feature'].tolist()
            id_vars = [sample_column]
            if 'phage' in full_feature_table.columns:
                id_vars.append('phage')
            # phenotype_column may be a single name or a list (multi-target);
            # keep ALL present target columns in the per-cutoff feature table.
            target_cols = phenotype_column if isinstance(phenotype_column, (list, tuple)) else [phenotype_column]
            for tcol in target_cols:
                if tcol in full_feature_table.columns:
                    id_vars.append(tcol)
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