import os
import pandas as pd
import numpy as np
import itertools
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef, roc_curve, precision_recall_curve, mean_squared_error, r2_score, average_precision_score
from plotnine import ggplot, aes, geom_line, geom_abline, labs, theme, guides, guide_legend, element_rect, element_blank, element_line, geom_vline, geom_jitter, geom_point, geom_smooth
from tqdm import tqdm
import time
import joblib
import re
from genophi.feature_selection import load_and_prepare_data, filter_data, train_and_evaluate, grid_search, save_feature_importances, grid_search_regressor, train_and_evaluate_regressor, grid_search_multilabel, grid_search_multitarget
from genophi.utils import detect_task_type, MULTILABEL_CLASSIFICATION, MULTITARGET_REGRESSION
import json
import shap
import matplotlib.pyplot as plt
import logging
import gc

# Set environment variables to control threading
os.environ['OMP_NUM_THREADS'] = '12'
os.environ['OPENBLAS_NUM_THREADS'] = '12'
os.environ['MKL_NUM_THREADS'] = '12'
os.environ['VECLIB_MAXIMUM_THREADS'] = '12'
os.environ['NUMEXPR_NUM_THREADS'] = '12'

def plot_custom_shap_beeswarm(shap_values_df, output_dir, prefix=None, binary_data=False):
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
    if binary_data:
        full_shap_values_df_top20['value'] = full_shap_values_df_top20['value'].astype(str)
    else:
        full_shap_values_df_top20['value'] = full_shap_values_df_top20['value'].astype(float)
        full_shap_values_df_top20['value'] = full_shap_values_df_top20.groupby('feature')['value'].transform(lambda x: x / x.max())

    # Custom beeswarm plot using plotnine
    shap_plot = (
        ggplot(full_shap_values_df_top20, aes(y='reorder(feature, shap_importance)', x='shap_value', fill='value', group='feature')) +
        geom_vline(xintercept=0, color='black', size=0.5) +
        geom_jitter(height=0.2, alpha=0.5, stroke=0, size=2) +
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

    del shap_plot, top_20_shap_df, full_shap_values_df_top20
    gc.collect()

def model_testing_select_MCC(
    input, 
    output_dir, 
    threads, 
    random_state, 
    task_type='classification', 
    set_filter='none', 
    sample_column=None, 
    phenotype_column=None, 
    phage_column='phage',
    use_dynamic_weights=False,
    weights_method='log10',
    binary_data=False, 
    max_ram=8, 
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    use_shap=False,
    target_mode='auto',
    strategy='joint'
):
    """
    Runs a single experiment for feature table, training a CatBoost model with grid search and saving results.

    Args:
        input (str): Path to input feature table.
        output_dir (str): Directory to store results.
        threads (int): Number of threads for training.
        random_state (int): Seed for reproducibility.
        task_type (str): Specifies 'classification' or 'regression' task.
        set_filter (str): Filter type for the dataset ('none', 'strain', 'phage', 'dataset').
        sample_column (str): Name of the sample column (optional).
        phenotype_column (str): Name of the phenotype column (optional).
        binary_data (bool): If True, plot SHAP jitter plot with binary data.
        use_shap (bool): If True, calculate and save SHAP values.
    """
    start_time = time.time()

    # Set global random seed for full reproducibility
    np.random.seed(random_state)

    # Load and prepare data
    X, y, full_feature_table = load_and_prepare_data(input, sample_column, phenotype_column, filter_type=set_filter, task_type=task_type)
    X_train, X_test, y_train, y_test, X_test_sample_ids, X_train_sample_ids = filter_data(
        X, y, 
        full_feature_table, 
        set_filter, 
        sample_column=sample_column if sample_column else 'strain', 
        random_state=random_state,
        output_dir=output_dir,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon
    )
    print(f"Training data shape: {X_train.shape}, Testing data shape: {X_test.shape}")

    if X_train is None:
        logging.info("Skipping this run due to insufficient training data.")
        return  # Exit this function, skipping further processing for this run

    # Resolve the specific task (binary / multiclass / multilabel / regression)
    # and strategy. For a DataFrame y with strategy='joint', this routes to the
    # joint multi-output trainer; single-target y is unchanged.
    specific_task, resolved_strategy = detect_task_type(
        y_train, task_type, target_mode=target_mode, strategy=strategy
    )

    metadata = None

    # Define task-specific grid search parameters
    if specific_task == MULTILABEL_CLASSIFICATION and resolved_strategy == 'joint':
        param_grid = {
            'iterations': [500, 1000],
            'learning_rate': [0.05, 0.1],
            'depth': [4, 6],
            'loss_function': ['MultiLogloss'],
            'thread_count': [threads]
        }
        best_model, best_params, best_macro_f1, metadata = grid_search_multilabel(
            X_train, y_train, X_test, y_test,
            X_train_sample_ids, X_test_sample_ids,
            param_grid, output_dir,
            phenotype_columns=list(y_train.columns),
            max_ram=max_ram
        )
        best_metric = best_macro_f1

    elif specific_task == MULTITARGET_REGRESSION and resolved_strategy == 'joint':
        param_grid = {
            'iterations': [500, 1000],
            'learning_rate': [0.05, 0.1],
            'depth': [4, 6],
            'loss_function': ['MultiRMSE'],
            'thread_count': [threads]
        }
        best_model, best_params, best_mean_r2, metadata = grid_search_multitarget(
            X_train, y_train, X_test, y_test,
            X_train_sample_ids, X_test_sample_ids,
            param_grid, output_dir,
            phenotype_columns=list(y_train.columns),
            max_ram=max_ram
        )
        best_metric = best_mean_r2

    elif task_type == 'classification':
        param_grid = {
            'iterations': [500, 1000],
            'learning_rate': [0.05, 0.1],
            'depth': [4, 6],
            'loss_function': ['Logloss'],
            'thread_count': [threads]
        }
        best_model, best_params, best_mcc = grid_search(
            X_train, 
            y_train, 
            X_test, 
            y_test, 
            X_train_sample_ids,
            X_test_sample_ids, 
            param_grid, 
            output_dir, 
            phenotype_column, 
            phage_column=phage_column,
            use_dynamic_weights=use_dynamic_weights,
            weights_method=weights_method,
            max_ram=max_ram
        )
        best_metric = best_mcc

    elif task_type == 'regression':
        param_grid = {
            'iterations': [500, 1000],
            'learning_rate': [0.05, 0.1],
            'depth': [4, 6],
            'loss_function': ['RMSE'],
            'thread_count': [threads]
        }
        best_model, best_params, best_r2 = grid_search_regressor(
            X_train, y_train, X_test, y_test, X_test_sample_ids, param_grid, output_dir, phenotype_column, max_ram=max_ram
        )
        best_metric = best_r2

    else:
        raise ValueError("task_type must be 'classification' or 'regression'")

    if best_model is None:
        print(f"Skipping iteration: No model found with the best metric value: {best_metric}")
        return

    print(f"Best Model Parameters: {best_params}, Best Metric ({specific_task}): {best_metric}")

    is_multilabel = (specific_task == MULTILABEL_CLASSIFICATION and resolved_strategy == 'joint')
    is_multitarget = (specific_task == MULTITARGET_REGRESSION and resolved_strategy == 'joint')
    is_multi_output = is_multilabel or is_multitarget

    # Save feature importances (single-output models only; multi-output SHAP/
    # importances are per-target tensors handled separately and deferred).
    if not is_multi_output:
        feature_importances_path = os.path.join(output_dir, "feature_importances.csv")
        save_feature_importances(best_model, X_train, feature_importances_path)

    # Save model
    best_model_path = os.path.join(output_dir, "best_model.pkl")
    with open(best_model_path, 'wb') as f:
        joblib.dump(best_model, f)
    print(f"Best model saved to {best_model_path}")

    # Save metadata sidecar (joint multi-output): target names, thresholds, etc.
    if metadata is not None:
        metadata_path = os.path.join(output_dir, "best_model_metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"Model metadata saved to {metadata_path}")

    # Calculate and save SHAP values if required (single-output only for now).
    if use_shap and not is_multi_output:
        process_shap_values(best_model, X_train, X_train_sample_ids, output_dir, binary_data)

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")

def process_shap_values(best_model, X_train, X_train_sample_ids, output_dir, binary_data):
    """
    Calculates and saves SHAP values, plots, and related files.

    Args:
        best_model: Trained CatBoost model.
        X_train (DataFrame): Training feature set.
        X_train_sample_ids (DataFrame): Sample IDs for the training set.
        output_dir (str): Directory to store SHAP results.
        binary_data (bool): If True, plot SHAP jitter plot with binary data.
    """
    explainer = shap.TreeExplainer(best_model, approximate=True)
    shap_values = explainer.shap_values(X_train)

    # Save SHAP values as CSV
    X_train_sample_ids = X_train_sample_ids.reset_index(drop=True)
    X_train_sample_ids = X_train_sample_ids.loc[:, ~X_train_sample_ids.columns.duplicated()]
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
    plot_custom_shap_beeswarm(shap_values_df, output_dir, binary_data=binary_data)

    del explainer, shap_values, shap_values_df, X_train_df
    gc.collect()


def process_multilabel_shap_values(joint_model, X, X_sample_ids, target_names,
                                   output_dir, binary_data=False):
    """Per-receptor SHAP for a joint MultiLogloss model.

    A joint model predicts all receptors from one shared ensemble, so a feature's
    attribution is decomposed PER receptor: receptor k's SHAP reflects what the
    joint model uses to predict k *given it also models the others*, which
    partials out shared/co-occurrence signal and surfaces receptor-SPECIFIC
    features (vs. independent SHAP, which includes correlated/confounded signal).

    Writes one ``<receptor>/shap_importances.csv`` per receptor, in the SAME
    long format as process_shap_values (columns: <sample_ids>, feature,
    shap_value, value), so the existing kmer-analysis / annotation pipeline runs
    per receptor unchanged. Also writes a combined mean|SHAP| ranking.

    Args:
        joint_model: trained CatBoostClassifier(loss_function='MultiLogloss').
        X (DataFrame): feature matrix (the model's feature set, correct column order).
        X_sample_ids (DataFrame): sample identifier columns to carry through.
        target_names (list[str]): receptor names, in the model's label order.
        output_dir (str): base dir; per-receptor subdirs are created.
    """
    from catboost import Pool

    X = X.reset_index(drop=True)
    X_sample_ids = X_sample_ids.reset_index(drop=True)
    X_sample_ids = X_sample_ids.loc[:, ~X_sample_ids.columns.duplicated()]
    id_cols = list(X_sample_ids.columns)
    feature_cols = list(X.columns)

    # CatBoost native per-output SHAP: shape (n_samples, n_labels, n_features + 1)
    # (last column per label is the base/expected-value term).
    shap = np.asarray(joint_model.get_feature_importance(type='ShapValues', data=Pool(X)))
    if shap.ndim != 3:
        raise ValueError(
            f"Expected a 3-D multi-output SHAP tensor (samples, labels, features+1); "
            f"got shape {shap.shape}. Is this a joint MultiLogloss model?")
    n_labels = shap.shape[1]
    if n_labels != len(target_names):
        logging.warning(
            f"SHAP has {n_labels} labels but {len(target_names)} target_names; "
            f"using the model's label order.")

    # Long-format feature values, shared across receptors (melt once).
    X_val = X.copy()
    X_val[id_cols] = X_sample_ids
    val_long = X_val.melt(id_vars=id_cols, var_name='feature', value_name='value')

    combined = []
    for k, receptor in enumerate(target_names[:n_labels]):
        rec_dir = os.path.join(output_dir, receptor)
        os.makedirs(rec_dir, exist_ok=True)

        sv = shap[:, k, :-1]   # (n_samples, n_features) -- drop base term
        sv_df = pd.DataFrame(sv, columns=feature_cols)
        sv_df[id_cols] = X_sample_ids
        sv_long = sv_df.melt(id_vars=id_cols, var_name='feature', value_name='shap_value')
        sv_long = sv_long.merge(val_long, on=id_cols + ['feature'], how='left')

        out_csv = os.path.join(rec_dir, 'shap_importances.csv')
        sv_long.to_csv(out_csv, index=False)
        logging.info(f"[SHAP] {receptor}: per-sample SHAP -> {out_csv}")

        # mean|SHAP| ranking for this receptor (the interpretable summary).
        rank = (pd.Series(np.abs(sv).mean(axis=0), index=feature_cols)
                .sort_values(ascending=False))
        rank.to_csv(os.path.join(rec_dir, 'mean_abs_shap_ranking.csv'),
                    header=['mean_abs_shap'])
        top = rank.head(1)
        combined.append(pd.DataFrame({'receptor': receptor,
                                      'feature': rank.index,
                                      'mean_abs_shap': rank.values}))

    # One combined long table: receptor x feature x mean|SHAP| (for cross-receptor views).
    combined_df = pd.concat(combined, ignore_index=True)
    combined_df.to_csv(os.path.join(output_dir, 'per_receptor_mean_abs_shap.csv'), index=False)
    logging.info(f"[SHAP] combined per-receptor ranking -> "
                 f"{os.path.join(output_dir, 'per_receptor_mean_abs_shap.csv')}")
    gc.collect()
    return combined_df

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

def parse_model_predictions_and_performance(model_dir, task_type='classification'):
    """
    Parses the prediction and performance data from the run directories within the specified model directory.
    
    Args:
        model_dir (str): Base directory for the select feature modeling containing cutoff subdirectories.
        task_type (str): Either 'classification' or 'regression', specifies the type of model.

    Returns:
        model_predictions_df (DataFrame): Combined DataFrame of all model predictions.
        model_performance_df (DataFrame): Combined DataFrame of model performance (top models summary).
    """
    model_predictions_df = pd.DataFrame()
    model_performance_df = pd.DataFrame()

    # Define the performance metric column based on task type. The actual
    # metric is detected per top_models_summary below (joint multi-label runs
    # write macro_f1 instead of mcc), so this is just the default expectation.
    if task_type == 'classification':
        metric_column = 'mcc'  # For binary classification, MCC.
    elif task_type == 'regression':
        metric_column = 'r2'  # For regression, R².
    else:
        raise ValueError("Invalid task_type. Must be 'classification' or 'regression'.")

    # Get all cutoff subdirectories (excluding CSV files)
    cut_offs = sorted([x for x in os.listdir(model_dir) if '.csv' not in x])

    # Loop through each cutoff directory and parse run directories
    for cut_off in cut_offs:
        print(f"Parsing cut-off: {cut_off}")
        cut_off_dir = os.path.join(model_dir, cut_off)
        run_dirs = sorted([x for x in os.listdir(cut_off_dir) if 'run' in x])

        # Parse predictions from each run
        predictions_temp = None
        for run in run_dirs:
            predictions_file = os.path.join(cut_off_dir, run, 'best_model_predictions.csv')
            if os.path.exists(predictions_file):
                predictions_temp = pd.read_csv(predictions_file)
                predictions_temp['run'] = run
                predictions_temp['cut_off'] = cut_off
                model_predictions_df = pd.concat([model_predictions_df, predictions_temp])

        # Parse top models summary for this cutoff
        top_models_temp = None
        top_models_file = os.path.join(cut_off_dir, 'top_models_summary.csv')
        if os.path.exists(top_models_file):
            top_models_temp = pd.read_csv(top_models_file)
            # Joint multi-label -> macro_f1; joint multi-target -> mean_r2;
            # prefer those when present.
            if 'macro_f1' in top_models_temp.columns:
                this_metric = 'macro_f1'
            elif 'mean_r2' in top_models_temp.columns:
                this_metric = 'mean_r2'
            else:
                this_metric = metric_column
            if this_metric not in top_models_temp.columns:
                raise ValueError(f"Expected metric '{this_metric}' not found in top_models_summary.csv")
            top_models_temp = top_models_temp[[this_metric]].reset_index()
            top_models_temp['index'] = ['_'.join(['run', str(x)]) for x in top_models_temp['index']]
            top_models_temp = top_models_temp.rename(columns={'index': 'run', this_metric: 'performance_metric'})
            top_models_temp['cut_off'] = cut_off
            model_performance_df = pd.concat([model_performance_df, top_models_temp])

    # Clean up temporary variables if they were created
    if 'predictions_temp' in locals() and predictions_temp is not None:
        del predictions_temp
    if 'top_models_temp' in locals() and top_models_temp is not None:
        del top_models_temp
    gc.collect()

    return model_predictions_df, model_performance_df

def evaluate_model_performance(predictions_file, output_dir, sample_column='strain', phenotype_column='interaction', task_type='classification'):
    """
    Evaluates model performances and generates performance plots grouped by cutoff for classification or regression.

    Args:
        predictions_file (str): Path to the CSV file containing model predictions.
        output_dir (str): Directory to save performance plots and evaluation metrics.
        sample_column (str): Column name for the sample identifier.
        phenotype_column (str): Column name for the phenotype.
        task_type (str): Either 'classification' or 'regression'.
    """
    model_performance_dir = os.path.join(output_dir, 'model_performance')
    os.makedirs(model_performance_dir, exist_ok=True)

    model_predictions_df_full = pd.read_csv(predictions_file)
    model_predictions_df_full['cut_off'] = model_predictions_df_full['cut_off'].astype(str)

    # Joint multi-output: predictions carry per-target columns rather than a
    # single trio. Multi-label (classification) has Confidence_<t>; multi-target
    # (regression) has Prediction_<t>/True_<t> only.
    if isinstance(phenotype_column, (list, tuple)):
        target_list = list(phenotype_column)
        if task_type == 'regression':
            evaluate_multitarget_performance(
                model_predictions_df_full, model_performance_dir,
                sample_column, target_list
            )
        else:
            evaluate_multilabel_performance(
                model_predictions_df_full, model_performance_dir,
                sample_column, target_list
            )
        del model_predictions_df_full
        gc.collect()
        return

    grouping_columns = ['cut_off', sample_column, phenotype_column]
    if 'phage' in model_predictions_df_full.columns and 'phage' not in grouping_columns:
        grouping_columns.insert(2, 'phage')

    if task_type == 'classification':
        evaluate_classifier_performance(model_predictions_df_full, model_performance_dir, grouping_columns, phenotype_column)
    elif task_type == 'regression':
        evaluate_regressor_performance(model_predictions_df_full, model_performance_dir, grouping_columns, phenotype_column)
    else:
        raise ValueError("Invalid task_type. Must be 'classification' or 'regression'.")

    del model_predictions_df_full
    gc.collect()

def evaluate_multilabel_performance(df, model_performance_dir, sample_column, target_names):
    """
    Evaluate joint multi-label predictions across runs and cutoffs.

    Predictions carry Prediction_<label>, Confidence_<label>, and True_<label>
    columns. For each (cut_off, sample), confidence is medianed across runs;
    per-label predictions are then thresholded at 0.5 of the median confidence.
    Per-label and aggregate (micro/macro-F1, Hamming, subset-accuracy) metrics
    are written, one row per cutoff, sorted by macro_f1.

    Args:
        df (DataFrame): combined predictions across runs/cutoffs.
        model_performance_dir (str): output directory for metrics.
        sample_column (str): sample identifier column.
        target_names (list[str]): the label column base names.
    """
    group_keys = ['cut_off', sample_column]
    if 'phage' in df.columns and 'phage' not in group_keys:
        group_keys.append('phage')

    # Only score targets whose columns are actually present (a target dropped in
    # every fold -- e.g. too sparse to ever model -- won't have columns here).
    present = [t for t in target_names
               if f'Confidence_{t}' in df.columns and f'True_{t}' in df.columns]
    missing = [t for t in target_names if t not in present]
    if missing:
        logging.warning(f"Targets absent from all pooled predictions; not scored: {missing}")
    conf_cols = [f'Confidence_{t}' for t in present]
    true_cols = [f'True_{t}' for t in present]

    # Median confidence per (cut_off, sample); true labels are constant per sample.
    agg_map = {c: 'median' for c in conf_cols}
    agg_map.update({c: 'first' for c in true_cols})
    df_calcs = df.groupby(group_keys).agg(agg_map).reset_index()

    metrics_list = []
    for cut_off, group in df_calcs.groupby('cut_off'):
        row = {'cut_off': cut_off}
        # Per-receptor metrics, each computed over only the samples where that
        # receptor was predicted (a fold that dropped it as single-class leaves
        # NaN for its held-out rows -> excluded from that receptor's metric, and
        # n_folds_<t> records the coverage). AUC/AUPR/MCC need both classes.
        valid_masks = {}
        for t in present:
            pr = group[f'Confidence_{t}'].to_numpy()
            yt_raw = group[f'True_{t}'].to_numpy()
            mask = ~(np.isnan(pr) | pd.isna(yt_raw))
            valid_masks[t] = mask
            if mask.sum() == 0:
                continue
            yk = yt_raw[mask].astype(int)
            pk = (pr[mask] >= 0.5).astype(int)
            prk = pr[mask]
            both = len(np.unique(yk)) > 1
            row[f'AUC_{t}'] = roc_auc_score(yk, prk) if both else float('nan')
            row[f'AUPR_{t}'] = average_precision_score(yk, prk) if both else float('nan')
            row[f'Accuracy_{t}'] = accuracy_score(yk, pk)
            row[f'Precision_{t}'] = precision_score(yk, pk, zero_division=0)
            row[f'Recall_{t}'] = recall_score(yk, pk, zero_division=0)
            row[f'F1_{t}'] = f1_score(yk, pk, zero_division=0)
            row[f'MCC_{t}'] = matthews_corrcoef(yk, pk) if both else float('nan')
            row[f'n_samples_{t}'] = int(mask.sum())
        # Aggregate (whole-matrix) metrics over samples where ALL present targets
        # were predicted (complete rows), so the matrix is well-defined.
        complete = np.all(np.column_stack([valid_masks[t] for t in present]), axis=1) if present else np.array([])
        if complete.sum() > 0:
            yt_m = group[true_cols].to_numpy()[complete].astype(int)
            pr_m = group[conf_cols].to_numpy()[complete]
            yp_m = (pr_m >= 0.5).astype(int)
            row['micro_f1'] = f1_score(yt_m, yp_m, average='micro', zero_division=0)
            row['macro_f1'] = f1_score(yt_m, yp_m, average='macro', zero_division=0)
            row['hamming_loss'] = float((yp_m != yt_m).mean())
            row['subset_accuracy'] = float((yp_m == yt_m).all(axis=1).mean())
            row['n_complete_samples'] = int(complete.sum())
        metrics_list.append(row)

    metrics_df = pd.DataFrame(metrics_list)
    try:
        metrics_df['cutoff_int'] = metrics_df['cut_off'].str.split('_').str[-1].astype(int)
        metrics_df = metrics_df.sort_values(['macro_f1', 'cutoff_int'], ascending=[False, True])
        metrics_df = metrics_df.drop('cutoff_int', axis=1)
    except (ValueError, TypeError):
        metrics_df = metrics_df.sort_values(['macro_f1'], ascending=[False])

    metrics_df.to_csv(os.path.join(model_performance_dir, 'model_performance_metrics.csv'), index=False)
    print(f"Multi-label performance metrics saved to "
          f"{os.path.join(model_performance_dir, 'model_performance_metrics.csv')}")

    del df_calcs, metrics_list, metrics_df
    gc.collect()


def _strong_responder_detection(yt, yp, top_frac=0.2):
    """Detection metrics for 'does the model rank strong-responder samples high'.

    For zero-inflated / heavy-tailed continuous targets (e.g. CRISPRi fitness),
    R2/Spearman are dominated by the un-rankable near-zero bulk and understate
    the model. The biological question is usually 'which samples have a STRONG
    effect', i.e. a top-K retrieval / detection problem. We label the top
    ``top_frac`` of samples by TRUE value as 'strong' and score how well the
    PREDICTED value ranks them.

    Args:
        yt (np.ndarray): true continuous values.
        yp (np.ndarray): predicted continuous values.
        top_frac (float): fraction of samples (by true value) treated as strong.

    Returns:
        dict with keys: n_strong, detect_AUPR, detect_AUROC, precision_at_k,
        recall_at_k (k = number of strong samples). NaN where undefined.
    """
    n = len(yt)
    k = max(1, int(round(n * top_frac)))
    # 'strong' = the k samples with the highest TRUE value.
    thresh = np.sort(yt)[-k]
    strong = (yt >= thresh).astype(int)
    n_strong = int(strong.sum())
    out = {'n_strong': n_strong, 'detect_AUPR': float('nan'),
           'detect_AUROC': float('nan'), 'precision_at_k': float('nan'),
           'recall_at_k': float('nan')}
    # Need both classes present and >1 strong for meaningful detection scoring.
    if n_strong == 0 or n_strong == n:
        return out
    out['detect_AUPR'] = float(average_precision_score(strong, yp))
    out['detect_AUROC'] = float(roc_auc_score(strong, yp))
    # Top-k by PREDICTED value; how many are truly strong.
    top_pred = np.argsort(yp)[::-1][:n_strong]
    n_hit = int(strong[top_pred].sum())
    out['precision_at_k'] = n_hit / n_strong
    out['recall_at_k'] = n_hit / n_strong  # equal when k == n_strong
    return out


def evaluate_multitarget_performance(df, model_performance_dir, sample_column,
                                     target_names, strong_top_frac=0.2):
    """
    Evaluate joint multi-target (MultiRMSE) predictions across runs and cutoffs.

    Predictions carry Prediction_<target> and True_<target> columns. For each
    (cut_off, sample), predictions are medianed across runs. Two metric families
    are written per target, one row per cutoff:

      - Regression fit: R2, RMSE, MAE, normalized-RMSE (+ aggregate mean_r2,
        mean_normalized_rmse). Sorted by mean_r2 descending.
      - Strong-responder DETECTION: detect_AUPR, detect_AUROC, precision@k for
        the top ``strong_top_frac`` of samples by true value (+ aggregate
        mean_detect_AUPR, mean_detect_AUROC). Appropriate for zero-inflated /
        heavy-tailed fitness targets where R2 understates the model.

    Args:
        df (DataFrame): combined predictions across runs/cutoffs.
        model_performance_dir (str): output directory for metrics.
        sample_column (str): sample identifier column.
        target_names (list[str]): the target column base names.
        strong_top_frac (float): fraction of samples (by true value) treated as
            'strong responders' for the detection metrics (default 0.2).
    """
    group_keys = ['cut_off', sample_column]
    if 'phage' in df.columns and 'phage' not in group_keys:
        group_keys.append('phage')

    # Only score targets present in the pooled predictions (a target dropped in
    # every fold won't have columns here).
    present = [t for t in target_names
               if f'Prediction_{t}' in df.columns and f'True_{t}' in df.columns]
    missing = [t for t in target_names if t not in present]
    if missing:
        logging.warning(f"Targets absent from all pooled predictions; not scored: {missing}")
    target_names = present
    pred_cols = [f'Prediction_{t}' for t in target_names]
    true_cols = [f'True_{t}' for t in target_names]

    agg_map = {c: 'median' for c in pred_cols}
    agg_map.update({c: 'first' for c in true_cols})
    df_calcs = df.groupby(group_keys).agg(agg_map).reset_index()

    metrics_list = []
    for cut_off, group in df_calcs.groupby('cut_off'):
        row = {'cut_off': cut_off}
        r2s, nrmses, dauprs, daurocs = [], [], [], []
        for k, t in enumerate(target_names):
            yt = group[true_cols[k]].to_numpy()
            yp = group[pred_cols[k]].to_numpy()
            rmse = float(np.sqrt(mean_squared_error(yt, yp)))
            mae = float(np.mean(np.abs(yt - yp)))
            r2 = float(r2_score(yt, yp)) if len(np.unique(yt)) > 1 else float('nan')
            spread = float(yt.max() - yt.min())
            nrmse = rmse / spread if spread > 0 else float('nan')
            row[f'R2_{t}'] = r2
            row[f'RMSE_{t}'] = rmse
            row[f'MAE_{t}'] = mae
            row[f'normalized_RMSE_{t}'] = nrmse
            if not np.isnan(r2):
                r2s.append(r2)
            if not np.isnan(nrmse):
                nrmses.append(nrmse)
            # Strong-responder detection (top-K retrieval) metrics.
            det = _strong_responder_detection(yt, yp, top_frac=strong_top_frac)
            row[f'detect_AUPR_{t}'] = det['detect_AUPR']
            row[f'detect_AUROC_{t}'] = det['detect_AUROC']
            row[f'precision_at_k_{t}'] = det['precision_at_k']
            row[f'n_strong_{t}'] = det['n_strong']
            if not np.isnan(det['detect_AUPR']):
                dauprs.append(det['detect_AUPR'])
            if not np.isnan(det['detect_AUROC']):
                daurocs.append(det['detect_AUROC'])
        row['mean_r2'] = float(np.mean(r2s)) if r2s else float('nan')
        row['mean_normalized_rmse'] = float(np.mean(nrmses)) if nrmses else float('nan')
        row['mean_detect_AUPR'] = float(np.mean(dauprs)) if dauprs else float('nan')
        row['mean_detect_AUROC'] = float(np.mean(daurocs)) if daurocs else float('nan')
        metrics_list.append(row)

    metrics_df = pd.DataFrame(metrics_list)
    try:
        metrics_df['cutoff_int'] = metrics_df['cut_off'].str.split('_').str[-1].astype(int)
        metrics_df = metrics_df.sort_values(['mean_r2', 'cutoff_int'], ascending=[False, True])
        metrics_df = metrics_df.drop('cutoff_int', axis=1)
    except (ValueError, TypeError):
        metrics_df = metrics_df.sort_values(['mean_r2'], ascending=[False])

    metrics_df.to_csv(os.path.join(model_performance_dir, 'model_performance_metrics.csv'), index=False)
    print(f"Multi-target performance metrics saved to "
          f"{os.path.join(model_performance_dir, 'model_performance_metrics.csv')}")

    del df_calcs, metrics_list, metrics_df
    gc.collect()


def evaluate_classifier_performance(df, model_performance_dir, grouping_columns, phenotype_column):
    """
    Evaluates classifier performance and generates performance plots including ROC, PR, hit rate, and hit ratio.

    Args:
        df (DataFrame): DataFrame containing model predictions and performance data.
        model_performance_dir (str): Directory to save performance plots.
        grouping_columns (list): Columns for grouping the data.
        phenotype_column (str): Column name for the phenotype.
    """
    # Calculate average prediction and confidence per group
    df_calcs = df.groupby(grouping_columns).agg({'Prediction': 'median', 'Confidence': 'median'}).reset_index()
    df_calcs['Prediction'] = (df_calcs['Confidence'] > 0.5).astype(int)

    def calculate_metrics(df):
        y_true = df[phenotype_column]
        y_pred_prob = df['Confidence']
        y_pred = df['Prediction']
        
        metrics = {
            'AUC': roc_auc_score(y_true, y_pred_prob),
            'AUPR': average_precision_score(y_true, y_pred_prob),
            'Accuracy': accuracy_score(y_true, y_pred),
            'Precision': precision_score(y_true, y_pred),
            'Recall': recall_score(y_true, y_pred),
            'F1': f1_score(y_true, y_pred),
            'MCC': matthews_corrcoef(y_true, y_pred)
        }
        
        fpr, tpr, roc_thresholds = roc_curve(y_true, y_pred_prob)
        roc_df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'cut_off': df['cut_off'].iloc[0]})
        
        precision, recall, pr_thresholds = precision_recall_curve(y_true, y_pred_prob)
        pr_df = pd.DataFrame({'precision': precision, 'recall': recall, 'cut_off': df['cut_off'].iloc[0]})
        
        return metrics, roc_df, pr_df

    metrics_list, roc_list, pr_list = [], [], []
    for cut_off, group in df_calcs.groupby('cut_off'):
        metrics, roc_df, pr_df = calculate_metrics(group)
        metrics['cut_off'] = cut_off
        metrics_list.append(metrics)
        roc_list.append(roc_df)
        pr_list.append(pr_df)

    metrics_df = pd.DataFrame(metrics_list)
    try:
        metrics_df['cutoff_int'] = metrics_df['cut_off'].str.split('_').str[-1].astype(int)
        metrics_df = metrics_df.sort_values(['MCC', 'cutoff_int'], ascending=[False, True])
        metrics_df = metrics_df.drop('cutoff_int', axis=1)
    except (ValueError, TypeError):
        # No numeric cutoff available, just sort by MCC
        metrics_df = metrics_df.sort_values(['MCC'], ascending=[False])
    roc_df = pd.concat(roc_list).reset_index(drop=True)
    pr_df = pd.concat(pr_list).reset_index(drop=True)

    plot_classifier_performance(roc_df, pr_df, model_performance_dir)
    metrics_df.to_csv(os.path.join(model_performance_dir, 'model_performance_metrics.csv'), index=False)

    del metrics_list, roc_list, pr_list, metrics_df, roc_df, pr_df
    gc.collect()

def evaluate_regressor_performance(df, model_performance_dir, grouping_columns, phenotype_column):
    """
    Evaluates regressor performance and generates comparison plots for predicted vs. actual values.

    Args:
        df (DataFrame): DataFrame containing model predictions.
        model_performance_dir (str): Directory to save performance plots.
        grouping_columns (list): Columns for grouping the data.
        phenotype_column (str): Column name for the phenotype.
    """
    df_calcs = df.groupby(grouping_columns)[['Prediction']].mean().reset_index()
    metrics_list, comparison_list = [], []

    for cut_off, group in df_calcs.groupby('cut_off'):
        y_true = group[phenotype_column]
        y_pred = group['Prediction']
        
        metrics = {
            'MSE': mean_squared_error(y_true, y_pred),
            'r2': r2_score(y_true, y_pred),
            'cut_off': cut_off
        }
        metrics_list.append(metrics)
        
        comparison_df = pd.DataFrame({'y_true': y_true, 'y_pred': y_pred, 'cut_off': cut_off})
        comparison_list.append(comparison_df)

    metrics_df = pd.DataFrame(metrics_list)
    try:
        metrics_df['cutoff_int'] = metrics_df['cut_off'].str.split('_').str[-1].astype(int)
        metrics_df = metrics_df.sort_values(['r2', 'cutoff_int'], ascending=[False, True])
        metrics_df = metrics_df.drop('cutoff_int', axis=1)
    except (ValueError, TypeError):
        metrics_df = metrics_df.sort_values(['r2'], ascending=[False])
    comparison_df = pd.concat(comparison_list).reset_index(drop=True)

    plot_regressor_performance(comparison_df, model_performance_dir)
    metrics_df.to_csv(os.path.join(model_performance_dir, 'model_performance_metrics.csv'), index=False)

    del metrics_list, comparison_list, metrics_df, comparison_df
    gc.collect()

def plot_classifier_performance(roc_df, pr_df, model_performance_dir):
    """
    Plots and saves classifier performance graphs including ROC and PR curves.
    
    Args:
        roc_df (DataFrame): DataFrame with ROC data for each cutoff.
        pr_df (DataFrame): DataFrame with PR data for each cutoff.
        model_performance_dir (str): Directory to save plots.
    """
    roc_curve_plot = (
        ggplot(roc_df, aes(x='fpr', y='tpr', color='cut_off')) +
        geom_line() +
        geom_abline(linetype='dashed') +
        labs(x='False Positive Rate', y='True Positive Rate', color='Cut Off') +
        theme(figure_size=(5, 4))
    )
    roc_curve_plot.save(os.path.join(model_performance_dir, 'roc_curve.png'))

    pr_curve_plot = (
        ggplot(pr_df, aes(x='recall', y='precision', color='cut_off')) +
        geom_line() +
        labs(x='Recall', y='Precision', color='Cut Off') +
        theme(figure_size=(5, 4))
    )
    pr_curve_plot.save(os.path.join(model_performance_dir, 'pr_curve.png'))
    del roc_curve_plot, pr_curve_plot
    gc.collect()

def plot_regressor_performance(comparison_df, model_performance_dir):
    """
    Plots and saves predicted vs. actual comparison plots for each cutoff.

    Args:
        comparison_df (DataFrame): DataFrame containing true and predicted values for each cutoff.
        model_performance_dir (str): Directory to save plots.
    """
    regressor_plot = (
        ggplot(comparison_df, aes(x='y_true', y='y_pred', color='cut_off')) +
        geom_point(alpha=0.5) +
        geom_abline(slope=1, intercept=0, linetype='dashed', color='black') +
        geom_smooth(method='lm', se=False) +
        labs(x='Actual Value', y='Predicted Value', color='Cut Off') +
        theme(figure_size=(6, 6))
    )
    regressor_plot.save(os.path.join(model_performance_dir, 'predicted_vs_actual_comparison.png'))
    del regressor_plot
    gc.collect()

def run_experiments(
    input_dir, 
    base_output_dir, 
    threads, 
    num_runs, 
    task_type='classification', 
    set_filter='none', 
    sample_column='strain', 
    phenotype_column='interaction', 
    phage_column='phage',
    use_dynamic_weights=False,
    weights_method='log10',
    binary_data=False, 
    max_ram=8, 
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    use_shap=False,
    target_mode='auto',
    strategy='joint'
):
    """
    Iterates through feature tables in a directory, running the model testing process for each.

        Args:
        input_dir (str): Directory containing feature tables.
        base_output_dir (str): Base directory to store results.
        threads (int): Number of threads for training.
        num_runs (int): Number of runs to perform per table.
        task_type (str): Specifies 'classification' or 'regression' task.
        set_filter (str): Filter type for the dataset ('none', 'strain', 'phage', 'dataset').
        sample_column (str): Name of the sample column (optional).
        phenotype_column (str): Name of the phenotype column (optional).
        phage_column (str): Name of the phage column (default: 'phage').
        use_dynamic_weights (bool): Whether to use dynamic weights for phage-based samples.
        weights_method (str): Method for calculating weights ('log10', 'inverse_frequency', or 'balanced').
        binary_data (bool): If True, plot SHAP jitter plot with binary data.
        max_ram (int): Maximum RAM to use for CatBoost training.
        use_clustering (bool): Whether to use clustering for filtering.
        cluster_method (str): Clustering method to use ('hdbscan' or 'hierarchical').
        n_clusters (int): Number of clusters for hierarchical clustering (default: 20).
        min_cluster_size (int): Minimum cluster size for HDBSCAN.
        min_samples (int): Minimum number of samples for HDBSCAN (default: None for same as min_cluster_size).
        cluster_selection_epsilon (float): Epsilon value for HDBSCAN.
        use_shap (bool): If True, calculate and save SHAP values.
    """
    start_total_time = time.time()
    if os.path.isdir(input_dir):
        feature_tables = sorted(os.listdir(input_dir))
    else:
        feature_tables = [os.path.basename(input_dir)]

    for feature_table in feature_tables:
        if os.path.isdir(input_dir):
            feature_table_path = os.path.join(input_dir, feature_table)
        else:
            feature_table_path = input_dir
        
        # Use the extract_cutoff_from_filename to generate the directory name
        if 'cutoff' in feature_table:
            cutoff_value = extract_cutoff_from_filename(feature_table)
            model_output_dir = os.path.join(base_output_dir, f'cutoff_{cutoff_value}')
        else:
            cutoff_value = feature_table.split('/')[-1].split('.')[0]
            model_output_dir = os.path.join(base_output_dir, cutoff_value)
        
        print(cutoff_value)
        top_models_summary_path = os.path.join(model_output_dir, 'top_models_summary.csv')

        if not os.path.exists(top_models_summary_path):
            if not os.path.exists(model_output_dir):
                os.makedirs(model_output_dir)

            top_models_df = pd.DataFrame()
            top_models_shap_df = pd.DataFrame()

            for i in tqdm(range(num_runs), desc=f"Running Experiments for cutoff {cutoff_value}"):
                output_dir = os.path.join(model_output_dir, f'run_{i}')
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                random_state = i
                # Set global random seeds for full reproducibility
                np.random.seed(random_state)

                model_performance_path = os.path.join(output_dir, 'model_performance.csv')
                if not os.path.exists(model_performance_path):
                    model_testing_select_MCC(
                        input=feature_table_path,
                        output_dir=output_dir,
                        threads=threads,
                        random_state=random_state,
                        task_type=task_type,
                        set_filter=set_filter,
                        sample_column=sample_column,
                        phenotype_column=phenotype_column,
                        phage_column=phage_column,
                        use_dynamic_weights=use_dynamic_weights,
                        weights_method=weights_method,
                        binary_data=binary_data,
                        max_ram=max_ram,
                        use_clustering=use_clustering,
                        cluster_method=cluster_method,
                        n_clusters=n_clusters,
                        min_cluster_size=min_cluster_size,
                        min_samples=min_samples,
                        cluster_selection_epsilon=cluster_selection_epsilon,
                        use_shap=use_shap,
                        target_mode=target_mode,
                        strategy=strategy
                    )
                else:
                    logging.info(f"Model performance already saved to {model_performance_path}")

                # Load model performance data for the best model selection
                if os.path.exists(model_performance_path):
                    run_results = pd.read_csv(model_performance_path)
                    # Select the best metric based on what the run actually wrote.
                    # Joint multi-label -> macro_f1; joint multi-target -> mean_r2;
                    # binary -> mcc; single regression -> r2.
                    if 'macro_f1' in run_results.columns:
                        top_model = run_results.nlargest(1, 'macro_f1')
                    elif 'mean_r2' in run_results.columns:
                        top_model = run_results.nlargest(1, 'mean_r2')
                    elif task_type == 'classification':
                        top_model = run_results.nlargest(1, 'mcc')
                    elif task_type == 'regression':
                        top_model = run_results.nlargest(1, 'r2')
                    top_models_df = pd.concat([top_models_df, top_model])

                # Load SHAP values for each run if available and use_shap is True
                if use_shap:
                    shap_values_csv_path = os.path.join(output_dir, "shap_importances.csv")
                    if os.path.exists(shap_values_csv_path):
                        shap_values_temp = pd.read_csv(shap_values_csv_path)
                        shap_values_temp = shap_values_temp.groupby(['feature', 'value']).agg({'shap_value': 'median'}).reset_index()
                        top_models_shap_df = pd.concat([top_models_shap_df, shap_values_temp])

            # Generate SHAP summary plot for the top models if use_shap is True
            if use_shap:
                model_performance_dir = os.path.join(base_output_dir, 'model_performance')
                os.makedirs(model_performance_dir, exist_ok=True)

                plot_custom_shap_beeswarm(top_models_shap_df, model_performance_dir, prefix=f'cutoff_{cutoff_value}', binary_data=binary_data)

            # Save top models summary for each cutoff
            top_models_df.to_csv(top_models_summary_path, index=False)
            print(f"Top models saved to {top_models_summary_path}")

            del top_models_df, top_models_shap_df
            gc.collect()
        else:
            logging.info(f"Top models summary already saved to {top_models_summary_path}")

    # After running all experiments, parse the predictions and performance
    model_predictions_output = os.path.join(base_output_dir, 'select_features_model_predictions.csv')
    model_performance_output = os.path.join(base_output_dir, 'select_features_model_performance.csv')

    if not os.path.exists(model_performance_output):
        logging.info("Parsing model predictions and performance data...")
        # Parse all predictions from the run directories
        model_predictions_df, model_performance_df = parse_model_predictions_and_performance(base_output_dir, task_type)

        # Save parsed predictions and performance data
        model_predictions_df.to_csv(model_predictions_output, index=False)
        model_performance_df.to_csv(model_performance_output, index=False)
        
        logging.info(f"Model predictions saved to {model_predictions_output}")
        logging.info(f"Model performance saved to {model_performance_output}")

        del model_predictions_df, model_performance_df
        gc.collect()
    else:
        logging.info(f"Model predictions and performance already saved to {model_predictions_output} and {model_performance_output}")

    model_metrics_path = os.path.join(base_output_dir, 'model_performance', 'model_performance_metrics.csv')
    if not os.path.exists(model_metrics_path):
        # Check if we have any predictions to evaluate (predictions file may be empty if all models failed)
        try:
            pred_check_df = pd.read_csv(model_predictions_output)
            if len(pred_check_df) == 0:
                logging.warning("No model predictions available. Skipping performance evaluation.")
            else:
                # Now evaluate model performance and generate performance plots, depending on task type
                evaluate_model_performance(
                    predictions_file=model_predictions_output,
                    output_dir=base_output_dir,
                    sample_column=sample_column,
                    phenotype_column=phenotype_column,
                    task_type=task_type
                )
        except pd.errors.EmptyDataError:
            logging.warning("No model predictions available. Skipping performance evaluation.")
    else:
        logging.info(f"Model performance metrics already saved to {model_metrics_path}")

    end_total_time = time.time()
    print(f"All experiments completed in {end_total_time - start_total_time:.2f} seconds.")
