"""
Nested k-fold cross-validation workflow for GenoPHI.

Runs deterministic (repeated) k-fold cross-validation entirely in-process:
for each fold, the strains are split into a modeling set and a held-out
validation set, a model is trained on the modeling set via
``run_protein_family_workflow``, the best cutoff (highest MCC) is selected, and
interactions for the held-out strains are predicted via
``assign_predict_workflow``. Per-fold predictions are aggregated into a single
predictions table plus summary statistics.

This is the portable, scheduler-agnostic equivalent of the SLURM job-array
submission script: folds run sequentially in this process rather than as
parallel SLURM array tasks. Only deterministic k-fold splitting is supported
(no random/bootstrap resampling).
"""

import os
import random
import logging

import pandas as pd

from genophi.workflows.protein_family_workflow import run_protein_family_workflow
from genophi.workflows.assign_predict_workflow import assign_predict_workflow

logger = logging.getLogger(__name__)

# Column names written by the prediction workflow into
# strain_median_predictions.csv (see prediction_workflow.calculate_median_predictions).
_SCORE_COL = 'Confidence'        # median predicted probability of the positive class
_PRED_COL = 'Final_Prediction'   # thresholded 0/1 call


def _get_full_strain_list(interaction_matrix, input_strain_dir, strain_column, suffix='faa'):
    """Return strains present in BOTH the interaction matrix and the strain directory.

    Mirrors the intersection logic used by the SLURM bootstrap workflow so that
    only strains we have both a phenotype and a FASTA file for are used.
    """
    logger.info(f"Reading interaction matrix: {interaction_matrix}")
    interaction_df = pd.read_csv(interaction_matrix)

    if strain_column not in interaction_df.columns:
        raise ValueError(
            f"Column '{strain_column}' not found in interaction matrix. "
            f"Available columns: {list(interaction_df.columns)}"
        )

    strains_in_matrix = {str(s) for s in interaction_df[strain_column].unique()}
    logger.info(f"Found {len(strains_in_matrix)} unique strains in interaction matrix")

    dot_suffix = f".{suffix}"
    strain_files = [f for f in os.listdir(input_strain_dir) if f.endswith(dot_suffix)]
    strains_in_dir = {f[: -len(dot_suffix)] for f in strain_files}
    logger.info(f"Found {len(strains_in_dir)} strain files in {input_strain_dir}")

    full_strain_list = sorted(strains_in_matrix & strains_in_dir)
    logger.info(
        f"Intersection: {len(full_strain_list)} strains found in both matrix and directory"
    )
    if not full_strain_list:
        logger.error(
            "No strains found in both the interaction matrix and the input directory. "
            "This is often caused by naming mismatches between the files and the matrix."
        )
    return full_strain_list


def split_strains_kfold(full_strain_list, iteration, n_folds=10):
    """Deterministically split strains into modeling/validation sets for one fold.

    ``iteration`` is 1-based. Iterations 1..n_folds form round 0, the next
    n_folds form round 1, and so on. Within a round the strains are shuffled
    with a per-round seed (the round index), then partitioned into ``n_folds``
    roughly equal folds with the remainder distributed across the leading
    folds. The fold for this iteration becomes the validation set; the rest is
    the modeling set.

    Returns ``(modeling_strains, validation_strains)``.
    """
    if not full_strain_list:
        raise ValueError("Cannot split an empty strain list.")

    cv_round = (iteration - 1) // n_folds
    fold_index = (iteration - 1) % n_folds

    logger.info(
        f"K-fold split: iteration {iteration} -> round {cv_round}, fold {fold_index} "
        f"({n_folds} folds)"
    )

    # Deterministic shuffle per round. Sort first so the input order does not
    # affect the result.
    random.seed(cv_round)
    strain_list_copy = sorted(full_strain_list)
    random.shuffle(strain_list_copy)

    fold_size = len(strain_list_copy) // n_folds
    remainder = len(strain_list_copy) % n_folds

    folds = []
    start = 0
    for i in range(n_folds):
        # Distribute the remainder across the leading folds.
        end = start + fold_size + (1 if i < remainder else 0)
        folds.append(strain_list_copy[start:end])
        start = end

    validation_strains = folds[fold_index]
    modeling_strains = [s for j, fold in enumerate(folds) for s in fold if j != fold_index]

    logger.info(
        f"Fold {fold_index}: {len(modeling_strains)} modeling strains, "
        f"{len(validation_strains)} validation strains"
    )
    return modeling_strains, validation_strains


def _select_best_cutoff(iteration_output_dir):
    """Return the cutoff with the highest MCC (tie-break on higher cut_off)."""
    metrics_file = os.path.join(
        iteration_output_dir,
        'modeling_results', 'model_performance', 'model_performance_metrics.csv',
    )
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    return metrics_df['cut_off'].values[0]


def _run_single_fold(
    iteration,
    full_strain_list,
    input_strain_dir,
    input_phage_dir,
    interaction_matrix,
    output_dir,
    n_folds,
    min_seq_id,
    coverage,
    sensitivity,
    threads,
    num_runs_fs,
    num_runs_modeling,
    use_dynamic_weights,
    weights_method,
    use_clustering,
    cluster_method,
    n_clusters,
    min_cluster_size,
    min_samples,
    cluster_selection_epsilon,
    check_feature_presence,
    filter_by_cluster_presence,
    min_cluster_presence,
    max_ram,
    use_feature_clustering,
    feature_cluster_method,
    feature_n_clusters,
    feature_min_cluster_presence,
    duplicate_all,
):
    """Run training + held-out prediction for a single fold.

    Returns the path to the fold's median-predictions CSV. Idempotent: if that
    file already exists the fold is skipped and the path is returned as-is.
    """
    iteration_output_dir = os.path.join(output_dir, f'iteration_{iteration}')
    median_predictions_file = os.path.join(
        iteration_output_dir, 'model_validation', 'predict_results',
        'strain_median_predictions.csv',
    )

    if os.path.exists(median_predictions_file):
        logger.info(f"Iteration {iteration} already complete, skipping.")
        return median_predictions_file

    logger.info(f"Starting iteration {iteration}...")
    os.makedirs(iteration_output_dir, exist_ok=True)
    modeling_tmp_dir = os.path.join(iteration_output_dir, 'tmp')

    # --- Strain split (deterministic k-fold) ---
    modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')

    if not os.path.exists(modeling_strains_path) or not os.path.exists(validation_strains_path):
        modeling_strains, validation_strains = split_strains_kfold(
            full_strain_list, iteration=iteration, n_folds=n_folds
        )
        if not modeling_strains or not validation_strains:
            raise RuntimeError(
                f"Empty strain split for iteration {iteration}: "
                f"{len(modeling_strains)} modeling, {len(validation_strains)} validation"
            )
        pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
        pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
    else:
        logger.info(f"Strain lists already exist for iteration {iteration}; loading.")

    # --- Step 1: train protein-family models on the modeling strains ---
    metrics_file = os.path.join(
        iteration_output_dir,
        'modeling_results', 'model_performance', 'model_performance_metrics.csv',
    )
    if not os.path.exists(metrics_file):
        logger.info(f"Iteration {iteration}: running protein family workflow...")
        run_protein_family_workflow(
            input_path_strain=input_strain_dir,
            input_path_phage=input_phage_dir,
            phenotype_matrix=interaction_matrix,
            output_dir=iteration_output_dir,
            tmp_dir=modeling_tmp_dir,
            min_seq_id=min_seq_id,
            coverage=coverage,
            sensitivity=sensitivity,
            threads=threads,
            phenotype_column='interaction',
            phage_column='phage',
            strain_list=modeling_strains_path,
            phage_list=interaction_matrix,
            filter_type='strain',
            num_runs_fs=num_runs_fs,
            num_runs_modeling=num_runs_modeling,
            use_dynamic_weights=use_dynamic_weights,
            weights_method=weights_method,
            use_clustering=use_clustering,
            cluster_method=cluster_method,
            n_clusters=n_clusters,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_epsilon=cluster_selection_epsilon,
            check_feature_presence=check_feature_presence,
            filter_by_cluster_presence=filter_by_cluster_presence,
            min_cluster_presence=min_cluster_presence,
            # Force bootstrapping=True so the clustering step resolves duplicate
            # protein IDs across ALL strains (modeling AND held-out validation),
            # writing every strain's "strain::protein"-prefixed FASTA into
            # modified_AAs/. Held-out prediction then assigns validation strains
            # from that directory against the matching cluster DB. With the
            # default (False), only modeling strains are prefixed and validation
            # assignment finds no matching files.
            bootstrapping=True,
            max_ram=max_ram,
            use_feature_clustering=use_feature_clustering,
            feature_cluster_method=feature_cluster_method,
            feature_n_clusters=feature_n_clusters,
            feature_min_cluster_presence=feature_min_cluster_presence,
        )
    else:
        logger.info(f"Iteration {iteration}: modeling results already exist; skipping training.")

    # --- Step 2: select the cutoff with the highest MCC ---
    best_cutoff = _select_best_cutoff(iteration_output_dir)
    model_dir = os.path.join(iteration_output_dir, 'modeling_results', f'{best_cutoff}')

    # --- Step 3: predict interactions for the held-out validation strains ---
    validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
    os.makedirs(validation_output_dir, exist_ok=True)
    validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')

    # Assign features to the held-out validation strains. Because training ran
    # with bootstrapping=True, the clustering step prefixed EVERY strain's
    # protein IDs to "strain::protein" and wrote them under modified_AAs/strain/
    # -- so the held-out validation strains are present there with IDs that
    # match the cluster database. Prefer that directory; fall back to the
    # original input directory only if it is missing (e.g. no duplicate IDs
    # existed, so no prefixing was needed and original IDs already match).
    modified_aa_dir = os.path.join(iteration_output_dir, 'strain', 'modified_AAs', 'strain')
    if os.path.isdir(modified_aa_dir) and os.listdir(modified_aa_dir):
        strain_input_dir = modified_aa_dir
    else:
        strain_input_dir = input_strain_dir
    logger.info(f"Iteration {iteration}: assigning validation features from {strain_input_dir}")

    select_feature_table = os.path.join(
        iteration_output_dir, 'feature_selection', 'filtered_feature_tables',
        f'select_feature_table_{best_cutoff}.csv',
    )

    logger.info(f"Iteration {iteration}: running prediction workflow...")
    assign_predict_workflow(
        input_dir=strain_input_dir,
        genome_list=validation_strains_path,
        mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
        clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
        feature_map=os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv'),
        tmp_dir=validation_tmp_dir,
        suffix='faa',
        model_dir=model_dir,
        feature_table=select_feature_table,
        phage_feature_table_path=os.path.join(
            iteration_output_dir, 'phage', 'features', 'feature_table.csv'
        ),
        output_dir=validation_output_dir,
        threads=threads,
        genome_type='strain',
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        duplicate_all=duplicate_all,
    )

    logger.info(f"Iteration {iteration} completed successfully.")
    return median_predictions_file


def _aggregate_results(output_dir, total_iterations):
    """Concatenate per-fold predictions and write summary statistics.

    Tolerant of missing folds: any iteration without a predictions file is
    logged and skipped. Returns the number of folds aggregated.
    """
    final_predictions = []
    completed = 0

    for i in range(1, total_iterations + 1):
        median_predictions_file = os.path.join(
            output_dir, f'iteration_{i}', 'model_validation', 'predict_results',
            'strain_median_predictions.csv',
        )
        if os.path.exists(median_predictions_file):
            df = pd.read_csv(median_predictions_file)
            df['iteration'] = i
            final_predictions.append(df)
            completed += 1
        else:
            logger.warning(f"Results missing for iteration {i}.")

    if not final_predictions:
        logger.error("No fold results found to aggregate.")
        return 0

    final_df = pd.concat(final_predictions, ignore_index=True)
    final_df.to_csv(os.path.join(output_dir, 'final_predictions.csv'), index=False)
    logger.info(
        f"Saved final_predictions.csv with {len(final_df)} predictions "
        f"from {completed} folds."
    )

    if {'strain', 'phage', 'prediction'}.issubset(final_df.columns):
        summary = (
            final_df.groupby(['strain', 'phage'])['prediction']
            .agg(['mean', 'std', 'count'])
            .round(4)
            .reset_index()
        )
        summary.columns = ['strain', 'phage', 'mean_prediction', 'std_prediction', 'n_iterations']
        summary.to_csv(os.path.join(output_dir, 'prediction_summary.csv'), index=False)
        logger.info(
            f"Saved prediction_summary.csv with {len(summary)} unique strain-phage pairs."
        )
    else:
        logger.warning(
            "Could not build prediction summary; expected columns "
            "'strain', 'phage', 'prediction' not all present."
        )

    return completed


def _compute_global_metrics(
    output_dir,
    interaction_matrix,
    strain_column='strain',
    phage_column='phage',
    phenotype_column='interaction',
):
    """Compute global performance metrics and PR/ROC curves over pooled outer-test folds.

    The per-fold prediction files contain a predicted score but not the ground
    truth, so the pooled predictions (final_predictions.csv) are merged back
    against the interaction matrix on (strain, phage) to recover true labels.
    Every held-out prediction across all folds is treated as one pooled test
    set (each strain is held out exactly once per round).

    Writes, into output_dir/performance/:
        - global_metrics.csv   (AUC, AP, Accuracy, Precision, Recall, F1, MCC, n)
        - roc_curve.png, pr_curve.png

    Returns the metrics dict, or None if metrics could not be computed.
    """
    from sklearn.metrics import (
        roc_curve, precision_recall_curve, roc_auc_score, average_precision_score,
        accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef,
    )
    import matplotlib
    matplotlib.use('Agg')  # headless / HPC-safe backend
    import matplotlib.pyplot as plt

    final_path = os.path.join(output_dir, 'final_predictions.csv')
    if not os.path.exists(final_path):
        logger.warning("No final_predictions.csv found; skipping global metrics.")
        return None

    preds = pd.read_csv(final_path)
    if _SCORE_COL not in preds.columns:
        logger.warning(
            f"Pooled predictions missing '{_SCORE_COL}' column; cannot compute "
            f"global metrics. Found columns: {list(preds.columns)}"
        )
        return None

    # Recover ground-truth labels by merging on (strain, phage).
    truth = pd.read_csv(interaction_matrix)
    required = {strain_column, phage_column, phenotype_column}
    if not required.issubset(truth.columns):
        logger.warning(
            f"Interaction matrix missing one of {required}; cannot recover true "
            f"labels for global metrics. Found: {list(truth.columns)}"
        )
        return None

    truth = truth[[strain_column, phage_column, phenotype_column]].copy()
    # Align identifier dtypes so the merge keys match (strain IDs are read as
    # strings from FASTA filenames elsewhere in the pipeline).
    for col in (strain_column, phage_column):
        preds[col] = preds[col].astype(str)
        truth[col] = truth[col].astype(str)

    merged = preds.merge(truth, on=[strain_column, phage_column], how='inner')
    merged = merged.dropna(subset=[phenotype_column, _SCORE_COL])
    if merged.empty:
        logger.warning(
            "No pooled predictions could be matched to ground-truth labels "
            "(check strain/phage naming between predictions and interaction matrix)."
        )
        return None

    n_matched, n_total = len(merged), len(preds)
    if n_matched < n_total:
        logger.warning(
            f"Matched {n_matched}/{n_total} pooled predictions to ground-truth "
            f"labels; unmatched pairs are excluded from global metrics."
        )

    y_true = merged[phenotype_column].astype(int)
    y_score = merged[_SCORE_COL].astype(float)
    if _PRED_COL in merged.columns:
        y_pred = merged[_PRED_COL].astype(int)
    else:
        y_pred = (y_score > 0.5).astype(int)

    perf_dir = os.path.join(output_dir, 'performance')
    os.makedirs(perf_dir, exist_ok=True)

    # AUC/AP need both classes present; guard so a degenerate pool doesn't crash.
    both_classes = y_true.nunique() > 1
    metrics = {
        'n_predictions': n_matched,
        'n_positive': int(y_true.sum()),
        'n_negative': int((1 - y_true).sum()),
        'AUC': roc_auc_score(y_true, y_score) if both_classes else float('nan'),
        'AP': average_precision_score(y_true, y_score) if both_classes else float('nan'),
        'Accuracy': accuracy_score(y_true, y_pred),
        'Precision': precision_score(y_true, y_pred, zero_division=0),
        'Recall': recall_score(y_true, y_pred, zero_division=0),
        'F1': f1_score(y_true, y_pred, zero_division=0),
        'MCC': matthews_corrcoef(y_true, y_pred),
    }
    pd.DataFrame([metrics]).to_csv(os.path.join(perf_dir, 'global_metrics.csv'), index=False)
    logger.info(
        f"Global metrics over {n_matched} pooled outer-test predictions: "
        f"AUC={metrics['AUC']:.3f}, AP={metrics['AP']:.3f}, MCC={metrics['MCC']:.3f}"
    )

    if not both_classes:
        logger.warning(
            "Pooled outer-test labels contain a single class; ROC/PR curves "
            "are not meaningful and were skipped."
        )
        return metrics

    # ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_score)
    plt.figure()
    plt.plot(fpr, tpr, label=f"ROC (AUC = {metrics['AUC']:.3f})")
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Outer-test ROC (pooled folds)')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(perf_dir, 'roc_curve.png'), dpi=150, bbox_inches='tight')
    plt.close()

    # Precision-recall curve
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    plt.figure()
    plt.step(recall, precision, where='post', label=f"PR (AP = {metrics['AP']:.3f})")
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Outer-test Precision-Recall (pooled folds)')
    plt.legend(loc='lower left')
    plt.savefig(os.path.join(perf_dir, 'pr_curve.png'), dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved ROC and PR curves to {perf_dir}")
    return metrics


def run_nested_cv_workflow(
    input_strain_dir,
    input_phage_dir,
    interaction_matrix,
    output_dir,
    n_folds=10,
    cv_rounds=1,
    strain_column='strain',
    suffix='faa',
    min_seq_id=0.4,
    coverage=0.8,
    sensitivity=7.5,
    threads=4,
    num_runs_fs=25,
    num_runs_modeling=50,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=False,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=2,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
    max_ram=40,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    duplicate_all=False,
    clear_tmp=False,
):
    """Run deterministic (repeated) k-fold nested cross-validation.

    Folds run sequentially in this process. The total number of folds trained
    is ``n_folds * cv_rounds`` (each round reshuffles the fold assignment with a
    new seed). A failing fold is logged and skipped; the run continues and
    aggregates whatever folds completed.

    After aggregation, global performance metrics and pooled outer-test PR/ROC
    curves are written to ``output_dir/performance/`` (global_metrics.csv,
    roc_curve.png, pr_curve.png). The pooled predictions are merged back against
    the interaction matrix to recover true labels.

    Args:
        input_strain_dir (str): Directory of strain FASTA files.
        input_phage_dir (str): Directory of phage FASTA files.
        interaction_matrix (str): Path to the interaction/phenotype matrix CSV.
        output_dir (str): Directory for all results (per-fold under iteration_N/).
        n_folds (int): Number of folds per round (default: 10).
        cv_rounds (int): Number of repeated k-fold rounds (default: 1).
        strain_column (str): Strain identifier column in the matrix (default: 'strain').
        suffix (str): FASTA file suffix for strain files (default: 'faa').
        min_seq_id (float): Minimum sequence identity for MMseqs2 (default: 0.4).
        coverage (float): Minimum coverage for MMseqs2 (default: 0.8).
        sensitivity (float): MMseqs2 sensitivity (default: 7.5).
        threads (int): Threads per fold (default: 4).
        clear_tmp (bool): Delete per-fold tmp dirs on success (default: False).
        (remaining args are passed through to run_protein_family_workflow /
        assign_predict_workflow; see those functions for details.)

    Returns:
        int: Number of folds successfully aggregated.
    """
    total_iterations = n_folds * cv_rounds

    logger.info("=== GenoPHI nested k-fold cross-validation ===")
    logger.info(f"Output directory: {output_dir}")
    logger.info(
        f"Folds: {n_folds} x {cv_rounds} round(s) = {total_iterations} total fits"
    )

    os.makedirs(output_dir, exist_ok=True)

    full_strain_list = _get_full_strain_list(
        interaction_matrix, input_strain_dir, strain_column, suffix=suffix
    )
    if not full_strain_list:
        raise RuntimeError(
            "No valid strains found (intersection of interaction matrix and "
            "strain directory is empty). Cannot run cross-validation."
        )
    if len(full_strain_list) < n_folds:
        raise ValueError(
            f"Cannot create {n_folds} folds from {len(full_strain_list)} strains."
        )

    succeeded, failed = [], []
    for iteration in range(1, total_iterations + 1):
        try:
            _run_single_fold(
                iteration=iteration,
                full_strain_list=full_strain_list,
                input_strain_dir=input_strain_dir,
                input_phage_dir=input_phage_dir,
                interaction_matrix=interaction_matrix,
                output_dir=output_dir,
                n_folds=n_folds,
                min_seq_id=min_seq_id,
                coverage=coverage,
                sensitivity=sensitivity,
                threads=threads,
                num_runs_fs=num_runs_fs,
                num_runs_modeling=num_runs_modeling,
                use_dynamic_weights=use_dynamic_weights,
                weights_method=weights_method,
                use_clustering=use_clustering,
                cluster_method=cluster_method,
                n_clusters=n_clusters,
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                cluster_selection_epsilon=cluster_selection_epsilon,
                check_feature_presence=check_feature_presence,
                filter_by_cluster_presence=filter_by_cluster_presence,
                min_cluster_presence=min_cluster_presence,
                max_ram=max_ram,
                use_feature_clustering=use_feature_clustering,
                feature_cluster_method=feature_cluster_method,
                feature_n_clusters=feature_n_clusters,
                feature_min_cluster_presence=feature_min_cluster_presence,
                duplicate_all=duplicate_all,
            )
            succeeded.append(iteration)
            if clear_tmp:
                _clear_fold_tmp(output_dir, iteration)
        except Exception as e:  # noqa: BLE001 - keep going across folds
            failed.append(iteration)
            logger.error(f"Iteration {iteration} failed: {e}", exc_info=True)

    logger.info(
        f"Folds complete: {len(succeeded)} succeeded, {len(failed)} failed."
    )
    if failed:
        logger.warning(f"Failed folds: {failed}")

    aggregated = _aggregate_results(output_dir, total_iterations)
    logger.info(f"Aggregated {aggregated} fold(s) into final predictions.")

    if aggregated > 0:
        try:
            _compute_global_metrics(
                output_dir=output_dir,
                interaction_matrix=interaction_matrix,
                strain_column=strain_column,
                phage_column='phage',
                phenotype_column='interaction',
            )
        except Exception as e:  # noqa: BLE001 - metrics are non-fatal
            logger.error(f"Failed to compute global metrics: {e}", exc_info=True)

    return aggregated


def _clear_fold_tmp(output_dir, iteration):
    """Remove the per-fold tmp directories after a successful fold.

    Only call this AFTER the fold's validation prediction has completed. The
    fold's MMseqs2 database (iteration_N/tmp/strain/mmseqs_db) is required to
    assign features to the held-out validation strains during prediction, so it
    must not be removed before _run_single_fold returns. Each fold builds its
    own mmseqs_db and nothing reuses it across folds, so clearing here is safe.
    """
    import shutil

    iteration_output_dir = os.path.join(output_dir, f'iteration_{iteration}')
    for tmp_path in (
        os.path.join(iteration_output_dir, 'tmp'),
        os.path.join(iteration_output_dir, 'model_validation', 'tmp'),
    ):
        if os.path.isdir(tmp_path):
            shutil.rmtree(tmp_path, ignore_errors=True)
            logger.info(f"Cleared tmp directory: {tmp_path}")
