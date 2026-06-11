"""
Cross-validation for GenoPHI multi-output models (phage -> target vector).

Thin orchestration over existing workflows. Per fold:
  1. split genomes into modeling / held-out (k-fold over genomes).
  2. run the FULL modeling workflow (``run_kmer_workflow``) on the MODELING
     genomes only -- so the feature table (k-mer co-occurrence collapsing),
     feature selection, ensemble modeling, and internal cutoff selection are all
     fit without ever seeing the held-out genomes. No leakage.
  3. pick the best cutoff (metric-aware), assign that fold's feature space to the
     held-out genomes (``assign_kmer_features``), and ensemble-predict their
     target vectors (``predict_multioutput``).
  4. aggregate held-out predictions across folds and score per target.

Mirrors the user's strain x phage nested-CV (which holds the fold out of the
whole modeling workflow + ensemble-predicts), but for the single-genome
multi-output case: genomes are the "strain" side, there is no phage side, and
held-out genomes are assigned from the fold's feature table rather than via
MMseqs re-assignment.
"""

import os
import json
import logging

import numpy as np
import pandas as pd

from genophi.workflows.kmer_full_workflow import run_kmer_workflow
from genophi.workflows.multioutput_prediction_workflow import (
    load_genome_sequences, assign_kmer_features, predict_multioutput,
)
from genophi.select_feature_modeling import (
    evaluate_multilabel_performance, evaluate_multitarget_performance,
)

logger = logging.getLogger(__name__)


def _kfold_splits(genomes, n_folds, cv_rounds=1):
    """Deterministic repeated k-fold over a sorted genome list.

    Returns a list of (iteration, modeling_genomes, heldout_genomes). No RNG
    (avoids reproducibility issues); folds are contiguous slices of a
    round-shifted ordering.
    """
    genomes = sorted(genomes)
    n = len(genomes)
    splits = []
    it = 0
    for r in range(cv_rounds):
        order = genomes[r:] + genomes[:r]   # rotate per round for variety
        for k in range(n_folds):
            heldout = [order[i] for i in range(n) if i % n_folds == k]
            modeling = [g for g in order if g not in set(heldout)]
            if not heldout or not modeling:
                continue
            it += 1
            splits.append((it, modeling, heldout))
    return splits


def _best_cutoff_dir(fold_modeling_dir, task_type):
    """Pick the best cutoff's model dir from a completed fold modeling run.

    Reads modeling_results/model_performance/model_performance_metrics.csv,
    which the evaluators already sort best-first (macro_f1 / mean_r2). Returns
    (model_dir, best_cutoff_label, metrics_row).
    """
    mr = os.path.join(fold_modeling_dir, 'modeling_results')
    metrics_file = os.path.join(mr, 'model_performance', 'model_performance_metrics.csv')
    if not os.path.exists(metrics_file):
        raise FileNotFoundError(f"No fold metrics at {metrics_file}")
    m = pd.read_csv(metrics_file)
    if m.empty:
        raise ValueError(f"Empty fold metrics at {metrics_file}")
    best = m.iloc[0]   # evaluators sort best-first
    cutoff = str(best['cut_off'])
    return os.path.join(mr, cutoff), cutoff, best


def _predict_heldout_joint(fold_dir, targets, task_type, seqs, feature_map,
                           sample_column, threads):
    """Assign + ensemble-predict held-out genomes for a JOINT fold.

    One shared model dir (best cutoff), one assignment, multi-output predict.
    Returns a DataFrame with sample_column + per-target Prediction_/Confidence_.
    """
    model_dir, cutoff, _ = _best_cutoff_dir(os.path.join(fold_dir, 'modeling'), task_type)
    predictive_table = os.path.join(
        fold_dir, 'modeling', 'feature_selection', 'filtered_feature_tables',
        f'select_feature_table_{cutoff}.csv')
    assigned = assign_kmer_features(
        seqs, feature_map, predictive_table,
        sample_column=sample_column, target_columns=targets, threads=threads)
    return predict_multioutput(model_dir, assigned, sample_column=sample_column,
                               threads=threads)


def _predict_heldout_independent(fold_dir, targets, task_type, seqs, feature_map,
                                 sample_column, threads):
    """Assign + ensemble-predict held-out genomes for an INDEPENDENT fold.

    Each target has its own model dir, best cutoff, and predictive feature table
    under modeling/<target>/. Predict each target separately (binary ensemble
    across runs, median confidence) and assemble per-target columns. Targets that
    failed to model in this fold (no metrics) are filled with NaN.
    """
    from genophi.workflows.prediction_workflow import (
        predict_interactions, calculate_median_predictions,
    )

    per_target_frames = []
    sample_ids = None
    for t in targets:
        tdir = os.path.join(fold_dir, 'modeling', t)
        try:
            model_dir, cutoff, _ = _best_cutoff_dir(tdir, task_type)
        except (FileNotFoundError, ValueError):
            logger.warning(f"[independent CV] target '{t}' has no fold metrics; "
                           f"held-out predictions set to NaN.")
            continue
        predictive_table = os.path.join(
            tdir, 'feature_selection', 'filtered_feature_tables',
            f'select_feature_table_{cutoff}.csv')
        assigned = assign_kmer_features(
            seqs, feature_map, predictive_table,
            sample_column=sample_column, target_columns=targets, threads=threads)
        # binary single-genome ensemble prediction across the modeling runs.
        raw = predict_interactions(model_dir, assigned, single_strain_mode=True,
                                   strain_source=sample_column, threads=threads)
        med = calculate_median_predictions(raw, single_strain_mode=True,
                                           strain_source=sample_column)
        frame = med[[sample_column, 'Confidence', 'Final_Prediction']].rename(
            columns={'Confidence': f'Confidence_{t}', 'Final_Prediction': f'Prediction_{t}'})
        per_target_frames.append(frame)
        if sample_ids is None:
            sample_ids = frame[[sample_column]]

    if not per_target_frames:
        # no target modeled this fold; return just the held-out ids.
        ids = sorted(seqs.keys())
        return pd.DataFrame({sample_column: ids})

    out = per_target_frames[0]
    for f in per_target_frames[1:]:
        out = out.merge(f, on=sample_column, how='outer')
    return out


def run_multioutput_cv_workflow(
    input_strain_dir,
    phenotype_matrix,
    output_dir,
    phenotype_column,
    task_type='classification',
    target_mode='auto',
    strategy='joint',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    suffix='faa',
    k=5,
    num_runs_fs=25,
    num_runs_modeling=50,
    method='rfe',
    max_features='none',
    threads=4,
    max_ram=8,
    strong_top_frac=0.2,
    use_clustering=False,
    cluster_method='hierarchical',
    n_clusters=20,
    min_cluster_size=5,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
):
    """Run k-fold CV for a multi-output model (in-process loop over all folds).

    Convenience entry point that runs every fold sequentially then aggregates.
    For cluster execution (one job per fold), call ``run_one_cv_fold`` per fold
    (job array) then ``aggregate_cv_results`` once (dependent job) instead.

    Returns:
        dict with paths to the pooled predictions CSV and the per-target metrics CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    targets = _as_target_list(phenotype_column)
    genomes, splits = _cv_genomes_and_splits(
        input_strain_dir, phenotype_matrix, sample_column, suffix, n_folds, cv_rounds)
    logger.info(f"CV over {len(genomes)} genomes, {n_folds} folds x {cv_rounds} rounds")

    for it, _, _ in splits:
        run_one_cv_fold(
            fold_idx=it,
            input_strain_dir=input_strain_dir, phenotype_matrix=phenotype_matrix,
            output_dir=output_dir, phenotype_column=targets, task_type=task_type,
            target_mode=target_mode, strategy=strategy, n_folds=n_folds,
            cv_rounds=cv_rounds, sample_column=sample_column, suffix=suffix, k=k,
            num_runs_fs=num_runs_fs, num_runs_modeling=num_runs_modeling,
            method=method, max_features=max_features, threads=threads, max_ram=max_ram,
            use_clustering=use_clustering, cluster_method=cluster_method,
            n_clusters=n_clusters, min_cluster_size=min_cluster_size,
            use_feature_clustering=use_feature_clustering,
            feature_cluster_method=feature_cluster_method,
            feature_n_clusters=feature_n_clusters,
            feature_min_cluster_presence=feature_min_cluster_presence,
            filter_by_cluster_presence=filter_by_cluster_presence,
            min_cluster_presence=min_cluster_presence,
        )

    return aggregate_cv_results(
        output_dir=output_dir, phenotype_column=targets, task_type=task_type,
        n_folds=n_folds, cv_rounds=cv_rounds, sample_column=sample_column,
        strong_top_frac=strong_top_frac,
    )


def _as_target_list(phenotype_column):
    """Normalize phenotype_column to a list of target names."""
    return list(phenotype_column) if isinstance(phenotype_column, (list, tuple)) else [phenotype_column]


def _cv_genomes_and_splits(input_strain_dir, phenotype_matrix, sample_column,
                           suffix, n_folds, cv_rounds):
    """Deterministically derive the genome list and fold splits.

    Pure function of the inputs (no RNG, no shared state), so every SLURM job
    that calls it with the same args computes the IDENTICAL splits -- each job
    can then select its own fold by index with no coordination.
    """
    pheno = pd.read_csv(phenotype_matrix)
    pheno[sample_column] = pheno[sample_column].astype(str)
    genomes_matrix = set(pheno[sample_column])
    genomes_dir = {f[:-len(suffix) - 1] for f in os.listdir(input_strain_dir)
                   if f.endswith(f'.{suffix}')}
    genomes = sorted(genomes_matrix & genomes_dir)
    if len(genomes) < n_folds:
        raise ValueError(f"Only {len(genomes)} usable genomes for {n_folds} folds.")
    splits = _kfold_splits(genomes, n_folds, cv_rounds)
    return genomes, splits


def _write_splits_manifest(output_dir, splits):
    """Write an auditable per-fold split manifest (idempotent).

    output_dir/cv_splits.csv: one row per fold with the held-out and modeling
    genomes (semicolon-joined, sorted). Two CV runs that share inputs/n_folds/
    cv_rounds produce byte-identical manifests -- diff them to PROVE joint and
    independent used the same train/test splits.
    """
    path = os.path.join(output_dir, 'cv_splits.csv')
    rows = []
    for it, modeling, heldout in splits:
        rows.append({
            'fold': it,
            'n_modeling': len(modeling),
            'n_heldout': len(heldout),
            'heldout_genomes': ';'.join(sorted(heldout)),
            'modeling_genomes': ';'.join(sorted(modeling)),
        })
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def run_one_cv_fold(
    fold_idx,
    input_strain_dir,
    phenotype_matrix,
    output_dir,
    phenotype_column,
    task_type='classification',
    target_mode='auto',
    strategy='joint',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    suffix='faa',
    k=5,
    num_runs_fs=25,
    num_runs_modeling=50,
    method='rfe',
    max_features='none',
    threads=4,
    max_ram=8,
    use_clustering=False,
    cluster_method='hierarchical',
    n_clusters=20,
    min_cluster_size=5,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
):
    """Run ONE CV fold end-to-end and write its held-out predictions.

    Computes the deterministic splits, selects fold ``fold_idx`` (1-based, the
    iteration index from _kfold_splits), trains the full modeling workflow on
    that fold's modeling genomes only, then assigns + ensemble-predicts the
    held-out genomes. Writes ``output_dir/fold_<fold_idx>/fold_predictions.csv``.

    Designed to be the body of a SLURM job-array task (fold_idx =
    SLURM_ARRAY_TASK_ID). Idempotent: skips modeling if the fold's metrics marker
    already exists.

    Returns:
        str: path to this fold's fold_predictions.csv.
    """
    targets = _as_target_list(phenotype_column)
    pheno = pd.read_csv(phenotype_matrix)
    pheno[sample_column] = pheno[sample_column].astype(str)
    _, splits = _cv_genomes_and_splits(
        input_strain_dir, phenotype_matrix, sample_column, suffix, n_folds, cv_rounds)
    os.makedirs(output_dir, exist_ok=True)
    _write_splits_manifest(output_dir, splits)  # auditable, identical across runs

    match = [s for s in splits if s[0] == fold_idx]
    if not match:
        raise ValueError(
            f"fold_idx {fold_idx} not in 1..{len(splits)} "
            f"(n_folds={n_folds} x cv_rounds={cv_rounds}).")
    it, modeling, heldout = match[0]

    fold_dir = os.path.join(output_dir, f'fold_{it}')
    os.makedirs(fold_dir, exist_ok=True)
    logger.info(f"[fold {it}/{len(splits)}] {len(modeling)} modeling, {len(heldout)} held-out")

    # Modeling-genome list restricts the run; load_genome_list reads it as a CSV
    # keyed by the sample/id column, so the header must match sample_column.
    modeling_list_path = os.path.join(fold_dir, 'modeling_genomes.csv')
    pd.DataFrame({sample_column: modeling}).to_csv(modeling_list_path, index=False)

    # --- Full modeling workflow on MODELING genomes only (no leakage) ---
    if strategy == 'independent':
        metrics_marker = os.path.join(fold_dir, 'modeling', 'independent_summary.csv')
    else:
        metrics_marker = os.path.join(
            fold_dir, 'modeling', 'modeling_results', 'model_performance',
            'model_performance_metrics.csv')
    if not os.path.exists(metrics_marker):
        run_kmer_workflow(
            input_strain_dir=input_strain_dir,
            output_dir=fold_dir,
            phenotype_matrix=phenotype_matrix,
            k=k, suffix=suffix,
            strain_list=modeling_list_path,
            sample_column=sample_column,
            strain_column=sample_column,
            phenotype_column=targets,
            target_mode=target_mode,
            strategy=strategy,
            task_type=task_type,
            modeling=True,
            num_features='none',
            filter_type='none',
            num_runs_fs=num_runs_fs,
            num_runs_modeling=num_runs_modeling,
            method=method,
            max_features=max_features,
            threads=threads,
            max_ram=max_ram,
            use_clustering=use_clustering,
            cluster_method=cluster_method,
            n_clusters=n_clusters,
            min_cluster_size=min_cluster_size,
            use_feature_clustering=use_feature_clustering,
            feature_cluster_method=feature_cluster_method,
            feature_n_clusters=feature_n_clusters,
            feature_min_cluster_presence=feature_min_cluster_presence,
            filter_by_cluster_presence=filter_by_cluster_presence,
            min_cluster_presence=min_cluster_presence,
        )

    # --- Held-out assignment + ensemble prediction ---
    feature_map = os.path.join(
        fold_dir, 'kmer_tables', 'feature_tables', 'selected_features.csv')
    if not os.path.exists(feature_map):
        hits = [os.path.join(r, 'selected_features.csv')
                for r, _, fs in os.walk(fold_dir) if 'selected_features.csv' in fs]
        feature_map = hits[0] if hits else feature_map

    seqs = load_genome_sequences(input_strain_dir, suffix=suffix, genome_list=heldout)
    if strategy == 'independent':
        preds = _predict_heldout_independent(
            fold_dir, targets, task_type, seqs, feature_map, sample_column, threads)
    else:
        preds = _predict_heldout_joint(
            fold_dir, targets, task_type, seqs, feature_map, sample_column, threads)

    truth = pheno[pheno[sample_column].isin(heldout)][[sample_column] + targets]
    merged = preds.merge(truth, on=sample_column, how='left')
    merged['fold'] = it
    fold_preds_path = os.path.join(fold_dir, 'fold_predictions.csv')
    merged.to_csv(fold_preds_path, index=False)
    logger.info(f"[fold {it}] held-out predictions -> {fold_preds_path}")
    return fold_preds_path


def aggregate_cv_results(
    output_dir,
    phenotype_column,
    task_type='classification',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    strong_top_frac=0.2,
):
    """Pool per-fold held-out predictions and score per target.

    ALL-OR-NOTHING: requires every expected fold's fold_predictions.csv to be
    present. If any are missing it raises, listing them, so the missing folds can
    be re-run before a (clean) CV estimate is produced.

    Returns:
        dict with paths to the pooled predictions CSV and the per-target metrics CSV.
    """
    targets = _as_target_list(phenotype_column)
    total_folds = n_folds * cv_rounds

    frames, missing = [], []
    for it in range(1, total_folds + 1):
        fp = os.path.join(output_dir, f'fold_{it}', 'fold_predictions.csv')
        if os.path.exists(fp):
            frames.append(pd.read_csv(fp))
        else:
            missing.append(it)
    if missing:
        raise RuntimeError(
            f"Aggregation requires all {total_folds} folds; missing fold(s): "
            f"{missing}. Re-run those folds, then aggregate.")

    final = pd.concat(frames, ignore_index=True)
    preds_path = os.path.join(output_dir, 'cv_predictions.csv')
    final.to_csv(preds_path, index=False)
    logger.info(f"CV predictions ({len(final)} held-out rows, {total_folds} folds) -> {preds_path}")

    # Reshape to the layout the evaluators expect: True_<t> + a single pooled cut_off.
    scored = final.copy()
    for t in targets:
        if t in scored.columns:
            scored[f'True_{t}'] = scored[t]
    scored['cut_off'] = 'cv_pooled'

    metrics_dir = os.path.join(output_dir, 'cv_performance')
    os.makedirs(metrics_dir, exist_ok=True)
    if task_type == 'regression':
        evaluate_multitarget_performance(scored, metrics_dir, sample_column,
                                         targets, strong_top_frac=strong_top_frac)
    else:
        evaluate_multilabel_performance(scored, metrics_dir, sample_column, targets)
    metrics_path = os.path.join(metrics_dir, 'model_performance_metrics.csv')

    logger.info(f"CV per-target metrics -> {metrics_path}")
    return {'predictions': preds_path, 'metrics': metrics_path}
