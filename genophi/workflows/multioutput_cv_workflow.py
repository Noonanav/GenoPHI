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
):
    """Run k-fold CV for a multi-output model, returning aggregated metrics.

    See module docstring for the per-fold logic. ``phenotype_column`` is the
    target list (single name also accepted). ``strategy``/``target_mode`` select
    joint vs. independent and the task mode, threaded into the per-fold modeling.

    Returns:
        dict with paths to the per-fold predictions CSV and the aggregated
        per-target metrics CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    targets = phenotype_column if isinstance(phenotype_column, (list, tuple)) else [phenotype_column]
    targets = list(targets)

    # Genomes present in BOTH the phenotype matrix and the FASTA dir.
    pheno = pd.read_csv(phenotype_matrix)
    pheno[sample_column] = pheno[sample_column].astype(str)
    genomes_matrix = set(pheno[sample_column])
    genomes_dir = {f[:-len(suffix) - 1] for f in os.listdir(input_strain_dir)
                   if f.endswith(f'.{suffix}')}
    genomes = sorted(genomes_matrix & genomes_dir)
    if len(genomes) < n_folds:
        raise ValueError(f"Only {len(genomes)} usable genomes for {n_folds} folds.")
    logger.info(f"CV over {len(genomes)} genomes, {n_folds} folds x {cv_rounds} rounds")

    splits = _kfold_splits(genomes, n_folds, cv_rounds)
    all_preds = []

    for it, modeling, heldout in splits:
        fold_dir = os.path.join(output_dir, f'fold_{it}')
        os.makedirs(fold_dir, exist_ok=True)
        logger.info(f"[fold {it}/{len(splits)}] {len(modeling)} modeling, "
                    f"{len(heldout)} held-out")

        # Write this fold's modeling-genome list (strain_list restricts the run).
        # load_genome_list reads it as a CSV keyed by the sample/id column, so
        # the header must match sample_column.
        modeling_list_path = os.path.join(fold_dir, 'modeling_genomes.csv')
        pd.DataFrame({sample_column: modeling}).to_csv(modeling_list_path, index=False)

        # --- Step 1+2: full modeling workflow on MODELING genomes only ---
        # Idempotency marker differs by strategy: joint writes one shared
        # modeling_results; independent writes a per-target tree + summary.
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
                use_clustering=False,
            )

        # --- Step 3: held-out assignment + ensemble prediction ---
        # The feature_map (k-mer -> feature) is shared across the fold regardless
        # of strategy; it is written by the k-mer table step.
        feature_map = os.path.join(
            fold_dir, 'kmer_tables', 'feature_tables', 'selected_features.csv')
        if not os.path.exists(feature_map):
            hits = [os.path.join(r, 'selected_features.csv')
                    for r, _, fs in os.walk(fold_dir) if 'selected_features.csv' in fs]
            feature_map = hits[0] if hits else feature_map

        seqs = load_genome_sequences(input_strain_dir, suffix=suffix, genome_list=heldout)

        if strategy == 'independent':
            preds = _predict_heldout_independent(
                fold_dir, targets, task_type, seqs, feature_map,
                sample_column, threads)
        else:
            preds = _predict_heldout_joint(
                fold_dir, targets, task_type, seqs, feature_map,
                sample_column, threads)

        # attach the true target values for the held-out genomes.
        truth = pheno[pheno[sample_column].isin(heldout)][[sample_column] + targets]
        merged = preds.merge(truth, on=sample_column, how='left')
        merged['fold'] = it
        all_preds.append(merged)

    final = pd.concat(all_preds, ignore_index=True)
    preds_path = os.path.join(output_dir, 'cv_predictions.csv')
    final.to_csv(preds_path, index=False)
    logger.info(f"CV predictions ({len(final)} held-out rows) -> {preds_path}")

    # --- Step 4: score the pooled held-out predictions per target ---
    # Reshape pooled predictions into the layout the evaluators expect:
    # Prediction_<t>/Confidence_<t>/True_<t> + cut_off (single pooled "cutoff").
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
