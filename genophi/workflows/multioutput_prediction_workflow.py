"""
Multi-output prediction for GenoPHI joint models.

Applies a trained joint multi-output model (or an ensemble of per-run models)
to a feature table for NEW samples and emits per-target predictions. Handles
multilabel classification (MultiLogloss), multiclass classification
(MultiClass), and multitarget regression (MultiRMSE), using the
``best_model_metadata.json`` sidecar written at training time to recover target
names, the task type, and per-label thresholds.

This is the inference counterpart to the training-time evaluators
(evaluate_multilabel_performance / evaluate_multitarget_performance): training
self-evaluates on its internal test split; this module applies the saved models
to held-out / new samples. It is what the multi-output cross-validation
orchestrator uses for the held-out fold, and what you would use to score brand-
new phages against a trained joint model.
"""

import os
import json
import logging
import multiprocessing
from functools import partial

import numpy as np
import pandas as pd
from Bio import SeqIO

from genophi.workflows.prediction_workflow import load_model, align_feature_names

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------
# Held-out k-mer feature assignment (ported from the RBP wrapper's proven path)
# --------------------------------------------------------------------------

def load_genome_sequences(input_dir, suffix='faa', genome_list=None):
    """Load AA sequences per genome from a directory of FASTA files.

    Returns {genome_name: [seq, ...]} where genome_name is the file basename
    (minus suffix). If genome_list is given, only those genomes are loaded.
    """
    sequences = {}
    for fn in os.listdir(input_dir):
        if not fn.endswith(f'.{suffix}'):
            continue
        name = fn[:-len(suffix) - 1]
        if genome_list and name not in set(genome_list):
            continue
        seqs = [str(rec.seq) for rec in SeqIO.parse(os.path.join(input_dir, fn), 'fasta')]
        sequences[name] = seqs
    logger.info(f"Loaded sequences for {len(sequences)} genomes from {input_dir}")
    return sequences


def _assign_one_genome(item, feature_to_kmers):
    """Presence/absence of each feature in one genome (any k-mer in any seq)."""
    name, seqs = item
    presence = {}
    for feature, kmers in feature_to_kmers.items():
        presence[feature] = 1 if any(any(km in s for s in seqs) for km in kmers) else 0
    return name, presence


def assign_kmer_features(genome_sequences, feature_map_file, predictive_feature_table,
                         sample_column='phage', target_columns=None, threads=4):
    """Build a presence/absence feature table for new genomes by k-mer matching.

    Mirrors the held-out assignment used in the RBP wrapper: each predictive
    feature is 1 for a genome if any of that feature's k-mers occurs in any of
    the genome's sequences. The feature set is restricted to the predictive
    features of the chosen cutoff (the model's feature space), so the result
    aligns with what the fold's model expects.

    Args:
        genome_sequences (dict): {genome: [seqs]} from load_genome_sequences.
        feature_map_file (str): selected_features.csv with Feature, Cluster_Label
            (k-mer) columns from the modeling fold.
        predictive_feature_table (str): the cutoff's select_feature_table CSV;
            its non-id / non-target columns are the predictive features.
        sample_column (str): id column name to emit.
        target_columns (list): target column names to exclude when identifying
            predictive features (so targets are never treated as features).
        threads (int): parallel workers.

    Returns:
        pd.DataFrame: sample_column + one 0/1 column per predictive feature.
    """
    target_columns = target_columns or []
    feature_map = pd.read_csv(feature_map_file)
    pred_df = pd.read_csv(predictive_feature_table, nrows=1)
    exclude = set([sample_column, 'phage', 'strain', 'interaction'] + list(target_columns))
    predictive_features = [c for c in pred_df.columns if c not in exclude]

    fmap = feature_map[feature_map['Feature'].isin(predictive_features)]
    feature_to_kmers = fmap.groupby('Feature')['Cluster_Label'].apply(list).to_dict()
    logger.info(f"Assigning {len(feature_to_kmers)} predictive features to "
                f"{len(genome_sequences)} genomes")

    func = partial(_assign_one_genome, feature_to_kmers=feature_to_kmers)
    if threads and threads > 1 and len(genome_sequences) > 1:
        with multiprocessing.Pool(processes=threads) as pool:
            results = pool.map(func, genome_sequences.items())
    else:
        results = [func(it) for it in genome_sequences.items()]

    presence = {name: feats for name, feats in results}
    df = pd.DataFrame.from_dict(presence, orient='index')
    # Ensure every predictive feature is present (0 if a genome never hit it).
    for feat in predictive_features:
        if feat not in df.columns:
            df[feat] = 0
    df.index.name = sample_column
    df = df.reset_index()[[sample_column] + predictive_features]
    return df


def load_model_metadata(model_run_dir):
    """Load the best_model_metadata.json sidecar from a run directory.

    Returns the metadata dict, or None if absent (e.g. a single-target/binary
    model that predates multi-output metadata).
    """
    meta_path = os.path.join(model_run_dir, 'best_model_metadata.json')
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            return json.load(f)
    return None


def _multilabel_proba_matrix(model, X):
    """N x K positive-class probability matrix from a MultiLogloss model."""
    proba = np.asarray(model.predict_proba(X))
    if proba.ndim == 3:          # (K, N, 2) -> (N, K)
        proba = proba[:, :, 1].T
    return proba


def _predict_one_model(model, metadata, X):
    """Predict per-target with a single model, returning a tidy long dict.

    Returns a dict of {column_name: array} for the per-sample outputs, keyed by
    the actual target names from metadata. Shapes:
      - multilabel:  Confidence_<t> (prob), Prediction_<t> (thresholded)
      - multitarget: Prediction_<t> (continuous)
      - multiclass:  Prediction (label), Confidence_<class> (prob per class)
    """
    task = metadata['specific_task']
    targets = metadata['target_names']
    out = {}

    if task == 'multilabel_classification':
        proba = _multilabel_proba_matrix(model, X)
        thresholds = metadata.get('per_label_thresholds') or {}
        for k, t in enumerate(targets):
            thr = thresholds.get(t, 0.5)
            out[f'Confidence_{t}'] = proba[:, k]
            out[f'Prediction_{t}'] = (proba[:, k] >= thr).astype(int)

    elif task == 'multitarget_regression':
        pred = np.asarray(model.predict(X))
        for k, t in enumerate(targets):
            out[f'Prediction_{t}'] = pred[:, k]

    elif task == 'multiclass_classification':
        proba = np.asarray(model.predict_proba(X))   # N x C
        classes = metadata.get('class_labels') or list(range(proba.shape[1]))
        out['Prediction'] = model.predict(X).ravel()
        for ci, cls in enumerate(classes):
            out[f'Confidence_{cls}'] = proba[:, ci]

    else:
        raise ValueError(
            f"multioutput prediction does not handle specific_task='{task}'. "
            f"Use the standard prediction workflow for binary/single-target."
        )
    return out


def predict_multioutput(model_dir, feature_table, sample_column='phage', threads=4):
    """Ensemble multi-output prediction across all per-run models in model_dir.

    Args:
        model_dir (str): directory containing run* subdirs, each with
            best_model.pkl + best_model_metadata.json.
        feature_table (str or pd.DataFrame): feature table for the samples to
            predict (sample_column + feature columns; no target columns needed).
        sample_column (str): identifier column to carry through.
        threads (int): CatBoost prediction threads.

    Returns:
        pd.DataFrame: one row per sample, with the sample_column plus per-target
        Prediction_<t>/Confidence_<t> columns. For classification, confidences
        are medianed across runs and predictions re-thresholded; for regression,
        predictions are medianed across runs.
    """
    if isinstance(feature_table, str):
        feature_table = pd.read_csv(feature_table)
    if sample_column not in feature_table.columns:
        raise ValueError(f"sample_column '{sample_column}' not in feature table.")

    sample_ids = feature_table[sample_column].reset_index(drop=True)
    # Drop any non-feature columns we know about; keep numeric features.
    X_full = feature_table.drop(columns=[sample_column], errors='ignore')
    X_full = X_full.select_dtypes(include=[np.number]).reset_index(drop=True)

    run_dirs = sorted(d for d in os.listdir(model_dir) if d.startswith('run'))
    if not run_dirs:
        raise ValueError(f"No run* model subdirectories found in {model_dir}.")

    metadata = None
    per_run = []   # list of dicts of column->array, one per run
    for rd in run_dirs:
        run_path = os.path.join(model_dir, rd)
        model_file = os.path.join(run_path, 'best_model.pkl')
        if not os.path.exists(model_file):
            model_file = os.path.join(run_path, 'best_model.cbm')
            if not os.path.exists(model_file):
                logger.warning(f"No model in {rd}; skipping.")
                continue
        md = load_model_metadata(run_path)
        if md is None:
            raise ValueError(
                f"No best_model_metadata.json in {run_path}; this model was not "
                f"trained as multi-output. Use the standard prediction workflow."
            )
        metadata = md
        model = load_model(model_file)
        aligned = align_feature_names(model, X_full)
        per_run.append(_predict_one_model(model, metadata, aligned))

    if not per_run:
        raise ValueError(f"No usable models found in {model_dir}.")

    task = metadata['specific_task']
    targets = metadata['target_names']
    result = pd.DataFrame({sample_column: sample_ids})

    if task == 'multilabel_classification':
        thresholds = metadata.get('per_label_thresholds') or {}
        for t in targets:
            # median confidence across runs, then re-threshold.
            conf = np.median(np.column_stack([r[f'Confidence_{t}'] for r in per_run]), axis=1)
            thr = thresholds.get(t, 0.5)
            result[f'Confidence_{t}'] = conf
            result[f'Prediction_{t}'] = (conf >= thr).astype(int)

    elif task == 'multitarget_regression':
        for t in targets:
            pred = np.median(np.column_stack([r[f'Prediction_{t}'] for r in per_run]), axis=1)
            result[f'Prediction_{t}'] = pred

    elif task == 'multiclass_classification':
        # Average class probabilities across runs, argmax for the label.
        classes = metadata.get('class_labels') or sorted(
            {c.split('Confidence_')[1] for r in per_run for c in r if c.startswith('Confidence_')}
        )
        prob_stack = {c: np.median(np.column_stack([r[f'Confidence_{c}'] for r in per_run]), axis=1)
                      for c in classes}
        prob_mat = np.column_stack([prob_stack[c] for c in classes])
        for c in classes:
            result[f'Confidence_{c}'] = prob_stack[c]
        result['Prediction'] = [classes[i] for i in prob_mat.argmax(axis=1)]

    return result


def run_multioutput_prediction_workflow(model_dir, feature_table, output_dir,
                                        sample_column='phage', threads=4):
    """Predict per-target outputs for new samples and write a CSV.

    Returns the path to the written predictions CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    preds = predict_multioutput(model_dir, feature_table,
                                sample_column=sample_column, threads=threads)
    out_path = os.path.join(output_dir, 'multioutput_predictions.csv')
    preds.to_csv(out_path, index=False)
    logger.info(f"Multi-output predictions ({len(preds)} samples) -> {out_path}")
    return out_path
