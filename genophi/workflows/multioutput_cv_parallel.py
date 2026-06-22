"""
Parallelized INDEPENDENT multi-output CV (one job per fold x receptor).

The in-process independent CV (multioutput_cv_workflow.run_one_cv_fold with
strategy='independent') runs all targets sequentially within a fold, which does
not fit a cluster walltime when there are many targets x folds. This module
splits the same computation into three granular, independently-callable stages
so each can be a separate SLURM job:

  Stage A  build_fold_table(fold_idx)          -- 1 job per fold
           Build that fold's k-mer feature table on the MODELING genomes only
           (no leakage). The expensive k-mer step happens ONCE per fold and is
           shared by all targets.

  Stage B  train_fold_target(fold_idx, target) -- 1 job per (fold, target)
           Load the fold's feature table, strip sibling target columns (no label
           leakage), and run feature selection + ensemble modeling for ONE
           target into fold_<f>/modeling/<target>/.

  Stage C  aggregate_independent_cv(...)        -- 1 job
           For each fold, assign + ensemble-predict the held-out genomes using
           the per-target models, pool across folds, and score per target.

Splits are the same deterministic _kfold_splits used everywhere, so the folds
match the joint run exactly (diff cv_splits.csv to prove it). The joint CV path
(multioutput_cv_workflow) is untouched.
"""

import os
import logging

import pandas as pd

from genophi.workflows.kmer_full_workflow import run_kmer_workflow
from genophi.workflows.select_and_model_workflow import (
    run_modeling_workflow_from_feature_table,
)
from genophi.workflows.multioutput_prediction_workflow import load_genome_sequences
from genophi.workflows.multioutput_cv_workflow import (
    _as_target_list, _cv_genomes_and_splits, _write_splits_manifest,
    _predict_heldout_independent, aggregate_cv_results,
)

logger = logging.getLogger(__name__)


def _fold_split(input_strain_dir, phenotype_matrix, sample_column, suffix,
                n_folds, cv_rounds, fold_idx, output_dir,
                group_metadata=None, group_strain_column='strain',
                group_column='group'):
    """Resolve (modeling, heldout) for a fold; write the manifest once.

    With ``group_metadata`` the folds are leave-one-group-out (phylo clades),
    identical to the joint path, so joint and independent share the same splits.
    """
    _, splits = _cv_genomes_and_splits(
        input_strain_dir, phenotype_matrix, sample_column, suffix, n_folds, cv_rounds,
        group_metadata=group_metadata, group_strain_column=group_strain_column,
        group_column=group_column)
    os.makedirs(output_dir, exist_ok=True)
    _write_splits_manifest(output_dir, splits)
    match = [s for s in splits if s[0] == fold_idx]
    if not match:
        raise ValueError(
            f"fold_idx {fold_idx} not in 1..{len(splits)} "
            f"(n_folds={n_folds} x cv_rounds={cv_rounds}, "
            f"group_metadata={'set' if group_metadata else 'none'}).")
    return match[0]   # (it, modeling, heldout)


def build_fold_table(
    fold_idx,
    input_strain_dir,
    phenotype_matrix,
    output_dir,
    phenotype_column,
    task_type='classification',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    suffix='faa',
    k=5,
    threads=4,
    max_ram=8,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    group_metadata=None,
    group_strain_column='strain',
    group_column='group',
):
    """Stage A: build fold ``fold_idx``'s k-mer feature table (no modeling).

    Writes the merged feature table (features + all target columns) for the
    fold's MODELING genomes only. Returns the fold dir.
    """
    targets = _as_target_list(phenotype_column)
    it, modeling, _ = _fold_split(
        input_strain_dir, phenotype_matrix, sample_column, suffix,
        n_folds, cv_rounds, fold_idx, output_dir,
        group_metadata=group_metadata, group_strain_column=group_strain_column,
        group_column=group_column)

    fold_dir = os.path.join(output_dir, f'fold_{it}')
    os.makedirs(fold_dir, exist_ok=True)
    modeling_list_path = os.path.join(fold_dir, 'modeling_genomes.csv')
    pd.DataFrame({sample_column: modeling}).to_csv(modeling_list_path, index=False)

    # Build the table only (modeling=False). Shared by all targets in Stage B.
    table_marker = os.path.join(fold_dir, 'full_feature_table.csv')
    if not os.path.exists(table_marker):
        run_kmer_workflow(
            input_strain_dir=input_strain_dir,
            output_dir=fold_dir,
            phenotype_matrix=phenotype_matrix,
            k=k, suffix=suffix,
            strain_list=modeling_list_path,
            sample_column=sample_column,
            strain_column=sample_column,
            phenotype_column=targets,
            strategy='independent',
            task_type=task_type,
            modeling=False,
            filter_type='none',
            threads=threads,
            max_ram=max_ram,
            use_clustering=False,
            use_feature_clustering=use_feature_clustering,
            feature_cluster_method=feature_cluster_method,
            feature_n_clusters=feature_n_clusters,
            feature_min_cluster_presence=feature_min_cluster_presence,
        )
    logger.info(f"[fold {it}] feature table built -> {fold_dir}")
    return fold_dir


def _fold_table_path(fold_dir):
    """Locate the fold's merged feature table (features + targets)."""
    p = os.path.join(fold_dir, 'full_feature_table.csv')
    if os.path.exists(p):
        return p
    hits = [os.path.join(r, f) for r, _, fs in os.walk(fold_dir)
            for f in fs if f == 'full_feature_table.csv']
    if hits:
        return hits[0]
    raise FileNotFoundError(f"No full_feature_table.csv under {fold_dir}")


def train_fold_target(
    fold_idx,
    target,
    output_dir,
    phenotype_column,
    task_type='classification',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    num_runs_fs=25,
    num_runs_modeling=50,
    method='rfe',
    max_features='none',
    min_features='none',
    threads=4,
    max_ram=8,
    use_clustering=False,
    cluster_method='hierarchical',
    n_clusters=20,
    min_cluster_size=5,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
):
    """Stage B: FS + ensemble modeling for ONE target on a fold's table.

    Loads fold ``fold_idx``'s feature table (from Stage A), strips the sibling
    target columns (so only ``target`` remains as the label -- no leakage), and
    trains into fold_<f>/modeling/<target>/. Idempotent.

    Requires build_fold_table(fold_idx) to have run first.
    """
    all_targets = _as_target_list(phenotype_column)
    if target not in all_targets:
        raise ValueError(f"target '{target}' not in {all_targets}")

    fold_dir = os.path.join(output_dir, f'fold_{fold_idx}')
    full_table = _fold_table_path(fold_dir)

    target_out = os.path.join(fold_dir, 'modeling', target)
    os.makedirs(target_out, exist_ok=True)

    # Idempotency: skip if this target already produced metrics.
    metrics = os.path.join(target_out, 'modeling_results', 'model_performance',
                           'model_performance_metrics.csv')
    if os.path.exists(metrics):
        logger.info(f"[fold {fold_idx}/{target}] already complete; skipping.")
        return target_out

    # Strip sibling targets so only `target` is in the table (as the label).
    siblings = [t for t in all_targets if t != target]
    full_df = pd.read_csv(full_table)
    target_table = os.path.join(target_out, f'{target}_feature_table.csv')
    full_df.drop(columns=siblings, errors='ignore').to_csv(target_table, index=False)

    logger.info(f"[fold {fold_idx}/{target}] training (siblings stripped)...")
    run_modeling_workflow_from_feature_table(
        full_feature_table=target_table,
        output_dir=target_out,
        phenotype_column=target,
        task_type=task_type,
        threads=threads,
        num_features='none',
        filter_type='none',
        num_runs_fs=num_runs_fs,
        num_runs_modeling=num_runs_modeling,
        sample_column=sample_column,
        method=method,
        max_features=max_features,
        min_features=min_features,
        binary_data=True,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        filter_by_cluster_presence=filter_by_cluster_presence,
        min_cluster_presence=min_cluster_presence,
        max_ram=max_ram,
    )
    logger.info(f"[fold {fold_idx}/{target}] done -> {target_out}")
    return target_out


def predict_fold_heldout(
    fold_idx,
    input_strain_dir,
    phenotype_matrix,
    output_dir,
    phenotype_column,
    task_type='classification',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    suffix='faa',
    threads=4,
    group_metadata=None,
    group_strain_column='strain',
    group_column='group',
):
    """Assign + ensemble-predict a fold's held-out genomes (independent models).

    Run after all of the fold's Stage-B target jobs complete. Writes
    fold_<f>/fold_predictions.csv (per-target Prediction_/Confidence_ + truth).
    Returns the predictions path.
    """
    targets = _as_target_list(phenotype_column)
    pheno = pd.read_csv(phenotype_matrix)
    pheno[sample_column] = pheno[sample_column].astype(str)
    it, _, heldout = _fold_split(
        input_strain_dir, phenotype_matrix, sample_column, suffix,
        n_folds, cv_rounds, fold_idx, output_dir,
        group_metadata=group_metadata, group_strain_column=group_strain_column,
        group_column=group_column)

    fold_dir = os.path.join(output_dir, f'fold_{it}')
    feature_map = os.path.join(
        fold_dir, 'kmer_tables', 'feature_tables', 'selected_features.csv')
    if not os.path.exists(feature_map):
        hits = [os.path.join(r, 'selected_features.csv')
                for r, _, fs in os.walk(fold_dir) if 'selected_features.csv' in fs]
        feature_map = hits[0] if hits else feature_map

    seqs = load_genome_sequences(input_strain_dir, suffix=suffix, genome_list=heldout)
    preds = _predict_heldout_independent(
        fold_dir, targets, task_type, seqs, feature_map, sample_column, threads)

    truth = pheno[pheno[sample_column].isin(heldout)][[sample_column] + targets]
    merged = preds.merge(truth, on=sample_column, how='left')
    merged['fold'] = it
    fold_preds_path = os.path.join(fold_dir, 'fold_predictions.csv')
    merged.to_csv(fold_preds_path, index=False)
    logger.info(f"[fold {it}] held-out predictions -> {fold_preds_path}")
    return fold_preds_path


def aggregate_independent_cv(
    input_strain_dir,
    phenotype_matrix,
    output_dir,
    phenotype_column,
    task_type='classification',
    n_folds=5,
    cv_rounds=1,
    sample_column='phage',
    suffix='faa',
    threads=4,
    strong_top_frac=0.2,
    group_metadata=None,
    group_strain_column='strain',
    group_column='group',
):
    """Stage C: predict each fold's held-out genomes, pool, and score per target.

    Runs predict_fold_heldout for every fold (cheap; the heavy modeling is done),
    then aggregate_cv_results (all-or-nothing) to produce per-target CV metrics.
    """
    _, splits = _cv_genomes_and_splits(
        input_strain_dir, phenotype_matrix, sample_column, suffix, n_folds, cv_rounds,
        group_metadata=group_metadata, group_strain_column=group_strain_column,
        group_column=group_column)
    for it, _, _ in splits:
        predict_fold_heldout(
            fold_idx=it, input_strain_dir=input_strain_dir,
            phenotype_matrix=phenotype_matrix, output_dir=output_dir,
            phenotype_column=phenotype_column, task_type=task_type,
            n_folds=n_folds, cv_rounds=cv_rounds, sample_column=sample_column,
            suffix=suffix, threads=threads,
            group_metadata=group_metadata, group_strain_column=group_strain_column,
            group_column=group_column)

    return aggregate_cv_results(
        output_dir=output_dir, phenotype_column=phenotype_column,
        task_type=task_type, n_folds=n_folds, cv_rounds=cv_rounds,
        sample_column=sample_column, strong_top_frac=strong_top_frac,
        group_metadata=group_metadata, input_strain_dir=input_strain_dir,
        phenotype_matrix=phenotype_matrix, suffix=suffix,
        group_strain_column=group_strain_column, group_column=group_column)
