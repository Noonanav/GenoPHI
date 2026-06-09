"""
Independent multi-output workflow for GenoPHI.

Trains one independent model per target column by looping the existing
single-target modeling workflow (``run_modeling_workflow_from_feature_table``)
over each target. Each target gets its own feature-selection pass and its own
model under ``output_dir/<target>/`` -- faithful to the per-target SLURM wrapper
behavior, just pulled into the package so data loading is shared and the
per-target SLURM fan-out is no longer required.

This is the ``--strategy independent`` path. It is orthogonal to the joint
multi-output path (a single CatBoost MultiLogloss/MultiRMSE model); the two are
benchmarkable against each other because they consume the same inputs.
"""

import os
import logging

import pandas as pd

from genophi.workflows.select_and_model_workflow import (
    run_modeling_workflow_from_feature_table,
)

logger = logging.getLogger(__name__)


def _resolve_targets(full_feature_table, phenotype_columns):
    """Normalize the requested targets to a validated list of column names."""
    if isinstance(phenotype_columns, str):
        targets = [phenotype_columns]
    else:
        targets = list(phenotype_columns)

    if not targets:
        raise ValueError("No phenotype target columns provided.")

    header = pd.read_csv(full_feature_table, nrows=0)
    missing = [t for t in targets if t not in header.columns]
    if missing:
        raise ValueError(
            f"Target column(s) not found in feature table {full_feature_table}: "
            f"{missing}. Available columns include: {list(header.columns)[:20]}..."
        )
    return targets


def _collect_target_summary(target, target_output_dir, task_type):
    """Read a single target's best metric row from its modeling results.

    Returns a one-row dict, or None if the metrics file is absent (e.g. the
    target's run failed). The metric column mirrors the single-target contract:
    classification ranks on MCC, regression on R2.
    """
    metrics_file = os.path.join(
        target_output_dir, 'modeling_results', 'model_performance',
        'model_performance_metrics.csv',
    )
    if not os.path.exists(metrics_file):
        logger.warning(
            f"No metrics file for target '{target}' at {metrics_file}; "
            f"skipping it in the combined summary."
        )
        return None

    metrics_df = pd.read_csv(metrics_file)
    if metrics_df.empty:
        logger.warning(f"Empty metrics file for target '{target}'.")
        return None

    metric_column = 'MCC' if task_type == 'classification' else 'r2'
    sort_ascending = False
    if metric_column not in metrics_df.columns:
        # Fall back to the first numeric column so the summary is still emitted.
        logger.warning(
            f"Metric column '{metric_column}' not in metrics for target "
            f"'{target}'; columns are {list(metrics_df.columns)}."
        )
        best = metrics_df.iloc[0].to_dict()
    else:
        best = metrics_df.sort_values(
            metric_column, ascending=sort_ascending
        ).iloc[0].to_dict()

    best = {'target': target, **best}
    return best


def run_independent_multioutput_workflow(
    full_feature_table,
    output_dir,
    phenotype_columns,
    task_type='classification',
    **modeling_kwargs,
):
    """Train one independent model per target by looping the single-target workflow.

    Args:
        full_feature_table (str): Path to the full feature table (contains all
            target columns plus features).
        output_dir (str): Base output directory. Each target's results are
            written under ``output_dir/<target>/``.
        phenotype_columns (str or list[str]): Target column name(s). A single
            string trains one target (and is equivalent to calling the
            single-target workflow directly).
        task_type (str): 'classification' or 'regression'.
        **modeling_kwargs: Passed through verbatim to
            ``run_modeling_workflow_from_feature_table`` for each target
            (threads, num_runs_fs, num_runs_modeling, method, clustering args,
            etc.). ``phenotype_column``/``output_dir``/``task_type`` are managed
            per-target and must not be supplied here.

    Returns:
        str: Path to the combined per-target summary CSV
        (``output_dir/independent_summary.csv``), or None if no target produced
        usable metrics.
    """
    for reserved in ('phenotype_column', 'task_type'):
        modeling_kwargs.pop(reserved, None)

    targets = _resolve_targets(full_feature_table, phenotype_columns)
    os.makedirs(output_dir, exist_ok=True)

    logger.info(
        f"Independent multi-output: training {len(targets)} independent "
        f"{task_type} model(s) for targets {targets}"
    )

    # When there are sibling targets, the OTHER target columns must be removed
    # from the feature table before training each target -- otherwise they remain
    # in X as features and the model both (a) leaks co-target labels and (b)
    # expects them at prediction time (held-out samples don't have them). Write a
    # per-target feature table once, with siblings stripped.
    needs_sibling_strip = len(targets) > 1
    if needs_sibling_strip:
        full_df = pd.read_csv(full_feature_table)

    summary_rows = []
    failed_targets = []
    for idx, target in enumerate(targets, start=1):
        target_output_dir = os.path.join(output_dir, target)
        os.makedirs(target_output_dir, exist_ok=True)

        logger.info(
            f"[{idx}/{len(targets)}] Training independent model for target "
            f"'{target}' -> {target_output_dir}"
        )

        # Drop sibling target columns so only `target` remains (as the label).
        if needs_sibling_strip:
            siblings = [t for t in targets if t != target]
            target_table = os.path.join(target_output_dir, f'{target}_feature_table.csv')
            full_df.drop(columns=siblings, errors='ignore').to_csv(target_table, index=False)
            target_input = target_table
        else:
            target_input = full_feature_table

        # A single hard target (e.g. too few positives -> 0 features survive the
        # cutoff filter -> modeling raises) must NOT kill the whole run. Isolate
        # each target: log the failure, record it, and continue with the rest so
        # the other targets and the combined summary are still produced.
        try:
            run_modeling_workflow_from_feature_table(
                full_feature_table=target_input,
                output_dir=target_output_dir,
                phenotype_column=target,
                task_type=task_type,
                **modeling_kwargs,
            )
        except Exception as e:
            logger.error(
                f"Target '{target}' failed during modeling and was skipped: "
                f"{type(e).__name__}: {e}"
            )
            failed_targets.append(target)
            continue

        row = _collect_target_summary(target, target_output_dir, task_type)
        if row is not None:
            summary_rows.append(row)
        else:
            # Modeling completed but produced no usable metrics (e.g. 0 features
            # at every cutoff -> no model trained).
            failed_targets.append(target)

    if failed_targets:
        logger.warning(
            f"{len(failed_targets)}/{len(targets)} target(s) failed and were "
            f"skipped: {failed_targets}"
        )

    if not summary_rows:
        logger.warning(
            "No per-target metrics collected; combined summary not written."
        )
        return None

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, 'independent_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    logger.info(
        f"Independent multi-output complete. Combined summary "
        f"({len(summary_rows)}/{len(targets)} targets, "
        f"{len(failed_targets)} failed) -> {summary_path}"
    )
    return summary_path
