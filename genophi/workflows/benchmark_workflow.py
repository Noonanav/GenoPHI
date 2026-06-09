"""
Joint vs. independent benchmark reconciliation for GenoPHI multi-output runs.

Reconciles the two multi-output strategies' outputs into a single per-target
comparison table:

  - independent: ``<independent_out>/independent_summary.csv``
      one row per target, columns: target, AUC, AUPR, Accuracy, Precision,
      Recall, F1, MCC, cut_off
  - joint: ``<joint_out>/modeling_results/model_performance/model_performance_metrics.csv``
      one row per cutoff (sorted best-first), with per-target columns
      ``<METRIC>_<target>`` plus aggregate columns.

The joint table's best-cutoff row is unpivoted to one row per target, then
stacked with the independent summary into a long-format ``strategy x target``
table, plus a wide pivot per metric for quick side-by-side reading.
"""

import os
import logging

import pandas as pd

logger = logging.getLogger(__name__)

# Per-target metrics shared by both strategies (classification).
CLASSIFICATION_METRICS = ['AUC', 'AUPR', 'Accuracy', 'Precision', 'Recall', 'F1', 'MCC']
# Per-target metrics for multi-target regression.
REGRESSION_METRICS = ['R2', 'RMSE', 'MAE', 'normalized_RMSE']


def _find_joint_metrics(joint_out):
    """Locate the joint model_performance_metrics.csv under a joint output dir."""
    candidate = os.path.join(
        joint_out, 'modeling_results', 'model_performance',
        'model_performance_metrics.csv',
    )
    if os.path.exists(candidate):
        return candidate
    # Fall back to a recursive search (in case of nested layout).
    for root, _, files in os.walk(joint_out):
        if 'model_performance_metrics.csv' in files:
            return os.path.join(root, 'model_performance_metrics.csv')
    raise FileNotFoundError(
        f"Could not find model_performance_metrics.csv under {joint_out}"
    )


def _unpivot_joint(joint_metrics_path, metrics):
    """Unpivot the joint best-cutoff row into one row per target.

    Returns a long DataFrame: target, <metric...>, cut_off, strategy='joint'.
    Targets are inferred from the ``<metric>_<target>`` column names.
    """
    df = pd.read_csv(joint_metrics_path)
    if df.empty:
        raise ValueError(f"Joint metrics file is empty: {joint_metrics_path}")
    best = df.iloc[0]  # already sorted best-first by the evaluator

    # Infer target names from the first metric's columns: "<metric>_<target>".
    prefix = None
    for m in metrics:
        if any(c.startswith(f'{m}_') for c in df.columns):
            prefix = m
            break
    if prefix is None:
        raise ValueError(
            f"No per-target columns found for metrics {metrics} in {joint_metrics_path}"
        )
    targets = [c[len(prefix) + 1:] for c in df.columns if c.startswith(f'{prefix}_')]

    rows = []
    for t in targets:
        row = {'target': t, 'strategy': 'joint', 'cut_off': best.get('cut_off')}
        for m in metrics:
            col = f'{m}_{t}'
            row[m] = best[col] if col in df.columns else float('nan')
        rows.append(row)
    return pd.DataFrame(rows)


def _load_independent(independent_out, metrics):
    """Load the independent summary into a long DataFrame with strategy='independent'."""
    path = os.path.join(independent_out, 'independent_summary.csv')
    if not os.path.exists(path):
        raise FileNotFoundError(f"independent_summary.csv not found in {independent_out}")
    df = pd.read_csv(path)
    df = df.rename(columns={'target': 'target'})
    keep = ['target'] + [m for m in metrics if m in df.columns]
    if 'cut_off' in df.columns:
        keep.append('cut_off')
    out = df[keep].copy()
    out['strategy'] = 'independent'
    return out


def run_benchmark_workflow(joint_out, independent_out, output_path, task_type='classification'):
    """Reconcile joint vs. independent outputs into a comparison table.

    Args:
        joint_out (str): joint run output dir (contains modeling_results/...).
        independent_out (str): independent run output dir (contains independent_summary.csv).
        output_path (str): where to write the long-format benchmark CSV. A wide
            per-metric pivot is written alongside as ``<stem>_wide.csv``.
        task_type (str): 'classification' (default) or 'regression' — selects the
            per-target metric set.

    Returns:
        pd.DataFrame: the long-format benchmark table.
    """
    metrics = REGRESSION_METRICS if task_type == 'regression' else CLASSIFICATION_METRICS

    joint_metrics_path = _find_joint_metrics(joint_out)
    joint_long = _unpivot_joint(joint_metrics_path, metrics)
    indep_long = _load_independent(independent_out, metrics)

    # Align columns and stack.
    common_cols = ['target', 'strategy'] + metrics + ['cut_off']
    for df in (joint_long, indep_long):
        for c in common_cols:
            if c not in df.columns:
                df[c] = float('nan')
    long_df = pd.concat([indep_long[common_cols], joint_long[common_cols]], ignore_index=True)
    long_df = long_df.sort_values(['target', 'strategy']).reset_index(drop=True)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    long_df.to_csv(output_path, index=False)
    logger.info(f"Benchmark (long) written to {output_path}")

    # Wide per-metric pivot: rows = target, columns = (metric, strategy).
    stem, ext = os.path.splitext(output_path)
    wide_path = f"{stem}_wide{ext or '.csv'}"
    wide = long_df.pivot(index='target', columns='strategy', values=metrics)
    # Flatten the MultiIndex columns to "<metric>_<strategy>".
    wide.columns = [f'{metric}_{strat}' for metric, strat in wide.columns]
    # Add joint-minus-independent deltas per metric for quick reading.
    for m in metrics:
        ji, jj = f'{m}_independent', f'{m}_joint'
        if ji in wide.columns and jj in wide.columns:
            wide[f'{m}_delta(joint-indep)'] = wide[jj] - wide[ji]
    wide = wide.reset_index()
    wide.to_csv(wide_path, index=False)
    logger.info(f"Benchmark (wide, with deltas) written to {wide_path}")

    return long_df
