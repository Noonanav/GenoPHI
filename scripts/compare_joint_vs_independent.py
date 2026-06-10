#!/usr/bin/env python3
"""
Compare joint vs. independent multi-output CV results, per target.

Reads the two cv_performance/model_performance_metrics.csv files (one row,
cut_off='cv_pooled', with per-target <METRIC>_<target> columns + n_samples_<t>),
and emits a tidy per-target comparison with joint, independent, and delta
(joint - independent) for each metric.

Usage:
  compare_joint_vs_independent.py \
      --joint   .../cv_joint_multilabel/cv_performance/model_performance_metrics.csv \
      --independent .../cv_independent/cv_performance/model_performance_metrics.csv \
      --task classification \
      -o joint_vs_independent.csv \
      [--joint_splits .../cv_joint_multilabel/cv_splits.csv \
       --independent_splits .../cv_independent/cv_splits.csv]

If both --*_splits are given, it checks they are identical (matched folds) and
refuses to compare otherwise unless --allow_split_mismatch is set.
"""
import os
import sys
import argparse

import numpy as np
import pandas as pd

CLASS_METRICS = ['AUC', 'AUPR', 'Accuracy', 'Precision', 'Recall', 'F1', 'MCC']
REG_METRICS = ['R2', 'RMSE', 'MAE', 'normalized_RMSE', 'detect_AUPR', 'detect_AUROC']
# For these, higher is better (used to color the delta interpretation).
HIGHER_BETTER = set(CLASS_METRICS + ['R2', 'detect_AUPR', 'detect_AUROC'])


def _targets_from_columns(cols, metrics):
    """Infer target names from '<metric>_<target>' columns."""
    for m in metrics:
        pref = f'{m}_'
        hits = [c[len(pref):] for c in cols if c.startswith(pref)]
        if hits:
            return hits
    return []


def _best_row(path):
    """Load the metrics CSV and return its best/only row as a Series."""
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"Empty metrics file: {path}")
    return df.iloc[0]  # cv_pooled (or evaluators' best-first sort)


def main():
    p = argparse.ArgumentParser(description="Compare joint vs independent CV per target")
    p.add_argument('--joint', required=True, help='joint cv_performance metrics CSV')
    p.add_argument('--independent', required=True, help='independent cv_performance metrics CSV')
    p.add_argument('--task', default='classification', choices=['classification', 'regression'])
    p.add_argument('-o', '--output', required=True, help='output comparison CSV')
    p.add_argument('--joint_splits')
    p.add_argument('--independent_splits')
    p.add_argument('--allow_split_mismatch', action='store_true')
    args = p.parse_args()

    metrics = REG_METRICS if args.task == 'regression' else CLASS_METRICS

    # Optional: verify the two runs used identical folds.
    if args.joint_splits and args.independent_splits:
        js = open(args.joint_splits).read()
        isp = open(args.independent_splits).read()
        if js != isp:
            msg = ("cv_splits.csv differ between joint and independent -- the folds "
                   "are NOT matched, so per-target deltas are confounded.")
            if not args.allow_split_mismatch:
                print(f"ERROR: {msg}\nPass --allow_split_mismatch to compare anyway.")
                return 1
            print(f"WARNING: {msg} (proceeding due to --allow_split_mismatch)")
        else:
            print("Splits verified identical between joint and independent. ✓")

    j = _best_row(args.joint)
    i = _best_row(args.independent)

    j_targets = _targets_from_columns(j.index, metrics)
    i_targets = _targets_from_columns(i.index, metrics)
    targets = [t for t in j_targets if t in i_targets]
    only_j = [t for t in j_targets if t not in i_targets]
    only_i = [t for t in i_targets if t not in j_targets]
    if only_j:
        print(f"NOTE: targets only in joint (skipped): {only_j}")
    if only_i:
        print(f"NOTE: targets only in independent (skipped): {only_i}")

    rows = []
    for t in targets:
        row = {'target': t}
        for m in metrics:
            jv = j.get(f'{m}_{t}', np.nan)
            iv = i.get(f'{m}_{t}', np.nan)
            row[f'{m}_joint'] = jv
            row[f'{m}_indep'] = iv
            row[f'{m}_delta'] = (jv - iv)
        # coverage (how many held-out samples each receptor was scored over)
        row['n_joint'] = j.get(f'n_samples_{t}', np.nan)
        row['n_indep'] = i.get(f'n_samples_{t}', np.nan)
        rows.append(row)

    comp = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    comp.to_csv(args.output, index=False)

    # Console summary: the headline metric (MCC for classification, R2 for regression).
    headline = 'MCC' if args.task == 'classification' else 'R2'
    print(f"\n{'='*72}\nJOINT vs INDEPENDENT per target ({headline}, delta = joint - independent)\n{'='*72}")
    print(f"{'target':12} {headline+'_joint':>10} {headline+'_indep':>10} {'delta':>9}  {'n_j':>4} {'n_i':>4}")
    for _, r in comp.iterrows():
        jv, iv, dv = r[f'{headline}_joint'], r[f'{headline}_indep'], r[f'{headline}_delta']
        print(f"{r['target']:12} {jv:10.3f} {iv:10.3f} {dv:+9.3f}  {int(r['n_joint']) if not pd.isna(r['n_joint']) else '-':>4} {int(r['n_indep']) if not pd.isna(r['n_indep']) else '-':>4}")

    # Mean deltas across targets for the key metrics.
    print(f"\nMean delta (joint - independent) across {len(comp)} targets:")
    for m in metrics:
        d = comp[f'{m}_delta'].dropna()
        if len(d):
            arrow = '(higher=better)' if m in HIGHER_BETTER else '(lower=better)'
            print(f"  {m:16} {d.mean():+.3f}  {arrow}")

    # Aggregate (whole-matrix) metrics if present in both.
    print("\nAggregate metrics:")
    for agg in ['micro_f1', 'macro_f1', 'hamming_loss', 'subset_accuracy',
                'mean_r2', 'mean_normalized_rmse', 'mean_detect_AUROC', 'mean_detect_AUPR']:
        if agg in j.index and agg in i.index:
            print(f"  {agg:20} joint={j[agg]:.3f}  indep={i[agg]:.3f}  delta={j[agg]-i[agg]:+.3f}")

    print(f"\nWrote comparison table -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
