#!/usr/bin/env python3
"""
Random vs phylogenetic (PEQ leave-one-clade-out) CV for the E. coli receptor
multi-output models, JOINT and INDEPENDENT.

Datasets are keyed (split, strategy):
  split:    random | phylo075 | phylo090
  strategy: joint  | indep
All are 255 phages, 14 receptors, identical modeling params; only the fold
definition (random vs leave-one-clade-out at PEQ 0.75/0.90) and the modeling
strategy differ. Within a split, joint and indep use IDENTICAL folds
(matched-split: diff cv_splits.csv).

Figures (phylo_plots/):
  <strategy>_<METRIC>_bars.png   per-receptor bars across the 3 splits, for each strategy
  jointVindep_<split>_<METRIC>.png   joint vs indep within a split, per receptor
  random_vs_phylo_delta_<strategy>.png   per-receptor MCC/AUC drop random->phylo
  summary.csv   aggregate macro/micro/subset for every (split,strategy) present
"""
import os
import re
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# Data + outputs live under cv_plots/ (gitignored). Resolve relative to repo root
# so this works whether run from scripts/ or cv_plots/.
_REPO = os.path.dirname(HERE) if os.path.basename(HERE) == 'scripts' else HERE
DATA = os.environ.get('CV_DATA', os.path.join(_REPO, 'cv_plots', 'data'))
OUT = os.environ.get('CV_OUT', os.path.join(_REPO, 'cv_plots', 'phylo_plots'))
os.makedirs(OUT, exist_ok=True)

# (split, strategy) -> run-dir under DATA
DIRS = {
    ('random', 'joint'):   'EDGE/target_modeling/cv_joint_multilabel',
    ('random', 'indep'):   'EDGE/target_modeling/cv_independent',
    ('phylo075', 'joint'): 'thr075_joint',
    ('phylo075', 'indep'): 'thr075_indep',
    ('phylo090', 'joint'): 'thr09_joint',
    ('phylo090', 'indep'): 'thr09_indep',
}
SPLITS = ['random', 'phylo075', 'phylo090']
STRATS = ['joint', 'indep']
SPLIT_LABEL = {'random': 'Random 10-fold', 'phylo075': 'Phylo LOCO (0.75)',
               'phylo090': 'Phylo LOCO (0.90)'}
STRAT_LABEL = {'joint': 'Joint', 'indep': 'Independent'}
SPLIT_COLOR = {'random': '#2c7fb8', 'phylo075': '#fb6a4a', 'phylo090': '#de2d26'}
STRAT_COLOR = {'joint': '#2c7fb8', 'indep': '#de2d26'}
RECEPTORS = ['Kdo', 'btuB', 'GluI', 'tsx', 'NGR', 'HepI', 'fhuA', 'ompC',
             'lptD', 'ompA', 'ompF', 'lamB', 'HepII', 'yncD']

METRICS = ('AUC', 'AUPR', 'MCC', 'Precision', 'Recall', 'F1')


def _metrics_csv(split, strat):
    d = DIRS.get((split, strat))
    if d is None:
        return None
    p = os.path.join(DATA, d, 'cv_performance', 'model_performance_metrics.csv')
    return p if os.path.exists(p) else None


def have(split, strat):
    return _metrics_csv(split, strat) is not None


def load_metric(split, strat, metric):
    p = _metrics_csv(split, strat)
    if p is None:
        return {}
    row = pd.read_csv(p).iloc[0]
    out = {}
    for col in row.index:
        m = re.match(rf'^{metric}_(.+)$', col)
        if m:
            out[m.group(1)] = pd.to_numeric(row[col], errors='coerce')
    return out


def _style_axis(ax, metric):
    if metric == 'AUC':
        ax.axhline(0.5, ls='--', c='grey', lw=1); ax.set_ylim(0, 1)
    elif metric == 'MCC':
        ax.axhline(0.0, ls='--', c='grey', lw=1)
    elif metric in ('AUPR', 'Precision', 'Recall', 'F1'):
        ax.set_ylim(0, 1)


def bars_across_splits(strat, metric):
    """Per-receptor bars: the 3 splits side by side, for ONE strategy."""
    present = [s for s in SPLITS if have(s, strat)]
    if not present:
        return
    vals = {s: load_metric(s, strat, metric) for s in present}
    base = vals.get('random') or vals[present[0]]
    recs = [r for r in RECEPTORS if r in base]
    x = np.arange(len(recs)); w = 0.8 / len(present)
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(recs) + 2), 4.4))
    for i, s in enumerate(present):
        ax.bar(x + (i - (len(present) - 1) / 2) * w,
               [vals[s].get(r, np.nan) for r in recs], w,
               label=SPLIT_LABEL[s], color=SPLIT_COLOR[s])
    _style_axis(ax, metric)
    ax.set_xticks(x); ax.set_xticklabels(recs, rotation=45, ha='right')
    ax.set_ylabel(metric)
    ax.set_title(f'{STRAT_LABEL[strat]} — {metric}: random vs phylogenetic CV')
    ax.legend(frameon=False)
    fig.tight_layout()
    f = os.path.join(OUT, f'{strat}_{metric}_bars.png')
    fig.savefig(f, dpi=200); plt.close(fig); print('wrote', f)


def joint_vs_indep(split, metric):
    """Per-receptor bars: joint vs indep within ONE split."""
    if not (have(split, 'joint') and have(split, 'indep')):
        return
    j = load_metric(split, 'joint', metric)
    i = load_metric(split, 'indep', metric)
    recs = [r for r in RECEPTORS if r in j or r in i]
    x = np.arange(len(recs)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(recs) + 2), 4.4))
    ax.bar(x - w / 2, [j.get(r, np.nan) for r in recs], w, label='Joint', color=STRAT_COLOR['joint'])
    ax.bar(x + w / 2, [i.get(r, np.nan) for r in recs], w, label='Independent', color=STRAT_COLOR['indep'])
    _style_axis(ax, metric)
    ax.set_xticks(x); ax.set_xticklabels(recs, rotation=45, ha='right')
    ax.set_ylabel(metric)
    ax.set_title(f'{SPLIT_LABEL[split]} — {metric}: joint vs independent')
    ax.legend(frameon=False)
    fig.tight_layout()
    f = os.path.join(OUT, f'jointVindep_{split}_{metric}.png')
    fig.savefig(f, dpi=200); plt.close(fig); print('wrote', f)


def delta_plot(strat):
    """Per-receptor drop (random - phylo) for MCC and AUC, one strategy."""
    if not have('random', strat):
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, metric in zip(axes, ('MCC', 'AUC')):
        rnd = load_metric('random', strat, metric)
        recs = [r for r in RECEPTORS if r in rnd]
        x = np.arange(len(recs)); w = 0.38
        for k, s in enumerate(['phylo075', 'phylo090']):
            if not have(s, strat):
                continue
            ph = load_metric(s, strat, metric)
            drop = [rnd.get(r, np.nan) - ph.get(r, np.nan) for r in recs]
            ax.bar(x + (k - 0.5) * w, drop, w, label=SPLIT_LABEL[s], color=SPLIT_COLOR[s])
        ax.axhline(0, c='k', lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(recs, rotation=45, ha='right')
        ax.set_ylabel(f'{metric} drop (random − phylo)')
        ax.set_title(f'{metric}')
        ax.legend(frameon=False)
    fig.suptitle(f'{STRAT_LABEL[strat]}: per-receptor degradation random → leave-one-clade-out', y=1.02)
    fig.tight_layout()
    f = os.path.join(OUT, f'random_vs_phylo_delta_{strat}.png')
    fig.savefig(f, dpi=200, bbox_inches='tight'); plt.close(fig); print('wrote', f)


def calibration_plot(split):
    """Predicted-positive vs true-positive count per receptor (joint vs indep).

    Shows the calibration story: under phylo shift, indep often ranks well
    (good AUC) but predicts ~0 positives -> its thresholded MCC/F1 collapse.
    """
    if not (have(split, 'joint') and have(split, 'indep')):
        return
    j = pd.read_csv(os.path.join(DATA, DIRS[(split, 'joint')], 'cv_predictions.csv'))
    i = pd.read_csv(os.path.join(DATA, DIRS[(split, 'indep')], 'cv_predictions.csv'))
    recs = [r for r in RECEPTORS if r in j.columns and int(j[r].fillna(0).sum()) > 0]
    true = [int(j[r].fillna(0).sum()) for r in recs]
    jp = [int(j[f'Prediction_{r}'].fillna(0).sum()) if f'Prediction_{r}' in j.columns else 0 for r in recs]
    ip = [int(i[f'Prediction_{r}'].fillna(0).sum()) if f'Prediction_{r}' in i.columns else 0 for r in recs]
    x = np.arange(len(recs)); w = 0.27
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(recs) + 2), 4.4))
    ax.bar(x - w, true, w, label='True positives', color='#555555')
    ax.bar(x, jp, w, label='Joint predicted', color=STRAT_COLOR['joint'])
    ax.bar(x + w, ip, w, label='Indep predicted', color=STRAT_COLOR['indep'])
    ax.set_xticks(x); ax.set_xticklabels(recs, rotation=45, ha='right')
    ax.set_ylabel('positive count (of 255 held-out)')
    ax.set_title(f'{SPLIT_LABEL[split]} — calibration: predicted vs true positives\n'
                 f'(indep predicting ~0 despite high AUC = threshold collapse)')
    ax.legend(frameon=False)
    fig.tight_layout()
    f = os.path.join(OUT, f'calibration_{split}.png')
    fig.savefig(f, dpi=200); plt.close(fig); print('wrote', f)


def summary():
    rows = {}
    for split, strat in itertools.product(SPLITS, STRATS):
        p = _metrics_csv(split, strat)
        if p is None:
            continue
        r = pd.read_csv(p).iloc[0]
        rows[f'{split}_{strat}'] = {k: pd.to_numeric(r.get(k), errors='coerce')
                                    for k in ('macro_f1', 'micro_f1', 'hamming_loss',
                                              'subset_accuracy', 'n_complete_samples')}
    df = pd.DataFrame(rows).T
    df.to_csv(os.path.join(OUT, 'summary.csv'))
    print('\n=== aggregate summary ==='); print(df.to_string())


def main():
    for strat in STRATS:
        for m in METRICS:
            bars_across_splits(strat, m)
        delta_plot(strat)
    for split in SPLITS:
        for m in METRICS:
            joint_vs_indep(split, m)
        calibration_plot(split)
    summary()


if __name__ == '__main__':
    main()
