#!/usr/bin/env python3
"""
Plot joint-vs-independent CV performance for the multi-output receptor models.

Reads, per dataset x strategy:
  data/<dataset>_<strategy>/model_performance_metrics.csv   (per-receptor bar metrics)
  data/<dataset>_<strategy>/cv_predictions.csv              (per-sample, for ROC/PR)

Produces, under plots/:
  <dataset>_AUC_bars.png      joint vs indep AUC per receptor
  <dataset>_AUPR_bars.png     joint vs indep AUPR per receptor
  <dataset>_ROC.png           ROC curve per receptor (joint + indep overlaid)
  <dataset>_PR.png            PR curve per receptor (joint + indep overlaid)

The metrics CSV has one row (cv_pooled) with columns AUC_<receptor>, AUPR_<receptor>,
... per receptor. cv_predictions.csv has one row per held-out phage with the true
label column <receptor> and the predicted probability Confidence_<receptor>.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score

HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(HERE) if os.path.basename(HERE) == 'scripts' else HERE
DATA = os.environ.get('CV_DATA', os.path.join(_REPO, 'cv_plots', 'data'))
OUT = os.environ.get('CV_OUT', os.path.join(_REPO, 'cv_plots', 'plots'))
os.makedirs(OUT, exist_ok=True)

DATASETS = ['ecoli', 'myco']
STRATS = ['joint', 'indep']
COLORS = {'joint': '#2c7fb8', 'indep': '#de2d26'}

# Run-dir roots as extracted from the LRC tar (data/<root>/cv_<strategy>/...).
ROOTS = {
    'ecoli': 'EDGE/target_modeling',
    'myco': 'BRaVE/Myco',
}
STRAT_DIR = {'joint': 'cv_joint_multilabel', 'indep': 'cv_independent'}


def _metrics_path(ds, st):
    return os.path.join(DATA, ROOTS[ds], STRAT_DIR[st],
                        'cv_performance', 'model_performance_metrics.csv')


def _preds_path(ds, st):
    return os.path.join(DATA, ROOTS[ds], STRAT_DIR[st], 'cv_predictions.csv')


def load_metric(ds, st, metric):
    """Return {receptor: value} for AUC_/AUPR_<receptor> columns of the pooled row."""
    p = _metrics_path(ds, st)
    if not os.path.exists(p):
        return {}
    df = pd.read_csv(p)
    row = df.iloc[0]
    out = {}
    for col in df.columns:
        m = re.match(rf'^{metric}_(.+)$', col)
        if m:
            out[m.group(1)] = pd.to_numeric(row[col], errors='coerce')
    return out


def receptor_order(ds):
    """Union of receptors across strategies, ordered by joint AUC (desc)."""
    j = load_metric(ds, 'joint', 'AUC')
    i = load_metric(ds, 'indep', 'AUC')
    recs = list(dict.fromkeys(list(j) + list(i)))
    recs.sort(key=lambda r: (-(j.get(r) if j.get(r) == j.get(r) else -1)))
    return recs


def bar_plot(ds, metric):
    recs = receptor_order(ds)
    if not recs:
        print(f'[{ds}] no metrics for {metric}; skipping bars')
        return
    j = load_metric(ds, 'joint', metric)
    i = load_metric(ds, 'indep', metric)
    x = np.arange(len(recs))
    w = 0.38
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * len(recs) + 2), 4.2))
    ax.bar(x - w / 2, [j.get(r, np.nan) for r in recs], w, label='Joint', color=COLORS['joint'])
    ax.bar(x + w / 2, [i.get(r, np.nan) for r in recs], w, label='Independent', color=COLORS['indep'])
    if metric == 'AUC':
        ax.axhline(0.5, ls='--', c='grey', lw=1, label='chance (0.5)')
        ax.set_ylim(0, 1)
    ax.set_xticks(x)
    ax.set_xticklabels(recs, rotation=45, ha='right')
    ax.set_ylabel(metric)
    ax.set_title(f'{ds.capitalize()}: {metric} per receptor (joint vs independent)')
    ax.legend(frameon=False)
    fig.tight_layout()
    f = os.path.join(OUT, f'{ds}_{metric}_bars.png')
    fig.savefig(f, dpi=200)
    plt.close(fig)
    print('wrote', f)


def _truth_score(preds, rec):
    """Extract (y_true, y_score) for a receptor from a cv_predictions frame."""
    score_col = f'Confidence_{rec}'
    if rec not in preds.columns or score_col not in preds.columns:
        return None
    sub = preds[[rec, score_col]].apply(pd.to_numeric, errors='coerce').dropna()
    y = sub[rec].astype(int).values
    s = sub[score_col].values
    if len(np.unique(y)) < 2:        # need both classes to draw a curve
        return None
    return y, s


def curve_plot(ds, kind):
    """kind in {'ROC','PR'}: one subplot grid, each receptor a panel, joint+indep overlaid."""
    jp, ip = _preds_path(ds, 'joint'), _preds_path(ds, 'indep')
    frames = {st: (pd.read_csv(p) if os.path.exists(p) else None)
              for st, p in [('joint', jp), ('indep', ip)]}
    if all(v is None for v in frames.values()):
        print(f'[{ds}] no cv_predictions; skipping {kind}')
        return
    recs = receptor_order(ds)
    n = len(recs)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.0 * nrow), squeeze=False)
    for idx, rec in enumerate(recs):
        ax = axes[idx // ncol][idx % ncol]
        for st in STRATS:
            preds = frames[st]
            if preds is None:
                continue
            ts = _truth_score(preds, rec)
            if ts is None:
                continue
            y, s = ts
            if kind == 'ROC':
                fpr, tpr, _ = roc_curve(y, s)
                auc = roc_auc_score(y, s)
                ax.plot(fpr, tpr, color=COLORS[st], lw=1.6,
                        label=f'{st} (AUC={auc:.2f})')
            else:
                prec, rsc, _ = precision_recall_curve(y, s)
                ap = average_precision_score(y, s)
                ax.plot(rsc, prec, color=COLORS[st], lw=1.6,
                        label=f'{st} (AP={ap:.2f})')
        if kind == 'ROC':
            ax.plot([0, 1], [0, 1], ls='--', c='grey', lw=0.8)
            ax.set_xlabel('FPR'); ax.set_ylabel('TPR')
        else:
            ax.set_xlabel('Recall'); ax.set_ylabel('Precision')
        ax.set_title(rec, fontsize=10)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
        ax.legend(fontsize=7, frameon=False, loc='lower left' if kind == 'PR' else 'lower right')
    # blank unused panels
    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis('off')
    fig.suptitle(f'{ds.capitalize()}: {kind} curves per receptor', y=1.0)
    fig.tight_layout()
    f = os.path.join(OUT, f'{ds}_{kind}.png')
    fig.savefig(f, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote', f)


def main():
    for ds in DATASETS:
        for metric in ('AUC', 'AUPR'):
            bar_plot(ds, metric)
        for kind in ('ROC', 'PR'):
            curve_plot(ds, kind)


if __name__ == '__main__':
    main()
