#!/usr/bin/env python3
"""
O157 / G4C small-n (38 phage) receptor CV: random 5-fold vs phylo 9-fold (LOCO).
Single-target (G4C; O157 is unscorable), so joint==independent.

Figures (cv_plots/o157/):
  o157_metrics.png          grouped bars (random vs phylo) per metric, G4C (+O157 where scorable)
  o157_calibration.png      predicted vs true positive counts per split
  o157_phylo_foldbreakdown  per-fold G4C: held-out positives, recovered, false positives
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
DATA = os.path.join(_REPO, 'cv_plots', 'data')
OUT = os.path.join(_REPO, 'cv_plots', 'o157')
os.makedirs(OUT, exist_ok=True)

RUNS = {  # split -> dir under DATA
    'random': 'cv/indep',
    'phylo':  'cv_phylo/indep',
}
LABEL = {'random': 'Random 5-fold', 'phylo': 'Phylo LOCO (9 folds)'}
COLOR = {'random': '#2c7fb8', 'phylo': '#de2d26'}
METRICS = ('AUC', 'AUPR', 'MCC', 'Precision', 'Recall', 'F1')


def metrics_row(split):
    p = os.path.join(DATA, RUNS[split], 'cv_performance', 'model_performance_metrics.csv')
    return pd.read_csv(p).iloc[0] if os.path.exists(p) else None


def preds(split):
    p = os.path.join(DATA, RUNS[split], 'cv_predictions.csv')
    return pd.read_csv(p) if os.path.exists(p) else None


def receptors_in(row):
    return [m.group(1) for c in row.index for m in [re.match(r'^AUC_(.+)$', c)] if m]


def metric_bars():
    rr = {s: metrics_row(s) for s in RUNS}
    recs = receptors_in(rr['random'])              # G4C, O157
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, metric in zip(axes.flat, METRICS):
        x = np.arange(len(recs)); w = 0.38
        for i, s in enumerate(['random', 'phylo']):
            vals = [pd.to_numeric(rr[s].get(f'{metric}_{r}'), errors='coerce')
                    if rr[s] is not None else np.nan for r in recs]
            ax.bar(x + (i - 0.5) * w, vals, w, label=LABEL[s], color=COLOR[s])
        if metric == 'AUC':
            ax.axhline(0.5, ls='--', c='grey', lw=1)
        if metric == 'MCC':
            ax.axhline(0.0, ls='--', c='grey', lw=1)
        ax.set_xticks(x); ax.set_xticklabels(recs)
        ax.set_title(metric); ax.set_ylim(min(0, ax.get_ylim()[0]), 1)
    axes.flat[0].legend(frameon=False, fontsize=9)
    fig.suptitle('O157 receptor CV: random vs phylogenetic holdout', y=1.0)
    fig.tight_layout()
    f = os.path.join(OUT, 'o157_metrics.png'); fig.savefig(f, dpi=200, bbox_inches='tight')
    plt.close(fig); print('wrote', f)


def calibration():
    fig, ax = plt.subplots(figsize=(7, 4.4))
    recs = receptors_in(metrics_row('random'))
    x = np.arange(len(recs)); w = 0.25
    pr = {s: preds(s) for s in RUNS}
    true = [int(pr['random'][r].fillna(0).sum()) for r in recs]
    ax.bar(x - w, true, w, label='True positives', color='#555555')
    for i, s in enumerate(['random', 'phylo']):
        pp = [int(pr[s][f'Prediction_{r}'].fillna(0).sum())
              if f'Prediction_{r}' in pr[s].columns else 0 for r in recs]
        ax.bar(x + i * w, pp, w, label=f'{LABEL[s]} predicted', color=COLOR[s])
    ax.set_xticks(x); ax.set_xticklabels(recs)
    ax.set_ylabel('positive count (of 38)')
    ax.set_title('O157: predicted vs true positives')
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    f = os.path.join(OUT, 'o157_calibration.png'); fig.savefig(f, dpi=200, bbox_inches='tight')
    plt.close(fig); print('wrote', f)


def phylo_foldbreakdown():
    """Per-fold G4C: held-out positives, true-pos recovered, false positives."""
    df = preds('phylo')
    if df is None or 'fold' not in df.columns or 'G4C' not in df.columns:
        print('skip foldbreakdown (no phylo preds/fold col)'); return
    pred_col = 'Prediction_G4C'
    rows = []
    for fold, sub in df.groupby('fold'):
        true = sub['G4C'].fillna(0).astype(int)
        pred = sub[pred_col].fillna(0).astype(int) if pred_col in sub.columns else pd.Series(0, index=sub.index)
        tp = int(((true == 1) & (pred == 1)).sum())
        fp = int(((true == 0) & (pred == 1)).sum())
        npos = int((true == 1).sum())
        rows.append((str(fold), npos, tp, fp))
    rows.sort(key=lambda r: (-r[1], r[0]))
    labels = [r[0] for r in rows]
    x = np.arange(len(labels)); w = 0.27
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(labels) + 2), 4.4))
    ax.bar(x - w, [r[1] for r in rows], w, label='Held-out G4C+', color='#555555')
    ax.bar(x, [r[2] for r in rows], w, label='Recovered (TP)', color='#2ca25f')
    ax.bar(x + w, [r[3] for r in rows], w, label='False positives', color='#de2d26')
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_xlabel('phylo fold'); ax.set_ylabel('phage count')
    ax.set_title('Phylo LOCO per-fold G4C: positives held out, recovered, false positives')
    ax.legend(frameon=False)
    fig.tight_layout()
    f = os.path.join(OUT, 'o157_phylo_foldbreakdown.png'); fig.savefig(f, dpi=200, bbox_inches='tight')
    plt.close(fig); print('wrote', f)


def roc(receptor='G4C'):
    """Solo square ROC for one receptor, random vs phylo overlaid."""
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    for s in ['random', 'phylo']:
        df = preds(s)
        sc = f'Confidence_{receptor}'
        if df is None or receptor not in df.columns or sc not in df.columns:
            continue
        sub = df[[receptor, sc]].apply(pd.to_numeric, errors='coerce').dropna()
        y = sub[receptor].astype(int).values
        p = sub[sc].values
        if len(np.unique(y)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y, p)
        ax.plot(fpr, tpr, color=COLOR[s], lw=2,
                label=f'{LABEL[s]} (AUC={roc_auc_score(y, p):.2f})')
    ax.plot([0, 1], [0, 1], ls='--', c='grey', lw=0.8)
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title(f'O157 {receptor}: ROC — random vs phylogenetic CV')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.set_aspect('equal')
    ax.legend(frameon=False, fontsize=9, loc='lower right')
    fig.tight_layout()
    f = os.path.join(OUT, f'o157_{receptor}_roc.png')
    fig.savefig(f, dpi=200, bbox_inches='tight'); plt.close(fig); print('wrote', f)


def metrics_table():
    """Tidy random-vs-phylo per-receptor metrics -> CSV + rendered PNG."""
    rr = {s: metrics_row(s) for s in RUNS}
    recs = receptors_in(rr['random'])
    rows = []
    for r in recs:
        for s in ['random', 'phylo']:
            row = rr[s]
            rec = {'receptor': r, 'split': LABEL[s]}
            for m in METRICS + ('n_samples',):
                v = pd.to_numeric(row.get(f'{m}_{r}'), errors='coerce') if row is not None else np.nan
                rec[m] = round(v, 3) if pd.notna(v) else np.nan
            rows.append(rec)
    tbl = pd.DataFrame(rows)
    csv = os.path.join(OUT, 'o157_metrics_table.csv')
    tbl.to_csv(csv, index=False); print('wrote', csv)

    # rendered table image
    fig, ax = plt.subplots(figsize=(0.95 * len(tbl.columns) + 1, 0.5 * len(tbl) + 1))
    ax.axis('off')
    disp = tbl.copy().fillna('—').astype(str)
    t = ax.table(cellText=disp.values, colLabels=disp.columns, loc='center', cellLoc='center')
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.4)
    for j in range(len(disp.columns)):
        t[(0, j)].set_facecolor('#333333'); t[(0, j)].set_text_props(color='white', weight='bold')
    for i in range(1, len(disp) + 1):
        c = '#e8f0f7' if 'Random' in disp.iloc[i - 1]['split'] else '#fde6e3'
        for j in range(len(disp.columns)):
            t[(i, j)].set_facecolor(c)
    ax.set_title('O157 receptor CV: random vs phylogenetic metrics', pad=12)
    f = os.path.join(OUT, 'o157_metrics_table.png')
    fig.savefig(f, dpi=200, bbox_inches='tight'); plt.close(fig); print('wrote', f)


if __name__ == '__main__':
    metric_bars()
    calibration()
    phylo_foldbreakdown()
    roc('G4C')
    metrics_table()
