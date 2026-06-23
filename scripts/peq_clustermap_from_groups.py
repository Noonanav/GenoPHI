#!/usr/bin/env python3
"""
PEQ clustermap with leave-one-group-out fold boxes, where folds come from a
PREDEFINED strain->group CSV (not a threshold cut). Subsets the PEQ matrix to the
phages in the group file. Each fold draws as ONE contiguous box (folds reordered
by mean dendrogram position, members kept in dendrogram order within).

Use for plotting the exact folds a CV run used (e.g. the O157 38-phage
0.75/min2 phylo folds).
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

sys.path.insert(0, '/usr2/people/anoonan/BRaVE/machine_learning/genophi/scripts')
from peq_phylo_groups import load_peq


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--peq', required=True)
    p.add_argument('--groups', required=True, help='strain,group CSV (the CV folds)')
    p.add_argument('--strain_column', default='strain')
    p.add_argument('--group_column', default='group')
    p.add_argument('--out', required=True)
    p.add_argument('--linkage', default='average')
    args = p.parse_args()

    names, sim = load_peq(args.peq)
    g = pd.read_csv(args.groups)
    gmap = dict(zip(g[args.strain_column].astype(str), g[args.group_column].astype(str)))
    keep = [i for i, n in enumerate(names) if n in gmap]
    sub_names = [names[i] for i in keep]
    S = sim[np.ix_(keep, keep)]
    print(f'{len(keep)} phages in {g[args.group_column].nunique()} folds', file=sys.stderr)

    dist = 1.0 - S
    np.fill_diagonal(dist, 0.0); dist = np.clip(dist, 0.0, None)
    Z = linkage(squareform(dist, checks=False), method=args.linkage)
    leaf_order = leaves_list(Z)

    # fold label per sub-index
    fold = np.array([gmap[sub_names[i]] for i in range(len(sub_names))])
    leaf_pos = {i: rank for rank, i in enumerate(leaf_order)}
    fold_mean = {}
    for i in range(len(sub_names)):
        fold_mean.setdefault(fold[i], []).append(leaf_pos[i])
    fold_mean = {f: np.mean(v) for f, v in fold_mean.items()}
    order = np.array(sorted(range(len(sub_names)),
                            key=lambda i: (fold_mean[fold[i]], leaf_pos[i])))
    M = S[np.ix_(order, order)]
    lbls = fold[order]

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(M, cmap='viridis', vmin=0, vmax=1, aspect='equal', interpolation='nearest')
    start = 0
    for i in range(1, len(lbls) + 1):
        if i == len(lbls) or lbls[i] != lbls[start]:
            size = i - start
            ax.add_patch(Rectangle((start - 0.5, start - 0.5), size, size,
                                   fill=False, edgecolor='red', lw=1.4))
            # label the fold at the box's diagonal centre
            ax.text(start + size / 2 - 0.5, start + size / 2 - 0.5, lbls[start].replace('clade_', ''),
                    color='white', fontsize=8, ha='center', va='center', weight='bold')
            start = i
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('phages (UPGMA dendrogram order)')
    nfold = len(set(lbls))
    ax.set_title(f'PEQ clustermap: {len(sub_names)} phages, {nfold} leave-one-clade-out folds')
    cbar = fig.colorbar(ax.images[0], ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('PEQ similarity')
    fig.tight_layout()
    fig.savefig(args.out, dpi=200, bbox_inches='tight')
    print('wrote', args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
