#!/usr/bin/env python3
"""
Clustermap of the PEQ similarity matrix with cluster boxes at several thresholds.

For each PEQ threshold, UPGMA-cluster (1-PEQ distance), cut, merge clusters below
--min_size into nearest neighbour (same logic as peq_phylo_groups.py), then draw
the similarity heatmap reordered by the dendrogram with boxes around the merged
clusters. One panel per threshold so the clade structure is directly comparable.
"""
import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform

# reuse the exact grouping logic
sys.path.insert(0, '/usr2/people/anoonan/BRaVE/machine_learning/genophi/scripts')
from peq_phylo_groups import load_peq, merge_small


def ordered_boxes(sim, peq_threshold, min_size, method='average'):
    """Return (leaf_order, boxed_labels_in_leaf_order, n_raw, n_merged)."""
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    dist = np.clip(dist, 0.0, None)
    Z = linkage(squareform(dist, checks=False), method=method)
    raw = fcluster(Z, t=1.0 - peq_threshold, criterion='distance')
    merged = merge_small(raw, dist, min_size)
    order = leaves_list(Z)                      # dendrogram leaf order
    return order, merged[order], len(np.unique(raw)), len(np.unique(merged))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--peq', required=True)
    p.add_argument('--out', required=True, help='output PNG')
    p.add_argument('--thresholds', default='0.75,0.9',
                   help='comma-separated PEQ thresholds to panel')
    p.add_argument('--min_size', type=int, default=10)
    p.add_argument('--linkage', default='average')
    args = p.parse_args()

    names, sim = load_peq(args.peq)
    thrs = [float(t) for t in args.thresholds.split(',')]
    n = len(thrs)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6.5), squeeze=False)

    for j, thr in enumerate(thrs):
        ax = axes[0][j]
        order, lbls, n_raw, n_merged = ordered_boxes(sim, thr, args.min_size, args.linkage)
        M = sim[np.ix_(order, order)]
        ax.imshow(M, cmap='viridis', vmin=0, vmax=1, aspect='equal', interpolation='nearest')
        # boxes around contiguous runs of the same cluster label (leaf order)
        start = 0
        for i in range(1, len(lbls) + 1):
            if i == len(lbls) or lbls[i] != lbls[start]:
                size = i - start
                ax.add_patch(Rectangle((start - 0.5, start - 0.5), size, size,
                                       fill=False, edgecolor='red', lw=1.3))
                start = i
        ax.set_title(f'PEQ thr={thr}  (min_size={args.min_size})\n'
                     f'{n_raw} raw clades -> {n_merged} folds', fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel('phages (dendrogram order)')

    cbar = fig.colorbar(axes[0][-1].images[0], ax=axes[0].tolist(),
                        fraction=0.025, pad=0.02)
    cbar.set_label('PEQ similarity')
    fig.suptitle('PEQ clustermap with leave-one-group-out fold boxes', y=1.02, fontsize=13)
    fig.savefig(args.out, dpi=180, bbox_inches='tight')
    print(f'wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
