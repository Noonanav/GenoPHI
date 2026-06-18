#!/usr/bin/env python3
"""
Plain PEQ clustermap (UPGMA dendrogram order, no fold boxes) with the held-out
phages marked via a colored side strip. Shows how distinct the left-out set is
from the training set in PEQ similarity space.

Train/held-out membership: --trained_table is a CSV with a phage column; every
phage in the PEQ matrix NOT in that table is "held-out".
"""
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

sys.path.insert(0, '/usr2/people/anoonan/BRaVE/machine_learning/genophi/scripts')
from peq_phylo_groups import load_peq

HELDOUT_C = '#d62728'   # red
TRAINED_C = '#bbbbbb'   # grey


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--peq', required=True)
    p.add_argument('--trained_table', required=True,
                   help='CSV with a phage column listing the TRAINING-set phages')
    p.add_argument('--phage_column', default='phage')
    p.add_argument('--out', required=True)
    p.add_argument('--linkage', default='average')
    p.add_argument('--label_heldout', action='store_true',
                   help='also print held-out phage names along the axis')
    args = p.parse_args()

    names, sim = load_peq(args.peq)
    trained = set(pd.read_csv(args.trained_table)[args.phage_column].astype(str))
    heldout = [n for n in names if n not in trained]
    print(f'{len(names)} phages: {len(names)-len(heldout)} trained, {len(heldout)} held-out',
          file=sys.stderr)

    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    dist = np.clip(dist, 0.0, None)
    Z = linkage(squareform(dist, checks=False), method=args.linkage)
    order = leaves_list(Z)
    M = sim[np.ix_(order, order)]
    is_held = np.array([names[i] not in trained for i in order])
    n = len(names)
    held_pos = np.where(is_held)[0]

    # gridspec: heatmap + a dedicated colorbar column, so the colorbar does NOT
    # steal width from the heatmap axis (that was misaligning the marker strip).
    fig = plt.figure(figsize=(9.5, 9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.03], wspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(M, cmap='viridis', vmin=0, vmax=1, aspect='equal',
                   interpolation='nearest')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('phages (UPGMA dendrogram order)')

    # Mark the held-out phages ON the diagonal (black squares) -- they live in the
    # same axis as the heatmap, so they cannot drift out of alignment.
    ax.scatter(held_pos, held_pos, s=14, c='k', marker='s',
               linewidths=0, zorder=3)
    # plus thin red ticks just outside the top/left edges, on the SAME axis extent
    for hp in held_pos:
        ax.plot([hp, hp], [-2.5, -0.5], c=HELDOUT_C, lw=0.8, clip_on=False, zorder=4)
        ax.plot([-2.5, -0.5], [hp, hp], c=HELDOUT_C, lw=0.8, clip_on=False, zorder=4)

    if args.label_heldout:
        ax.set_yticks(held_pos)
        ax.set_yticklabels([names[order[i]] for i in held_pos], fontsize=5)

    ax.legend(handles=[
        Patch(color='k', label=f'Held-out / predicted (n={len(heldout)})'),
        Patch(color=HELDOUT_C, label='held-out edge ticks'),
    ], loc='lower left', bbox_to_anchor=(0, 1.01), ncol=2, frameon=False, fontsize=9)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('PEQ similarity')
    fig.suptitle('PEQ clustermap: held-out (predicted) phages vs training set', y=0.95)
    fig.savefig(args.out, dpi=200, bbox_inches='tight')
    print('wrote', args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
