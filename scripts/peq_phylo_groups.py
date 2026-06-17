#!/usr/bin/env python3
"""
Build phylogenetic CV groups from a phamClust PEQ pairwise-similarity matrix.

Input: a square TSV of pairwise PEQ similarities (1.0 = identical). The first
header token is the genome count, followed by the N genome names; each data row
is "<name>\t<sim_1>...<sim_N>".

Approach (per the chosen design):
  1. distance = 1 - PEQ (symmetrized), condensed.
  2. UPGMA (average linkage) hierarchical clustering.
  3. Cut at a PEQ-similarity threshold -> natural clusters (clades).
  4. Merge clusters below --min_size into their nearest neighbouring cluster
     (by smallest average inter-cluster distance) until every group is >= min_size,
     so every leave-one-group-out fold has enough held-out phages to be trainable.

Output: a strain,group CSV consumable by nested_cv_workflow's --group_metadata
(leave-one-group-out). Column names default to strain/group.
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def load_peq(path):
    """Return (names, sim) where sim is an NxN symmetric similarity matrix."""
    # First row: count + names. Data rows: name + N floats.
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    names = header[1:]                      # drop the leading count token
    df = pd.read_csv(path, sep='\t', skiprows=1, header=None,
                     index_col=0, names=['__row__'] + names)
    # Reorder columns to match row order; coerce to float.
    sim = df.reindex(index=names)[names].astype(float).values
    # Symmetrize (phamClust is symmetric up to float noise) and pin diagonal.
    sim = (sim + sim.T) / 2.0
    np.fill_diagonal(sim, 1.0)
    return names, sim


def cluster_threshold(sim, peq_threshold, method='average'):
    """UPGMA on distance=1-sim; cut so members are within (1-peq_threshold)."""
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    dist = np.clip(dist, 0.0, None)
    Z = linkage(squareform(dist, checks=False), method=method)
    height = 1.0 - peq_threshold            # distance cut corresponding to PEQ thr
    labels = fcluster(Z, t=height, criterion='distance')
    return labels, dist


def merge_small(labels, dist, min_size):
    """Merge clusters smaller than min_size into their nearest cluster.

    Nearest = smallest average inter-cluster distance. Repeats until no cluster
    below min_size remains (or only one cluster is left).
    """
    labels = labels.copy()
    while True:
        uniq, counts = np.unique(labels, return_counts=True)
        if len(uniq) == 1:
            break
        small = uniq[counts < min_size]
        if len(small) == 0:
            break
        # Merge the single smallest cluster first (stable, deterministic).
        order = sorted(small, key=lambda c: (counts[list(uniq).index(c)], c))
        src = order[0]
        src_idx = np.where(labels == src)[0]
        # nearest OTHER cluster by mean pairwise distance
        best, best_d = None, np.inf
        for c in uniq:
            if c == src:
                continue
            tgt_idx = np.where(labels == c)[0]
            d = dist[np.ix_(src_idx, tgt_idx)].mean()
            if d < best_d:
                best_d, best = d, c
        labels[labels == src] = best
    # Relabel 1..K contiguous, ordered by descending size for readability.
    uniq, counts = np.unique(labels, return_counts=True)
    order = [u for _, u in sorted(zip(-counts, uniq))]
    remap = {old: i + 1 for i, old in enumerate(order)}
    return np.array([remap[l] for l in labels])


def summarize(names, labels):
    s = pd.Series(labels, index=names)
    sizes = s.value_counts().sort_index()
    return sizes


def main():
    p = argparse.ArgumentParser(description="PEQ -> phylo CV groups (UPGMA, threshold+merge)")
    p.add_argument('--peq', required=True, help='pairwise_peq_similarities.tsv')
    p.add_argument('--out', help='output strain,group CSV (omit for --scan)')
    p.add_argument('--peq_threshold', type=float, default=0.75,
                   help='PEQ similarity cut (members within this similarity cluster together)')
    p.add_argument('--min_size', type=int, default=5,
                   help='merge clusters smaller than this into nearest neighbour')
    p.add_argument('--linkage', default='average',
                   help="linkage method (default average=UPGMA)")
    p.add_argument('--strain_column', default='strain')
    p.add_argument('--group_column', default='group')
    p.add_argument('--scan', action='store_true',
                   help='print cluster counts across a range of thresholds and exit')
    args = p.parse_args()

    names, sim = load_peq(args.peq)
    print(f"Loaded {len(names)} genomes; PEQ range "
          f"[{sim[~np.eye(len(names),dtype=bool)].min():.3f}, "
          f"{sim[~np.eye(len(names),dtype=bool)].max():.3f}]", file=sys.stderr)

    if args.scan:
        print("peq_thr  raw_clusters  after_merge(min={})  fold_sizes".format(args.min_size))
        for thr in [0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]:
            raw, dist = cluster_threshold(sim, thr, args.linkage)
            merged = merge_small(raw, dist, args.min_size)
            sizes = summarize(names, merged)
            print(f"{thr:5.2f}   {len(np.unique(raw)):4d}          {len(sizes):4d}"
                  f"               {list(sizes.values)}")
        return 0

    labels, dist = cluster_threshold(sim, args.peq_threshold, args.linkage)
    raw_k = len(np.unique(labels))
    merged = merge_small(labels, dist, args.min_size)
    sizes = summarize(names, merged)
    print(f"PEQ threshold {args.peq_threshold}: {raw_k} raw clusters -> "
          f"{len(sizes)} after merging (<{args.min_size}).", file=sys.stderr)
    print(f"Fold sizes: {dict(sizes)}", file=sys.stderr)

    out = pd.DataFrame({
        args.strain_column: names,
        args.group_column: [f"clade_{l:02d}" for l in merged],
    })
    if args.out:
        out.to_csv(args.out, index=False)
        print(f"Wrote {args.out} ({len(out)} strains, {len(sizes)} groups).", file=sys.stderr)
    else:
        print(out.to_csv(index=False))
    return 0


if __name__ == '__main__':
    sys.exit(main())
