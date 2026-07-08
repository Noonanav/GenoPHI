#!/usr/bin/env python3
"""
Build phylogeny-stratified (leave-one-clade-out) CV folds for the Gaborieau
E. coli strain x phage dataset -- the phylo analogue of the leave-one-serotype-out
(LOSO-O) run. Each fold holds out one phylogenetic clade of strains as
``validation``; all other strains are ``modeling`` (shared phages, novel strains
-- the same strain-only pattern as LOSO).

Grouping method (mirrors the PEQ-clade approach, scripts/peq_phylo_groups.py on
the multioutput-support branch, commit f16ca87, but starting from a TREE instead
of a similarity matrix):

  1. Read the PPanGGOLiN persistent-genome newick tree.
  2. Build the pairwise PATRISTIC distance matrix over the strains that are
     actually modeled (the 402 in the interaction data; tree tips not in that set
     -- e.g. H27 -- are dropped).
  3. Ward hierarchical clustering on those distances (ward BALANCES cluster sizes;
     average/complete linkage leave one 150-200 strain mega-clade that dominates
     the pooled metric -- unusable on E. coli's tight persistent genome).
  4. Cut at a target number of clusters K (criterion='maxclust').
  5. MERGE any clade smaller than ``min_size`` into its nearest neighbour (by mean
     inter-cluster patristic distance) so every leave-one-clade-out fold is
     trainable -- the PEQ merge step.

Ward clusters patristic-distance vectors, so groups are phylogenetically coherent
though not strictly monophyletic (accepted design choice; balanced folds matter
more than strict monophyly for a fair cross-lineage generalization test).

Emits, per cutoff K:
  * phylo_groups_K<K>.csv        -- strain,phylogroup  (the new metadata artifact)
  * folds_phylo_K<K>.csv         -- fold_label,strain,role  (LOSO-format; feeds
                                    run_predefined_fold, identical execution path)
  * folds_phylo_K<K>_summary.csv -- per-clade sizes

Everything is deterministic (no RNG). Run in the genophi conda env:
    python scripts/build_phylo_folds.py
"""
import os
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ============================ CONFIG (edit me) ============================
TREE = ("/usr2/people/anoonan/BRaVE/resources/genome_data/e_coli/ecoli_genomes/"
        "Gaborieau/ppanggolin_update/persistent.proteins.newick")
# Strains that are actually modeled -- pulled from the existing LOSO folds file so
# the phylo strain pool is byte-identical to the serotype run's pool.
LOSO_FOLDS = "/usr2/people/anoonan/BRaVE/resources/genome_data/e_coli/folds_loso_O.csv"
STRAIN_COL = "strain"

OUT_DIR = "/usr2/people/anoonan/BRaVE/resources/genome_data/e_coli/phylo_cv"
LINKAGE_METHOD = "ward"
# The two cutoffs (robustness pair). K14 = coarse floor (largest clade ~88);
# K20 = where the mega-clades resolve (largest ~64), better balanced.
K_VALUES = [14, 20]
MIN_SIZE = 5   # merge any clade < MIN_SIZE into its nearest neighbour
# =========================================================================


def patristic_matrix(tree, tips):
    """Full pairwise patristic distance matrix over ``tips`` (in that order)."""
    name2cl = {cl.name: cl for cl in tree.get_terminals()}
    n = len(tips)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(name2cl[tips[i]], name2cl[tips[j]])
            D[i, j] = D[j, i] = d
    return D


def merge_small(D, labels, min_size):
    """Merge clusters smaller than ``min_size`` into the nearest cluster (by mean
    inter-cluster patristic distance), smallest first. Returns labels relabelled
    1..K by descending size. Deterministic."""
    labels = labels.copy()
    while True:
        ids, counts = np.unique(labels, return_counts=True)
        if len(ids) <= 2:
            break
        small = [(c, i) for i, c in zip(ids, counts) if c < min_size]
        if not small:
            break
        # smallest cluster first; tie-break on lowest label id for determinism
        small.sort()
        s = small[0][1]
        smask = labels == s
        others = [i for i in ids if i != s]
        # nearest other cluster by mean inter-cluster distance; tie-break on id
        best = min(others, key=lambda o: (D[np.ix_(smask, labels == o)].mean(), o))
        labels[smask] = best
    # relabel 1..K by descending size (deterministic; ties by original id)
    ids, counts = np.unique(labels, return_counts=True)
    order = sorted(zip(ids, counts), key=lambda t: (-t[1], t[0]))
    remap = {old: new + 1 for new, (old, _) in enumerate(order)}
    return np.array([remap[x] for x in labels])


def build_for_k(D, strains, Z, K, min_size, out_dir):
    """Cut the ward tree at K, merge small clades, write group + fold + summary CSVs."""
    raw = fcluster(Z, t=K, criterion="maxclust")
    labels = merge_small(D, raw, min_size)
    n_groups = labels.max()

    groups = pd.DataFrame({STRAIN_COL: strains, "phylogroup":
                           [f"clade_{l:02d}" for l in labels]})
    grp_path = os.path.join(out_dir, f"phylo_groups_K{K}.csv")
    groups.to_csv(grp_path, index=False)

    # LOSO-format folds: one fold per clade; validation=clade, modeling=all others.
    all_strains = sorted(strains)
    rows = []
    for g in sorted(groups["phylogroup"].unique()):
        val = set(groups.loc[groups["phylogroup"] == g, STRAIN_COL])
        for s in all_strains:
            rows.append({"fold_label": g, "strain": s,
                         "role": "validation" if s in val else "modeling"})
    folds = pd.DataFrame(rows, columns=["fold_label", "strain", "role"])
    folds_path = os.path.join(out_dir, f"folds_phylo_K{K}.csv")
    folds.to_csv(folds_path, index=False)

    # Summary + integrity checks.
    summary = (groups.groupby("phylogroup").size()
               .rename("n_validation").reset_index()
               .sort_values("n_validation", ascending=False))
    summary["n_modeling"] = len(all_strains) - summary["n_validation"]
    summary.to_csv(os.path.join(out_dir, f"folds_phylo_K{K}_summary.csv"), index=False)

    # Integrity: every strain validated exactly once, no role overlap per fold.
    val_counts = folds[folds.role == "validation"].strain.value_counts()
    assert (val_counts == 1).all(), "some strain held out != once"
    assert set(folds.strain) == set(all_strains), "strain set mismatch"
    for g in folds.fold_label.unique():
        sub = folds[folds.fold_label == g]
        m = set(sub[sub.role == "modeling"].strain)
        v = set(sub[sub.role == "validation"].strain)
        assert not (m & v), f"{g}: role overlap"
        assert m and v, f"{g}: empty role"

    print(f"\n=== K={K} ({LINKAGE_METHOD}, min_size={min_size}) -> "
          f"{n_groups} clades over {len(all_strains)} strains ===")
    for _, r in summary.iterrows():
        print(f"  {r.phylogroup:10s} validation={r.n_validation:>4}  "
              f"modeling={r.n_modeling:>4}")
    print(f"  groups: {grp_path}")
    print(f"  folds:  {folds_path}")
    return folds_path


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    modeled = set(pd.read_csv(LOSO_FOLDS)[STRAIN_COL].astype(str))
    tree = Phylo.read(TREE, "newick")
    tips = [cl.name for cl in tree.get_terminals()]
    strains = [t for t in tips if t in modeled]
    dropped = sorted(set(tips) - modeled)
    missing = sorted(modeled - set(tips))
    print(f"tree tips: {len(tips)}; modeled strains: {len(modeled)}; "
          f"used (intersection): {len(strains)}")
    if dropped:
        print(f"  tree tips dropped (not in interaction data): {dropped}")
    if missing:
        # These strains would be silently excluded from CV -- flag loudly.
        print(f"  WARNING: {len(missing)} modeled strain(s) NOT in tree, "
              f"excluded from phylo folds: {missing}")

    D = patristic_matrix(tree, strains)
    Z = linkage(squareform(D, checks=False), method=LINKAGE_METHOD)

    for K in K_VALUES:
        build_for_k(D, strains, Z, K, MIN_SIZE, OUT_DIR)


if __name__ == "__main__":
    main()
