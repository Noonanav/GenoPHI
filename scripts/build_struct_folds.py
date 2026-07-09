#!/usr/bin/env python3
"""
Build structure-stratified (leave-one-structural-cluster-out) CV folds for the
Gaborieau E. coli strain x phage dataset -- the O-antigen-CHEMISTRY analogue of the
leave-one-serotype-out (LOSO-O) and leave-one-clade-out (phylo) runs.

Motivation. LOSO-O never actually tests structural generalization: every held-out
serogroup still has a structurally-near training serogroup. Measured on this panel:

    LOSO-O   median min-GED-to-training = 0.167   (max 0.350)
    panel    median PAIRWISE GED        = 0.500   (max 0.733)

i.e. the held-out serogroups sit in the bottom ~8% of the structural-distance
distribution. Holding out whole *structural clusters* of serogroups roughly doubles
the held-out chemistry distance (median min-to-training GED ~0.30-0.32 at K=6-8),
which is the only way to probe the far end of the surface-distance axis.

Grouping method (mirrors build_phylo_folds.py, but on O-antigen GED instead of
patristic distance):

  1. Read the serogroup x serogroup Graph Edit Distance matrix (EDGE_2/tables).
  2. Restrict to serogroups carried by the modeled strains AND having an ECODAB
     structure (unstructured / untypeable serogroups can't be clustered).
  3. WARD hierarchical clustering on the GED matrix (ward balances cluster sizes;
     average linkage leaves one mega-cluster -- same failure mode as the phylo run).
  4. Cut at target K (criterion='maxclust').
  5. MERGE any structural cluster whose GENOME count is < MIN_SIZE into its nearest
     neighbouring cluster (by mean inter-cluster GED), so every fold is trainable.

Strain assignment. A strain belongs to the cluster of its O-type. Multi-O calls
('O2/O50') are assigned to the cluster of their FIRST resolvable structured antigen,
so every strain is validated exactly once (the phylo builder's invariant). Strains
whose O-type is untypeable or lacks a structure form the ALWAYS-TRAIN pool: they are
in every fold's modeling set and are never held out.

Emits, per cutoff K:
  * struct_groups_K<K>.csv        -- strain,structgroup
  * folds_struct_K<K>.csv         -- fold_label,strain,role  (LOSO-format)
  * folds_struct_K<K>_summary.csv -- per-cluster sizes + min-GED-to-training

Deterministic (no RNG). Run in the genophi conda env:
    python scripts/build_struct_folds.py
"""
import csv
import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ============================ CONFIG (edit me) ============================
BASE = "/usr2/people/anoonan/BRaVE/resources/genome_data/e_coli"
GED = f"{BASE}/ecoli_genomes/EDGE_2/tables/ged_dist_ged.csv"
ECTYPER = f"{BASE}/ecoli_genomes/Gaborieau_genomes/ectyper_out/output.tsv"
# Strain pool taken from the existing LOSO folds so it is byte-identical to the
# serotype run's pool (same convention as build_phylo_folds.py).
LOSO_FOLDS = f"{BASE}/folds_loso_O.csv"
STRAIN_COL = "strain"

OUT_DIR = f"{BASE}/struct_cv"
# AVERAGE (not ward). Ward balances cluster sizes, which is what the PHYLO builder
# wanted -- but here it defeats the purpose: ward's clusters are structurally
# heterogeneous, so min-GED-to-training stays at ~0.167 (i.e. no better than plain
# LOSO). Average linkage yields structurally TIGHT, isolated clusters -- separation
# is the entire point of this design, so we accept less-balanced folds.
LINKAGE_METHOD = "average"
# AVERAGE, not ward: ward balances cluster SIZES (right for the phylo builder) but
# produces structurally heterogeneous clusters, leaving min-GED-to-training at ~0.167
# -- no better than plain LOSO, which defeats the purpose. Average linkage yields
# structurally tight, isolated clusters; we accept less-balanced folds because
# SEPARATION is the entire point of this design.
#
# K=20 / MIN_SIZE=3 chosen from a (K x MIN_SIZE) scan:
#   * Beyond K~18 the biggest fold stops shrinking -- floor = 63 strains (16% of the
#     panel), one genuinely tight structural cluster average linkage will not split.
#   * MERGE (not drop) undersized clusters, so every structured strain is held out
#     exactly once (329/329). Dropping tiny folds instead would leave 6-21 strains
#     never validated.
#   * minsz=3 recovers the separation that minsz=4 loses (median 0.264 vs 0.250) at
#     zero coverage cost; smallest fold is 3 strains (usable).
#
#   K=20,minsz=3 -> 16 folds, biggest 63 (16%), covered 329/329,
#                   median min-GED-to-training 0.264, max 0.385
#   [LOSO baseline: 19 folds, biggest 29 (7%),  median 0.167, max 0.350]
K_VALUES = [20]
MIN_SIZE = 3          # merge any cluster with < MIN_SIZE genomes into its neighbour
# =========================================================================


def parse_o_groups(o_type):
    """'O2' -> ['O2']; 'O2/O50' -> ['O2','O50'] (order preserved); '-' -> []."""
    if not isinstance(o_type, str):
        return []
    s = o_type.strip()
    if s in ("", "-", "O-", "ND", "NA"):
        return []
    out = []
    for part in s.split("/"):
        part = part.strip()
        if not part or part == "-":
            continue
        if not part.startswith("O"):
            part = "O" + part
        out.append(part)
    return out


def load_ged():
    d = {}
    with open(GED) as fh:
        for r in csv.DictReader(fh):
            if r["distance"] != "":
                a, b, v = r["serogroup_a"], r["serogroup_b"], float(r["distance"])
                d.setdefault(a, {})[b] = v
                d.setdefault(b, {})[a] = v
    return d


def ged_matrix(ged, serogroups):
    """Pairwise GED matrix over serogroups (dropping any with missing entries)."""
    labs = sorted(s for s in serogroups if s in ged)
    n = len(labs)
    M = np.zeros((n, n))
    for i, a in enumerate(labs):
        for j, b in enumerate(labs):
            if i != j:
                M[i, j] = ged[a].get(b, np.nan)
    ok = [i for i in range(n) if not np.isnan(M[i]).any()]
    return M[np.ix_(ok, ok)], [labs[i] for i in ok]


def merge_small_by_genomes(M, labels, sero_labels, sero_ngenomes, min_size):
    """Merge clusters whose total GENOME count < min_size into the nearest cluster
    (mean inter-cluster GED), smallest first. Deterministic; mirrors the phylo
    builder's merge_small but sizes clusters by strains, not serogroups."""
    labels = labels.copy()
    while True:
        ids = np.unique(labels)
        if len(ids) <= 2:
            break
        counts = {c: sum(sero_ngenomes.get(s, 0)
                         for s, l in zip(sero_labels, labels) if l == c) for c in ids}
        small = sorted([(counts[c], c) for c in ids if counts[c] < min_size])
        if not small:
            break
        s = small[0][1]
        smask = labels == s
        others = [i for i in ids if i != s]
        best = min(others, key=lambda o: (M[np.ix_(smask, labels == o)].mean(), o))
        labels[smask] = best
    # relabel 1..K by descending genome count (deterministic)
    ids = np.unique(labels)
    counts = {c: sum(sero_ngenomes.get(s, 0)
                     for s, l in zip(sero_labels, labels) if l == c) for c in ids}
    order = sorted(ids, key=lambda c: (-counts[c], c))
    remap = {old: new + 1 for new, old in enumerate(order)}
    return np.array([remap[x] for x in labels])


def build_for_k(M, sero_labels, Z, K, strain_sero, all_strains, min_size, out_dir):
    sero_ngenomes = {}
    for s, sg in strain_sero.items():
        if sg:
            sero_ngenomes[sg] = sero_ngenomes.get(sg, 0) + 1

    raw = fcluster(Z, t=K, criterion="maxclust")
    labels = merge_small_by_genomes(M, raw, sero_labels, sero_ngenomes, min_size)
    clmap = dict(zip(sero_labels, labels))
    n_groups = int(labels.max())

    # strain -> structgroup (always-train strains get no group)
    rows_g = []
    for s in all_strains:
        sg = strain_sero.get(s)
        c = clmap.get(sg) if sg else None
        rows_g.append({STRAIN_COL: s,
                       "structgroup": f"struct_{c:02d}" if c else "always_train"})
    groups = pd.DataFrame(rows_g)
    groups.to_csv(os.path.join(out_dir, f"struct_groups_K{K}.csv"), index=False)

    # folds: one per structural cluster; validation = its strains, modeling = rest
    fold_labels = sorted(g for g in groups.structgroup.unique() if g != "always_train")
    rows = []
    for g in fold_labels:
        val = set(groups.loc[groups.structgroup == g, STRAIN_COL])
        for s in all_strains:
            rows.append({"fold_label": g, "strain": s,
                         "role": "validation" if s in val else "modeling"})
    folds = pd.DataFrame(rows, columns=["fold_label", "strain", "role"])
    folds_path = os.path.join(out_dir, f"folds_struct_K{K}.csv")
    folds.to_csv(folds_path, index=False)

    # summary incl. the thing we care about: min GED from each cluster to training
    idx = {s: i for i, s in enumerate(sero_labels)}
    srows = []
    for g in fold_labels:
        c = int(g.split("_")[1])
        sgs = [s for s in sero_labels if clmap[s] == c]
        out = [s for s in sero_labels if clmap[s] != c]
        mn = min((M[idx[a], idx[b]] for a in sgs for b in out), default=np.nan)
        nval = int((groups.structgroup == g).sum())
        srows.append({"fold_label": g, "n_serogroups": len(sgs), "n_validation": nval,
                      "n_modeling": len(all_strains) - nval,
                      "min_ged_to_training": round(float(mn), 4),
                      "serogroups": ";".join(sgs)})
    summary = pd.DataFrame(srows).sort_values("n_validation", ascending=False)
    summary.to_csv(os.path.join(out_dir, f"folds_struct_K{K}_summary.csv"), index=False)

    # integrity: every held-out strain validated exactly once; no role overlap
    vc = folds[folds.role == "validation"].strain.value_counts()
    assert (vc == 1).all(), "some strain held out != once"
    assert set(folds.strain) == set(all_strains), "strain set mismatch"
    for g in fold_labels:
        sub = folds[folds.fold_label == g]
        v = set(sub[sub.role == "validation"].strain)
        m = set(sub[sub.role == "modeling"].strain)
        assert not (v & m), f"{g}: role overlap"
        assert v and m, f"{g}: empty role"

    n_always = int((groups.structgroup == "always_train").sum())
    print(f"\n=== K={K} ({LINKAGE_METHOD}, min_size={min_size}) -> {n_groups} "
          f"structural clusters over {len(all_strains)} strains ===")
    for _, r in summary.iterrows():
        print(f"  {r.fold_label:11s} val={r.n_validation:>3}  model={r.n_modeling:>3}  "
              f"minGED->train={r.min_ged_to_training:.3f}  ({r.n_serogroups} serogroups)")
    print(f"  always-train (untypeable / no structure): {n_always} strains")
    print(f"  folds: {folds_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    ged = load_ged()

    modeled = sorted(set(pd.read_csv(LOSO_FOLDS)[STRAIN_COL].astype(str)))
    ect = pd.read_csv(ECTYPER, sep="\t", dtype=str)
    ect = ect[ect["Name"].isin(modeled)]

    # strain -> first resolvable STRUCTURED O-group (one cluster per strain)
    strain_sero = {}
    for _, r in ect.iterrows():
        if r["Species"] != "Escherichia coli":
            continue
        for g in parse_o_groups(r["O-type"]):
            if g in ged:
                strain_sero[r["Name"]] = g
                break

    panel = sorted(set(strain_sero.values()))
    print(f"modeled strains: {len(modeled)}; with a structured O-type: {len(strain_sero)}")
    print(f"structured serogroups in panel: {len(panel)}")
    print(f"always-train pool: {len(modeled) - len(strain_sero)} strains")

    M, sero_labels = ged_matrix(ged, panel)
    Z = linkage(squareform(M, checks=False), method=LINKAGE_METHOD)

    for K in K_VALUES:
        build_for_k(M, sero_labels, Z, K, strain_sero, modeled, MIN_SIZE, OUT_DIR)


if __name__ == "__main__":
    main()
