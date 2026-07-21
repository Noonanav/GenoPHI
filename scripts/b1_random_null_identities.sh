#!/bin/bash
#SBATCH --job-name=b1_random_null
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# B1 / A1a -- RANDOM-CONTROL (the "gate").
#
# Question: is the model's SELECTED-feature identity-to-training a REAL discovered
# gradient, or would ANY same-size random set of protein clusters give the same
# per-strain identity-to-training profile? We build a NULL: for each fold, draw M
# random cluster-sets, each the SAME SIZE as that fold's selected-feature set, from
# ALL clusters in the pool (selected + not). For each draw we recompute each held-out
# strain's mean identity-to-training over that random cluster set. Downstream (on the
# analysis host) we ask whether the OBSERVED (selected) metric's correlation with
# per-strain AUROC beats the null distribution of that same correlation.
#
# This job NEVER moves the big search file. It streams each fold's ~9 GB val search
# ONCE, and emits only a tiny per-strain x per-draw summary:
#     fold  draw_id  strain  mean_identity  n_hits
# plus the observed (selected-cluster) row as draw_id = -1.
#
# ---- Inputs (LRC scratch) ----
#   VAL search  : $STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv
#                 col1 = held-out(test) protein, col2 = training protein, col4 = identity
#   full protein->cluster map for the POOL: shared clustering assigned_clusters
#                 (representative mapping) -- every training protein -> its cluster rep.
#   selected feature -> cluster : $STRUCT/folds/$f/strain/features/selected_features.csv
#                 (Feature, Cluster_Label) ; Cluster_Label is the cluster representative.
#
# ---- Method (one streaming pass per fold) ----
#   1. Read selected_features.csv -> SET S of selected cluster reps; k = |S|.
#   2. Read the full protein->cluster-rep map -> universe U of clusters (the draw pool).
#   3. Draw M random size-k subsets of U (seeded). For each cluster c, precompute the
#      list of draws that include c (a bitmap column); the observed set S is "draw -1".
#   4. Stream the val search once. For each row (test_protein, train_protein, ident):
#        - look up train_protein's cluster c (skip if unmapped)
#        - test_strain = prefix of col1 before '::'
#        - for observed (if c in S) and for every draw d that sampled c:
#            acc[d][strain].sum += ident ; acc[d][strain].n += 1
#   5. Emit fold,draw_id,strain,mean_identity(=sum/n),n_hits.
#
# NOTE the leak filter is already handled by the val search itself: col1 is held-out,
# col2 is training (fold-external) -- self-hits and test-test hits are absent by
# construction (same guarantee we relied on for the observed metric).

set -euo pipefail

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
PHYLO=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results
ANNOT=/global/scratch/users/anoonan/BRaVE/struct_cv/annotate
OUT=/global/scratch/users/anoonan/BRaVE/struct_cv/b1_random_null
mkdir -p "$OUT"

# full protein -> cluster-representative map for the entire pool (label-independent
# shared clustering). Format: clusters.tsv <rep>\t<member> (rep-first, MMseqs createtsv).
# (This is the SAME shared clustering used in filter_pred_identities.sh.)
POOL_MAP="$PHYLO/shared_clustering/strain/clusters.tsv"

# TRUE predictive-cluster set = distinct 'cluster' values in the annotate overview
# (strain_predictive_feature_overview.csv: strain,Feature,cluster,protein_ID,...).
# This is the SAME source as the observed feature-identity metric -- the RFE-selected
# predictive features, NOT selected_features.csv (which is the full collapsed feature
# space that GOES INTO selection, ~26k clusters -> a meaningless null size).

M=1000                 # random draws per fold
SEED=20260721
# Orientation of POOL_MAP columns. MMseqs `createtsv` emits <rep>\t<member> (rep
# first) -> REP_COL=1, PROT_COL=2. If the file is <protein>\t<rep> instead, set
# REP_COL=2 PROT_COL=1. Preflight check #1 tells you which; default = MMseqs.
REP_COL=${REP_COL:-1}
PROT_COL=${PROT_COL:-2}

for i in $(seq -w 1 16); do
  f=struct_${i}
  val="$STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv"
  self="$ANNOT/$f/strain_predictive_feature_overview.csv"   # predictive-cluster source
  dst="$OUT/${f}_null_summary.csv"
  if [ -s "$dst" ]; then echo "[$f] done, skip"; continue; fi
  if [ ! -s "$val" ] || [ ! -s "$self" ] || [ ! -s "$POOL_MAP" ]; then
    echo "[$f] MISSING input (val:$([ -s "$val" ]&&echo y||echo n) self:$([ -s "$self" ]&&echo y||echo n) pool:$([ -s "$POOL_MAP" ]&&echo y||echo n))"
    continue
  fi

  echo "[$(date)] [$f] building draws + streaming $(du -h "$val" | cut -f1) ..."

  python3 - "$f" "$val" "$self" "$POOL_MAP" "$dst" "$M" "$SEED" "$REP_COL" "$PROT_COL" <<'PY'
import sys, csv, random

fold, val, self_csv, pool_map, dst, M, seed, rep_col, prot_col = sys.argv[1:10]
M = int(M); seed = int(seed)
REP = int(rep_col) - 1        # 0-based index of the cluster-rep column in POOL_MAP
PROT = int(prot_col) - 1      # 0-based index of the protein column in POOL_MAP

# 1. FULL protein -> cluster-rep map for the pool (needed both to resolve selected
#    proteins to their clusters AND to map streamed training-protein hits). The pool
#    map is long-format <rep>\t<member>; universe = distinct reps (the draw pool).
prot2clu = {}
universe = set()
with open(pool_map) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) > max(REP, PROT):
            rep = parts[REP]
            prot2clu[parts[PROT]] = rep
            universe.add(rep)
universe = list(universe)
print(f"[{fold}] protein->cluster entries={len(prot2clu)}  "
      f"pool clusters={len(universe)}", file=sys.stderr)

# 2. TRUE predictive-cluster set for this fold, from the annotate overview
#    (strain_predictive_feature_overview.csv: strain,Feature,cluster,protein_ID,...).
#    The 'cluster' column is a cluster REP (same namespace as clusters.tsv col1, verified),
#    so match directly. k = number of distinct predictive clusters = the null-draw size,
#    the SAME set the observed feature-identity metric is built on.
selected = set()
n_unmapped = 0
with open(self_csv) as fh:
    r = csv.DictReader(fh)
    if "cluster" not in r.fieldnames:
        sys.exit(f"[{fold}] ERROR: no 'cluster' column in {self_csv}; "
                 f"cols={r.fieldnames}")
    for row in r:
        c = row["cluster"]
        if c:
            selected.add(c)
# sanity: how many selected clusters are actually reps in the pool universe
uni_set = set(universe)
in_pool = sum(1 for c in selected if c in uni_set)
n_unmapped = len(selected) - in_pool
k = len(selected)
print(f"[{fold}] predictive clusters k={k}  in_pool={in_pool} unmapped={n_unmapped}",
      file=sys.stderr)
if k == 0:
    sys.exit(f"[{fold}] ERROR: k=0 -- no predictive clusters read from {self_csv}")

# 3. draws: for each cluster rep, which draws include it. Observed set = draw id -1.
rng = random.Random(seed + int(fold.split('_')[-1]))   # stable per-fold seed
draw_members = {}   # cluster -> list[int]
for c in selected:
    draw_members.setdefault(c, []).append(-1)
for d in range(M):
    for c in rng.sample(universe, k):
        draw_members.setdefault(c, []).append(d)

# 4. stream the val search once; accumulate per (draw, strain)
# acc[draw][strain] = [sum_ident, n]
acc = {}
def bump(d, strain, ident):
    dd = acc.setdefault(d, {})
    e = dd.get(strain)
    if e is None:
        dd[strain] = [ident, 1]
    else:
        e[0] += ident; e[1] += 1

nrow = 0
with open(val) as fh:
    for line in fh:
        p = line.split("\t")
        if len(p) < 4:
            continue
        train = p[1]
        c = prot2clu.get(train)
        if c is None:
            continue
        draws = draw_members.get(c)
        if not draws:
            continue
        strain = p[0].split("::", 1)[0]
        try:
            ident = float(p[3])
        except ValueError:
            continue
        for d in draws:
            bump(d, strain, ident)
        nrow += 1
print(f"[{fold}] kept-hit rows={nrow}", file=sys.stderr)

# 5. emit summary
with open(dst, "w") as out:
    w = csv.writer(out)
    w.writerow(["fold", "draw_id", "strain", "mean_identity", "n_hits"])
    for d in sorted(acc):
        for strain, (s, n) in acc[d].items():
            w.writerow([fold, d, strain, f"{s/n:.5f}", n])
print(f"[{fold}] wrote {dst}", file=sys.stderr)
PY

  echo "[$(date)] [$f] -> $(wc -l < "$dst") rows"
done

echo "=== DONE ==="
du -sh "$OUT"
wc -l "$OUT"/*_null_summary.csv | tail -1
