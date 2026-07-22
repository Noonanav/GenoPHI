#!/bin/bash
#SBATCH --job-name=b_null_feature
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=256G
#SBATCH --time=16:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# RANDOM-FEATURE null for the best-per-feature identity metric (the -0.382 axis).
#
# The pooled-mean null draws random CLUSTERS. But the strongest metric,
# best-per-(feature, training-strain), needs FEATURES (a feature = group of clusters
# with identical presence/absence). So this null draws random FEATURES matched on the
# selected FEATURE COUNT (Option 1), reconstructing the FULL feature space from the
# presence_absence matrix exactly as the model does (hash identical P/A columns).
#
# Emits per held-out strain, per draw:
#   best_per_feature  best hit per (feature, training-strain), then mean  [the -0.382 axis]
#   pooled_mean       mean over all hits                                  [continuity]
# Leak-filtered (training strain must be fold-external).
#
# Output: struct_NN_featnull.csv  (draw_id -1 = observed selected features; 0..M-1 = null)

set -euo pipefail

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
PHYLO=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results
ANNOT=/global/scratch/users/anoonan/BRaVE/struct_cv/annotate
SC=$PHYLO/shared_clustering/strain
POOL_MAP="$SC/clusters.tsv"                     # <rep>\t<member>
PA="$SC/presence_absence_matrix.csv"           # rows=strains, cols=clusters (confirm!)
HELDOUT_DIR=/global/scratch/users/anoonan/BRaVE/struct_cv/fold_heldout
OUT=/global/scratch/users/anoonan/BRaVE/struct_cv/featnull
mkdir -p "$OUT"

M="${M:-1000}"; SEED="${SEED:-20260721}"
# P/A matrix orientation: PA_ORIENT=strain_rows (rows=strains, cols=clusters, col0=strain)
# is the default from feature_selection_optimized. If transposed, set PA_ORIENT=cluster_rows.
PA_ORIENT="${PA_ORIENT:-strain_rows}"

# ---- Step 0: reconstruct FULL feature space from the P/A matrix (once) ----
FULL_FEATURES="$OUT/full_features.tsv"          # feature_id <TAB> cluster
if [ ! -s "$FULL_FEATURES" ]; then
  command -v python3 >/dev/null || export PATH="$HOME/.conda/envs/genophi/bin:$PATH"
  echo "[$(date)] collapsing full P/A matrix -> feature space ..."
  python3 - "$PA" "$FULL_FEATURES" "$PA_ORIENT" <<'PY'
import sys, csv
pa, out, orient = sys.argv[1:4]
# Read P/A; hash each CLUSTER's presence/absence vector across strains; group identical.
# feature_selection_optimized: columns = clusters, rows = strains, binary (>0 -> 1).
import numpy as np
with open(pa) as fh:
    r = csv.reader(fh)
    header = next(r)
    def present(v):
        # model uses (presence_absence > 0). counts are ints; be robust to floats/blanks.
        try: return '1' if float(v) > 0 else '0'
        except ValueError: return '0'
    if orient == "strain_rows":
        clusters = header[1:]                       # col0 = strain id
        cols = {c: [] for c in clusters}
        for row in r:
            for c, v in zip(clusters, row[1:]):
                cols[c].append(present(v))
        patterns = {c: "".join(cols[c]) for c in clusters}
    else:  # cluster_rows: col0 = cluster id, rest = strains
        patterns = {}
        for row in r:
            cl = row[0]
            patterns[cl] = "".join(present(v) for v in row[1:])
# group clusters by identical pattern -> feature
hash2cl = {}
for cl, pat in patterns.items():
    hash2cl.setdefault(pat, []).append(cl)
with open(out, "w") as o:
    for idx, group in enumerate(hash2cl.values()):
        for cl in group:
            o.write(f"sc_{idx}\t{cl}\n")
print(f"full features: {len(hash2cl)}  clusters: {len(patterns)}", file=sys.stderr)
PY
fi
nfeat_full=$(cut -f1 "$FULL_FEATURES" | sort -u | wc -l)
echo "[$(date)] full feature space: $nfeat_full features"

# ---- per fold ----
for i in $(seq -w 1 16); do
  f=struct_${i}
  val="$STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv"
  ov="$ANNOT/$f/strain_predictive_feature_overview.csv"
  heldout="$HELDOUT_DIR/${f}_heldout.txt"
  dst="$OUT/${f}_featnull.csv"
  if [ -s "$dst" ]; then echo "[$f] done, skip"; continue; fi
  for x in "$val" "$ov" "$heldout" "$FULL_FEATURES" "$POOL_MAP"; do
    [ -s "$x" ] || { echo "[$f] MISSING $x"; continue 2; }
  done
  echo "[$(date)] [$f] streaming $(du -h "$val"|cut -f1) ..."
  python3 - "$f" "$val" "$POOL_MAP" "$ov" "$FULL_FEATURES" "$heldout" "$dst" "$M" "$SEED" <<'PY'
import sys, csv, random
fold, val, pool_map, ov, full_feat, heldout_file, dst, M, seed = sys.argv[1:10]
M = int(M); seed = int(seed)

held = set(l.strip() for l in open(heldout_file) if l.strip())
prot2clu = {}
with open(pool_map) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            prot2clu[p[1]] = p[0]

# full feature space: feature -> clusters ; and the list of features to draw from
feat2cl = {}
with open(full_feat) as fh:
    for line in fh:
        a = line.rstrip("\n").split("\t")
        if len(a) >= 2:
            feat2cl.setdefault(a[0], []).append(a[1])
all_features = list(feat2cl)

# observed selected features for this fold (from annotate overview): feature -> clusters
sel_feat2cl = {}
with open(ov) as fh:
    for r in csv.DictReader(fh):
        if r.get("Feature") and r.get("cluster"):
            sel_feat2cl.setdefault(r["Feature"], set()).add(r["cluster"])
k = len(sel_feat2cl)
print(f"[{fold}] selected features k={k}  full pool={len(all_features)}", file=sys.stderr)

# Build, per draw, a map cluster -> feature_label (namespaced by draw so best-per-feature
# groups correctly). draw -1 = observed (uses selected feature grouping).
# cluster -> list of (draw_id, feature_key)
clu_draws = {}
def add(draw_id, feat_key, clusters):
    for c in clusters:
        clu_draws.setdefault(c, []).append((draw_id, feat_key))
for ft, cls in sel_feat2cl.items():
    add(-1, f"-1:{ft}", cls)
rng = random.Random(seed + int(fold.split('_')[-1]))
for d in range(M):
    for ft in rng.sample(all_features, k):
        add(d, f"{d}:{ft}", feat2cl[ft])

def strain_of(x): return x.split("::", 1)[0]

# accumulators
# best per (draw, strain, feature_key, train_strain) -> best ident   (for best_per_feature)
best = {}          # key -> best ident
pooled = {}        # (draw,strain) -> [sum,n]
nrow = 0
with open(val) as fh:
    for line in fh:
        p = line.split("\t")
        if len(p) < 4: continue
        c = prot2clu.get(p[1])
        if c is None: continue
        entries = clu_draws.get(c)
        if not entries: continue
        tstrain = strain_of(p[1])
        if tstrain in held: continue        # leak filter
        try: ident = float(p[3])
        except ValueError: continue
        hstrain = strain_of(p[0])
        for (d, fkey) in entries:
            pk = (d, hstrain)
            e = pooled.get(pk); pooled[pk] = [ident,1] if e is None else [e[0]+ident,e[1]+1]
            bk = (d, hstrain, fkey, tstrain)
            b = best.get(bk)
            if b is None or ident > b: best[bk] = ident
        nrow += 1
print(f"[{fold}] kept rows={nrow}", file=sys.stderr)

# best_per_feature: for each (draw,strain), mean over the best-per-(feature,train-strain)
from collections import defaultdict
bpf_vals = defaultdict(list)
for (d, s, fkey, t), v in best.items():
    bpf_vals[(d, s)].append(v)

with open(dst, "w") as o:
    w = csv.writer(o)
    w.writerow(["fold","draw_id","strain","best_per_feature","pooled_mean","n_hits"])
    for (d, s), (sm, n) in pooled.items():
        vals = bpf_vals.get((d, s), [])
        bpf = sum(vals)/len(vals) if vals else ""
        w.writerow([fold, d, s, f"{bpf:.5f}" if bpf!="" else "", f"{sm/n:.5f}", n])
print(f"[{fold}] wrote {dst}", file=sys.stderr)
PY
  echo "[$(date)] [$f] -> $(wc -l < "$dst") rows"
done
echo "=== DONE ==="; du -sh "$OUT"