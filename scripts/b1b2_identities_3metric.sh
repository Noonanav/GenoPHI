#!/bin/bash
#SBATCH --job-name=b1b2_3metric
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# 3-METRIC identity-to-training for the gate (B1 random null + observed) and, in a
# second mode, the curated set (B2). Replaces the single pooled-mean metric with THREE
# per-strain summaries, computed IDENTICALLY for observed / random / curated so all
# gate comparisons use the same aggregation:
#
#   pooled_mean       mean identity over ALL hits (protein->training) in the set
#                     -- overall typicality of the strain's proteins to training.
#   best_per_gene     for each of the strain's proteins, take its BEST hit to training,
#                     then average -- how close is each sampled gene to its NEAREST
#                     training match. (the "closest match per gene" metric)
#   nearest_training  mean identity per training STRAIN (over its hits), take the MAX
#                     across training strains -- is there any close ANALOG strain.
#
# All three come from ONE streaming pass per fold. Memory: pooled (draw,strain);
# best_per_gene (draw,strain,heldout_protein); nearest (draw,strain,train_strain).
#
# MODE:
#   MODE=null     -> observed(-1) + M random draws(0..M-1). Draws k clusters from the pool,
#                    k = |predictive clusters| for that fold (annotate overview).
#   MODE=curated  -> single curated set (draw -2). Clusters = CURATED_CLUSTERS file.
#
# Output: struct_NN_<MODE>_3metric.csv
#   fold, draw_id, strain, pooled_mean, best_per_gene, nearest_training, n_hits

set -euo pipefail

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
PHYLO=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results
ANNOT=/global/scratch/users/anoonan/BRaVE/struct_cv/annotate
POOL_MAP="$PHYLO/shared_clustering/strain/clusters.tsv"     # <rep>\t<member>
# Per-fold held-out strain lists (scp'd from niya). REQUIRED for the leak filter:
# clusters.tsv contains ALL 402 strains including held-out/val strains, so a held-out
# strain's proteins can hit OTHER held-out-fold strains (and itself). We drop any hit
# whose TRAINING strain is in the held-out fold -- matching val_identity_metric.py.
HELDOUT_DIR="${HELDOUT_DIR:-/global/scratch/users/anoonan/BRaVE/struct_cv/fold_heldout}"

MODE="${MODE:-null}"                 # null | curated
M="${M:-1000}"                       # random draws (null mode)
SEED="${SEED:-20260721}"
CURATED_CLUSTERS="${CURATED_CLUSTERS:-/global/scratch/users/anoonan/BRaVE/struct_cv/curated/curated_clusters.txt}"
OUT="${OUT:-/global/scratch/users/anoonan/BRaVE/struct_cv/gate_3metric}"
mkdir -p "$OUT"

for i in $(seq -w 1 16); do
  f=struct_${i}
  val="$STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv"
  ov="$ANNOT/$f/strain_predictive_feature_overview.csv"    # predictive clusters (null mode)
  dst="$OUT/${f}_${MODE}_3metric.csv"
  if [ -s "$dst" ]; then echo "[$f] done, skip"; continue; fi
  [ -s "$val" ] && [ -s "$POOL_MAP" ] || { echo "[$f] MISSING val/pool"; continue; }
  if [ "$MODE" = "null" ] && [ ! -s "$ov" ]; then echo "[$f] MISSING overview"; continue; fi
  if [ "$MODE" = "curated" ] && [ ! -s "$CURATED_CLUSTERS" ]; then echo "[$f] MISSING curated"; continue; fi

  heldout="$HELDOUT_DIR/${f}_heldout.txt"
  [ -s "$heldout" ] || { echo "[$f] MISSING heldout list $heldout"; continue; }
  echo "[$(date)] [$f] MODE=$MODE streaming $(du -h "$val"|cut -f1) ..."
  python3 - "$f" "$val" "$POOL_MAP" "$ov" "$dst" "$MODE" "$M" "$SEED" "$CURATED_CLUSTERS" "$heldout" <<'PY'
import sys, csv, random
fold, val, pool_map, ov, dst, mode, M, seed, cur_file, heldout_file = sys.argv[1:11]
M = int(M); seed = int(seed)

# held-out fold strains -- drop any hit whose TRAINING strain is one of these (leak:
# clusters.tsv includes all strains, so held-out strains hit each other / themselves).
held = set(l.strip() for l in open(heldout_file) if l.strip())

# protein -> cluster rep, and universe of reps
prot2clu = {}; universe = set()
with open(pool_map) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            prot2clu[p[1]] = p[0]; universe.add(p[0])
universe = list(universe)

# build the set of draws: draw_id -> set of clusters (as membership via draw_members)
draw_members = {}   # cluster -> list[draw_id]
if mode == "null":
    # predictive clusters (observed, draw -1)
    selected = set()
    with open(ov) as fh:
        r = csv.DictReader(fh)
        for row in r:
            if row.get("cluster"):
                selected.add(row["cluster"])
    k = len(selected)
    for c in selected:
        draw_members.setdefault(c, []).append(-1)
    rng = random.Random(seed + int(fold.split('_')[-1]))
    for d in range(M):
        for c in rng.sample(universe, k):
            draw_members.setdefault(c, []).append(d)
    print(f"[{fold}] null: k={k} draws={M}", file=sys.stderr)
else:  # curated
    curated = set(l.strip() for l in open(cur_file) if l.strip())
    for c in curated:
        draw_members.setdefault(c, []).append(-2)
    print(f"[{fold}] curated clusters={len(curated)}", file=sys.stderr)

def strain_of(pid):
    return pid.split("::", 1)[0]

# accumulators
pooled = {}          # (d,strain) -> [sum,n]
best_gene = {}       # (d,strain) -> {heldout_prot: best_ident}
near = {}            # (d,strain) -> {train_strain: [sum,n]}

nrow = 0
with open(val) as fh:
    for line in fh:
        p = line.split("\t")
        if len(p) < 4:
            continue
        c = prot2clu.get(p[1])
        if c is None:
            continue
        draws = draw_members.get(c)
        if not draws:
            continue
        tstrain = strain_of(p[1])
        if tstrain in held:          # LEAK FILTER: training strain must be fold-external
            continue
        try:
            ident = float(p[3])
        except ValueError:
            continue
        hprot = p[0]; hstrain = strain_of(hprot)
        for d in draws:
            key = (d, hstrain)
            e = pooled.get(key)
            if e is None: pooled[key] = [ident, 1]
            else: e[0] += ident; e[1] += 1
            bg = best_gene.get(key)
            if bg is None: best_gene[key] = {hprot: ident}
            else:
                if hprot not in bg or ident > bg[hprot]: bg[hprot] = ident
            nn = near.get(key)
            if nn is None: near[key] = {tstrain: [ident, 1]}
            else:
                t = nn.get(tstrain)
                if t is None: nn[tstrain] = [ident, 1]
                else: t[0] += ident; t[1] += 1
        nrow += 1
print(f"[{fold}] kept-hit rows={nrow}", file=sys.stderr)

with open(dst, "w") as o:
    w = csv.writer(o)
    w.writerow(["fold","draw_id","strain","pooled_mean","best_per_gene",
                "nearest_training","n_hits"])
    for (d, s), (sm, n) in pooled.items():
        pm = sm / n
        bg = best_gene[(d, s)]
        bpg = sum(bg.values()) / len(bg)
        nn = near[(d, s)]
        nearest = max(v[0]/v[1] for v in nn.values())
        w.writerow([fold, d, s, f"{pm:.5f}", f"{bpg:.5f}", f"{nearest:.5f}", n])
print(f"[{fold}] wrote {dst}", file=sys.stderr)
PY
  echo "[$(date)] [$f] -> $(wc -l < "$dst") rows"
done
echo "=== DONE ($MODE) ==="; du -sh "$OUT"