#!/bin/bash
#SBATCH --job-name=b2_curated
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# B2 / A1b -- CURATED determinant set (the receptor/known-determinant control).
#
# Curated set = 71 BW25113 RB-TnSeq phage-fitness hits (receptor + entry + surface
# genes). We expand them to the shared-clustering families via MMseqs2, then compute
# each held-out strain's identity-to-training over the CURATED clusters -- exactly the
# same per-strain metric as the model's SELECTED set (b1 observed) and the RANDOM null.
# Downstream we compare curated vs selected vs the 1000-draw random null, vs AUROC and
# vs phenotypic novelty. Does a hand-curated determinant set track difficulty better
# than the model's discovered features, or than random?
#
# ---- Steps ----
#   0. MMseqs2 search the 71 curated proteins vs the pool proteome; map hits -> the
#      shared-clustering cluster reps they fall in = CURATED CLUSTER SET.
#   1-4. For each fold, stream the ~9 GB val search once; for training-protein hits whose
#      cluster is in the curated set, accumulate per held-out strain -> mean identity.
#   5. Emit fold,draw_id(=-2 for curated),strain,mean_identity,n_hits.  (draw_id -2 keeps
#      it distinct from b1's observed -1 and null 0..999, so files can be concatenated.)
#
# Inputs mirror b1_random_null_identities.sh.

set -euo pipefail

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
PHYLO=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results
CURATED_FAA=/global/scratch/users/anoonan/BRaVE/struct_cv/curated/curated_71.faa   # scp'd from niya
WORK=/global/scratch/users/anoonan/BRaVE/struct_cv/curated
OUT=$WORK
mkdir -p "$OUT"

POOL_MAP="$PHYLO/shared_clustering/strain/clusters.tsv"          # <rep>\t<member>
# Ready-made MMseqs sequence DB of the pool proteins (the sequences that went INTO the
# clustering). Search the 71 curated proteins directly against this -- no createdb needed.
POOL_DB="$PHYLO/shared_clustering/tmp/strain/mmseqs_db"

MINID=0.4          # match the clustering threshold
COV=0.8            # coverage floor for a curated hit
THREADS=4

# ---- Step 0: MMseqs search 71 curated proteins -> pool DB, map to cluster reps ----
CUR_CLUSTERS="$OUT/curated_clusters.txt"
if [ ! -s "$CUR_CLUSTERS" ]; then
  source ~/miniforge3/etc/profile.d/conda.sh && conda activate genophi   # has mmseqs
  tmp="$OUT/mmtmp"; mkdir -p "$tmp"
  hits="$OUT/curated_search.m8"
  echo "[$(date)] MMseqs search 71 curated vs pool DB ($POOL_DB) ..."
  # build a query DB from the 71, search vs the existing pool sequence DB, convert to m8.
  mmseqs createdb "$CURATED_FAA" "$OUT/curated_qdb" 2>&1 | tail -2
  mmseqs search "$OUT/curated_qdb" "$POOL_DB" "$OUT/curated_resdb" "$tmp" \
    --min-seq-id "$MINID" -c "$COV" --cov-mode 1 --threads "$THREADS" 2>&1 | tail -5
  mmseqs convertalis "$OUT/curated_qdb" "$POOL_DB" "$OUT/curated_resdb" "$hits" \
    --format-output "query,target,pident,alnlen,evalue,bits" 2>&1 | tail -2

  # map each hit target protein -> its cluster rep (via POOL_MAP protein->rep), collect
  # the distinct curated cluster reps.
  echo "[$(date)] mapping hit proteins -> cluster reps ..."
  python3 - "$POOL_MAP" "$hits" "$CUR_CLUSTERS" <<'PY'
import sys
pool_map, hits, out = sys.argv[1:4]
prot2rep = {}
with open(pool_map) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            prot2rep[p[1]] = p[0]        # member -> rep
reps = set()
unmapped = 0
with open(hits) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            t = c[1]
            r = prot2rep.get(t)
            if r is None:
                unmapped += 1
            else:
                reps.add(r)
with open(out, "w") as o:
    for r in sorted(reps):
        o.write(r + "\n")
print(f"curated hit proteins -> {len(reps)} distinct clusters (unmapped {unmapped})",
      file=sys.stderr)
PY
fi
k=$(wc -l < "$CUR_CLUSTERS")
echo "[$(date)] curated cluster set size k=$k"

# ---- Steps 1-5: per-fold streaming (same as b1) ----
for i in $(seq -w 1 16); do
  f=struct_${i}
  val="$STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv"
  dst="$OUT/${f}_curated_summary.csv"
  if [ -s "$dst" ]; then echo "[$f] done, skip"; continue; fi
  if [ ! -s "$val" ] || [ ! -s "$POOL_MAP" ] || [ ! -s "$CUR_CLUSTERS" ]; then
    echo "[$f] MISSING input"; continue
  fi
  echo "[$(date)] [$f] streaming $(du -h "$val" | cut -f1) over curated clusters ..."
  python3 - "$f" "$val" "$POOL_MAP" "$CUR_CLUSTERS" "$dst" <<'PY'
import sys, csv
fold, val, pool_map, cur_file, dst = sys.argv[1:6]

# protein -> cluster rep (whole pool)
prot2clu = {}
with open(pool_map) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            prot2clu[p[1]] = p[0]
# curated cluster set
curated = set(l.strip() for l in open(cur_file) if l.strip())
print(f"[{fold}] curated clusters={len(curated)} pool prot={len(prot2clu)}", file=sys.stderr)

acc = {}   # strain -> [sum, n]
nrow = 0
with open(val) as fh:
    for line in fh:
        p = line.split("\t")
        if len(p) < 4:
            continue
        c = prot2clu.get(p[1])
        if c is None or c not in curated:
            continue
        s = p[0].split("::", 1)[0]
        try:
            ident = float(p[3])
        except ValueError:
            continue
        e = acc.get(s)
        if e is None:
            acc[s] = [ident, 1]
        else:
            e[0] += ident; e[1] += 1
        nrow += 1
print(f"[{fold}] curated kept-hit rows={nrow}", file=sys.stderr)

with open(dst, "w") as o:
    w = csv.writer(o)
    w.writerow(["fold", "draw_id", "strain", "mean_identity", "n_hits"])
    for s, (sm, n) in acc.items():
        w.writerow([fold, -2, s, f"{sm/n:.5f}", n])
print(f"[{fold}] wrote {dst}", file=sys.stderr)
PY
  echo "[$(date)] [$f] -> $(wc -l < "$dst") rows"
done

echo "=== DONE ==="
du -sh "$OUT"; wc -l "$OUT"/*_curated_summary.csv | tail -1
