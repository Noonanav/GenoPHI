#!/bin/bash
#SBATCH --job-name=b2c_curated_orth
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# EXPANDED curated set: 541 published phage-fitness determinant proteins (high+low,
# RB-TnSeq) across 5 backgrounds (BW25113/ECRC98/100/101/102). Supersedes the 71-gene
# BW25113-only set. Same pipeline: MMseqs search -> curated clusters -> per-strain
# identity-to-training (3 metrics, leak-filtered) for the curated gate.
#
# Step 0: search 541 curated proteins vs pool DB -> curated_published_clusters.txt
# Step 1: stream each fold's val search over those clusters -> pooled_mean/best_per_gene/
#         nearest_training per held-out strain (draw -2), leak-filtered.

set -euo pipefail

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
PHYLO=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results
POOL_MAP="$PHYLO/shared_clustering/strain/clusters.tsv"
POOL_DB="$PHYLO/shared_clustering/tmp/strain/mmseqs_db"
HELDOUT_DIR=/global/scratch/users/anoonan/BRaVE/struct_cv/fold_heldout
# SET selects which curated FASTA to ortholog-map: "pub" (541, 5 backgrounds) or "71"
# (original BW25113). Run the job twice, once per SET, to compare both to the all-hits
# versions.
SET="${SET:-pub}"
if [ "$SET" = "71" ]; then
  CURATED_FAA=/global/scratch/users/anoonan/BRaVE/struct_cv/curated/curated_71.faa
  TAG=orth71
else
  CURATED_FAA=/global/scratch/users/anoonan/BRaVE/struct_cv/curated_published/curated_published.faa
  TAG=orthpub
fi
OUT="/global/scratch/users/anoonan/BRaVE/struct_cv/curated_published/ortholog_${SET}"
mkdir -p "$OUT"
MINID=0.4; COV=0.8; THREADS=4

command -v mmseqs >/dev/null 2>&1 || export PATH="$HOME/.conda/envs/genophi/bin:$PATH"
command -v mmseqs >/dev/null 2>&1 || { echo "ERROR: mmseqs not found"; exit 1; }

# ---- Step 0: search 541 -> pool, map to cluster reps ----
CUR="$OUT/curated_published_clusters.txt"
if [ ! -s "$CUR" ]; then
  tmp="$OUT/mmtmp"; mkdir -p "$tmp"; hits="$OUT/search.m8"
  echo "[$(date)] MMseqs search 541 curated-published vs pool DB ..."
  mmseqs createdb "$CURATED_FAA" "$OUT/qdb" 2>&1 | tail -2
  mmseqs search "$OUT/qdb" "$POOL_DB" "$OUT/resdb" "$tmp" \
    --min-seq-id "$MINID" -c "$COV" --cov-mode 1 --threads "$THREADS" 2>&1 | tail -4
  mmseqs convertalis "$OUT/qdb" "$POOL_DB" "$OUT/resdb" "$hits" \
    --format-output "query,target,pident" 2>&1 | tail -2
  python3 - "$POOL_MAP" "$hits" "$CUR" <<'PY'
import sys
pool_map, hits, out = sys.argv[1:4]
prot2rep={}
for line in open(pool_map):
    p=line.rstrip("\n").split("\t")
    if len(p)>=2: prot2rep[p[1]]=p[0]

# ORTHOLOG variant: for each (determinant query, TARGET STRAIN), keep only the SINGLE
# BEST-pident target protein -> the determinant's ortholog in that strain, not all its
# paralogs/homologs. This gives a tighter, ortholog-focused curated set than taking every
# hit >= min-seq-id (which pulls in distant family members).
def strain_of(x): return x.split("::",1)[0]
best = {}   # (query, target_strain) -> (pident, target_protein)
for line in open(hits):
    c=line.rstrip("\n").split("\t")
    if len(c)>=3:
        q, t, pid = c[0], c[1], float(c[2])
        key=(q, strain_of(t))
        cur=best.get(key)
        if cur is None or pid > cur[0]:
            best[key]=(pid, t)
reps=set(); un=0
for (pid, t) in best.values():
    r=prot2rep.get(t)
    if r is None: un+=1
    else: reps.add(r)
open(out,"w").write("\n".join(sorted(reps))+"\n")
print(f"curated ORTHOLOG -> {len(reps)} clusters from {len(best)} best (query,strain) hits "
      f"(unmapped {un})", file=sys.stderr)
PY
fi
k=$(wc -l < "$CUR"); echo "[$(date)] curated-published clusters k=$k"

# ---- Step 1: per fold, 3-metric leak-filtered, draw -2 ----
for i in $(seq -w 1 16); do
  f=struct_${i}
  val="$STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv"
  heldout="$HELDOUT_DIR/${f}_heldout.txt"
  dst="$OUT/${f}_${TAG}_3metric.csv"
  [ -s "$dst" ] && { echo "[$f] skip"; continue; }
  for x in "$val" "$heldout" "$CUR" "$POOL_MAP"; do [ -s "$x" ] || { echo "[$f] MISSING $x"; continue 2; }; done
  echo "[$(date)] [$f] streaming $(du -h "$val"|cut -f1) ..."
  python3 - "$f" "$val" "$POOL_MAP" "$CUR" "$heldout" "$dst" <<'PY'
import sys, csv
fold, val, pool_map, cur_file, heldout_file, dst = sys.argv[1:7]
held=set(l.strip() for l in open(heldout_file) if l.strip())
prot2clu={}
for line in open(pool_map):
    p=line.rstrip("\n").split("\t")
    if len(p)>=2: prot2clu[p[1]]=p[0]
curated=set(l.strip() for l in open(cur_file) if l.strip())
def so(x): return x.split("::",1)[0]
pooled={}; best_gene={}; near={}; nrow=0
with open(val) as fh:
    for line in fh:
        p=line.split("\t")
        if len(p)<4: continue
        c=prot2clu.get(p[1])
        if c is None or c not in curated: continue
        ts=so(p[1])
        if ts in held: continue
        try: ident=float(p[3])
        except ValueError: continue
        hp=p[0]; s=so(hp); key=s
        e=pooled.get(key); pooled[key]=[ident,1] if e is None else [e[0]+ident,e[1]+1]
        bg=best_gene.setdefault(key,{})
        if hp not in bg or ident>bg[hp]: bg[hp]=ident
        nn=near.setdefault(key,{})
        t=nn.get(ts); nn[ts]=[ident,1] if t is None else [t[0]+ident,t[1]+1]
        nrow+=1
print(f"[{fold}] kept={nrow}", file=sys.stderr)
with open(dst,"w") as o:
    w=csv.writer(o); w.writerow(["fold","draw_id","strain","pooled_mean","best_per_gene","nearest_training","n_hits"])
    for s,(sm,n) in pooled.items():
        bg=best_gene[s]; bpg=sum(bg.values())/len(bg)
        nearest=max(v[0]/v[1] for v in near[s].values())
        w.writerow([fold,-2,s,f"{sm/n:.5f}",f"{bpg:.5f}",f"{nearest:.5f}",n])
print(f"[{fold}] wrote {dst}", file=sys.stderr)
PY
  echo "[$(date)] [$f] -> $(wc -l < "$dst") rows"
done
echo "=== DONE ==="; du -sh "$OUT"