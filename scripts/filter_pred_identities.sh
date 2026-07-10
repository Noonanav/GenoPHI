#!/bin/bash
#SBATCH --job-name=filter_pred_ident
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Extract within-family % identities among GenoPHI predictive-feature proteins, for
# a PEQ-style (sharing x identity) strain distance restricted to predictive features.
#
# assigned_clusters.tsv is the all-vs-all MMseqs2 search table (55 GB, 488M rows):
#   col1=query protein, col2=target protein, col4=%identity.
# We stream it ONCE and keep only rows where BOTH proteins are predictive-feature
# proteins AND assigned to the SAME predictive feature (sc_*). Output is tiny.
#
# Inputs are label-independent shared clustering (phylo run) + the struct folds'
# selected_features.csv (protein -> feature map).

set -euo pipefail

PHYLO=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results
STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
OUTDIR=/global/scratch/users/anoonan/BRaVE/struct_cv/pred_identities
ASSIGNED=$PHYLO/shared_clustering/strain/assigned_clusters.tsv

mkdir -p "$OUTDIR"

# 1. protein -> feature map, per fold (a protein may map to different sc_* in
#    different folds, so we keep fold in the key). Build one combined map:
#       protein <TAB> fold::feature
#    We only need same-feature-same-fold pairs, so encode fold+feature together.
echo "[$(date)] building protein->feature map from struct selected_features.csv ..."
: > "$OUTDIR/protein_feature.tsv"
for sf in "$STRUCT"/folds/struct_*/strain/features/selected_features.csv; do
  fold=$(echo "$sf" | sed -E 's#.*/(struct_[0-9]+)/.*#\1#')
  # selected_features.csv: Feature,Cluster_Label  (skip header)
  tail -n +2 "$sf" | awk -F, -v f="$fold" '{print $2"\t"f"::"$1}'
done | sort -u > "$OUTDIR/protein_feature.tsv"
echo "[$(date)] map rows: $(wc -l < "$OUTDIR/protein_feature.tsv")"

# 2. stream the 55 GB search table once.
#    Keep a row iff both query and target are predictive proteins that share at
#    least one (fold::feature). Emit: query, target, identity, fold::feature.
echo "[$(date)] streaming assigned_clusters.tsv ($(du -h "$ASSIGNED" | cut -f1)) ..."
awk -F'\t' '
  # pass 1: protein -> set of fold::feature labels (space-joined)
  NR==FNR { pf[$1] = (($1 in pf) ? pf[$1] " " $2 : $2); next }
  # pass 2: search table. col1=query col2=target col4=identity
  ($1 in pf) && ($2 in pf) {
    nq = split(pf[$1], q, " ")
    # build a set of query labels
    delete qs
    for (i = 1; i <= nq; i++) qs[q[i]] = 1
    nt = split(pf[$2], t, " ")
    for (j = 1; j <= nt; j++) {
      if (t[j] in qs) { print $1"\t"$2"\t"$4"\t"t[j]; break }
    }
  }
' "$OUTDIR/protein_feature.tsv" "$ASSIGNED" > "$OUTDIR/pred_within_family_identities.tsv"

echo "[$(date)] DONE"
echo "output: $OUTDIR/pred_within_family_identities.tsv"
wc -l "$OUTDIR/pred_within_family_identities.tsv"
du -h  "$OUTDIR/pred_within_family_identities.tsv"
