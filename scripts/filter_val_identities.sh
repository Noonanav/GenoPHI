#!/bin/bash
#SBATCH --job-name=filter_val_ident
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# For the PEQ-style (sharing x identity) feature distance we need TEST-vs-TRAINING
# protein identities among PREDICTIVE-FEATURE proteins only.
#
# Each fold's model_validation/tmp/assigned_clusters.tsv is a 9 GB search of
# held-out (test) proteins (col1) vs training proteins (col2), identity in col4.
# Stream it once per fold, keep rows where col2 (training protein) is a member of a
# predictive cluster, and emit the FULL protein-by-protein record + its feature:
#   held_out_protein  held_out_strain  feature  training_protein  identity
#
# This keeps ALL protein-by-protein distances so any aggregation is a downstream
# groupby (mean across all predictive-gene hits [primary], max-per-feature then
# mean, min, etc.) -- nothing is collapsed here.
#
# protein->feature map = strain_predictive_protein_info.csv (Feature col1,
# protein_ID col4), which spans the predictive clusters' member proteins.
# ~410k rows / fold from the sample rate (0.5% of 77M) -> ~200 MB total, transferable.

set -euo pipefail

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
ANNOT=/global/scratch/users/anoonan/BRaVE/struct_cv/annotate
OUT=/global/scratch/users/anoonan/BRaVE/struct_cv/val_identities
mkdir -p "$OUT"

for i in $(seq -w 1 16); do
  f=struct_${i}
  asn="$STRUCT/folds/$f/model_validation/tmp/assigned_clusters.tsv"
  info="$ANNOT/$f/strain_predictive_protein_info.csv"
  dst="$OUT/${f}_val_pred_identity.tsv"
  if [ -s "$dst" ]; then echo "[$f] done, skip"; continue; fi
  if [ ! -s "$asn" ] || [ ! -s "$info" ]; then
    echo "[$f] MISSING input (asn:$([ -s "$asn" ]&&echo y||echo n) info:$([ -s "$info" ]&&echo y||echo n))"
    continue
  fi

  # Map: predictive-cluster MEMBER training protein -> its predictive FEATURE (sc_N).
  # From protein_info (Feature,Importance,cluster,protein_ID,...): col1=Feature,
  # col4=protein_ID. A protein can belong to >1 feature; join keeps the first (rare).
  # This is the full protein->feature key so we can group by feature in the parser.
  tail -n +2 "$info" | awk -F, '{print $4"\t"$1}' | sort -u > "/tmp/${f}_prot2feat.tsv"
  echo "[$(date)] [$f] pred member proteins: $(cut -f1 /tmp/${f}_prot2feat.tsv | sort -u | wc -l); streaming $(du -h "$asn"|cut -f1) ..."

  # assigned_clusters.tsv: col1 = held-out (test) protein, col2 = training protein it
  # hit, col4 = identity. Keep rows where col2 (the TRAINING protein) is a predictive
  # member; emit the FULL protein-by-protein record so ANY aggregation is a groupby
  # downstream (mean-across-genes, max-per-feature, etc.):
  #   held_out_protein  held_out_strain  feature  training_protein  identity
  # (Do NOT filter col1 -- the held-out strain's own proteins aren't in the training
  #  predictive-member set.)
  awk -F'\t' '
    NR==FNR { feat[$1]=$2; next }                       # train protein -> feature
    ($2 in feat) {
      s=$1; sub(/::.*/,"",s)                            # held-out strain from col1
      print $1"\t"s"\t"feat[$2]"\t"$2"\t"$4
    }' "/tmp/${f}_prot2feat.tsv" "$asn" > "$dst"

  echo "[$(date)] [$f] kept $(wc -l < "$dst") rows -> $(du -h "$dst"|cut -f1)"
  rm -f "/tmp/${f}_prot2feat.tsv"
done

echo "=== DONE ==="
du -sh "$OUT"
wc -l "$OUT"/*_val_pred_identity.tsv | tail -1
