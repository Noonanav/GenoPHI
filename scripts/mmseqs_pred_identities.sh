#!/bin/bash
#SBATCH --job-name=mmseqs_pred_ident
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Targeted all-vs-all sequence identity among the GenoPHI PREDICTIVE proteins
# (struct-CV folds), for the PEQ-style sharing x identity feature distance.
#
# We only need identities among predictive proteins (~union of the 16 folds'
# predictive_AA_seqs.faa), NOT the 55 GB genome-wide all-vs-all. Concatenate +
# dedup those FASTAs and run one MMseqs search. Low --min-seq-id so diverged
# within-family pairs are captured (the sparse tail is exactly what matters).
#
# Output: pred_pairwise_identity.tsv  (query, target, pident, qcov, tcov, evalue)
# Small enough to transfer.

set -euo pipefail
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate genophi   # has mmseqs (used by shared clustering)

A=/global/scratch/users/anoonan/BRaVE/struct_cv/annotate
OUT=/global/scratch/users/anoonan/BRaVE/struct_cv/mmseqs_pred
mkdir -p "$OUT"
TMP="$OUT/tmp"; mkdir -p "$TMP"

# 1. union of predictive proteins across folds, dedup by header (keep first seq).
COMBINED="$OUT/predictive_union.faa"
echo "[$(date)] building deduplicated predictive FASTA ..."
awk '/^>/{h=$0; keep=!(h in seen); seen[h]=1} keep' \
  "$A"/struct_*/*_predictive_AA_seqs.faa > "$COMBINED"
echo "  distinct proteins: $(grep -c '^>' "$COMBINED")"

# 2. all-vs-all identity. easy-search does createdb+search+convertalis in one.
#    --min-seq-id 0.4 to MATCH the shared clustering threshold (min_seq_id=0.4):
#    within-cluster members are all >=40% to the representative, so 0.4 captures the
#    relevant within-family pairs; pairs below 0.4 have no hit and are treated as
#    maximally distant in the PEQ metric (the sensible "effectively different" case).
#    -c 0.8 coverage matches shared clustering (coverage=0.8).
echo "[$(date)] mmseqs easy-search (all-vs-all, min-seq-id 0.4) ..."
mmseqs easy-search "$COMBINED" "$COMBINED" \
  "$OUT/pred_pairwise_identity.tsv" "$TMP" \
  --min-seq-id 0.4 -c 0.8 --cov-mode 0 \
  --format-output "query,target,pident,qcov,tcov,evalue" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" -s 6.0

echo "[$(date)] DONE"
wc -l "$OUT/pred_pairwise_identity.tsv"
du -h  "$OUT/pred_pairwise_identity.tsv"
rm -rf "$TMP"
