#!/bin/bash
#SBATCH --job-name=annotate_struct
#SBATCH --account=pc_crispriart
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Run `genophi annotate` (run_predictive_proteins_workflow) for each structure-CV
# fold, to resolve selected features -> clusters -> member proteins -> strain, and
# extract the predictive-protein AA sequences. This produces, per fold:
#   strain_predictive_feature_overview.csv   (strain, Feature, cluster, protein_ID)
#   predictive_AA_seqs.faa                    (the selected proteins' sequences)
#
# The FASTA lets us do a SMALL targeted MMseqs all-vs-all on just the predictive
# proteins (a few thousand) for the PEQ-style sharing x identity metric -- no need
# to touch the 55 GB all-vs-all assigned_clusters.tsv.
#
# Uses the best-MCC cutoff's select_feature_table as the feature file (matches the
# model that made the predictions). clusters.tsv comes from the reused phylo
# shared_clustering (label-independent, valid for the struct folds).

set -euo pipefail
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate genophi

STRUCT=/global/scratch/users/anoonan/BRaVE/struct_cv/K20_results
SHARED=/global/scratch/users/anoonan/BRaVE/phylo_cv/K20_results/shared_clustering/strain
CLUSTERS=$SHARED/clusters.tsv
# MUST use the MODIFIED AAs (genome::protein namespace) -- the raw strain_AAs_update
# proteomes have un-namespaced headers that do NOT match the cluster/overview
# protein_IDs, so parse_and_filter_aa_sequences extracts 0 sequences from them.
FASTA_DIR=$SHARED/modified_AAs/strain
OUTROOT=/global/scratch/users/anoonan/BRaVE/struct_cv/annotate

mkdir -p "$OUTROOT"

for fdir in "$STRUCT"/folds/struct_*; do
  f=$(basename "$fdir")
  out="$OUTROOT/$f"
  if [ -s "$out/strain_predictive_feature_overview.csv" ]; then
    echo "[$f] already done, skip"; continue
  fi
  mkdir -p "$out"

  # best cutoff by MCC (tie-break higher cutoff) -- matches _select_best_cutoff
  best=$(tail -n +2 "$fdir/modeling_results/model_performance/model_performance_metrics.csv" \
    | awk -F, '{n=$NF; sub(/cutoff_/,"",n); print $6","n","$NF}' \
    | sort -t, -k1,1gr -k2,2nr | head -1 | cut -d, -f3)
  feat_tbl="$fdir/feature_selection/filtered_feature_tables/select_feature_table_${best}.csv"

  # modeling_dir is the cutoff subdir (production uses modeling_results/cutoff_<best>),
  # NOT the parent. $best is like 'cutoff_3', so this resolves correctly.
  model_dir="$fdir/modeling_results/${best}"

  echo "[$f] best=$best  feature_file=$feat_tbl  modeling_dir=$model_dir"
  genophi annotate \
    --feature_file_path        "$feat_tbl" \
    --feature2cluster_path     "$fdir/strain/features/selected_features.csv" \
    --cluster2protein_path     "$CLUSTERS" \
    --fasta_dir_or_file        "$FASTA_DIR" \
    --modeling_dir             "$model_dir" \
    --feature_assignments_path "$fdir/strain/features/feature_assignments.csv" \
    --output_dir               "$out" \
    --output_fasta             "${f}_predictive_AA_seqs.faa" \
    --strain_column            strain \
    --feature_type             strain \
    || { echo "[$f] annotate FAILED"; continue; }

  echo "[$f] done -> $out"
  ls -lh "$out" | head
done

echo "=== ALL FOLDS DONE ==="
echo "overview rows per fold:"
for f in "$OUTROOT"/struct_*/strain_predictive_feature_overview.csv; do
  echo "  $(wc -l < "$f")  $f"
done
echo "distinct predictive protein_IDs across folds:"
cut -d, -f4 "$OUTROOT"/struct_*/strain_predictive_feature_overview.csv 2>/dev/null \
  | grep -v protein_ID | sort -u | wc -l
