#!/usr/bin/env python3
"""
Predict an OOD dataset with a trained all-4-genera pure-k-mer model.

The second half of the train-then-predict pair: kmer_train_all4.py builds the
model once, this applies it to one out-of-distribution dataset. Predicts the
full host x phage grid; score afterwards by merging against the dataset's
<dataset>_pairs.csv on (strain, phage).

The model's predictive k-mers are assigned to the OOD genomes by substring
matching (_assign_kmer_features), which reads sequences per genome FILE and
never touches protein IDs -- so the OOD FASTAs are read from their own
directories, with no shared namespace or symlinking needed.

The cutoff is selected by MCC (_select_best_cutoff_kmer), matching every other
genophi prediction path, rather than being hardcoded.

  python kmer_predict_ood.py --model /path/to/kmer_all4 \
      --host-dir  /path/to/OOD/cellulophaga/host_faa \
      --phage-dir /path/to/OOD/cellulophaga/phage_faa \
      --output    /path/to/OOD/kmer_predictions/cellulophaga
"""
import os
import argparse
import logging

from genophi.workflows.prediction_workflow import run_prediction_workflow
from genophi.workflows.nested_cv_workflow import (
    _assign_kmer_features, _select_best_cutoff_kmer,
)


def genome_ids(directory, suffix='faa'):
    if not os.path.isdir(directory):
        raise SystemExit(f"Missing directory: {directory}")
    dot = f".{suffix}"
    ids = sorted(f[:-len(dot)] for f in os.listdir(directory) if f.endswith(dot))
    if not ids:
        raise SystemExit(f"No .{suffix} files in {directory}")
    return ids


def main():
    ap = argparse.ArgumentParser(
        description="Predict an OOD dataset with a trained all-4 pure-k-mer model.")
    ap.add_argument('--model', required=True,
                    help="kmer_train_all4.py output directory.")
    ap.add_argument('--host-dir', required=True, help="OOD host .faa directory.")
    ap.add_argument('--phage-dir', required=True, help="OOD phage .faa directory.")
    ap.add_argument('--output', required=True, help="Prediction output directory.")
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--suffix', default='faa')
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    hosts = genome_ids(args.host_dir, args.suffix)
    phages = genome_ids(args.phage_dir, args.suffix)
    logging.info(f"{len(hosts)} hosts x {len(phages)} phages "
                 f"= {len(hosts) * len(phages)} predictions.")

    # Model artefacts, in the k-mer workflow's layout (modeling/ subdir).
    best_cutoff = _select_best_cutoff_kmer(args.model)
    model_dir = os.path.join(args.model, 'modeling', 'modeling_results', str(best_cutoff))
    # Per-axis feature maps. kmer_table_workflow writes these as sibling FILES; passing the
    # strain map to the phage axis matches no pc_ features and yields an empty table.
    strain_feature_map = os.path.join(args.model, 'feature_tables', 'selected_features.csv')
    phage_feature_map = os.path.join(args.model, 'feature_tables', 'phage_selected_features.csv')
    select_feature_table = os.path.join(
        args.model, 'modeling', 'feature_selection', 'filtered_feature_tables',
        f'select_feature_table_{best_cutoff}.csv',
    )
    for path in (model_dir, strain_feature_map, phage_feature_map, select_feature_table):
        if not os.path.exists(path):
            raise SystemExit(f"Model artefact missing: {path}")
    logging.info(f"Using cutoff {best_cutoff} ({model_dir}).")

    os.makedirs(args.output, exist_ok=True)
    predict_output_dir = os.path.join(args.output, 'predict_results')
    os.makedirs(predict_output_dir, exist_ok=True)

    # Describe the OOD genomes in the model's k-mer feature vocabulary.
    logging.info("Assigning predictive k-mers to OOD phages...")
    phage_features = _assign_kmer_features(
        args.phage_dir, phages, phage_feature_map, select_feature_table,
        id_col='phage', source_prefix='p', suffix=args.suffix)
    phage_feature_path = os.path.join(args.output, 'phage_feature_table.csv')
    phage_features.to_csv(phage_feature_path, index=False)

    logging.info("Assigning predictive k-mers to OOD hosts...")
    strain_features = _assign_kmer_features(
        args.host_dir, hosts, strain_feature_map, select_feature_table,
        id_col='strain', source_prefix='s', suffix=args.suffix)
    # run_prediction_workflow reads strain feature tables from input_dir.
    strain_features.to_csv(
        os.path.join(args.output, 'strain_feature_table.csv'), index=False)

    logging.info("Predicting host x phage grid...")
    run_prediction_workflow(
        input_dir=args.output,
        phage_feature_table_path=phage_feature_path,
        model_dir=model_dir,
        output_dir=predict_output_dir,
        feature_table=select_feature_table,
        strain_source='strain',
        phage_source='phage',
        threads=args.threads,
    )

    logging.info(f"Predictions written to {predict_output_dir}")
    logging.info("Score by merging strain_median_predictions.csv against "
                 "<dataset>_pairs.csv on (strain, phage).")


if __name__ == "__main__":
    main()
