#!/usr/bin/env python3
"""
Train an all-4-genera model with PURE k-mer features (no held-out genus).

The k-mer counterpart of genophi_combined_full_model: one model over every
strain x phage pair in the combined 4-genera matrix, for predicting OOD datasets
afterwards with kmer_predict_ood.py.

Training and prediction are separate steps on purpose. The alternative -- one
corner fold per OOD dataset -- retrains this same model once per dataset, and
training is by far the expensive half.

Reuses the internals of run_corner_kmer_fold (nested_cv_workflow.py:1136) with
the held-out corner removed: resolve duplicate protein IDs, build the training
FASTAs and protein-mapping CSVs, then run_kmer_table_workflow with modeling
enabled. ignore_families=True, so no MMseqs2 and no clustering.

  python kmer_train_all4.py --output /path/to/kmer_all4 --threads 32
"""
import os
import argparse
import logging

import pandas as pd

from genophi.workflows.kmer_table_workflow import run_kmer_table_workflow
from genophi.workflows.kmer_full_workflow import detect_duplicate_ids, modify_duplicate_ids
from genophi.workflows.nested_cv_workflow import _write_genome_fasta, _write_protein_mapping

DEFAULT_MATRIX = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
DEFAULT_STRAIN_DIR = "/global/home/groups/pc_phiml/data/combined/strain_AAs"
DEFAULT_PHAGE_DIR = "/global/home/groups/pc_phiml/data/combined/phage_AAs"


def genome_ids(directory, suffix='faa'):
    dot = f".{suffix}"
    return {f[:-len(dot)] for f in os.listdir(directory) if f.endswith(dot)}


def main():
    ap = argparse.ArgumentParser(description="Train an all-4-genera pure-k-mer model.")
    ap.add_argument('--matrix', default=DEFAULT_MATRIX)
    ap.add_argument('--strain-dir', default=DEFAULT_STRAIN_DIR)
    ap.add_argument('--phage-dir', default=DEFAULT_PHAGE_DIR)
    ap.add_argument('--output', required=True, help="Model output directory.")
    ap.add_argument('--k', type=int, default=4)
    ap.add_argument('--k-range', action='store_true')
    ap.add_argument('--one-gene', action='store_true')
    ap.add_argument('--num-runs-fs', type=int, default=25)
    ap.add_argument('--num-runs-modeling', type=int, default=50)
    ap.add_argument('--num-features', type=int, default=100)
    ap.add_argument('--threads', type=int, default=32)
    ap.add_argument('--max-ram', type=int, default=40)
    ap.add_argument('--strain-column', default='strain')
    ap.add_argument('--phage-column', default='phage')
    ap.add_argument('--suffix', default='faa')
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    os.makedirs(args.output, exist_ok=True)

    # Genomes present in BOTH the matrix and the FASTA directories -- the same
    # intersection the fold workflows apply.
    df = pd.read_csv(args.matrix, dtype={args.strain_column: str, args.phage_column: str})
    strains = sorted({str(s) for s in df[args.strain_column].unique()}
                     & genome_ids(args.strain_dir, args.suffix))
    phages = sorted({str(p) for p in df[args.phage_column].unique()}
                    & genome_ids(args.phage_dir, args.suffix))
    if not strains or not phages:
        raise SystemExit(
            f"Empty training set: {len(strains)} strains, {len(phages)} phages. "
            "Check that the matrix IDs match the FASTA basenames."
        )
    logging.info(f"Training on {len(strains)} strains x {len(phages)} phages "
                 f"({len(df)} interactions).")

    # Prefix protein IDs (genome::protein) if any collide, so the protein_csv ->
    # k-mer genome attribution stays unambiguous. Prediction is unaffected by
    # this: _assign_kmer_features matches k-mers as substrings of sequences and
    # keys only on genome filename, never on protein ID.
    def _resolve_dir(input_dir, ids, column):
        if detect_duplicate_ids(input_dir, args.suffix, ids, 'directory'):
            logging.info(f"Duplicate {column} protein IDs -> prefixing genome::protein.")
            return modify_duplicate_ids(input_dir, args.output, args.suffix, ids, column)
        return input_dir

    strain_dir = _resolve_dir(args.strain_dir, strains, 'strain')
    phage_dir = _resolve_dir(args.phage_dir, phages, 'phage')

    strain_fasta, n_s = _write_genome_fasta(
        strain_dir, strains, os.path.join(args.output, 'train_strains.faa'), args.suffix)
    phage_fasta, n_p = _write_genome_fasta(
        phage_dir, phages, os.path.join(args.output, 'train_phages.faa'), args.suffix)
    logging.info(f"Wrote {n_s} strain and {n_p} phage sequences.")

    strain_csv = _write_protein_mapping(
        strain_dir, strains, os.path.join(args.output, 'train_strain_proteins.csv'),
        'strain', args.suffix)
    phage_csv = _write_protein_mapping(
        phage_dir, phages, os.path.join(args.output, 'train_phage_proteins.csv'),
        'phage', args.suffix)

    run_kmer_table_workflow(
        strain_fasta=strain_fasta,
        protein_csv=strain_csv,
        k=args.k,
        id_col='strain',
        one_gene=args.one_gene,
        output_dir=args.output,
        k_range=args.k_range,
        phenotype_matrix=args.matrix,
        phage_fasta=phage_fasta,
        protein_csv_phage=phage_csv,
        remove_suffix=False,
        sample_column='strain',
        phenotype_column='interaction',
        modeling=True,
        filter_type='strain',
        num_features=args.num_features,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        method='rfe',
        threads=args.threads,
        task_type='classification',
        ignore_families=True,
        max_ram=args.max_ram,
        use_clustering=True,
        cluster_method='hierarchical',
        use_dynamic_weights=True,
        weights_method='inverse_frequency',
    )

    logging.info(f"Model written to {args.output}")
    logging.info("Predict OOD datasets with kmer_predict_ood.py --model "
                 f"{args.output}")


if __name__ == "__main__":
    main()
