#!/usr/bin/env python3
"""
Build a strain->dataset metadata CSV for leave-one-dataset-out CV
(genophi nested-cv --cv_mode group --group_column dataset).

Collapses the per-interaction 'dataset' column in the combined interaction matrix
to one row per strain, and sanity-checks that each strain belongs to exactly one
dataset (group mode is one-label-per-strain; a strain in two datasets would be
silently assigned to one).

Edit the CONFIG paths, then run:  python build_dataset_metadata.py
"""
import pandas as pd

# ============================ CONFIG ============================
MATRIX = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
OUT = "/global/home/groups/pc_phiml/embeddings/combined/strain_dataset_metadata.csv"
STRAIN_COL = "strain"
DATASET_COL = "dataset"
# ===============================================================


def main():
    df = pd.read_csv(MATRIX, dtype=str)
    for c in (STRAIN_COL, DATASET_COL):
        if c not in df.columns:
            raise SystemExit(f"Column '{c}' not in matrix ({list(df.columns)}).")

    # One (strain, dataset) per strain.
    pairs = df[[STRAIN_COL, DATASET_COL]].drop_duplicates()

    # Sanity: a strain mapping to >1 dataset is ambiguous for group mode.
    counts = pairs.groupby(STRAIN_COL)[DATASET_COL].nunique()
    multi = counts[counts > 1]
    if len(multi):
        print(f"WARNING: {len(multi)} strain(s) appear in MULTIPLE datasets:")
        for s in multi.index[:20]:
            ds = sorted(pairs.loc[pairs[STRAIN_COL] == s, DATASET_COL])
            print(f"  {s}: {ds}")
        print("Group mode would assign each to ONE dataset arbitrarily. Resolve "
              "these (or use predefined folds) before running.\n")

    meta = pairs.drop_duplicates(subset=[STRAIN_COL]).sort_values(STRAIN_COL)
    meta.columns = ['strain', 'dataset']
    meta.to_csv(OUT, index=False)

    print(f"Strains: {meta['strain'].nunique()}  Datasets: {meta['dataset'].nunique()}")
    print("\nDataset sizes (strains per dataset = held-out fold size):")
    print(meta['dataset'].value_counts().to_string())
    print(f"\nWrote: {OUT}")


if __name__ == "__main__":
    main()
