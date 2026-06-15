#!/usr/bin/env python3
"""
Build a predefined leave-one-dataset-out folds CSV from a strain->dataset
metadata file (e.g. the 4-genera datasets). Output feeds genophi nested-cv via
run_predefined_fold (same pipeline as the LOSO-serotype run).

For each dataset D: validation = all strains in D; modeling = all other strains.
Datasets are clean single-label groups, so no cross-reactive/multi-call logic.

Edit CONFIG, then:  python build_dataset_folds.py
"""
import pandas as pd

# ============================ CONFIG ============================
METADATA = "/global/home/groups/pc_phiml/embeddings/combined/strain_dataset_metadata.csv"
OUT_FOLDS = "/global/home/groups/pc_phiml/embeddings/combined/folds_dataset.csv"
STRAIN_COL = "strain"
DATASET_COL = "dataset"
# ===============================================================


def main():
    meta = pd.read_csv(METADATA, dtype=str)
    meta = meta[[STRAIN_COL, DATASET_COL]].dropna().drop_duplicates(subset=[STRAIN_COL])
    all_strains = sorted(meta[STRAIN_COL])
    datasets = sorted(meta[DATASET_COL].unique())

    rows = []
    for d in datasets:
        val = set(meta.loc[meta[DATASET_COL] == d, STRAIN_COL])
        for s in all_strains:
            rows.append({
                'fold_label': d,
                'strain': s,
                'role': 'validation' if s in val else 'modeling',
            })

    folds = pd.DataFrame(rows, columns=['fold_label', 'strain', 'role'])
    folds.to_csv(OUT_FOLDS, index=False)

    print(f"{len(datasets)} dataset folds over {len(all_strains)} strains:")
    for d in datasets:
        nval = (meta[DATASET_COL] == d).sum()
        print(f"  {d:14s} validation={nval:>4}  modeling={len(all_strains)-nval:>4}")
    print(f"\nWrote: {OUT_FOLDS}")


if __name__ == "__main__":
    main()
