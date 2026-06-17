#!/usr/bin/env python3
"""
Build leave-one-genus-out CORNER folds from a full interaction matrix.

For each genus (value in the 'dataset' column):
  validation_strains = that genus's strains
  validation_phages  = that genus's phages
  training_strains   = all OTHER genera's strains
  training_phages    = all OTHER genera's phages

Writes, under OUT_ROOT, one directory per genus in the corner-fold convention
(the same format run_corner_fold_from_shared / submit_corner_logo_slurm.py
consume):

  OUT_ROOT/fold_<genus>/
    training_strains.csv  training_phages.csv
    validation_strains.csv validation_phages.csv
    held_out_group.txt

This is the full-dataset equivalent of cross_genus_splits/logo (uncapped).

Edit CONFIG, then:  python build_corner_folds.py
"""
import os
import re
import pandas as pd

# ============================ CONFIG ============================
MATRIX = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
OUT_ROOT = "/global/home/groups/pc_phiml/set_transformer/manuscript/genophi_corner_folds_FULL"
STRAIN_COL = "strain"
PHAGE_COL = "phage"
DATASET_COL = "dataset"
# ===============================================================


def safe_name(s):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s.strip())


def main():
    df = pd.read_csv(MATRIX, dtype={STRAIN_COL: str, PHAGE_COL: str})
    for c in (STRAIN_COL, PHAGE_COL, DATASET_COL):
        if c not in df.columns:
            raise SystemExit(f"Column '{c}' not in matrix ({list(df.columns)}).")
    df[STRAIN_COL] = df[STRAIN_COL].astype(str)
    df[PHAGE_COL] = df[PHAGE_COL].astype(str)

    # Sanity: each strain / phage should belong to exactly one genus.
    for axis in (STRAIN_COL, PHAGE_COL):
        multi = df.groupby(axis)[DATASET_COL].nunique()
        bad = multi[multi > 1]
        if len(bad):
            print(f"WARNING: {len(bad)} {axis}(s) span multiple datasets: {list(bad.index[:10])}")

    genera = sorted(df[DATASET_COL].unique())
    os.makedirs(OUT_ROOT, exist_ok=True)
    print(f"{len(genera)} genera: {genera}")

    for g in genera:
        in_g = df[df[DATASET_COL] == g]
        out_g = df[df[DATASET_COL] != g]
        val_strains = sorted(in_g[STRAIN_COL].unique())
        val_phages = sorted(in_g[PHAGE_COL].unique())
        train_strains = sorted(out_g[STRAIN_COL].unique())
        train_phages = sorted(out_g[PHAGE_COL].unique())

        fold_dir = os.path.join(OUT_ROOT, f"fold_{safe_name(g)}")
        os.makedirs(fold_dir, exist_ok=True)
        pd.DataFrame({'strain': train_strains}).to_csv(os.path.join(fold_dir, 'training_strains.csv'), index=False)
        pd.DataFrame({'phage': train_phages}).to_csv(os.path.join(fold_dir, 'training_phages.csv'), index=False)
        pd.DataFrame({'strain': val_strains}).to_csv(os.path.join(fold_dir, 'validation_strains.csv'), index=False)
        pd.DataFrame({'phage': val_phages}).to_csv(os.path.join(fold_dir, 'validation_phages.csv'), index=False)
        with open(os.path.join(fold_dir, 'held_out_group.txt'), 'w') as fh:
            fh.write(g + "\n")

        # Corner size + class balance (the scorable held-out pairs).
        corner = in_g  # rows where strain in val AND phage in val == this genus's rows
        bal = corner['interaction'].value_counts().to_dict()
        print(f"  {g:14s} val {len(val_strains)}s x {len(val_phages)}p | "
              f"train {len(train_strains)}s x {len(train_phages)}p | "
              f"corner rows={len(corner)} balance={bal}")

    print(f"\nWrote corner folds under: {OUT_ROOT}")


if __name__ == "__main__":
    main()
