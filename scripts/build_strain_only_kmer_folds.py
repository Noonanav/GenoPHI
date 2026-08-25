#!/usr/bin/env python3
"""
Convert strain_only splits into the 4-file corner format run_corner_kmer_fold reads.

The strain_only in-distribution folds hold out STRAINS while phages stay shared,
so each fold ships only training_strains.csv and validation_strains.csv. The
k-mer fold runner takes two genome lists per axis, and nothing requires the two
phage lists to be disjoint -- setting both to the full phage set turns the
"corner" into a full row block:

    training quadrant = training_strains x all_phages
    predicted block   = validation_strains x all_phages

which is a strain-level holdout. This is not a leak: phages are shared by
design in this arm, exactly as they are for the pLM strain_only models.

Reads the pLM's own strain_only splits and writes a parallel tree with the two
phage lists added. The source splits are NOT modified -- the protein-family and
pLM results were produced from them.

The phage list is per-fold: the phages appearing in that fold's TRAINING rows,
which is exactly the modeling_phages set run_fold_from_shared derives from the
modeling interaction matrix (nested_cv_workflow.py:731-735). This matters --
using every phage instead would pull the sequences of any validation-only phage
into training feature-building, which the protein-family strain_only run
excludes. Matching it keeps the two comparable.

The same list is used on both sides of the split, since the model's phage
vocabulary is what the held-out strains get predicted against. Predictions
therefore cover validation_strains x those phages, a superset of the assayed
pairs; scoring merges against the matrix on (strain, phage), so unassayed cells
drop out.

Edit CONFIG, then:  python build_strain_only_kmer_folds.py
"""
import os
import pandas as pd

# ============================ CONFIG ============================
SRC_SPLITS = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_only/splits/outer_all4_held"
OUT_SPLITS = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_only/splits_kmer_corner"

MATRIX = "/global/home/groups/pc_phiml/embeddings/combined/combined_interactions_full_4genera.csv"
PHAGE_DIR = "/global/home/groups/pc_phiml/data/combined/phage_AAs"

STRAIN_COL = "strain"
PHAGE_COL = "phage"
SUFFIX = "faa"
# ================================================================


def main():
    if not os.path.isdir(SRC_SPLITS):
        raise SystemExit(f"Missing source splits: {SRC_SPLITS}")

    df = pd.read_csv(MATRIX, dtype={STRAIN_COL: str, PHAGE_COL: str})
    df[STRAIN_COL] = df[STRAIN_COL].astype(str)
    df[PHAGE_COL] = df[PHAGE_COL].astype(str)
    dot = f".{SUFFIX}"
    on_disk = {f[:-len(dot)] for f in os.listdir(PHAGE_DIR) if f.endswith(dot)}
    all_phages = {str(p) for p in df[PHAGE_COL].unique()} & on_disk
    if not all_phages:
        raise SystemExit("No phages found in both the matrix and the FASTA directory.")
    print(f"{len(all_phages)} phages in the matrix with a FASTA")

    folds = sorted(
        d for d in os.listdir(SRC_SPLITS)
        if d.startswith("fold_") and os.path.isdir(os.path.join(SRC_SPLITS, d))
    )
    if not folds:
        raise SystemExit(f"No fold_* dirs under {SRC_SPLITS}")

    os.makedirs(OUT_SPLITS, exist_ok=True)
    for fold in folds:
        src = os.path.join(SRC_SPLITS, fold)
        dst = os.path.join(OUT_SPLITS, fold)
        os.makedirs(dst, exist_ok=True)

        strain_lists = {}
        for role in ("training", "validation"):
            src_file = os.path.join(src, f"{role}_strains.csv")
            if not os.path.exists(src_file):
                raise SystemExit(f"Missing {src_file}")
            strains = pd.read_csv(src_file)
            strains.to_csv(os.path.join(dst, f"{role}_strains.csv"), index=False)
            col = STRAIN_COL if STRAIN_COL in strains.columns else strains.columns[0]
            strain_lists[role] = [str(s) for s in strains[col]]

        # Phages appearing in the TRAINING rows -- exactly the modeling_phages
        # that run_fold_from_shared derives from the modeling interaction matrix
        # (nested_cv_workflow.py:731-735). Using the full phage set instead
        # would pull sequences of any validation-only phage into training
        # feature-building, which the protein-family strain_only run excludes.
        train_rows = df[df[STRAIN_COL].isin(set(strain_lists["training"]))]
        fold_phages = sorted({str(p) for p in train_rows[PHAGE_COL].unique()} & all_phages)
        if not fold_phages:
            raise SystemExit(f"{fold}: no phages in the training rows.")

        # Shared across the split: the model's phage vocabulary is what the
        # held-out strains are predicted against.
        for role in ("training", "validation"):
            pd.DataFrame({PHAGE_COL: fold_phages}).to_csv(
                os.path.join(dst, f"{role}_phages.csv"), index=False)

        dropped = len(all_phages) - len(fold_phages)
        note = f"  ({dropped} phage(s) absent from training rows, excluded)" if dropped else ""
        print(f"{fold}: {len(strain_lists['training'])} training / "
              f"{len(strain_lists['validation'])} validation strains "
              f"x {len(fold_phages)} phages{note}")

    print(f"\nfolds written to {OUT_SPLITS}")
    print("Point submit_indist_strain_only_kmer_slurm.py at these and submit.")


if __name__ == "__main__":
    main()
