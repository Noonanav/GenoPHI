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

Every phage in the interaction matrix (that has a FASTA) goes into both lists.
Predicting validation_strains x ALL phages yields a superset of the pairs that
were actually assayed; scoring merges against the matrix on (strain, phage), so
the extra cells drop out.

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

    # Phages present in BOTH the matrix and the FASTA directory -- the same
    # intersection the workflows apply.
    df = pd.read_csv(MATRIX, dtype={STRAIN_COL: str, PHAGE_COL: str})
    dot = f".{SUFFIX}"
    on_disk = {f[:-len(dot)] for f in os.listdir(PHAGE_DIR) if f.endswith(dot)}
    all_phages = sorted({str(p) for p in df[PHAGE_COL].unique()} & on_disk)
    if not all_phages:
        raise SystemExit("No phages found in both the matrix and the FASTA directory.")
    print(f"shared phage set: {len(all_phages)} phages")

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

        for role in ("training", "validation"):
            src_file = os.path.join(src, f"{role}_strains.csv")
            if not os.path.exists(src_file):
                raise SystemExit(f"Missing {src_file}")
            strains = pd.read_csv(src_file)
            strains.to_csv(os.path.join(dst, f"{role}_strains.csv"), index=False)
            # Phages are shared: the same full set on both sides.
            pd.DataFrame({PHAGE_COL: all_phages}).to_csv(
                os.path.join(dst, f"{role}_phages.csv"), index=False)

        n_train = len(pd.read_csv(os.path.join(dst, "training_strains.csv")))
        n_val = len(pd.read_csv(os.path.join(dst, "validation_strains.csv")))
        print(f"{fold}: {n_train} training / {n_val} validation strains "
              f"x {len(all_phages)} phages")

    print(f"\nfolds written to {OUT_SPLITS}")
    print("Point submit_indist_strain_only_kmer_slurm.py at these and submit.")


if __name__ == "__main__":
    main()
