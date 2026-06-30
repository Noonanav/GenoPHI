#!/usr/bin/env python3
"""
Build a genophi predefined-folds CSV from the in-distribution strain_only splits.

The pLM in-distribution comparison stores each outer fold as a directory of
per-fold CSVs:

  strain_only/splits/outer_all4_held/fold_N/{training_strains,validation_strains}.csv

genophi's predefined-fold entrypoint (run_predefined_fold / load_predefined_folds)
instead reads ONE long-format CSV with columns ``fold_label, strain, role`` where
role is 'modeling' or 'validation'. This converts the former into the latter so
the GenoPHI baseline runs on the EXACT same held-out strain sets as the pLM
models (apples-to-apples).

  python build_indist_strain_only_folds.py            # write folds_strain_only.csv
  python build_indist_strain_only_folds.py --check    # also sanity-report
"""
import os
import glob
import argparse
import pandas as pd

# ============================ CONFIG (edit me) ============================
SPLITS_DIR = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_only/splits/outer_all4_held"
OUT_CSV = "/global/scratch/users/anoonan/set_transformer/manuscript/indist_crossover/strain_only/folds_strain_only.csv"
STRAIN_COL = "strain"
# =========================================================================


def _read_ids(path, col):
    df = pd.read_csv(path)
    # named column if present, else first column (robust to header variants)
    if col in df.columns:
        return df[col].astype(str).tolist()
    return df.iloc[:, 0].astype(str).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help="Print per-fold sanity report.")
    args = ap.parse_args()

    fold_dirs = sorted(glob.glob(os.path.join(SPLITS_DIR, "fold_*")))
    if not fold_dirs:
        raise SystemExit(f"No fold_* dirs under {SPLITS_DIR}")

    rows = []
    for fd in fold_dirs:
        label = os.path.basename(fd)
        train_f = os.path.join(fd, "training_strains.csv")
        val_f = os.path.join(fd, "validation_strains.csv")
        for f in (train_f, val_f):
            if not os.path.exists(f):
                raise SystemExit(f"Missing {f}")
        train = _read_ids(train_f, STRAIN_COL)
        val = _read_ids(val_f, STRAIN_COL)

        overlap = set(train) & set(val)
        if overlap:
            raise SystemExit(f"{label}: {len(overlap)} strains are BOTH modeling and "
                             f"validation (e.g. {sorted(overlap)[:5]}) -- splits are not disjoint.")

        for s in train:
            rows.append({'fold_label': label, STRAIN_COL: s, 'role': 'modeling'})
        for s in val:
            rows.append({'fold_label': label, STRAIN_COL: s, 'role': 'validation'})

        if args.check:
            print(f"  {label}: modeling={len(train)} validation={len(val)} total={len(train)+len(val)}")

    out = pd.DataFrame(rows, columns=['fold_label', STRAIN_COL, 'role'])
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}  ({len(out)} rows, {out['fold_label'].nunique()} folds)")

    if args.check:
        # Validation sets should be disjoint across folds (proper outer CV).
        val = out[out['role'] == 'validation']
        dup = val[STRAIN_COL].duplicated().sum()
        print(f"  validation strains appearing in >1 fold: {dup} (expect 0 for outer CV)")


if __name__ == "__main__":
    main()
