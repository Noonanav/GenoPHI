#!/usr/bin/env python3
"""Step 0: turn a raw phenotype/fitness matrix into the layout GenoPHI expects.

Assay tables almost never arrive in the right shape. RB-TnSeq and CRISPRi
exports are usually **genes as rows, samples as columns** — the transpose of
what the modeling code wants — and they carry non-sample columns (media
controls, blanks), free-text target names that become directory names during
independent modeling, and gaps that must be resolved before a regression can
run at all.

GenoPHI wants, for every workflow:

    sample_id, target_1, target_2, ...
    Bxb1,      0.42,     7.31
    DMS3,      -0.11,    0.05

one row per genome, one column per thing to predict.

This script performs that conversion explicitly and reports what it did, so the
transform is reviewable rather than buried in a notebook cell. Its output feeds
straight into prepare_inputs.py.

Example (a genes-as-rows RB-TnSeq matrix with an LB control column)
-------------------------------------------------------------------
    python reshape_matrix.py \
        --input 20260710_PA14_s_pos.csv \
        --samples_in columns \
        --drop_samples LB \
        --output pseudomonas_fitness.csv
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

# Column/sample names that are usually assay controls rather than genomes.
CONTROL_HINTS = ('lb', 'control', 'ctrl', 'blank', 'media', 'medium', 'nc',
                 'pc', 'mock', 'untreated', 'no_phage', 'nophage', 'wt')


def sanitize(name):
    """Make a target name safe to use as a directory name.

    Independent modeling writes one directory per target, so a target called
    ``MSMEG5483 |  | MspA family porin`` becomes a path. Spaces, pipes and
    slashes have to go; the original is preserved in the name map.
    """
    s = str(name).strip()
    s = re.sub(r'[\\/|]+', '_', s)
    s = re.sub(r'\s+', '_', s)
    s = re.sub(r'[^A-Za-z0-9._-]', '', s)
    s = re.sub(r'_+', '_', s).strip('_.')
    return s or 'target'


def dedupe(names):
    """Ensure sanitized names stay unique."""
    seen, out = {}, []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f'{n}__{seen[n]}')
        else:
            seen[n] = 0
            out.append(n)
    return out


def main():
    ap = argparse.ArgumentParser(
        description='Reshape a raw phenotype/fitness matrix into the '
                    'samples-as-rows layout GenoPHI expects.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--input', required=True, help='Raw matrix CSV (or TSV)')
    ap.add_argument('--output', required=True, help='Where to write the reshaped CSV')
    ap.add_argument('--samples_in', required=True, choices=['rows', 'columns'],
                    help="Where the SAMPLES (genomes) are in the input. Assay "
                         "exports are usually 'columns' (genes as rows).")
    ap.add_argument('--sep', default=',', help="Input delimiter (use '\\t' for TSV)")
    ap.add_argument('--id_column', default=None,
                    help='Name of the id column/row header (default: first column)')
    ap.add_argument('--sample_column', default='phage',
                    help='Name for the sample id column in the output')
    ap.add_argument('--drop_samples', default='',
                    help='Comma-separated samples to remove (media controls, blanks)')
    ap.add_argument('--drop_targets', default='',
                    help='Comma-separated target columns to remove')
    ap.add_argument('--fillna', default='0',
                    help="Value for missing cells, or 'error' to refuse them. "
                         'Regression needs a dense matrix; 0 means "no confident '
                         'effect", which is the right reading for low-coverage cells.')
    ap.add_argument('--top_n_targets', type=int, default=0,
                    help='Keep only the N most responsive targets (0 = keep all)')
    ap.add_argument('--top_by', default='max', choices=['max', 'absmax', 'var', 'nonzero'],
                    help="Ranking for --top_n_targets. 'max' (per-target maximum) "
                         "surfaces strong positive responders; 'absmax' also "
                         'surfaces single extreme negatives, which is usually not '
                         'what you want.')
    ap.add_argument('--no_sanitize', action='store_true',
                    help='Keep original target names verbatim (they become '
                         'directory names during independent modeling)')
    args = ap.parse_args()

    sep = '\t' if args.sep in ('\\t', 'tab') else args.sep
    df = pd.read_csv(args.input, sep=sep)
    print(f'read {args.input}: {df.shape[0]} rows x {df.shape[1]} columns')

    id_col = args.id_column or df.columns[0]
    if id_col not in df.columns:
        sys.exit(f"ERROR: id column '{id_col}' not found. Columns: {list(df.columns)[:10]}")
    df = df.set_index(id_col)

    if args.samples_in == 'columns':
        print(f"  samples are COLUMNS -> transposing "
              f"({df.shape[1]} samples x {df.shape[0]} targets)")
        df = df.T
    else:
        print(f'  samples are ROWS ({df.shape[0]} samples x {df.shape[1]} targets)')

    df.index = df.index.astype(str).str.strip()
    df.columns = [str(c).strip() for c in df.columns]

    # --- controls -----------------------------------------------------------
    drop = [s.strip() for s in args.drop_samples.split(',') if s.strip()]
    unknown = [s for s in drop if s not in df.index]
    if unknown:
        sys.exit(f'ERROR: --drop_samples names not present: {unknown}')
    if drop:
        df = df.drop(index=drop)
        print(f'  dropped {len(drop)} sample(s): {drop}')

    suspects = [s for s in df.index if str(s).strip().lower() in CONTROL_HINTS]
    if suspects:
        print(f'  WARNING likely control sample(s) still present: {suspects}')
        print('          these are not genomes; drop them with --drop_samples')

    dt = [t.strip() for t in args.drop_targets.split(',') if t.strip()]
    if dt:
        missing = [t for t in dt if t not in df.columns]
        if missing:
            sys.exit(f'ERROR: --drop_targets names not present: {missing}')
        df = df.drop(columns=dt)
        print(f'  dropped {len(dt)} target(s)')

    # --- numeric + missing --------------------------------------------------
    df = df.apply(pd.to_numeric, errors='coerce')
    n_missing = int(df.isna().sum().sum())
    if n_missing:
        if args.fillna == 'error':
            sys.exit(f'ERROR: {n_missing} missing cell(s) and --fillna error. '
                     'Regression needs a dense matrix.')
        fill = float(args.fillna)
        print(f'  filled {n_missing} missing cell(s) '
              f'({100 * n_missing / df.size:.2f}%) with {fill}')
        df = df.fillna(fill)
    else:
        print('  no missing cells')

    # --- target selection ---------------------------------------------------
    if args.top_n_targets and args.top_n_targets < df.shape[1]:
        if args.top_by == 'max':
            score = df.max(axis=0)
        elif args.top_by == 'absmax':
            score = df.abs().max(axis=0)
        elif args.top_by == 'var':
            score = df.var(axis=0)
        else:
            score = (df != 0).sum(axis=0)
        keep = list(score.sort_values(ascending=False).head(args.top_n_targets).index)
        print(f'  kept top {len(keep)} of {df.shape[1]} targets by {args.top_by}')
        df = df[keep]

    # --- names --------------------------------------------------------------
    original = list(df.columns)
    if args.no_sanitize:
        clean = original
    else:
        clean = dedupe([sanitize(c) for c in original])
        changed = sum(a != b for a, b in zip(original, clean))
        if changed:
            print(f'  sanitized {changed} target name(s) for use as directory names')
    df.columns = clean

    dupes = df.index[df.index.duplicated()].tolist()
    if dupes:
        sys.exit(f'ERROR: duplicate sample ids after reshape: {sorted(set(dupes))}')

    # --- write --------------------------------------------------------------
    out = df.reset_index().rename(columns={df.index.name or 'index': args.sample_column})
    out.columns = [args.sample_column] + list(df.columns)
    out.to_csv(args.output, index=False)

    if not args.no_sanitize:
        map_path = os.path.splitext(args.output)[0] + '_target_names.csv'
        pd.DataFrame({'target': clean, 'original': original}).to_csv(map_path, index=False)
        print(f'  target name map -> {map_path}')

    vals = df.values.astype(float)
    print()
    print(f'wrote {args.output}')
    print(f'  {df.shape[0]} samples x {df.shape[1]} targets')
    print(f'  range {np.nanmin(vals):.3f} to {np.nanmax(vals):.3f}')
    print(f'  {100 * np.mean(np.abs(vals) < 1):.1f}% of cells within +/-1 of zero')
    strong = 100 * np.mean(np.abs(vals) > 4)
    print(f'  {strong:.1f}% of cells |value| > 4')
    if np.mean(np.abs(vals) < 1) > 0.4:
        print()
        print('  NOTE this matrix is zero-inflated. Report strong-responder')
        print('       detection (rank-AUROC, AUPR, precision@K), NOT R2 or')
        print('       Spearman, and do not log-transform -- the outliers are')
        print('       the signal. See README, "Regression".')
    print()
    print('Next: prepare_inputs.py --phenotype_table '
          f'{args.output} --sample_column {args.sample_column} '
          '--task_type regression --host_pattern <stem>')
    return 0


if __name__ == '__main__':
    sys.exit(main())
