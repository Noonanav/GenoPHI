#!/usr/bin/env python3
"""Step 1 of the GenoPHI multi-output workflow: validate inputs and build the genome set.

Nothing downstream is safe to run until the phenotype table and the genome FASTAs
agree on names. This script is the gate:

  1. Validate the phenotype table (id column, duplicate ids, target columns,
     class balance, degenerate/nested targets).
  2. Resolve every sample id to a genome file via the phage_genomes manifest.
  3. Extract CDS translations from each genome (.gbk) to <out>/genomes/<id>.faa.
  4. Re-read what actually landed on disk and assert the FASTA basenames match
     the phenotype table ids exactly.
  5. Write a filtered phenotype CSV containing only samples with a genome, plus
     a human-readable report and a machine-readable JSON summary.

Exit status is 0 only when the prepared inputs are internally consistent.

Example
-------
    python prepare_inputs.py \
        --phenotype_table Binary_table_only_Msp.csv \
        --sample_column Phages \
        --host_pattern yco \
        --output_dir prepared/myco_msp
"""

import argparse
import json
import os
import shutil
import sys
from collections import defaultdict

import pandas as pd

DEFAULT_MANIFEST = '/usr2/people/phages/phage_genomes/manifest.tsv'

# Manifest columns this script depends on.
MF_PATH = 'Genome_file'
MF_NAME = 'Phage'
MF_HOST = 'Host'
MF_LATEST = 'Latest'


class Report:
    """Accumulates findings so every problem is reported, not just the first."""

    def __init__(self):
        self.lines = []
        self.errors = []
        self.warnings = []
        self.data = {}

    def say(self, msg=''):
        self.lines.append(msg)
        print(msg, flush=True)

    def head(self, title):
        self.say()
        self.say(title)
        self.say('-' * len(title))

    def error(self, msg):
        self.errors.append(msg)
        self.say(f'  ERROR   {msg}')

    def warn(self, msg):
        self.warnings.append(msg)
        self.say(f'  WARNING {msg}')

    def ok(self, msg):
        self.say(f'  ok      {msg}')


def load_phenotype(path, sample_column, targets_arg, rep):
    """Read and validate the phenotype table. Returns (df, target_columns)."""
    rep.head(f'Phenotype table: {path}')
    df = pd.read_csv(path)
    rep.say(f'  {len(df)} rows x {len(df.columns)} columns')

    if sample_column not in df.columns:
        # A pandas-exported index shows up as 'Unnamed: 0'; a common case worth naming.
        hint = ''
        if 'Unnamed: 0' in df.columns:
            hint = " (found 'Unnamed: 0' -- pass --sample_column 'Unnamed: 0')"
        rep.error(f"sample column '{sample_column}' not in table{hint}. "
                  f'Columns: {list(df.columns)[:10]}')
        return df, []

    ids = df[sample_column].astype(str).str.strip()
    if (ids != df[sample_column].astype(str)).any():
        rep.warn('sample ids had leading/trailing whitespace; stripped')
    df[sample_column] = ids

    dupes = sorted(ids[ids.duplicated()].unique())
    if dupes:
        rep.error(f'{len(dupes)} duplicate sample id(s): {dupes[:10]}')
    else:
        rep.ok(f'{len(ids)} unique sample ids')

    if targets_arg and targets_arg != 'auto':
        targets = [t.strip() for t in targets_arg.split(',') if t.strip()]
        missing = [t for t in targets if t not in df.columns]
        if missing:
            rep.error(f'target column(s) not in table: {missing}')
            targets = [t for t in targets if t in df.columns]
    else:
        targets = [c for c in df.columns if c != sample_column]
        rep.ok(f'auto-detected {len(targets)} target column(s)')

    return df, targets


def describe_targets(df, targets, task_type, rep):
    """Report per-target balance and flag targets that cannot be modeled."""
    rep.head(f'Targets ({task_type})')
    summary = {}
    for t in targets:
        col = df[t]
        n_missing = int(col.isna().sum())
        if task_type == 'classification':
            vals = sorted(v for v in col.dropna().unique())
            n_pos = int((col == 1).sum())
            summary[t] = {'positives': n_pos, 'n': int(col.notna().sum()),
                          'values': [str(v) for v in vals], 'missing': n_missing}
            rep.say(f'  {t}: {n_pos}/{col.notna().sum()} positive, values={vals}')
            if len(vals) < 2:
                rep.error(f"target '{t}' is constant ({vals}) -- cannot be modeled")
            elif set(vals) - {0, 1}:
                rep.warn(f"target '{t}' is not strictly 0/1: {vals}")
            elif n_pos < 5 or (col.notna().sum() - n_pos) < 5:
                rep.warn(f"target '{t}' has <5 in one class -- folds may drop it "
                         f'(see min_features / single-class handling in the README)')
        else:
            desc = col.describe()
            n_zero = int((col == 0).sum())
            summary[t] = {'n': int(col.notna().sum()), 'zeros': n_zero,
                          'min': float(desc['min']), 'max': float(desc['max']),
                          'missing': n_missing}
            rep.say(f'  {t}: n={int(desc["count"])} zeros={n_zero} '
                    f'range=[{desc["min"]:.3g}, {desc["max"]:.3g}]')
            if n_missing:
                rep.error(f"target '{t}' has {n_missing} missing value(s) -- "
                          'MultiRMSE needs a dense matrix; fill with 0')
        if task_type == 'classification' and n_missing:
            rep.error(f"target '{t}' has {n_missing} missing value(s)")
    return summary


def check_target_independence(df, targets, rep):
    """Warn when two targets are identical or one is nested inside another.

    Not an error -- nested targets are modelable -- but per-target metrics and any
    joint-vs-independent comparison mean something different when the label sets
    overlap this heavily, so the user should know before interpreting results.
    """
    if len(targets) < 2:
        return []
    findings = []
    for i, a in enumerate(targets):
        for b in targets[i + 1:]:
            sa = set(df.index[df[a] == 1])
            sb = set(df.index[df[b] == 1])
            if not sa or not sb:
                continue
            if sa == sb:
                findings.append((a, b, 'identical'))
                rep.warn(f"targets '{a}' and '{b}' have identical positive sets")
            elif sb < sa:
                findings.append((b, a, 'nested'))
                rep.warn(f"target '{b}' positives are a strict subset of '{a}' "
                         f'({len(sb)} of {len(sa)})')
            elif sa < sb:
                findings.append((a, b, 'nested'))
                rep.warn(f"target '{a}' positives are a strict subset of '{b}' "
                         f'({len(sa)} of {len(sb)})')
    return findings


def resolve_from_dir(ids, faa_dir, suffix, rep, already):
    """Fall back to a curated directory of per-sample FASTAs.

    Not every genome is in the manifest -- curated or unreleased genomes often
    live in a project directory. Anything found here is used only for samples the
    manifest could not resolve, so the manifest stays authoritative.
    """
    rep.head(f'Local genome directory: {faa_dir}')
    if not os.path.isdir(faa_dir):
        rep.error(f'not a directory: {faa_dir}')
        return {}
    found = {}
    for sid in ids:
        if sid in already:
            continue
        for ext in (suffix, 'faa', 'fasta', 'fa', 'gbk'):
            cand = os.path.join(faa_dir, f'{sid}.{ext}')
            if os.path.exists(cand):
                found[sid] = cand
                break
    if found:
        rep.ok(f'{len(found)} sample(s) resolved here that the manifest lacked: '
               f'{sorted(found)}')
    else:
        rep.say('  nothing additional found here')
    return found


def resolve_genomes(ids, manifest_path, host_pattern, rep):
    """Map each sample id to one genome file using the manifest.

    The manifest's ``Latest`` column only disambiguates samples that have more
    than one row -- many single-row entries leave it blank, so filtering on
    ``Latest == 'latest'`` up front silently discards usable genomes.
    """
    rep.head(f'Manifest: {manifest_path}')
    if str(manifest_path).lower() == 'none':
        rep.say('  skipped (--manifest none)')
        return {}, {'unresolved': list(ids), 'missing_on_disk': []}
    if not os.path.exists(manifest_path):
        rep.error(f'manifest not found: {manifest_path}')
        return {}, {}
    mf = pd.read_csv(manifest_path, sep='\t', dtype=str).fillna('')
    rep.say(f'  {len(mf)} rows')

    for col in (MF_PATH, MF_NAME, MF_HOST):
        if col not in mf.columns:
            rep.error(f"manifest missing required column '{col}'")
            return {}, {}

    if host_pattern:
        sub = mf[mf[MF_HOST].str.contains(host_pattern, case=False, na=False)]
        hosts = sorted(set(sub[MF_HOST]))
        rep.say(f"  host filter '{host_pattern}' -> {len(sub)} rows, "
                f'{len(hosts)} distinct host string(s)')
        if len(hosts) > 1:
            rep.say(f'          host strings: {hosts}')
    else:
        sub = mf

    by_name = defaultdict(list)
    for _, row in sub.iterrows():
        by_name[row[MF_NAME]].append(row)

    resolved, unresolved = {}, []
    multi = 0
    for sid in ids:
        rows = by_name.get(sid)
        if not rows:
            unresolved.append(sid)
            continue
        if len(rows) > 1:
            multi += 1
            latest = [r for r in rows if r[MF_LATEST] == 'latest'] if MF_LATEST in sub.columns else []
            chosen = latest[0] if latest else rows[0]
            if not latest:
                rep.warn(f"'{sid}' has {len(rows)} manifest rows, none marked "
                         f'latest; using the first')
        else:
            chosen = rows[0]
        resolved[sid] = chosen[MF_PATH]

    rep.ok(f'{len(resolved)} of {len(ids)} sample(s) resolved to a genome file')
    if multi:
        rep.say(f'  {multi} sample(s) had multiple manifest rows (Latest used to pick)')
    if unresolved:
        rep.warn(f'{len(unresolved)} sample(s) absent from the manifest: {unresolved}')

    missing_on_disk = {s: p for s, p in resolved.items() if not os.path.exists(p)}
    if missing_on_disk:
        rep.error(f'{len(missing_on_disk)} genome path(s) in the manifest do not '
                  f'exist on disk: {list(missing_on_disk)[:5]}')
        for s in missing_on_disk:
            resolved.pop(s)
    else:
        rep.ok('every resolved genome path exists on disk')

    extras = sorted(set(by_name) - set(ids))
    if extras:
        rep.say(f'  {len(extras)} manifest sample(s) not in the phenotype table '
                f'(unused): {extras[:10]}')

    return resolved, {'unresolved': unresolved, 'missing_on_disk': list(missing_on_disk)}


def extract_faa(resolved, out_dir, rep, force=False):
    """Write <out_dir>/<id>.faa from each genome's CDS translations."""
    from Bio import SeqIO

    rep.head(f'Protein FASTA extraction -> {out_dir}')
    os.makedirs(out_dir, exist_ok=True)
    written, skipped, empty = [], [], []

    for sid, gpath in sorted(resolved.items()):
        dest = os.path.join(out_dir, f'{sid}.faa')
        if os.path.exists(dest) and not force:
            skipped.append(sid)
            continue
        if gpath.endswith(('.faa', '.fasta', '.fa')):
            # Already protein FASTA -- copy it under the sample's name.
            shutil.copyfile(gpath, dest)
            written.append(sid)
            continue
        if not gpath.endswith(('.gbk', '.gb', '.genbank', '.gbff')):
            rep.warn(f"'{sid}' genome is not GenBank or protein FASTA "
                     f'({os.path.basename(gpath)}); cannot use')
            empty.append(sid)
            continue
        records = []
        for rec in SeqIO.parse(gpath, 'genbank'):
            for feat in rec.features:
                if feat.type != 'CDS':
                    continue
                tr = feat.qualifiers.get('translation')
                if not tr:
                    continue
                tag = (feat.qualifiers.get('locus_tag') or
                       feat.qualifiers.get('protein_id') or
                       [f'{sid}_CDS_{len(records) + 1:04d}'])[0]
                records.append(f'>{tag}\n{tr[0]}\n')
        if not records:
            rep.warn(f"'{sid}' yielded 0 CDS translations from {gpath}")
            empty.append(sid)
            continue
        with open(dest, 'w') as fh:
            fh.writelines(records)
        written.append(sid)

    rep.ok(f'{len(written)} written, {len(skipped)} already present '
           f'(use --force to rebuild), {len(empty)} failed')
    return written, skipped, empty


def verify_names(df, sample_column, faa_dir, suffix, rep):
    """The gate: FASTA basenames on disk must equal the phenotype table ids."""
    rep.head('Name reconciliation (table ids vs FASTA files on disk)')
    on_disk = {f[: -(len(suffix) + 1)] for f in os.listdir(faa_dir)
               if f.endswith('.' + suffix)}
    in_table = set(df[sample_column])

    only_table = sorted(in_table - on_disk)
    only_disk = sorted(on_disk - in_table)

    rep.say(f'  {len(in_table)} ids in table, {len(on_disk)} .{suffix} files in {faa_dir}')
    if only_table:
        rep.error(f'{len(only_table)} id(s) with no FASTA: {only_table}')
    if only_disk:
        rep.error(f'{len(only_disk)} FASTA(s) with no table row: {only_disk}')
    if not only_table and not only_disk:
        rep.ok('exact match -- every sample has a genome and vice versa')
    return only_table, only_disk


def main():
    ap = argparse.ArgumentParser(
        description='Validate a phenotype table against the phage genome manifest '
                    'and build the matching protein FASTA set.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--phenotype_table', required=True,
                    help='CSV with a sample id column plus one column per target')
    ap.add_argument('--sample_column', default='phage',
                    help='Name of the sample id column in the phenotype table')
    ap.add_argument('--targets', default='auto',
                    help="Comma-separated target columns, or 'auto' for all "
                         'non-id columns')
    ap.add_argument('--task_type', default='classification',
                    choices=['classification', 'regression'])
    ap.add_argument('--output_dir', required=True,
                    help='Where prepared inputs and the report are written')
    ap.add_argument('--manifest', default=DEFAULT_MANIFEST,
                    help="phage_genomes manifest.tsv, or 'none' to use only "
                         '--faa_source_dir')
    ap.add_argument('--faa_source_dir', default=None,
                    help='Directory of curated per-sample genome files '
                         '(<sample>.faa / .fasta / .gbk), used for samples the '
                         'manifest cannot resolve. Many projects keep genomes '
                         'that are not in the manifest.')
    ap.add_argument('--host_pattern', default='',
                    help="Case-insensitive substring to filter manifest Host "
                         "(e.g. 'yco' for Mycobacterium + Mycolicibacterium). "
                         'Empty matches all hosts.')
    ap.add_argument('--faa_dir', default=None,
                    help='Where to write per-sample .faa (default: <output_dir>/genomes)')
    ap.add_argument('--suffix', default='faa', help='FASTA suffix')
    ap.add_argument('--out_sample_column', default='phage',
                    help='Sample column name in the written phenotype CSV; must '
                         'match what you pass downstream as --sample_column')
    ap.add_argument('--force', action='store_true',
                    help='Re-extract FASTAs that already exist')
    ap.add_argument('--allow_missing_genomes', action='store_true',
                    help='Drop samples that have no genome instead of failing')
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    faa_dir = args.faa_dir or os.path.join(args.output_dir, 'genomes')
    rep = Report()

    rep.say('=' * 72)
    rep.say('GenoPHI multi-output workflow -- input preparation')
    rep.say('=' * 72)

    df, targets = load_phenotype(args.phenotype_table, args.sample_column,
                                 args.targets, rep)
    if not targets:
        rep.error('no usable target columns; nothing to prepare')
        return finish(rep, args, 1)

    target_summary = describe_targets(df, targets, args.task_type, rep)
    nesting = check_target_independence(df, targets, rep) \
        if args.task_type == 'classification' else []

    ids = list(df[args.sample_column])
    resolved, gaps = resolve_genomes(ids, args.manifest, args.host_pattern, rep)
    n_from_manifest = len(resolved)
    if args.faa_source_dir:
        resolved.update(resolve_from_dir(ids, args.faa_source_dir, args.suffix,
                                         rep, resolved))
        rep.say(f'  sources: {n_from_manifest} from manifest, '
                f'{len(resolved) - n_from_manifest} from {args.faa_source_dir}')

    dropped = sorted(set(ids) - set(resolved))
    if dropped:
        if args.allow_missing_genomes:
            rep.warn(f'dropping {len(dropped)} sample(s) with no genome: {dropped}')
        else:
            rep.error(f'{len(dropped)} sample(s) have no genome: {dropped}. '
                      'Re-run with --allow_missing_genomes to drop them.')
            return finish(rep, args, 1, target_summary, nesting, gaps)

    extract_faa(resolved, faa_dir, rep, force=args.force)

    kept = df[df[args.sample_column].isin(resolved)].copy()
    kept = kept[[args.sample_column] + targets]
    kept = kept.rename(columns={args.sample_column: args.out_sample_column})

    only_table, only_disk = verify_names(kept, args.out_sample_column, faa_dir,
                                         args.suffix, rep)

    pheno_out = os.path.join(args.output_dir, 'phenotype.csv')
    kept.to_csv(pheno_out, index=False)
    targets_out = os.path.join(args.output_dir, 'targets.txt')
    with open(targets_out, 'w') as fh:
        fh.write('\n'.join(targets) + '\n')

    rep.head('Prepared inputs')
    rep.say(f'  phenotype matrix : {pheno_out}  ({len(kept)} samples)')
    rep.say(f'  genome dir       : {faa_dir}')
    rep.say(f'  targets file     : {targets_out}  ({len(targets)} targets)')
    rep.say(f'  sample column    : {args.out_sample_column}')

    status = 1 if (rep.errors or only_table or only_disk) else 0
    return finish(rep, args, status, target_summary, nesting, gaps,
                  extra={'n_samples': len(kept), 'targets': targets,
                         'phenotype_matrix': pheno_out, 'genome_dir': faa_dir,
                         'dropped_no_genome': dropped})


def finish(rep, args, status, target_summary=None, nesting=None, gaps=None, extra=None):
    rep.head('Summary')
    rep.say(f'  {len(rep.errors)} error(s), {len(rep.warnings)} warning(s)')
    rep.say('  RESULT: ' + ('READY' if status == 0 else 'NOT READY -- fix the errors above'))

    payload = {'status': 'ready' if status == 0 else 'not_ready',
               'errors': rep.errors, 'warnings': rep.warnings,
               'targets': target_summary or {}, 'nested_targets': nesting or [],
               'genome_gaps': gaps or {}}
    payload.update(extra or {})

    with open(os.path.join(args.output_dir, 'prepare_report.txt'), 'w') as fh:
        fh.write('\n'.join(rep.lines) + '\n')
    with open(os.path.join(args.output_dir, 'prepare_report.json'), 'w') as fh:
        json.dump(payload, fh, indent=2)
    return status


if __name__ == '__main__':
    sys.exit(main())
