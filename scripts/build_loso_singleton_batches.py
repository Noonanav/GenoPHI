#!/usr/bin/env python3
"""
Build batched held-out folds for the SMALL-serogroup tail of the LOSO-O E. coli set.

Background. The completed LOSO-O run only held out the 19 serogroups with n>=5
strains; all smaller serogroups (singletons + n<5) sat in the always-train pool
and were NEVER held out, so they have no held-out AUROC. To populate the
small-distance end of the AUROC-vs-surface-distance curve (boss's Q1), we hold
out this tail too -- in a few BATCHES (n=1 can't be a stable solo fold), scored
per-serogroup afterward from final_predictions.csv.

Ground truth + leakage rules (locked 2026-07-08):
  * The 228 LOSO strains and their fold labels (final_predictions.csv) are
    authoritative. A tail strain that also appears in the scored-19 set is a
    reconciliation collision (the ECtyper call disagrees with LOSO membership,
    only for O17/O153) -> it BELONGS to the 19, EXCLUDED from the tail.
  * LOCAL leakage rule: a strain is admissible to a batch iff none of its
    O-labels has a member left in that batch's TRAINING set. Concretely: exclude
    any tail strain whose serogroup is one of the scored-19 labels (those train
    in every tail-batch model). Single-call only; multi-call (slash) strains are
    excluded (ambiguous label + most collide with big serogroups).

Batching: whole serogroups are kept intact (never split across batches) and
greedy-packed into N_BATCHES balanced batches (~16-18 strains each -> matches the
prior LOSO fold-size range 5-29, mean 12), so each batch's held-out fraction is
comparable to the existing 19 folds.

Output: writes batch_1..N rows (fold_label, strain, role) to a folds CSV in the
same long format load_predefined_folds / run_predefined_fold consume, so the
existing submit_loso_slurm.py machinery runs them unchanged. load_predefined_folds
derives the modeling set ONLY from explicit role='modeling' rows (it does NOT
auto-assign the complement), so each batch lists BOTH: role='validation' for its
~16 held-out tail strains, and role='modeling' for ALL OTHER usable strains
(the scored-19 + the other batches' tail + the rest of the always-train pool) --
exactly "leave this batch out", the same rule the original 19 folds used.

The full usable-strain universe is taken from the LOSO interaction matrix's
strain set (every strain that had a phenotype + FASTA in the LOSO run), so the
tail batches train on precisely the strains the original run could use.

  python build_loso_singleton_batches.py            # write the folds CSV
  python build_loso_singleton_batches.py --check     # sanity report only, no write
"""
import os
import csv
import argparse
from collections import defaultdict

# ============================ CONFIG (edit me) ============================
# Paths AS ON LAWRENCIUM (build there so the interaction-matrix universe resolves).
# ectyper output.tsv is NOT natively on Lawrencium -- scp it from the workstation:
#   scp <workstation>/.../ectyper_out/output.tsv \
#       anoonan@lrc-login.lbl.gov:/global/scratch/users/anoonan/BRaVE/LOSO/tables/ectyper_output.tsv
ECTYPER_TSV = "/global/scratch/users/anoonan/BRaVE/LOSO/tables/ectyper_output.tsv"
# final_predictions.csv from the completed LOSO-O run = the scored-19 ground truth.
LOSO_PRED = "/global/scratch/users/anoonan/BRaVE/LOSO/loso_O_results/final_predictions.csv"
# Strain AA FASTA dir -- a tail strain is only runnable if its .faa exists here.
AA_DIR = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/strain_AAs_update"
# The LOSO interaction matrix -- its strain set = the usable-strain universe
# (every strain with a phenotype). Each batch's modeling set = this universe
# minus that batch's validation strains. Falls back to AA_DIR listing (with a
# warning) if unreadable.
INTERACTION_MATRIX = "/global/scratch/users/anoonan/BRaVE/manuscript/ecoli/Gaborieau_interaction_matrix_long_mod.csv"

OUT_CSV = "/global/scratch/users/anoonan/BRaVE/LOSO/tables/folds_loso_tail.csv"

N_BATCHES = 6
MAX_N = 5           # "small" = serogroup with fewer than MAX_N strains
STRAIN_COL = "strain"
# =========================================================================

# The 19 serogroups the LOSO-O run scored (n>=5); these TRAIN in every tail batch.
SCORED_19 = {"O104", "O1", "O153", "O15", "O16", "O17", "O18", "O2", "O25",
             "O45", "O4", "O6", "O78", "O7", "O81", "O8", "O83", "O86", "O9"}


def read_ectyper():
    """Return {strain: o_type} for single-call typed strains (skip '-' and slash)."""
    out = {}
    with open(ECTYPER_TSV) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        name_i = header.index("Name")
        o_i = header.index("O-type")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= o_i:
                continue
            s, o = f[name_i], f[o_i]
            if o == "-" or "/" in o:      # untypeable or ambiguous multi-call -> skip
                continue
            out[s] = o
    return out


def read_scored230():
    """Strains that were held out (scored) in the LOSO-19 run -> off-limits for tail."""
    scored = set()
    with open(LOSO_PRED) as fh:
        header = fh.readline().rstrip("\n").split(",")
        si = header.index(STRAIN_COL)
        for line in fh:
            f = line.rstrip("\n").split(",")
            scored.add(f[si])
    return scored


def has_faa(strain):
    return any(os.path.exists(os.path.join(AA_DIR, strain + ext))
               for ext in (".faa", ".fasta", ".fa"))


def read_usable_universe():
    """All strains with a phenotype (the usable-strain universe) from the LOSO matrix.

    Each batch's modeling set = this universe minus the batch's validation strains.
    Falls back to the AA_DIR basenames (with a warning) if the matrix is absent
    (e.g. building off-cluster) -- less authoritative but keeps the builder usable.
    """
    if os.path.exists(INTERACTION_MATRIX):
        strains = set()
        with open(INTERACTION_MATRIX) as fh:
            header = fh.readline().rstrip("\n").split(",")
            si = header.index(STRAIN_COL)
            for line in fh:
                f = line.rstrip("\n").split(",")
                if len(f) > si:
                    strains.add(f[si])
        return strains, "interaction_matrix"
    # fallback: AA FASTA basenames
    uni = set()
    for fn in os.listdir(AA_DIR):
        for ext in (".faa", ".fasta", ".fa"):
            if fn.endswith(ext):
                uni.add(fn[: -len(ext)])
    return uni, "AA_DIR(fallback)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help="Report only; do not write.")
    args = ap.parse_args()

    otype = read_ectyper()
    scored230 = read_scored230()
    universe, uni_src = read_usable_universe()
    print(f"usable-strain universe: {len(universe)} strains (source: {uni_src})")
    if uni_src.startswith("AA_DIR"):
        print("  WARNING: interaction matrix not found; using AA_DIR basenames as the "
              "training universe. Verify this matches the LOSO strain set on Lawrencium.")

    # serogroup -> members (single-call only)
    members = defaultdict(list)
    for s, o in otype.items():
        members[o].append(s)

    # small tail serogroups (n < MAX_N)
    small = {o: ss for o, ss in members.items() if len(ss) < MAX_N}

    # Admissibility: drop scored-19-labelled (leak vs training), drop collisions
    # (strain already scored in the 19), require an AA FASTA (runnable).
    clean = []          # (strain, serogroup)
    dropped = {'scored19_label': [], 'collision_230': [], 'no_faa': []}
    for o, ss in small.items():
        if o in SCORED_19:
            dropped['scored19_label'] += [(s, o) for s in ss]
            continue
        for s in ss:
            if s in scored230:
                dropped['collision_230'].append((s, o))
            elif not has_faa(s):
                dropped['no_faa'].append((s, o))
            else:
                clean.append((s, o))

    # Whole-serogroup greedy packing into N_BATCHES balanced batches.
    bys = defaultdict(list)
    for s, o in clean:
        bys[o].append(s)
    batches = [[] for _ in range(N_BATCHES)]
    for o in sorted(bys, key=lambda g: -len(bys[g])):       # largest serogroup first
        i = min(range(N_BATCHES), key=lambda j: len(batches[j]))
        batches[i].extend((s, o) for s in bys[o])

    # ---- report ----
    print(f"small single-call serogroups (n<{MAX_N}): {len(small)}")
    print(f"admissible clean tail: {len(clean)} strains, {len(bys)} serogroups")
    for k, v in dropped.items():
        print(f"  dropped [{k}]: {len(v)}" + (f"  e.g. {v[:4]}" if v else ""))
    for bi, b in enumerate(batches, 1):
        print(f"  batch_{bi}: {len(b)} strains, {len(set(o for _, o in b))} serogroups")

    # ---- invariants (fail loud) ----
    all_tail = [s for b in batches for s, _ in b]
    assert len(all_tail) == len(set(all_tail)), "tail strain assigned to >1 batch"
    assert not (set(all_tail) & scored230), "tail overlaps scored-19 strains (leak)"
    assert all(o not in SCORED_19 for b in batches for _, o in b), "scored-19 label in a batch"

    if args.check:
        print("\n--check: no file written. Invariants passed.")
        return

    # Each batch's modeling set = usable universe MINUS that batch's validation
    # strains (== "leave this batch out"; the same rule the 19 folds used).
    val_by_batch = {f'batch_{bi}': set(s for s, _ in b) for bi, b in enumerate(batches, 1)}

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    n_val = n_mod = 0
    with open(OUT_CSV, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['fold_label', STRAIN_COL, 'role', 'serogroup'])
        for bi, b in enumerate(batches, 1):
            label = f'batch_{bi}'
            val = val_by_batch[label]
            # serogroup only meaningful for the tail (validation) rows
            for s, o in b:
                w.writerow([label, s, 'validation', o])
                n_val += 1
            for s in sorted(universe - val):
                w.writerow([label, s, 'modeling', otype.get(s, '')])
                n_mod += 1

    print(f"\nwrote {OUT_CSV}")
    print(f"  {n_val} validation rows + {n_mod} modeling rows across {N_BATCHES} batches")
    print(f"  each batch: ~{len(all_tail)//N_BATCHES} validation, ~{len(universe)-len(all_tail)//N_BATCHES} modeling")
    print("Run with submit_loso_slurm.py pointed at this folds_file (batches are the fold labels).")
    print("Score per-serogroup afterward from final_predictions.csv (group on the serogroup col).")


if __name__ == "__main__":
    main()
