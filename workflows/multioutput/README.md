# GenoPHI multi-output workflow

Predict **several phenotypes at once** from phage (or bacterial) genomes — one
model per target, or one shared model over all targets.

This directory is the supported entry point. It wraps the multi-output code in
`genophi/workflows/` so you do not have to reconstruct the command lines, and it
sets the defaults that previous runs had to discover the hard way.

**Typical use:** you have a table of phages and, for each phage, several binary
labels (which receptor it uses, which host it kills). You want (a) an honest
estimate of how well those labels can be predicted from genome sequence, and
(b) a trained model to apply to new phages.

---

## Contents

| File | Purpose |
|---|---|
| `reshape_matrix.py` | Step 0 — turn a raw assay matrix into the samples-as-rows layout (skip if already in that shape) |
| `prepare_inputs.py` | Step 1 — validate the phenotype table, resolve genomes via the manifest, extract protein FASTAs |
| `submit_run.py` | Steps 2 and 3 — submit cross-validation (`--mode cv`) or final model training (`--mode final`) |
| `SKILL.md` | The same workflow written for an AI agent |
| `README.md` | This file |

---

## The four steps

```
    raw assay matrix (genes x samples, controls, gaps)
              │
        [0]   reshape_matrix.py     (skip if already samples-as-rows)
              │
    phenotype table (.csv)                genome manifest
              \                                 /
               \                               /
          [1]   prepare_inputs.py  ────────────
                        │
                        ├── phenotype.csv   (validated, filtered)
                        ├── genomes/*.faa   (one per sample)
                        └── targets.txt
                        │
          ┌─────────────┴─────────────┐
          │                           │
    [2] submit_run.py           [3] submit_run.py
        --mode cv                   --mode final
          │                           │
    per-target metrics           trained models
    (what you report)            (what you deploy)
```

**Step 2 and step 3 answer different questions and are both needed.** CV tells
you whether the labels are predictable at all; it produces no model you can
ship, because every model it builds has seen only 9/10 of the data. Final
training produces the deployable model but no honest metric, because it has no
held-out set. Run CV first — if CV says the signal is chance, a final model is
a memorisation of your table and should not be deployed.

---

## Before you start

You need three things.

**1. A GenoPHI checkout on the `multioutput-support` branch.** The multi-output
code does not exist on `main` or `nested-cv`.

```bash
git worktree add ~/software/GenoPHI-multioutput multioutput-support
```

Use a **worktree**, not `git checkout` in an existing clone. `pip install -e .`
binds an environment to a *directory path*, so switching that directory's branch
changes the code every running job imports. A worktree gives the branch its own
directory.

**2. A conda environment installed against that directory.**

```bash
conda create -n genophi-mo --clone genophi
conda activate genophi-mo
pip install -e ~/software/GenoPHI-multioutput
```

> **`PYTHONPATH` will not work as a shortcut.** Editable installs register an
> import finder that takes precedence over `PYTHONPATH`, so
> `PYTHONPATH=/path/to/worktree python -c 'import genophi'` still loads the
> *other* checkout. Verify with
> `python -c "import genophi; print(genophi.__file__)"` — it must print a path
> inside your worktree. If it does not, nothing below will behave as documented.

**3. Somewhere to run it.** You have two options, and the workflow supports
both through `--executor`:

| | `--executor local` | `--executor slurm` |
|---|---|---|
| Runs | here, in your shell, one fold at a time | as cluster jobs, folds in parallel |
| Needs | nothing extra | a SLURM allocation |
| Good for | small datasets, smoke tests, a first look | anything real |

Cross-validation is `n_folds × n_targets` independent model fits. A 10-fold,
2-target run is 20 fits; a 10-fold, 5-target run on 69 genomes previously took
about **17 hours** serially. Rough guidance: under ~50 samples and ~3 targets,
local is fine; above that, use a cluster.

Both executors call the **same functions** — local is not a reduced pipeline, it
is the same computation without the parallelism, so results are directly
comparable.

## Step 0 — reshape a raw assay matrix

Skip this if your table already has **one row per genome and one column per
target**. Assay exports usually do not.

RB-TnSeq and CRISPRi tables typically arrive as **genes in rows, samples in
columns** — the transpose of what the modeling code wants — plus non-sample
columns (media controls, blanks), free-text target names, and gaps.

```bash
python reshape_matrix.py \
    --input        20260710_PA14_s_pos.csv \
    --samples_in   columns \
    --drop_samples LB \
    --output       pa14_fitness.csv
```

`--samples_in` is required and has no default. It is the one thing the script
cannot safely infer, and getting it wrong silently produces a matrix of genomes
that do not exist.

| Flag | What it handles |
|---|---|
| `--samples_in {rows,columns}` | orientation; `columns` transposes |
| `--drop_samples LB,control` | non-genome columns. The script also warns about likely controls it recognises (`LB`, `blank`, `mock`, …) |
| `--fillna 0` | missing cells. Regression needs a dense matrix; `0` means "no confident effect", the right reading for low-coverage cells. `--fillna error` refuses instead |
| `--top_n_targets N --top_by max` | keep the N most responsive targets. Use `max`, not `absmax` — `absmax` surfaces single extreme *negative* outliers rather than the responders you care about |
| target name sanitisation | on by default. Independent modeling creates **one directory per target**, so `MSMEG5483 \| \| MspA family porin` would become a path. Originals are preserved in `<output>_target_names.csv` |

The script reports the reshaped dimensions, the value range, and how
zero-inflated the matrix is — and if it is, it tells you which metrics to use.

## Step 1 — prepare and validate inputs

```bash
python prepare_inputs.py \
    --phenotype_table  Binary_table_only_Msp.csv \
    --sample_column    Phages \
    --host_pattern     yco \
    --output_dir       prepared/msp_binary \
    --allow_missing_genomes
```

### What the phenotype table must look like

One row per sample, one column of sample ids, one column per target:

```csv
Phages,Msp_TRUE_dom,Msp_FALSE_dom
Bxb1,1,0
EagleEye,1,0
DS6A,0,0
```

- Sample ids must match the `Phage` column of the genome manifest **exactly**
  (case-sensitive).
- Classification targets must be `0`/`1` with no blanks.
- Regression targets must have **no missing values** — CatBoost's `MultiRMSE`
  needs a dense matrix. Fill unmeasured cells with `0` (see *Regression* below).

### What this step checks

| Check | Behaviour |
|---|---|
| Sample column exists, ids unique | error if not |
| Target columns exist and are non-constant | error if constant |
| Missing values in targets | error |
| Class balance per target | warns under 5 in either class |
| **Targets identical or nested** | warns (see below) |
| Every id resolves to a genome in the manifest | error, unless `--allow_missing_genomes` |
| Every manifest path exists on disk | error |
| **FASTA basenames == table ids, both directions** | error if not |

The last one is the point of the whole step. Everything downstream silently
produces garbage — or an all-zero feature vector — when a sample's genome is
missing under the name the table uses.

### Genome resolution

Genomes come from the manifest, by default
`/usr2/people/phages/phage_genomes/manifest.tsv`. Always use the live
`manifest.tsv`; the `manifest.tsv.<timestamp>` files beside it are older backups.

Two traps this step handles for you:

- **Do not filter on `Latest == 'latest'`.** That column only disambiguates
  samples with *several* manifest rows. Many single-row entries leave it blank —
  in the Mycobacterium set, 21 of 93 phages have a blank `Latest` and perfectly
  good genomes. Filtering on it up front would silently drop them.
- **`Host` strings are inconsistent.** The Mycobacterium phages are spread over
  seven host spellings (`Mycobacterium`, `Mycobacterium smegmatis mc#155`,
  `Mycolicibacterium smegmatis mc2155`, …). `--host_pattern` is a
  case-insensitive *substring*: use a short stem like `yco` or `oli`, not a full
  name. Omit it to search all hosts.

Protein FASTAs are written by extracting CDS `/translation` records from each
GenBank file, one `<sample>.faa` per sample. Re-running skips files that already
exist; `--force` rebuilds them.

### Outputs

```
prepared/msp_binary/
├── phenotype.csv          filtered to samples that have a genome; id column renamed to `phage`
├── genomes/<id>.faa       one protein FASTA per sample
├── targets.txt            one target name per line
├── prepare_report.txt     the full human-readable log
└── prepare_report.json    machine-readable; step 2 refuses to run unless status == "ready"
```

---

## Step 2 — cross-validation (the numbers you report)

```bash
python submit_run.py \
    --prepared_dir prepared/msp_binary \
    --mode cv \
    --strategy independent \
    --output_dir /global/scratch/users/$USER/msp_cv \
    --n_folds 10
```

Add `--dry_run` to write the SLURM scripts without submitting — always do this
once and read the generated script before committing a few hundred core-hours.

To run it here instead, with no cluster and no allocation:

```bash
python submit_run.py \
    --prepared_dir prepared/msp_binary \
    --mode cv --executor local \
    --output_dir runs/msp_cv
```

Local runs print per-fold progress and elapsed time. **Interrupting is safe** —
completed folds are skipped when you re-run into the same `--output_dir`.

Per fold, the pipeline rebuilds the k-mer feature table **from the training
genomes only**, runs feature selection and ensemble modelling, then assigns and
predicts the held-out genomes. Held-out samples never touch feature selection,
so there is no leakage.

### `independent` vs `joint`

`independent` is the default and the right choice unless you have a reason.

| | `independent` | `joint` |
|---|---|---|
| Models | one per target | one shared multi-output model |
| Loss | standard | `MultiLogloss` / `MultiRMSE` |
| Feature set | each target picks its own | one set shared by all targets |
| Cutoff | each target picks its own | one shared cutoff |
| SLURM jobs | `n_folds × n_targets` + 1 | `n_folds` + 1 |
| Sparse targets | can fail individually | borrow strength from co-targets |

On the 255-phage E. coli receptor set with proper CV the two are **statistically
a wash** (mean MCC difference +0.004 across 14 receptors). The difference is
per-target and depends on whether a target's signal correlates with its
co-targets, not simply on how sparse it is — of the two largest differences,
both on sparse receptors, one favoured joint and the other independent.

Use `joint` when you specifically want the robustness of one model that always
produces a value for every target, or when you are benchmarking the two.

> Aggregate metrics are **not comparable between the two strategies**. They are
> computed over different sets of complete samples (`n_complete_samples` in the
> output). Compare per-target metrics on equal coverage, or use
> `genophi benchmark`.

### Fold structure

Random k-fold is the default. For a harder, more honest test, hold out whole
phylogenetic clades:

```bash
--group_metadata clades.csv --group_strain_column strain --group_column group
```

with a CSV of `strain,group`. Folds become leave-one-group-out and the fold count
is the number of groups. Joint and independent use identical splits, so
`diff` the two runs' `cv_splits.csv` to prove it.

### Running on a cluster

**`--account` is required and has no default.** SLURM allocations are
per-project, and you can only charge one you belong to — the account in someone
else's example will be rejected. To find yours:

```bash
sacctmgr show associations user=$USER format=Account,Partition,QOS
# on Lawrencium:
lrc-associations -u $USER
```

Then pass it, along with a partition and QOS your allocation permits:

```bash
--account <your_project> --partition lr7 --qos lr_normal
```

Running without `--account` prints these instructions rather than guessing.

**Your inputs must be where the compute is.** The prepared directory is
self-contained (a phenotype CSV, a `genomes/` directory, `targets.txt`), so
copying it across is all that is needed:

```bash
PREP=<local>/prepared/<dataset>
REMOTE=<user>@<transfer-host>:/global/scratch/users/<user>/<project>/prepared/

ssh <user>@<transfer-host> "mkdir -p /global/scratch/users/<user>/<project>/prepared"
rsync -avP "$PREP" "$REMOTE"          # no trailing slash: creates prepared/<dataset>/
```

On Lawrencium the transfer host is `lrc-xfer.lbl.gov`; submit from the login
node. Then point `--prepared_dir` at the **remote** copy — the paths are baked
into the generated SLURM scripts, so a local path will not resolve on a compute
node.

Bring the results back the same way; the summaries are small:

```bash
rsync -avP <user>@<transfer-host>:<output_dir>/cv_performance/ ./cv_performance/
rsync -avP <user>@<transfer-host>:<output_dir>/cv_predictions.csv ./
```

The per-fold directories are large and rarely needed off-cluster — pull them
only if you intend to re-analyse folds individually.

### Reading the results

```
<output_dir>/
├── cv_splits.csv                             which sample was held out in which fold
├── cv_predictions.csv                        per-sample predictions + true labels + fold
├── cv_performance/model_performance_metrics.csv   ← the headline table
└── fold_1/ … fold_N/                         per-fold feature tables and models
```

`model_performance_metrics.csv` has one row (`cv_pooled`) with per-target
`AUC_<target>`, `AUPR_<target>`, `MCC_<target>`, `F1_<target>`,
`n_samples_<target>`, plus aggregate `macro_f1`, `micro_f1`, `hamming_loss`,
`subset_accuracy`, `n_complete_samples`.

**Read `n_samples_<target>` before you read anything else.** If it is below your
sample count, some folds dropped that target and its metrics rest on a subset.

**Judge by MCC and AUPR, not accuracy.** With a target that is 10% positive,
predicting all-negative scores 90% accuracy and is worthless; MCC ≈ 0 and AUPR
near the base rate reveal it. AUC ≈ 0.5 and MCC ≤ 0 mean chance.

**Sample size is the dominant factor.** The same pipeline gave AUC > 0.97 on 255
E. coli phages and chance-level performance (AUC 0.46–0.65, mostly negative MCC)
on 69 Mycobacterium phages — while fitting the Mycobacterium training data
perfectly (in-sample macro-F1 = 1.0). Below roughly 100 samples, expect
memorisation, and never quote an in-sample number.

---

## Step 3 — final model training (what you deploy)

Only after CV says the targets are predictable:

```bash
python submit_run.py \
    --prepared_dir prepared/msp_binary \
    --mode final \
    --strategy independent \
    --output_dir /global/scratch/users/$USER/msp_final
```

Add `--executor local` to run it here instead. Final training is a single fit
per target, so it is far cheaper than CV — often practical locally even when CV
was not.

This trains on **all** samples and writes one model directory per target:
`<output_dir>/<target>/modeling_results/cutoff_<n>/run_*/`, plus
`independent_summary.csv`.

> **The metrics in this run are in-sample. Never report them.** They will look
> excellent even when CV showed chance performance. The numbers you report come
> from step 2; the model you apply comes from step 3.

> **Clustering defaults differ between the CV path and `kmer-workflow`.** The CV
> workflows default `use_clustering=False`; the `genophi kmer-workflow` CLI
> defaults it to *True*. Calling `kmer-workflow` by hand with its own defaults
> therefore trains a model that does not match the CV you reported.
> `submit_run.py` passes `--no-clustering` in both modes so they agree. If you
> enable `--use_clustering`, enable it in **both** steps.

### Applying the model to new genomes

**This step is not yet wrapped by `submit_run.py`.** It works, but you must call
GenoPHI directly, and the command needs several artifacts from the training run:

```bash
genophi kmer-assign-predict \
    --input_dir        new_genomes/ \
    --model_dir        <output_dir>/<target>/modeling_results/cutoff_<n> \
    --mmseqs_db        <from the training run> \
    --clusters_tsv     <from the training run> \
    --feature_map      <from the training run> \
    --filtered_kmers   <from the training run> \
    --aa_sequence_file <the training set's combined .faa> \
    --genome_type      phage \
    --tmp_dir          <scratch tmp> \
    --output_dir       predictions/
```

New genomes must be mapped into the *same* feature space the model was trained
in — that is what the four artifact arguments are for. You cannot rebuild the
feature table independently and predict against it; the k-mer columns would not
line up.

The working reference implementation is `_predict_heldout_independent` in
`genophi/workflows/multioutput_cv_parallel.py`, which is exactly what CV uses to
score held-out genomes. For a **joint** model, `genophi multioutput-predict
--model_dir <cutoff dir> --feature_table <table> --output_dir <out>` takes an
already-assigned feature table and medians probabilities across the ensemble.

New genomes still need protein FASTAs — reuse `prepare_inputs.py` with a
phenotype table of dummy labels to get the extraction and name checking.

> Wrapping this as `submit_run.py --mode predict` is the obvious next addition;
> it needs the artifact paths inside a `kmer-workflow` output tree pinned down
> first.

## Regression (multi-target)

Set `--task_type regression` in step 1 and step 2. The joint path uses
`MultiRMSE`; targets must be dense (no NaN).

Fitness-style targets are **zero-inflated and heavy-tailed** — most values near
zero, a few strong responders — and that breaks the obvious metrics:

- **Do not headline R² or Spearman.** On the Mycobacterium fitness set, mean
  R² = 0.27 looked like failure but was a metric mismatch. Pearson is carried by
  the extremes; the unrankable near-zero bulk drags Spearman *below* Pearson.
- **Report strong-responder detection instead** — rank-AUROC, AUPR, precision@K.
  The same models scored mean rank-AUROC 0.844 and AUPR 3–7× the strong-fraction
  baseline. The outliers are the signal.
- **Do not transform the targets.** A signed-log compresses exactly the extremes
  you are trying to detect.
- Targets with a single strong responder are unmeasurable, not failures.
- When building the matrix, low-confidence cells (e.g. TnSeq coverage below
  threshold) become **0**, not NaN.

`--strong_top_frac` (default 0.2) sets what fraction counts as "strong" for the
detection metrics.

**Use `joint` for regression, not `independent`.** Fitness targets are usually
many and highly correlated, and both facts point the same way: joint shares
strength across correlated targets, and it is one model instead of
`n_folds × n_targets`. With 73 targets, independent 10-fold CV is **730** jobs;
joint is 10.

`--target_mode` is set to `multitarget` automatically when `--task_type
regression` is passed, so `MultiRMSE` is used rather than a classification loss.

### Worked example: Pseudomonas RB-TnSeq fitness

37 phages profiled against 73 PA14 genes — pilus, flagellar and LPS
biosynthesis loci. The raw export is genes-as-rows with an `LB` media control
column.

```bash
python reshape_matrix.py --input 20260710_PA14_s_pos.csv \
    --samples_in columns --drop_samples LB --output pa14_fitness.csv

python prepare_inputs.py --phenotype_table pa14_fitness.csv \
    --sample_column phage --task_type regression \
    --host_pattern seudomon --output_dir prepared/pa14_fitness \
    --allow_missing_genomes

python submit_run.py --prepared_dir prepared/pa14_fitness \
    --mode cv --strategy joint --task_type regression \
    --output_dir <scratch>/pa14_cv --account <your_project>
```

What the preparation steps report, and why each matters:

- **37 → 36 samples.** `JBD30` is not in the manifest under any spelling; there
  is no near-match. It is dropped and named.
- **7 phages had multiple manifest rows**, resolved by `Latest`. This is the case
  the `Latest` column exists for — unlike the Mycobacterium set, where every
  phage has exactly one row and blank `Latest` values must be ignored.
- **73 target names sanitised**, e.g. `IKLFDK_00745_pilus assembly protein PilE`
  → `IKLFDK_00745_pilus_assembly_protein_PilE`, with the originals kept in the
  name map.
- **59% of cells within ±1 of zero, 25% above |4|** — textbook zero-inflated.
  Report detection metrics, not R².

**Before running it, know what the data can support.** Two principal components
explain **92%** of the variance across all 73 targets (one explains 63%), and
31% of target pairs correlate above |r| = 0.8. So this is not 73 independent
phenotypes — it is roughly **two or three**, almost certainly
pilus-dependent versus flagellum/LPS-dependent entry. That is an argument for
joint modeling, and a caution against reading 73 per-target metrics as 73
independent results.

At n = 36 this is well below the n = 69 that produced chance-level performance
on the Mycobacterium set. The favourable factors are the low effective
dimensionality and a median of 8 strong responders per target (only 4 targets
have just one, and are therefore unmeasurable). Run it, but treat a weak result
as the expected outcome of the sample size rather than something to tune away.

---

## Parameters worth knowing

| Flag | Default here | Notes |
|---|---|---|
| `--min_features` | `3` | **Independent strategy only** — the joint fold runner has no such knob. GenoPHI's own default is 5, which **silently** drops folds on small datasets: feature selection yields 1–4 stable features, every cutoff table is rejected, no model is trained, and no error is raised. 3 is safe under ~150 samples. |
| `--k` | `5` | Protein k-mer size. 5 is what every validated run used. |
| `--num_runs_fs` | `25` | Feature-selection repeats. Lower only for smoke tests. |
| `--num_runs_modeling` | `50` | Ensemble size. Lower only for smoke tests. |
| `--max_features` | `100` | Cap on selected features. |
| `--n_folds` | `10` | |
| `--use_clustering` | off | Sample clustering. Off matches the published runs. |
| `--use_feature_clustering` | off | Collapses phylogenetically linked features. Untested as a rescue for small-n. |
| `--mem`, `--time_train` | 120 GB, 24 h | Per job. Feature-table jobs are the memory-hungry ones. |

**Runs are resumable.** Resubmit into the *same* `--output_dir` and completed
fold/target cells are skipped. This is the correct response to a walltime kill.

---

## Traps

Every one of these has cost someone a run.

1. **Wrong branch.** `multioutput-cv` does not exist on `main`/`nested-cv`. Check
   `genophi --help`.
2. **`PYTHONPATH` does not override an editable install.** Verify
   `genophi.__file__`.
3. **`min_features` fails silently.** "No feature tables to model" and no model,
   no error. Pass `--min_features 3` on small data.
4. **The two CV submitters take different flags.** The joint one accepts
   `--target_mode` and not `--min_features`; the independent one is the reverse.
   `submit_run.py` handles this — if you call the submitters directly, don't get
   it backwards.
5. **`kmer-workflow` defaults clustering ON, CV defaults it OFF.** See step 3.
6. **`Latest == 'latest'` drops valid genomes.** See step 1.
7. **Host strings vary.** Use a substring, not a full host name.
8. **Aggregate joint-vs-independent metrics are not comparable.** Different
   `n_complete_samples`.
9. **In-sample metrics look perfect and mean nothing.** A dataset that CV'd at
   chance still fits its training data at macro-F1 1.0.
10. **Nested or duplicate targets.** If one target's positives are a subset of
    another's, per-target metrics are not independent evidence and a
    joint-vs-independent comparison measures something other than what you
    think. Step 1 warns; the interpretation is yours.
11. **Output to scratch, never `$HOME`.**
12. **`--account` is not portable.** Allocations are per-project; an account
    copied from someone else's command will be rejected. There is deliberately
    no default.
13. **`--prepared_dir` must be a path that exists on the compute node.** It is
    written into the SLURM scripts verbatim. Transfer the prepared directory
    first and point at the remote copy.
14. **`genophi multioutput-cv` cannot express every run.** The CLI subcommand
    accepts neither `--min_features` nor `--group_metadata`, so it cannot do a
    small-n or a phylogenetic-fold run. `submit_run.py --executor local` calls
    the underlying functions directly and does not have that limitation.
15. **A rare target can be all-one-value inside a fold**, which crashes
    `MultiLogloss`. Handled by `drop_single_class_columns`, but it is why a
    target can vanish from some folds and why `n_samples_<target>` shrinks.

---

## Worked example: E. coli receptor prediction (published)

The reference run: **255 phages, 14 receptor targets**, 10-fold CV, both
strategies. This is the analysis behind the manuscript, and the numbers below are
what a run that worked looks like.

Inputs (paths as used on LRC, where the run executed):

| | |
|---|---|
| Phage proteomes | `LRC/TnSeq_phage_AAs/` — 255 per-phage `.faa` |
| Label matrix | `LRC/tables/255_phages/phage_target_summary_pivot_updated.csv` — 255 × 14 binary |
| Clade folds (optional) | `cv_plots/phylo/groups/peq_groups_thr075_min10.csv`, `…thr09_min10.csv` |

with `LRC/` = `/global/scratch/users/anoonan/EDGE/target_modeling/`.

```bash
# 1. prepare
python prepare_inputs.py \
    --phenotype_table  $LRC/tables/255_phages/phage_target_summary_pivot_updated.csv \
    --sample_column    phage \
    --host_pattern     oli \
    --output_dir       prepared/ecoli_receptors

# 2. cross-validate, both strategies into separate directories
python submit_run.py --prepared_dir prepared/ecoli_receptors \
    --mode cv --strategy independent --output_dir $LRC/cv_independent
python submit_run.py --prepared_dir prepared/ecoli_receptors \
    --mode cv --strategy joint       --output_dir $LRC/cv_joint_multilabel

# 3. compare them on equal coverage
genophi benchmark --joint_out $LRC/cv_joint_multilabel \
                  --independent_out $LRC/cv_independent \
                  --output benchmark_receptors.csv

# 4. final deployable models
python submit_run.py --prepared_dir prepared/ecoli_receptors \
    --mode final --strategy independent --output_dir $LRC/final_models
```

### What the results looked like

Per-receptor, pooled over the 10 held-out folds (`joint` strategy), alongside the
independent run's MCC:

| Receptor | AUC | AUPR | MCC (joint) | MCC (indep) | Δ (J−I) |
|---|---|---|---|---|---|
| tsx | 1.000 | 1.000 | 1.000 | 0.967 | +0.033 |
| NGR | 1.000 | 1.000 | 1.000 | 1.000 | 0.000 |
| ompC | 1.000 | 1.000 | 1.000 | 0.970 | +0.030 |
| Kdo | 0.996 | 0.967 | 0.942 | 0.942 | 0.000 |
| HepI | 0.978 | 0.949 | 0.923 | 0.941 | −0.018 |
| btuB | 0.980 | 0.949 | 0.911 | 0.927 | −0.015 |
| fhuA | 0.985 | 0.932 | 0.870 | 0.870 | 0.000 |
| ompF | 0.998 | 0.938 | 0.853 | 0.853 | 0.000 |
| GluI | 0.914 | 0.876 | 0.812 | 0.735 | +0.077 |
| ompA | 0.993 | 0.866 | 0.708 | 0.768 | −0.060 |
| yncD | 0.995 | 0.750 | 0.704 | 0.704 | 0.000 |
| lptD | 0.861 | 0.689 | 0.728 | 0.659 | +0.069 |
| **HepII** | 0.807 | 0.444 | 0.444 | 0.247 | **+0.197** |
| **lamB** | 0.783 | 0.229 | 0.279 | 0.538 | **−0.258** |

Aggregate: joint micro-F1 0.903 / macro-F1 0.791 / subset-accuracy 0.839 /
Hamming 0.015; independent 0.892 / 0.798 / 0.824 / 0.017. Every receptor has
`n_samples = 255` in both runs and `n_complete_samples = 255`, so these two
aggregates *are* comparable — which is not usually the case.

Three things to take from this table:

1. **Most receptors are highly predictable** from k-mer features alone — nine of
   fourteen above MCC 0.85.
2. **Joint and independent are a wash**: mean MCC difference **+0.004** across
   the fourteen. Neither strategy dominates.
3. **The two largest differences are both on sparse receptors and point in
   opposite directions** — HepII (5 positives) is rescued by joint modelling,
   lamB (6 positives) is diluted by it. So the joint/independent choice is
   *receptor-specific* and depends on whether a target's signal correlates with
   its co-targets, not simply on how rare it is. Do not assume joint rescues
   sparse targets.

Compare against a run that did **not** work — the same pipeline on 69
Mycobacterium phages returned AUC 0.46–0.65 with mostly negative MCC, while
fitting its training data at macro-F1 1.0. Same code, same parameters; the
difference is n = 255 versus n = 69.

### Second example: Mycobacterium Msp dominance

A smaller, live dataset — 93 phages, 2 binary targets:

```bash
DATA=/usr2/people/anoonan/BRaVE/resources/genome_data/mycobacterium/CRISPRi/data

python prepare_inputs.py \
    --phenotype_table $DATA/Binary_table_only_Msp.csv \
    --sample_column Phages --host_pattern yco \
    --output_dir $DATA/prepared/msp_binary \
    --allow_missing_genomes

python submit_run.py --prepared_dir $DATA/prepared/msp_binary \
    --mode cv --strategy independent \
    --output_dir /global/scratch/users/$USER/msp_cv
```

Step 1 reports what you want to know before spending anything: 93 ids, 92
resolved, `L5` absent from the manifest entirely, `FM` in the manifest but
unphenotyped, and a warning that `Msp_FALSE_dom`'s 25 positives are a strict
subset of `Msp_TRUE_dom`'s 38 — so those two targets are not independent
evidence of anything.

At n = 92 this sits in the range where the E. coli/Mycobacterium contrast says
to expect weak generalisation. Treat whatever CV returns as the finding.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `invalid choice: 'multioutput-cv'` | wrong branch or wrong env | check `genophi.__file__` |
| `No feature tables to model` | `min_features` too high | `--min_features 3` |
| `All train targets are equal` | a target is single-class in a fold | expected on rare targets; that target drops from that fold |
| `Missing required features` at predict time | model trained with sibling targets in X | old run predating the sibling-strip fix — retrain |
| `FileNotFoundError: filtered_feature_tables` | a (fold, target) yielded 0 features | expected on small-n; that cell is skipped |
| Stage C never runs | a Stage B task failed, `afterok` unsatisfied | check `logs/train_*.err`, resubmit into the same `--output_dir` |
| Metrics look perfect | you are reading an in-sample or final-training run | use `cv_performance/model_performance_metrics.csv` |
| `n_samples_<target>` < sample count | folds dropped that target | metrics rest on a subset; say so |
| Job killed at walltime | too few resources | resubmit to the same `--output_dir`; finished cells skip |

---

## Provenance

The defaults and warnings here come from runs on three datasets: E. coli
receptors (255 phages, 14 targets — random and phylogenetic CV, joint and
independent), Mycobacterium CRISPRi (69 phages, 5 receptors plus fitness
regression), and O157/G4C (38 phages). The sample-size finding, the
joint-vs-independent wash, and the zero-inflated regression metrics all come
from those runs.
