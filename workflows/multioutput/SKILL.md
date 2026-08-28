---
name: genophi-multioutput
description: Run the GenoPHI multi-output workflow — predict several phenotypes (receptors, host range, fitness) per genome, either one model per target or one shared model. Use when asked to cross-validate multi-target genome-to-phenotype models, train deployable multi-output models, prepare a phenotype table against the phage genome manifest, or interpret multi-output CV results. Covers multilabel classification and multi-target regression.
---

# GenoPHI multi-output workflow

You are running or assisting with a three-step pipeline: **prepare → cross-validate
→ train final**. Human documentation is in `README.md` beside this file; read it
if the user asks anything it covers in more depth.

Scripts live in `<repo>/workflows/multioutput/`:
- `prepare_inputs.py` — validation gate and genome extraction
- `submit_run.py` — `--mode cv` and `--mode final`

## Non-negotiables

1. **Never edit anything under `genophi/`** without describing the change and
   getting explicit user approval first. The package is extensively debugged.
   These wrapper scripts are yours to edit freely.
2. **Never report metrics from a `--mode final` run as performance.** They are
   in-sample and will look excellent regardless of whether the model generalises.
   Performance comes only from `--mode cv`.
3. **Never run step 2 or 3 if step 1 did not report `READY`.** `submit_run.py`
   enforces this; do not reach for `--skip_prepare_check` to get past it.
4. **Always `--dry_run` first** and show the user the generated SLURM script
   before a real submission. These runs cost hundreds of core-hours.
5. **Output to cluster scratch, never `$HOME`.**

## Step 0 — verify the environment

Before anything else:

```bash
python -c "import genophi; print(genophi.__file__)"
genophi --help | grep multioutput-cv
```

The path must be inside a checkout on branch `multioutput-support`, and
`multioutput-cv` must appear. If either fails, stop and tell the user — the
multi-output code is not on `main` or `nested-cv`.

**`PYTHONPATH` does not fix this.** Editable installs register an import finder
that outranks `PYTHONPATH`. The only fixes are activating an environment
installed against the right directory, or `pip install -e <worktree>` into a new
environment. Do not suggest `PYTHONPATH` as a workaround; it silently loads the
wrong branch.

## Step 1 — prepare inputs

```bash
python prepare_inputs.py \
    --phenotype_table <table.csv> \
    --sample_column <id column> \
    --host_pattern <short host substring> \
    --output_dir <prepared dir> \
    [--task_type regression] [--allow_missing_genomes]
```

Then **read `prepare_report.json`** and report to the user, in this order:

1. `status` — anything but `"ready"` blocks the next step.
2. `dropped_no_genome` — name every dropped sample explicitly. Silent drops are
   the most common way a run becomes wrong.
3. `nested_targets` — if non-empty, tell the user which targets overlap and that
   per-target metrics are not independent evidence. Do not silently proceed past
   this.
4. Per-target positive counts — flag any target with under 5 in either class as
   likely to be dropped by some folds.

Guidance for the flags:
- `--host_pattern` is a **case-insensitive substring** of the manifest `Host`
  column. Host strings are inconsistent (`Mycobacterium`,
  `Mycolicibacterium smegmatis mc2155`, …), so use a short stem (`yco`, `oli`),
  never a full host name. Confirm the reported row count is plausible.
- Never pre-filter the manifest on `Latest == 'latest'`. That column only
  disambiguates samples with multiple rows; many valid single-row entries leave
  it blank. `prepare_inputs.py` already handles this correctly.
- Use the live `manifest.tsv`. The `manifest.tsv.<timestamp>` siblings are older
  backups.
- Regression targets must be dense. Missing values are an error; unmeasured or
  low-confidence cells become `0`, not NaN.

## Step 2 — cross-validation

```bash
python submit_run.py --prepared_dir <prepared dir> --mode cv \
    --strategy independent --output_dir <scratch>/<name>_cv --dry_run
```

`independent` is the default; use `joint` only when the user asks for it or
wants a benchmark. Job count is `n_folds × n_targets` (+1) for independent and
`n_folds` (+1) for joint — state the number before submitting.

Pass `--min_features 3` (already the default here) for datasets under ~150
samples. GenoPHI's own default of 5 causes folds to produce **no model and no
error**.

For phylogenetic folds add `--group_metadata <csv> --group_column group`.

### Interpreting results

Read `<output_dir>/cv_performance/model_performance_metrics.csv`.

- Check `n_samples_<target>` **first**. Below the full sample count means folds
  dropped that target; its metrics rest on a subset and you must say so.
- Judge by **MCC and AUPR**, never accuracy — an all-negative prediction on a
  10%-positive target scores 90% accuracy.
- AUC ≈ 0.5 with MCC ≤ 0 is chance. Say so plainly; do not soften it.
- Do not compare aggregate metrics across strategies — `n_complete_samples`
  differs. Compare per-target on equal coverage, or use `genophi benchmark`.

Sample size dominates: the same pipeline reached AUC > 0.97 at n=255 and chance
at n=69, while fitting the n=69 training data perfectly. If a result is at
chance and n is small, the finding is the sample size — report it as a result,
not a failure to be tuned away.

### Regression results

For `--task_type regression` on fitness-style targets, **do not headline R² or
Spearman.** Zero-inflated heavy-tailed targets break both: Pearson rides the
extremes and the unrankable near-zero bulk drags Spearman below it. Report
rank-AUROC, AUPR and precision@K — strong-responder detection. Never propose a
signed-log or similar transform; it destroys exactly the outliers that carry the
signal. A target with one strong responder is unmeasurable, not a failure.

## Step 3 — final model

Only after CV shows real signal. If CV was at chance, say so and recommend
against deploying; proceed only if the user still wants the model.

```bash
python submit_run.py --prepared_dir <prepared dir> --mode final \
    --strategy independent --output_dir <scratch>/<name>_final
```

Models land in `<output_dir>/<target>/modeling_results/cutoff_<n>/run_*/`.

**Prediction on new genomes is not wrapped.** `submit_run.py` has no `predict`
mode. Use `genophi kmer-assign-predict`, which additionally requires
`--mmseqs_db`, `--clusters_tsv`, `--feature_map`, `--filtered_kmers` and
`--aa_sequence_file` from the training run — new genomes must be mapped into the
same feature space, so you cannot build a fresh feature table and predict against
it. The reference implementation is `_predict_heldout_independent` in
`genophi/workflows/multioutput_cv_parallel.py`. If a user asks for predictions,
say this step needs the artifact paths resolved rather than guessing a command.

Clustering settings **must match** between step 2 and step 3, or the deployed
model is not the one you measured. `submit_run.py` keeps them consistent by
passing `--no-clustering` in both modes; if the user enables `--use_clustering`,
enable it in both.

## Failure signatures

| Signature | Meaning | Action |
|---|---|---|
| `invalid choice: 'multioutput-cv'` | wrong branch or env | re-check step 0 |
| `No feature tables to model` | `min_features` too high | `--min_features 3` |
| `All train targets are equal` | target single-class in a fold | expected; that target drops from that fold |
| `FileNotFoundError: filtered_feature_tables` | 0 features for a (fold, target) | expected on small-n; the cell is skipped |
| `Missing required features` at predict | model predates the sibling-strip fix | retrain; do not patch the feature table |
| Stage C never starts | a Stage B task failed, `afterok` unsatisfied | read `logs/train_*.err`, resubmit to the **same** `--output_dir` |
| Job hit walltime | under-resourced | resubmit to the same `--output_dir`; completed cells skip |

Runs are resumable: the same `--output_dir` skips completed fold/target cells.
Resubmitting is almost always right; starting a fresh directory throws away work.

## Reporting

When you report a completed run, include: which strategy, fold scheme and
`n_folds`, per-target AUC/AUPR/MCC **with `n_samples`**, any target dropped from
any fold, and any sample dropped at step 1. If aggregate numbers are quoted,
state `n_complete_samples`. If the result is chance, say it is chance.
