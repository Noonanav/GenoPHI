"""
Nested k-fold cross-validation workflow for GenoPHI.

Runs deterministic (repeated) k-fold cross-validation entirely in-process:
for each fold, the strains are split into a modeling set and a held-out
validation set, a model is trained on the modeling set via
``run_protein_family_workflow``, the best cutoff (highest MCC) is selected, and
interactions for the held-out strains are predicted via
``assign_predict_workflow``. Per-fold predictions are aggregated into a single
predictions table plus summary statistics.

This is the portable, scheduler-agnostic equivalent of the SLURM job-array
submission script: folds run sequentially in this process rather than as
parallel SLURM array tasks.

Two split modes are supported: deterministic (repeated) k-fold, and
leave-one-group-out (e.g. leave-one-serotype-out) driven by a strain->group
metadata CSV.

Optionally (shared_clustering=True), MMseqs2 clustering is run ONCE over all
genomes and reused across folds. Clustering is sequence-only and
label-independent, so this is leakage-free as long as the hash-based feature
collapse stays per-fold on the training-only (filtered) presence/absence
matrix -- which it does. The functions run_shared_clustering() and
run_fold_from_shared() are exposed so a SLURM submitter can run one shared
clustering job and fan out one fold job per fold as dependencies.
"""

import os
import random
import logging

import pandas as pd

from genophi.workflows.protein_family_workflow import run_protein_family_workflow
from genophi.workflows.assign_predict_workflow import assign_predict_workflow
from genophi.workflows.assign_features_workflow import run_assign_features_workflow
from genophi.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table
from genophi.mmseqs2_clustering import (
    run_clustering_workflow,
    run_feature_assignment,
    merge_feature_tables,
)

logger = logging.getLogger(__name__)

# Column names written by the prediction workflow into
# strain_median_predictions.csv (see prediction_workflow.calculate_median_predictions).
_SCORE_COL = 'Confidence'        # median predicted probability of the positive class
_PRED_COL = 'Final_Prediction'   # thresholded 0/1 call


def _get_full_strain_list(interaction_matrix, input_strain_dir, strain_column, suffix='faa'):
    """Return strains present in BOTH the interaction matrix and the strain directory.

    Mirrors the intersection logic used by the SLURM bootstrap workflow so that
    only strains we have both a phenotype and a FASTA file for are used.
    """
    logger.info(f"Reading interaction matrix: {interaction_matrix}")
    interaction_df = pd.read_csv(interaction_matrix)

    if strain_column not in interaction_df.columns:
        raise ValueError(
            f"Column '{strain_column}' not found in interaction matrix. "
            f"Available columns: {list(interaction_df.columns)}"
        )

    strains_in_matrix = {str(s) for s in interaction_df[strain_column].unique()}
    logger.info(f"Found {len(strains_in_matrix)} unique strains in interaction matrix")

    dot_suffix = f".{suffix}"
    strain_files = [f for f in os.listdir(input_strain_dir) if f.endswith(dot_suffix)]
    strains_in_dir = {f[: -len(dot_suffix)] for f in strain_files}
    logger.info(f"Found {len(strains_in_dir)} strain files in {input_strain_dir}")

    full_strain_list = sorted(strains_in_matrix & strains_in_dir)
    logger.info(
        f"Intersection: {len(full_strain_list)} strains found in both matrix and directory"
    )
    if not full_strain_list:
        logger.error(
            "No strains found in both the interaction matrix and the input directory. "
            "This is often caused by naming mismatches between the files and the matrix."
        )
    return full_strain_list


def _write_modeling_matrix(interaction_matrix, modeling_strains, output_path,
                           strain_column='strain', phage_column='phage'):
    """Write an interaction matrix filtered to the modeling strains.

    Subsetting to the modeling strains' rows means any phage that only appeared
    with held-out strains drops out entirely. Passing this filtered matrix to
    training ensures phage features (and phage clustering) are built only from
    phages present in the remaining interactions -- important for group CV where
    a held-out group (e.g. a dataset) may take the only strains a phage was
    tested against.

    Returns (output_path, n_rows, n_phages, dropped_phages).
    """
    df = pd.read_csv(interaction_matrix)
    if strain_column not in df.columns:
        raise ValueError(
            f"Column '{strain_column}' not found in interaction matrix "
            f"({list(df.columns)})."
        )
    modeling_set = {str(s) for s in modeling_strains}
    df[strain_column] = df[strain_column].astype(str)
    filtered = df[df[strain_column].isin(modeling_set)].copy()

    dropped_phages = []
    if phage_column in df.columns:
        before = set(df[phage_column].astype(str).unique())
        after = set(filtered[phage_column].astype(str).unique())
        dropped_phages = sorted(before - after)

    filtered.to_csv(output_path, index=False)
    n_phages = filtered[phage_column].nunique() if phage_column in filtered.columns else 0
    return output_path, len(filtered), n_phages, dropped_phages


def _write_quadrant_matrix(interaction_matrix, train_strains, train_phages, output_path,
                           strain_column='strain', phage_column='phage'):
    """Write an interaction matrix filtered to the training quadrant (both axes).

    Keeps only rows where strain is in ``train_strains`` AND phage is in
    ``train_phages``. Used for corner-prediction CV (leave-one-genus-out where
    BOTH the held-out genus's strains and its phages are excluded from training).

    Returns (output_path, n_rows, n_strains, n_phages).
    """
    df = pd.read_csv(interaction_matrix)
    for col in (strain_column, phage_column):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in interaction matrix ({list(df.columns)}).")
    ts = {str(s) for s in train_strains}
    tp = {str(p) for p in train_phages}
    df[strain_column] = df[strain_column].astype(str)
    df[phage_column] = df[phage_column].astype(str)
    filtered = df[df[strain_column].isin(ts) & df[phage_column].isin(tp)].copy()
    filtered.to_csv(output_path, index=False)
    return (output_path, len(filtered),
            filtered[strain_column].nunique(), filtered[phage_column].nunique())


def split_strains_kfold(full_strain_list, iteration, n_folds=10):
    """Deterministically split strains into modeling/validation sets for one fold.

    ``iteration`` is 1-based. Iterations 1..n_folds form round 0, the next
    n_folds form round 1, and so on. Within a round the strains are shuffled
    with a per-round seed (the round index), then partitioned into ``n_folds``
    roughly equal folds with the remainder distributed across the leading
    folds. The fold for this iteration becomes the validation set; the rest is
    the modeling set.

    Returns ``(modeling_strains, validation_strains)``.
    """
    if not full_strain_list:
        raise ValueError("Cannot split an empty strain list.")

    cv_round = (iteration - 1) // n_folds
    fold_index = (iteration - 1) % n_folds

    logger.info(
        f"K-fold split: iteration {iteration} -> round {cv_round}, fold {fold_index} "
        f"({n_folds} folds)"
    )

    # Deterministic shuffle per round. Sort first so the input order does not
    # affect the result.
    random.seed(cv_round)
    strain_list_copy = sorted(full_strain_list)
    random.shuffle(strain_list_copy)

    fold_size = len(strain_list_copy) // n_folds
    remainder = len(strain_list_copy) % n_folds

    folds = []
    start = 0
    for i in range(n_folds):
        # Distribute the remainder across the leading folds.
        end = start + fold_size + (1 if i < remainder else 0)
        folds.append(strain_list_copy[start:end])
        start = end

    validation_strains = folds[fold_index]
    modeling_strains = [s for j, fold in enumerate(folds) for s in fold if j != fold_index]

    logger.info(
        f"Fold {fold_index}: {len(modeling_strains)} modeling strains, "
        f"{len(validation_strains)} validation strains"
    )
    return modeling_strains, validation_strains


def _load_group_map(group_metadata, full_strain_list, strain_column='strain', group_column='group'):
    """Read a strain -> group mapping CSV, restricted to the usable strains.

    Strains that are absent from the metadata, or whose group value is missing
    (NaN/empty), are logged and excluded from cross-validation (they cannot be
    held out as part of any group).

    Returns ``(group_map, ungrouped)`` where ``group_map`` is ``{strain: group}``
    for strains that have a valid group AND appear in ``full_strain_list``, and
    ``ungrouped`` is the sorted list of usable strains that have no valid group
    (missing from the metadata or with an empty/NaN group value).
    """
    logger.info(f"Reading group metadata: {group_metadata}")
    meta = pd.read_csv(group_metadata)

    for col in (strain_column, group_column):
        if col not in meta.columns:
            raise ValueError(
                f"Column '{col}' not found in group metadata. "
                f"Available columns: {list(meta.columns)}"
            )

    full_set = set(full_strain_list)
    group_map = {}
    missing_group = []
    for _, row in meta.iterrows():
        strain = str(row[strain_column])
        if strain not in full_set:
            continue  # not a usable strain (no FASTA/phenotype); ignore silently
        group = row[group_column]
        if pd.isna(group) or str(group).strip() == '':
            missing_group.append(strain)
            continue
        group_map[strain] = str(group)

    # Strains usable for CV but with no group entry at all.
    no_metadata = sorted(full_set - set(group_map) - set(missing_group))

    ungrouped = sorted(set(missing_group) | set(no_metadata))
    if ungrouped:
        logger.info(
            f"{len(ungrouped)} usable strain(s) have no group label: {ungrouped}"
        )
    logger.info(
        f"Group map covers {len(group_map)} strain(s) across "
        f"{len(set(group_map.values()))} group(s)."
    )
    return group_map, ungrouped


def split_strains_by_group(full_strain_list, group_map, ungrouped=None,
                           ungrouped_mode='drop'):
    """Build leave-one-group-out splits: one fold per unique group.

    For each group (sorted by name for determinism), the validation set is every
    strain in that group and the modeling set is every other grouped strain.

    ``ungrouped`` strains (no group label) are handled by ``ungrouped_mode``:
      - 'drop' (default): excluded from every fold (not trained, not validated).
        Clean leave-one-group-out: every strain in every fold has a known group.
      - 'train': added to the modeling set of every fold (never held out, since
        they belong to no group). Maximizes training data; the group defines
        only what is left OUT. These strains never receive a held-out prediction.

    Returns an ordered list of ``(group_name, modeling_strains, validation_strains)``.
    """
    if not group_map:
        raise ValueError("Cannot build group splits from an empty group map.")

    ungrouped = list(ungrouped or [])
    if ungrouped_mode not in ('drop', 'train'):
        raise ValueError(
            f"ungrouped_mode must be 'drop' or 'train', got '{ungrouped_mode}'."
        )

    grouped = sorted(set(full_strain_list) & set(group_map))
    groups = {}
    for strain in grouped:
        groups.setdefault(group_map[strain], []).append(strain)

    # Ungrouped strains added to every fold's modeling set when mode == 'train'.
    always_train = sorted(set(ungrouped) & set(full_strain_list)) if ungrouped_mode == 'train' else []
    if always_train:
        logger.info(
            f"{len(always_train)} ungrouped strain(s) added to every fold's "
            f"modeling set (ungrouped_mode='train'); they are never held out."
        )

    splits = []
    for group_name in sorted(groups):
        validation_strains = sorted(groups[group_name])
        modeling_grouped = [s for s in grouped if group_map[s] != group_name]
        modeling_strains = sorted(set(modeling_grouped) | set(always_train))
        if not modeling_strains:
            logger.warning(
                f"Group '{group_name}': no modeling strains remain (all strains "
                f"belong to this group); skipping this fold."
            )
            continue
        splits.append((group_name, modeling_strains, validation_strains))
        logger.info(
            f"Group '{group_name}': {len(modeling_strains)} modeling, "
            f"{len(validation_strains)} validation strain(s)."
        )

    if not splits:
        raise ValueError("No usable leave-one-group-out folds could be built.")
    return splits


def load_predefined_folds(folds_file, full_strain_list, strain_column='strain'):
    """Build an explicit, possibly-overlapping split plan from a fold-definition CSV.

    The caller supplies the exact fold membership (e.g. for leave-one-serotype-out
    with symmetric hold-out of cross-reactive strains), so the package applies no
    grouping/splitting logic of its own -- it runs precisely the folds given.

    Expected long-format CSV with columns ``fold_label``, ``<strain_column>``,
    ``role`` where role is 'validation' or 'modeling'. Folds may overlap (a strain
    may be validation in one fold and modeling/validation in another). Within a
    single fold a strain must not be both validation and modeling.

    Strains not present in ``full_strain_list`` (no FASTA/phenotype) are dropped
    from the fold with a warning. Folds are returned in first-seen order.

    Returns an ordered list of ``(fold_label, modeling_strains, validation_strains)``.
    """
    logger.info(f"Reading predefined folds: {folds_file}")
    df = pd.read_csv(folds_file)

    required = {'fold_label', strain_column, 'role'}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Folds file must contain columns {sorted(required)}; "
            f"found {list(df.columns)}."
        )

    df['role'] = df['role'].astype(str).str.strip().str.lower()
    bad_roles = sorted(set(df['role']) - {'validation', 'modeling'})
    if bad_roles:
        raise ValueError(
            f"Folds file 'role' must be 'validation' or 'modeling'; "
            f"found unexpected value(s): {bad_roles}."
        )

    df['fold_label'] = df['fold_label'].astype(str)
    df[strain_column] = df[strain_column].astype(str)
    full_set = set(full_strain_list)

    # Preserve first-seen fold order.
    fold_order = list(dict.fromkeys(df['fold_label'].tolist()))

    splits = []
    for label in fold_order:
        sub = df[df['fold_label'] == label]
        validation = sub.loc[sub['role'] == 'validation', strain_column].tolist()
        modeling = sub.loc[sub['role'] == 'modeling', strain_column].tolist()

        # A strain must not be both validation and modeling within one fold.
        overlap = set(validation) & set(modeling)
        if overlap:
            raise ValueError(
                f"Fold '{label}': {len(overlap)} strain(s) listed as BOTH validation "
                f"and modeling: {sorted(overlap)}. Each strain must have one role per fold."
            )

        # Drop strains we cannot actually use (no FASTA/phenotype).
        val_usable = sorted(s for s in dict.fromkeys(validation) if s in full_set)
        mod_usable = sorted(s for s in dict.fromkeys(modeling) if s in full_set)
        dropped = sorted((set(validation) | set(modeling)) - full_set)
        if dropped:
            logger.warning(
                f"Fold '{label}': {len(dropped)} strain(s) in the folds file are not "
                f"usable (missing FASTA/phenotype) and were dropped: {dropped}"
            )

        if not val_usable or not mod_usable:
            logger.warning(
                f"Fold '{label}': {len(mod_usable)} modeling / {len(val_usable)} "
                f"validation usable strain(s); skipping (a fold needs both)."
            )
            continue

        splits.append((label, mod_usable, val_usable))
        logger.info(
            f"Fold '{label}': {len(mod_usable)} modeling, {len(val_usable)} validation strain(s)."
        )

    if not splits:
        raise ValueError("No usable predefined folds could be built from the folds file.")
    return splits


def _select_best_cutoff(iteration_output_dir):
    """Return the cutoff with the highest MCC (tie-break on higher cut_off)."""
    metrics_file = os.path.join(
        iteration_output_dir,
        'modeling_results', 'model_performance', 'model_performance_metrics.csv',
    )
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    return metrics_df['cut_off'].values[0]


def _run_single_fold(
    iteration,
    modeling_strains,
    validation_strains,
    input_strain_dir,
    input_phage_dir,
    interaction_matrix,
    output_dir,
    min_seq_id,
    coverage,
    sensitivity,
    threads,
    num_runs_fs,
    num_runs_modeling,
    use_dynamic_weights,
    weights_method,
    use_clustering,
    cluster_method,
    n_clusters,
    min_cluster_size,
    min_samples,
    cluster_selection_epsilon,
    check_feature_presence,
    filter_by_cluster_presence,
    min_cluster_presence,
    max_ram,
    use_feature_clustering,
    feature_cluster_method,
    feature_n_clusters,
    feature_min_cluster_presence,
    duplicate_all,
):
    """Run training + held-out prediction for a single fold.

    Returns the path to the fold's median-predictions CSV. Idempotent: if that
    file already exists the fold is skipped and the path is returned as-is.
    """
    iteration_output_dir = os.path.join(output_dir, f'iteration_{iteration}')
    median_predictions_file = os.path.join(
        iteration_output_dir, 'model_validation', 'predict_results',
        'strain_median_predictions.csv',
    )

    if os.path.exists(median_predictions_file):
        logger.info(f"Iteration {iteration} already complete, skipping.")
        return median_predictions_file

    logger.info(f"Starting iteration {iteration}...")
    os.makedirs(iteration_output_dir, exist_ok=True)
    modeling_tmp_dir = os.path.join(iteration_output_dir, 'tmp')

    # --- Write the (precomputed) strain split for this fold ---
    modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')

    if not modeling_strains or not validation_strains:
        raise RuntimeError(
            f"Empty strain split for iteration {iteration}: "
            f"{len(modeling_strains)} modeling, {len(validation_strains)} validation"
        )

    if not os.path.exists(modeling_strains_path) or not os.path.exists(validation_strains_path):
        pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
        pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
    else:
        logger.info(f"Strain lists already exist for iteration {iteration}; reusing.")

    # --- Filter the interaction matrix to the modeling strains ---
    # Training builds phage features only from phages present in these remaining
    # interactions; phages tested solely against held-out strains drop out.
    modeling_matrix_path = os.path.join(iteration_output_dir, 'modeling_interaction_matrix.csv')
    _, n_rows, n_phages, dropped_phages = _write_modeling_matrix(
        interaction_matrix, modeling_strains, modeling_matrix_path,
    )
    logger.info(
        f"Iteration {iteration}: modeling interaction matrix has {n_rows} rows "
        f"across {n_phages} phage(s)."
    )
    if dropped_phages:
        logger.info(
            f"Iteration {iteration}: {len(dropped_phages)} phage(s) absent from "
            f"modeling interactions, excluded from feature building: {dropped_phages}"
        )

    # --- Step 1: train protein-family models on the modeling strains ---
    metrics_file = os.path.join(
        iteration_output_dir,
        'modeling_results', 'model_performance', 'model_performance_metrics.csv',
    )
    if not os.path.exists(metrics_file):
        logger.info(f"Iteration {iteration}: running protein family workflow...")
        run_protein_family_workflow(
            input_path_strain=input_strain_dir,
            input_path_phage=input_phage_dir,
            phenotype_matrix=modeling_matrix_path,
            output_dir=iteration_output_dir,
            tmp_dir=modeling_tmp_dir,
            min_seq_id=min_seq_id,
            coverage=coverage,
            sensitivity=sensitivity,
            threads=threads,
            phenotype_column='interaction',
            phage_column='phage',
            strain_list=modeling_strains_path,
            phage_list=modeling_matrix_path,
            filter_type='strain',
            num_runs_fs=num_runs_fs,
            num_runs_modeling=num_runs_modeling,
            use_dynamic_weights=use_dynamic_weights,
            weights_method=weights_method,
            use_clustering=use_clustering,
            cluster_method=cluster_method,
            n_clusters=n_clusters,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_epsilon=cluster_selection_epsilon,
            check_feature_presence=check_feature_presence,
            filter_by_cluster_presence=filter_by_cluster_presence,
            min_cluster_presence=min_cluster_presence,
            # Force bootstrapping=True so the clustering step resolves duplicate
            # protein IDs across ALL strains (modeling AND held-out validation),
            # writing every strain's "strain::protein"-prefixed FASTA into
            # modified_AAs/. Held-out prediction then assigns validation strains
            # from that directory against the matching cluster DB. With the
            # default (False), only modeling strains are prefixed and validation
            # assignment finds no matching files.
            bootstrapping=True,
            max_ram=max_ram,
            use_feature_clustering=use_feature_clustering,
            feature_cluster_method=feature_cluster_method,
            feature_n_clusters=feature_n_clusters,
            feature_min_cluster_presence=feature_min_cluster_presence,
        )
    else:
        logger.info(f"Iteration {iteration}: modeling results already exist; skipping training.")

    # --- Step 2: select the cutoff with the highest MCC ---
    best_cutoff = _select_best_cutoff(iteration_output_dir)
    model_dir = os.path.join(iteration_output_dir, 'modeling_results', f'{best_cutoff}')

    # --- Step 3: predict interactions for the held-out validation strains ---
    validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
    os.makedirs(validation_output_dir, exist_ok=True)
    validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')

    # Assign features to the held-out validation strains. Because training ran
    # with bootstrapping=True, the clustering step prefixed EVERY strain's
    # protein IDs to "strain::protein" and wrote them under modified_AAs/strain/
    # -- so the held-out validation strains are present there with IDs that
    # match the cluster database. Prefer that directory; fall back to the
    # original input directory only if it is missing (e.g. no duplicate IDs
    # existed, so no prefixing was needed and original IDs already match).
    modified_aa_dir = os.path.join(iteration_output_dir, 'strain', 'modified_AAs', 'strain')
    if os.path.isdir(modified_aa_dir) and os.listdir(modified_aa_dir):
        strain_input_dir = modified_aa_dir
    else:
        strain_input_dir = input_strain_dir
    logger.info(f"Iteration {iteration}: assigning validation features from {strain_input_dir}")

    select_feature_table = os.path.join(
        iteration_output_dir, 'feature_selection', 'filtered_feature_tables',
        f'select_feature_table_{best_cutoff}.csv',
    )

    logger.info(f"Iteration {iteration}: running prediction workflow...")
    assign_predict_workflow(
        input_dir=strain_input_dir,
        genome_list=validation_strains_path,
        mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
        clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
        feature_map=os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv'),
        tmp_dir=validation_tmp_dir,
        suffix='faa',
        model_dir=model_dir,
        feature_table=select_feature_table,
        phage_feature_table_path=os.path.join(
            iteration_output_dir, 'phage', 'features', 'feature_table.csv'
        ),
        output_dir=validation_output_dir,
        threads=threads,
        genome_type='strain',
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        duplicate_all=duplicate_all,
    )

    logger.info(f"Iteration {iteration} completed successfully.")
    return median_predictions_file


def run_shared_clustering(
    input_strain_dir,
    input_phage_dir,
    output_dir,
    min_seq_id=0.4,
    coverage=0.8,
    sensitivity=7.5,
    suffix='faa',
    threads=4,
    strain_column='strain',
    phage_column='phage',
    force_restart=False,
):
    """Cluster ALL strains and phages once; produce shared artifacts for every fold.

    MMseqs2 clustering is sequence-only and label-independent, so it can be run
    a single time over the full dataset and reused across folds. This writes the
    raw presence/absence matrices, cluster TSVs, and MMseqs2 databases that each
    fold then filters (per-fold) into training-only feature tables.

    Runs with bootstrapping=True so duplicate protein IDs are resolved across ALL
    genomes (so held-out strains carry the same 'strain::protein' namespace the
    cluster DB uses, which validation feature-assignment depends on).

    Returns a dict of shared artifact paths::

        {
          'strain': {'presence_absence': ..., 'clusters_tsv': ..., 'mmseqs_db': ...},
          'phage':  {'presence_absence': ..., 'clusters_tsv': ..., 'mmseqs_db': ...},
        }
    """
    os.makedirs(output_dir, exist_ok=True)
    artifacts = {}

    for genome_type, input_dir, column in (
        ('strain', input_strain_dir, strain_column),
        ('phage', input_phage_dir, phage_column),
    ):
        logger.info(f"Shared clustering: {genome_type}s (once, over all genomes)...")
        type_output_dir = os.path.join(output_dir, genome_type)
        type_tmp_dir = os.path.join(output_dir, 'tmp', genome_type)

        run_clustering_workflow(
            input_dir, type_output_dir, type_tmp_dir,
            min_seq_id, coverage, sensitivity, suffix, threads,
            'none', column, False, bootstrapping=True, clear_tmp=False,
            force_restart=force_restart,
        )

        artifacts[genome_type] = {
            'presence_absence': os.path.join(type_output_dir, 'presence_absence_matrix.csv'),
            'clusters_tsv': os.path.join(type_output_dir, 'clusters.tsv'),
            'mmseqs_db': os.path.join(type_tmp_dir, 'mmseqs_db'),
        }

    logger.info(f"Shared clustering complete. Artifacts under {output_dir}/")
    return artifacts


def run_fold_from_shared(
    shared_dir,
    iteration_output_dir,
    modeling_strains_path,
    validation_strains_path,
    modeling_matrix_path,
    input_strain_dir,
    threads=4,
    num_runs_fs=25,
    num_runs_modeling=50,
    min_seq_id=0.4,
    coverage=0.8,
    sensitivity=7.5,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=False,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=2,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
    max_ram=40,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    duplicate_all=True,
):
    """Run one fold using shared clustering artifacts (no per-fold MMseqs2).

    Given the shared presence/absence matrices + MMseqs2 DBs from
    ``run_shared_clustering``, this performs, for a single fold:

      B. filter the shared presence/absence matrix to the modeling strains and
         re-derive the hash-collapsed features (per-fold, training rows only --
         this is where leakage would occur if shared, so it is NOT shared);
      C. merge strain + phage features with the fold's modeling interaction matrix;
      D. feature selection + modeling on that merged table;
      E. assign + predict the held-out strains against the SHARED MMseqs2 DB,
         mapped through the fold's (per-fold) selected_features.csv.

    This is a standalone function so a SLURM submitter can call it as one job per
    fold, depending on a single shared-clustering job.

    Returns the path to the fold's median-predictions CSV.
    """
    shared = {
        gt: {
            'presence_absence': os.path.join(shared_dir, gt, 'presence_absence_matrix.csv'),
            'clusters_tsv': os.path.join(shared_dir, gt, 'clusters.tsv'),
            'mmseqs_db': os.path.join(shared_dir, 'tmp', gt, 'mmseqs_db'),
        }
        for gt in ('strain', 'phage')
    }

    os.makedirs(iteration_output_dir, exist_ok=True)

    # --- B. Per-fold filter -> hash-collapse -> assign (strain and phage) ---
    # Strain features are filtered to the modeling strains; phage features are
    # filtered to the phages present in the modeling interaction matrix.
    strain_features_dir = os.path.join(iteration_output_dir, 'strain', 'features')
    run_feature_assignment(
        shared['strain']['presence_absence'], strain_features_dir,
        source='strain', select=modeling_strains_path, select_column='strain',
        max_ram=max_ram,
    )

    # Build a phage list (phages present in the modeling interaction matrix).
    modeling_phages_path = os.path.join(iteration_output_dir, 'modeling_phages.csv')
    mm = pd.read_csv(modeling_matrix_path)
    pd.DataFrame({'phage': sorted(mm['phage'].astype(str).unique())}).to_csv(
        modeling_phages_path, index=False
    )
    phage_features_dir = os.path.join(iteration_output_dir, 'phage', 'features')
    run_feature_assignment(
        shared['phage']['presence_absence'], phage_features_dir,
        source='phage', select=modeling_phages_path, select_column='phage',
        max_ram=max_ram,
    )

    # --- C. Merge features with the fold's modeling interaction matrix ---
    merged_dir = os.path.join(iteration_output_dir, 'merged')
    os.makedirs(merged_dir, exist_ok=True)
    merged_table = merge_feature_tables(
        strain_features=os.path.join(strain_features_dir, 'feature_table.csv'),
        phenotype_matrix=modeling_matrix_path,
        output_dir=merged_dir,
        sample_column='strain',
        phage_features=os.path.join(phage_features_dir, 'feature_table.csv'),
        remove_suffix=False,
        use_feature_clustering=use_feature_clustering,
        feature_cluster_method=feature_cluster_method,
        feature_n_clusters=feature_n_clusters,
        feature_min_cluster_presence=feature_min_cluster_presence,
        phenotype_column='interaction',
    )

    # --- D. Feature selection + modeling (no clustering) ---
    run_modeling_workflow_from_feature_table(
        full_feature_table=merged_table,
        output_dir=iteration_output_dir,
        threads=threads,
        filter_type='strain',
        num_runs_fs=num_runs_fs,
        num_runs_modeling=num_runs_modeling,
        sample_column='strain',
        phage_column='phage',
        phenotype_column='interaction',
        task_type='classification',
        max_ram=max_ram,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        check_feature_presence=check_feature_presence,
        filter_by_cluster_presence=filter_by_cluster_presence,
        min_cluster_presence=min_cluster_presence,
    )

    # --- E. Assign + predict held-out strains against the SHARED DB ---
    best_cutoff = _select_best_cutoff(iteration_output_dir)
    model_dir = os.path.join(iteration_output_dir, 'modeling_results', f'{best_cutoff}')

    validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
    os.makedirs(validation_output_dir, exist_ok=True)
    validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')

    select_feature_table = os.path.join(
        iteration_output_dir, 'feature_selection', 'filtered_feature_tables',
        f'select_feature_table_{best_cutoff}.csv',
    )

    logger.info("Fold (shared): assigning + predicting held-out strains against shared DB...")
    assign_predict_workflow(
        input_dir=input_strain_dir,
        genome_list=validation_strains_path,
        mmseqs_db=shared['strain']['mmseqs_db'],
        clusters_tsv=shared['strain']['clusters_tsv'],
        feature_map=os.path.join(strain_features_dir, 'selected_features.csv'),
        tmp_dir=validation_tmp_dir,
        suffix='faa',
        model_dir=model_dir,
        feature_table=select_feature_table,
        phage_feature_table_path=os.path.join(phage_features_dir, 'feature_table.csv'),
        output_dir=validation_output_dir,
        threads=threads,
        genome_type='strain',
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        duplicate_all=duplicate_all,
    )

    median_predictions_file = os.path.join(
        validation_output_dir, 'predict_results', 'strain_median_predictions.csv',
    )
    return median_predictions_file


def run_corner_fold_from_shared(
    fold_label,
    training_strains_file,
    training_phages_file,
    validation_strains_file,
    validation_phages_file,
    shared_dir,
    interaction_matrix,
    input_strain_dir,
    input_phage_dir,
    output_dir,
    threads=4,
    num_runs_fs=25,
    num_runs_modeling=50,
    min_seq_id=0.4,
    coverage=0.8,
    sensitivity=7.5,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=False,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=2,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
    max_ram=40,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    duplicate_all=True,
    strain_column='strain',
    phage_column='phage',
):
    """Run ONE leave-one-genus-out CORNER fold against shared clustering artifacts.

    The four ``*_file`` arguments are CSVs each listing one genome ID per row (a
    header is allowed; the first column is used). The caller controls the fold
    membership; the package applies no grouping logic and is agnostic to any
    fold-directory naming convention. BOTH the held-out genus's strains AND its
    phages are excluded from training; the model is trained on the training
    quadrant (training_strains x training_phages) and predicts the held-out
    CORNER (validation_strains x validation_phages -- pairs where both are unseen).

    Steps:
      B. Build training-only features for BOTH axes from the shared P/A matrices
         (filter to training strains / training phages, then hash-collapse). This
         is the leak-prevention: held-out strains AND phages are absent from the
         feature collapse.
      C. Merge into the training-quadrant interaction matrix.
      D. Feature selection + modeling on that quadrant.
      E. Corner prediction (composed two-step):
         E1. Assign validation PHAGES against the shared phage DB, mapped through
             the fold's phage selected_features.csv -> a validation-phage feature
             table in the model's phage feature vocabulary.
         E2. Assign validation STRAINS against the shared strain DB and predict
             them paired with that validation-phage feature table -> predictions
             for exactly validation_strains x validation_phages.

    Returns the path to the corner median-predictions CSV.
    """
    def _read_ids(path):
        s = pd.read_csv(path)
        return [str(x) for x in s.iloc[:, 0].tolist()]

    train_strains = _read_ids(training_strains_file)
    train_phages = _read_ids(training_phages_file)
    val_strains = _read_ids(validation_strains_file)
    val_phages = _read_ids(validation_phages_file)

    median_predictions_file = os.path.join(
        output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv',
    )
    if os.path.exists(median_predictions_file):
        logger.info(f"Corner fold '{fold_label}' already complete, skipping.")
        return median_predictions_file

    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, 'fold_group.txt'), 'w') as fh:
        fh.write(f"{fold_label}\n")

    logger.info(
        f"Corner fold '{fold_label}': train {len(train_strains)} strains x "
        f"{len(train_phages)} phages -> predict corner {len(val_strains)} x "
        f"{len(val_phages)}."
    )

    shared = {
        gt: {
            'presence_absence': os.path.join(shared_dir, gt, 'presence_absence_matrix.csv'),
            'clusters_tsv': os.path.join(shared_dir, gt, 'clusters.tsv'),
            'mmseqs_db': os.path.join(shared_dir, 'tmp', gt, 'mmseqs_db'),
        }
        for gt in ('strain', 'phage')
    }

    # Write the per-axis training lists genophi's select-filtering expects.
    train_strains_path = os.path.join(output_dir, 'training_strains.csv')
    train_phages_path = os.path.join(output_dir, 'training_phages.csv')
    val_strains_path = os.path.join(output_dir, 'validation_strains.csv')
    val_phages_path = os.path.join(output_dir, 'validation_phages.csv')
    pd.DataFrame({'strain': train_strains}).to_csv(train_strains_path, index=False)
    pd.DataFrame({'phage': train_phages}).to_csv(train_phages_path, index=False)
    pd.DataFrame({'strain': val_strains}).to_csv(val_strains_path, index=False)
    pd.DataFrame({'phage': val_phages}).to_csv(val_phages_path, index=False)

    # Training-quadrant interaction matrix (both axes filtered).
    train_matrix_path = os.path.join(output_dir, 'training_quadrant_matrix.csv')
    _, n_rows, n_s, n_p = _write_quadrant_matrix(
        interaction_matrix, train_strains, train_phages, train_matrix_path,
        strain_column=strain_column, phage_column=phage_column,
    )
    logger.info(
        f"Corner fold '{fold_label}': training quadrant {n_rows} rows "
        f"({n_s} strains x {n_p} phages)."
    )

    # --- B. Training-only features for both axes (leak-free collapse) ---
    strain_features_dir = os.path.join(output_dir, 'strain', 'features')
    run_feature_assignment(
        shared['strain']['presence_absence'], strain_features_dir,
        source='strain', select=train_strains_path, select_column='strain',
        max_ram=max_ram,
    )
    phage_features_dir = os.path.join(output_dir, 'phage', 'features')
    run_feature_assignment(
        shared['phage']['presence_absence'], phage_features_dir,
        source='phage', select=train_phages_path, select_column='phage',
        max_ram=max_ram,
    )

    # --- C. Merge into the training-quadrant matrix ---
    merged_dir = os.path.join(output_dir, 'merged')
    os.makedirs(merged_dir, exist_ok=True)
    merged_table = merge_feature_tables(
        strain_features=os.path.join(strain_features_dir, 'feature_table.csv'),
        phenotype_matrix=train_matrix_path,
        output_dir=merged_dir,
        sample_column='strain',
        phage_features=os.path.join(phage_features_dir, 'feature_table.csv'),
        remove_suffix=False,
        use_feature_clustering=use_feature_clustering,
        feature_cluster_method=feature_cluster_method,
        feature_n_clusters=feature_n_clusters,
        feature_min_cluster_presence=feature_min_cluster_presence,
        phenotype_column='interaction',
    )

    # --- D. Feature selection + modeling ---
    run_modeling_workflow_from_feature_table(
        full_feature_table=merged_table,
        output_dir=output_dir,
        threads=threads,
        filter_type='strain',
        num_runs_fs=num_runs_fs,
        num_runs_modeling=num_runs_modeling,
        sample_column='strain',
        phage_column='phage',
        phenotype_column='interaction',
        task_type='classification',
        max_ram=max_ram,
        use_dynamic_weights=use_dynamic_weights,
        weights_method=weights_method,
        use_clustering=use_clustering,
        cluster_method=cluster_method,
        n_clusters=n_clusters,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=cluster_selection_epsilon,
        check_feature_presence=check_feature_presence,
        filter_by_cluster_presence=filter_by_cluster_presence,
        min_cluster_presence=min_cluster_presence,
    )

    best_cutoff = _select_best_cutoff(output_dir)
    model_dir = os.path.join(output_dir, 'modeling_results', f'{best_cutoff}')
    select_feature_table = os.path.join(
        output_dir, 'feature_selection', 'filtered_feature_tables',
        f'select_feature_table_{best_cutoff}.csv',
    )

    # --- E1. Assign validation PHAGES -> validation-phage feature table ---
    # Mapped through the fold's phage selected_features.csv so held-out phages
    # are described in the MODEL's phage feature vocabulary.
    valphage_dir = os.path.join(output_dir, 'validation_phage_features')
    valphage_tmp = os.path.join(valphage_dir, 'tmp')
    os.makedirs(valphage_dir, exist_ok=True)
    logger.info(f"Corner fold '{fold_label}': assigning validation phages to model phage features...")
    run_assign_features_workflow(
        input_dir=input_phage_dir,
        mmseqs_db=shared['phage']['mmseqs_db'],
        tmp_dir=valphage_tmp,
        output_dir=valphage_dir,
        feature_map=os.path.join(phage_features_dir, 'selected_features.csv'),
        clusters_tsv=shared['phage']['clusters_tsv'],
        genome_type='phage',
        genome_list=val_phages_path,
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        threads=threads,
        duplicate_all=duplicate_all,
    )
    val_phage_feature_table = os.path.join(valphage_dir, 'phage_combined_feature_table.csv')

    # --- E2. Assign validation STRAINS, predict paired with validation phages ---
    validation_output_dir = os.path.join(output_dir, 'model_validation')
    os.makedirs(validation_output_dir, exist_ok=True)
    validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')
    logger.info(f"Corner fold '{fold_label}': predicting held-out corner (val strains x val phages)...")
    assign_predict_workflow(
        input_dir=input_strain_dir,
        genome_list=val_strains_path,
        mmseqs_db=shared['strain']['mmseqs_db'],
        clusters_tsv=shared['strain']['clusters_tsv'],
        feature_map=os.path.join(strain_features_dir, 'selected_features.csv'),
        tmp_dir=validation_tmp_dir,
        suffix='faa',
        model_dir=model_dir,
        feature_table=select_feature_table,
        phage_feature_table_path=val_phage_feature_table,
        output_dir=validation_output_dir,
        threads=threads,
        genome_type='strain',
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        duplicate_all=duplicate_all,
    )

    return median_predictions_file


def _run_fold_shared_wrapper(
    iteration, modeling_strains, validation_strains, shared_dir,
    input_strain_dir, interaction_matrix, output_dir, **fold_kwargs
):
    """Prepare per-fold split/matrix files, then run the shared-clustering fold.

    Mirrors the pre-steps of _run_single_fold (write split CSVs, write the
    modeling-strain-filtered interaction matrix, idempotency skip) so the
    shared-clustering path produces the same per-fold layout.
    """
    iteration_output_dir = os.path.join(output_dir, f'iteration_{iteration}')
    median_predictions_file = os.path.join(
        iteration_output_dir, 'model_validation', 'predict_results',
        'strain_median_predictions.csv',
    )
    if os.path.exists(median_predictions_file):
        logger.info(f"Iteration {iteration} already complete, skipping.")
        return median_predictions_file

    os.makedirs(iteration_output_dir, exist_ok=True)

    if not modeling_strains or not validation_strains:
        raise RuntimeError(
            f"Empty strain split for iteration {iteration}: "
            f"{len(modeling_strains)} modeling, {len(validation_strains)} validation"
        )

    modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')
    pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
    pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)

    modeling_matrix_path = os.path.join(iteration_output_dir, 'modeling_interaction_matrix.csv')
    _, n_rows, n_phages, dropped_phages = _write_modeling_matrix(
        interaction_matrix, modeling_strains, modeling_matrix_path,
    )
    logger.info(
        f"Iteration {iteration}: modeling matrix {n_rows} rows / {n_phages} phage(s)."
    )
    if dropped_phages:
        logger.info(
            f"Iteration {iteration}: {len(dropped_phages)} phage(s) excluded "
            f"(no modeling interactions): {dropped_phages}"
        )

    return run_fold_from_shared(
        shared_dir=shared_dir,
        iteration_output_dir=iteration_output_dir,
        modeling_strains_path=modeling_strains_path,
        validation_strains_path=validation_strains_path,
        modeling_matrix_path=modeling_matrix_path,
        input_strain_dir=input_strain_dir,
        **fold_kwargs,
    )


def _read_fold_label_dir(fold_dir, default):
    """Return the fold label written in fold_group.txt, else the default."""
    label_file = os.path.join(fold_dir, 'fold_group.txt')
    if os.path.exists(label_file):
        with open(label_file) as fh:
            t = fh.read().strip()
            if t:
                return t
    return default


def run_predefined_fold(
    fold_label,
    folds_file,
    shared_dir,
    interaction_matrix,
    input_strain_dir,
    fold_output_dir,
    strain_column='strain',
    **fold_kwargs
):
    """Run ONE predefined fold (by label) against shared clustering artifacts.

    Standalone entrypoint for SLURM fan-out: a single fold job reads its fold's
    modeling/validation strains from ``folds_file``, writes the per-fold split +
    filtered interaction matrix into ``fold_output_dir``, and runs
    ``run_fold_from_shared`` against the already-built ``shared_dir``.

    ``fold_output_dir`` is keyed by the caller (e.g. .../folds/<label>/) so each
    job is independent and resumable. Returns the median-predictions path, or the
    existing path if the fold is already complete.
    """
    median_predictions_file = os.path.join(
        fold_output_dir, 'model_validation', 'predict_results',
        'strain_median_predictions.csv',
    )
    if os.path.exists(median_predictions_file):
        logger.info(f"Fold '{fold_label}' already complete, skipping.")
        return median_predictions_file

    # Pull this fold's split from the folds file (validates roles/disjointness).
    df = pd.read_csv(folds_file)
    df['fold_label'] = df['fold_label'].astype(str)
    df['role'] = df['role'].astype(str).str.strip().str.lower()
    df[strain_column] = df[strain_column].astype(str)
    sub = df[df['fold_label'] == str(fold_label)]
    if sub.empty:
        raise ValueError(
            f"Fold label '{fold_label}' not found in {folds_file}. "
            f"Available: {sorted(df['fold_label'].unique())}"
        )
    modeling = sorted(set(sub.loc[sub['role'] == 'modeling', strain_column]))
    validation = sorted(set(sub.loc[sub['role'] == 'validation', strain_column]))
    overlap = set(modeling) & set(validation)
    if overlap:
        raise ValueError(
            f"Fold '{fold_label}': strain(s) in both roles: {sorted(overlap)}."
        )
    if not modeling or not validation:
        raise ValueError(
            f"Fold '{fold_label}': {len(modeling)} modeling / {len(validation)} "
            f"validation strain(s); both required."
        )

    os.makedirs(fold_output_dir, exist_ok=True)
    with open(os.path.join(fold_output_dir, 'fold_group.txt'), 'w') as fh:
        fh.write(f"{fold_label}\n")

    modeling_strains_path = os.path.join(fold_output_dir, 'modeling_strains.csv')
    validation_strains_path = os.path.join(fold_output_dir, 'validation_strains.csv')
    pd.DataFrame(modeling, columns=['strain']).to_csv(modeling_strains_path, index=False)
    pd.DataFrame(validation, columns=['strain']).to_csv(validation_strains_path, index=False)

    modeling_matrix_path = os.path.join(fold_output_dir, 'modeling_interaction_matrix.csv')
    _, n_rows, n_phages, dropped_phages = _write_modeling_matrix(
        interaction_matrix, modeling, modeling_matrix_path,
    )
    logger.info(
        f"Fold '{fold_label}': {len(modeling)} modeling / {len(validation)} validation "
        f"strain(s); matrix {n_rows} rows / {n_phages} phage(s)."
    )
    if dropped_phages:
        logger.info(
            f"Fold '{fold_label}': {len(dropped_phages)} phage(s) excluded "
            f"(no modeling interactions)."
        )

    return run_fold_from_shared(
        shared_dir=shared_dir,
        iteration_output_dir=fold_output_dir,
        modeling_strains_path=modeling_strains_path,
        validation_strains_path=validation_strains_path,
        modeling_matrix_path=modeling_matrix_path,
        input_strain_dir=input_strain_dir,
        **fold_kwargs,
    )


def aggregate_predefined_folds(folds_dir, interaction_matrix, output_dir=None,
                               strain_column='strain', phage_column='phage',
                               phenotype_column='interaction'):
    """Aggregate per-fold predictions from a SLURM fan-out + compute metrics.

    Scans ``folds_dir`` for subdirectories each containing a completed fold
    (model_validation/predict_results/strain_median_predictions.csv + a
    fold_group.txt label), pools them, and writes final_predictions.csv,
    prediction_summary.csv, and performance/ (global + per-fold metrics, curves).

    ``output_dir`` defaults to ``folds_dir``. Returns the number of folds pooled.
    """
    output_dir = output_dir or folds_dir
    os.makedirs(output_dir, exist_ok=True)

    final_predictions = []
    pooled = 0
    for name in sorted(os.listdir(folds_dir)):
        fold_dir = os.path.join(folds_dir, name)
        if not os.path.isdir(fold_dir):
            continue
        pred_file = os.path.join(
            fold_dir, 'model_validation', 'predict_results',
            'strain_median_predictions.csv',
        )
        if not os.path.exists(pred_file):
            continue
        label = _read_fold_label_dir(fold_dir, name)
        df = pd.read_csv(pred_file)
        df['fold'] = label
        df['iteration'] = name
        final_predictions.append(df)
        pooled += 1

    if not final_predictions:
        logger.error(f"No completed folds found under {folds_dir}.")
        return 0

    final_df = pd.concat(final_predictions, ignore_index=True)
    final_df.to_csv(os.path.join(output_dir, 'final_predictions.csv'), index=False)
    logger.info(f"Pooled {pooled} fold(s) into final_predictions.csv.")

    # Reuse the metrics path (it reads final_predictions.csv + merges truth).
    try:
        _compute_global_metrics(
            output_dir=output_dir,
            interaction_matrix=interaction_matrix,
            strain_column=strain_column,
            phage_column=phage_column,
            phenotype_column=phenotype_column,
        )
    except Exception as e:  # noqa: BLE001
        logger.error(f"Failed to compute global metrics: {e}", exc_info=True)

    return pooled


def _aggregate_results(output_dir, total_iterations):
    """Concatenate per-fold predictions and write summary statistics.

    Tolerant of missing folds: any iteration without a predictions file is
    logged and skipped. Returns the number of folds aggregated.
    """
    final_predictions = []
    completed = 0

    for i in range(1, total_iterations + 1):
        median_predictions_file = os.path.join(
            output_dir, f'iteration_{i}', 'model_validation', 'predict_results',
            'strain_median_predictions.csv',
        )
        if os.path.exists(median_predictions_file):
            df = pd.read_csv(median_predictions_file)
            df['iteration'] = i
            final_predictions.append(df)
            completed += 1
        else:
            logger.warning(f"Results missing for iteration {i}.")

    if not final_predictions:
        logger.error("No fold results found to aggregate.")
        return 0

    final_df = pd.concat(final_predictions, ignore_index=True)
    final_df.to_csv(os.path.join(output_dir, 'final_predictions.csv'), index=False)
    logger.info(
        f"Saved final_predictions.csv with {len(final_df)} predictions "
        f"from {completed} folds."
    )

    if {'strain', 'phage', 'prediction'}.issubset(final_df.columns):
        summary = (
            final_df.groupby(['strain', 'phage'])['prediction']
            .agg(['mean', 'std', 'count'])
            .round(4)
            .reset_index()
        )
        summary.columns = ['strain', 'phage', 'mean_prediction', 'std_prediction', 'n_iterations']
        summary.to_csv(os.path.join(output_dir, 'prediction_summary.csv'), index=False)
        logger.info(
            f"Saved prediction_summary.csv with {len(summary)} unique strain-phage pairs."
        )
    else:
        logger.warning(
            "Could not build prediction summary; expected columns "
            "'strain', 'phage', 'prediction' not all present."
        )

    return completed


def _read_fold_label(output_dir, iteration):
    """Return the held-out fold/group label for an iteration, or a default."""
    label_file = os.path.join(output_dir, f'iteration_{iteration}', 'fold_group.txt')
    if os.path.exists(label_file):
        with open(label_file) as fh:
            label = fh.read().strip()
            if label:
                return label
    return f'fold_{iteration}'


def _compute_global_metrics(
    output_dir,
    interaction_matrix,
    strain_column='strain',
    phage_column='phage',
    phenotype_column='interaction',
):
    """Compute global performance metrics and PR/ROC curves over pooled outer-test folds.

    The per-fold prediction files contain a predicted score but not the ground
    truth, so the pooled predictions (final_predictions.csv) are merged back
    against the interaction matrix on (strain, phage) to recover true labels.
    Every held-out prediction across all folds is treated as one pooled test
    set (each strain is held out exactly once per round).

    Writes, into output_dir/performance/:
        - global_metrics.csv   (AUC, AP, Accuracy, Precision, Recall, F1, MCC, n)
        - roc_curve.png, pr_curve.png

    Returns the metrics dict, or None if metrics could not be computed.
    """
    from sklearn.metrics import (
        roc_curve, precision_recall_curve, roc_auc_score, average_precision_score,
        accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef,
    )
    import matplotlib
    matplotlib.use('Agg')  # headless / HPC-safe backend
    import matplotlib.pyplot as plt

    final_path = os.path.join(output_dir, 'final_predictions.csv')
    if not os.path.exists(final_path):
        logger.warning("No final_predictions.csv found; skipping global metrics.")
        return None

    preds = pd.read_csv(final_path)
    if _SCORE_COL not in preds.columns:
        logger.warning(
            f"Pooled predictions missing '{_SCORE_COL}' column; cannot compute "
            f"global metrics. Found columns: {list(preds.columns)}"
        )
        return None

    # Recover ground-truth labels by merging on (strain, phage).
    truth = pd.read_csv(interaction_matrix)
    required = {strain_column, phage_column, phenotype_column}
    if not required.issubset(truth.columns):
        logger.warning(
            f"Interaction matrix missing one of {required}; cannot recover true "
            f"labels for global metrics. Found: {list(truth.columns)}"
        )
        return None

    truth = truth[[strain_column, phage_column, phenotype_column]].copy()
    # Align identifier dtypes so the merge keys match (strain IDs are read as
    # strings from FASTA filenames elsewhere in the pipeline).
    for col in (strain_column, phage_column):
        preds[col] = preds[col].astype(str)
        truth[col] = truth[col].astype(str)

    merged = preds.merge(truth, on=[strain_column, phage_column], how='inner')
    merged = merged.dropna(subset=[phenotype_column, _SCORE_COL])
    if merged.empty:
        logger.warning(
            "No pooled predictions could be matched to ground-truth labels "
            "(check strain/phage naming between predictions and interaction matrix)."
        )
        return None

    n_matched, n_total = len(merged), len(preds)
    if n_matched < n_total:
        logger.warning(
            f"Matched {n_matched}/{n_total} pooled predictions to ground-truth "
            f"labels; unmatched pairs are excluded from global metrics."
        )

    def _metrics_for(sub):
        """Classification metrics for one prediction subset (a DataFrame)."""
        yt = sub[phenotype_column].astype(int)
        ys = sub[_SCORE_COL].astype(float)
        yp = sub[_PRED_COL].astype(int) if _PRED_COL in sub.columns else (ys > 0.5).astype(int)
        both = yt.nunique() > 1
        return {
            'n_predictions': len(sub),
            'n_positive': int(yt.sum()),
            'n_negative': int((1 - yt).sum()),
            'AUC': roc_auc_score(yt, ys) if both else float('nan'),
            'AP': average_precision_score(yt, ys) if both else float('nan'),
            'Accuracy': accuracy_score(yt, yp),
            'Precision': precision_score(yt, yp, zero_division=0),
            'Recall': recall_score(yt, yp, zero_division=0),
            'F1': f1_score(yt, yp, zero_division=0),
            'MCC': matthews_corrcoef(yt, yp),
        }

    y_true = merged[phenotype_column].astype(int)
    y_score = merged[_SCORE_COL].astype(float)

    perf_dir = os.path.join(output_dir, 'performance')
    os.makedirs(perf_dir, exist_ok=True)

    # --- Global (pooled outer-test) metrics ---
    both_classes = y_true.nunique() > 1
    metrics = _metrics_for(merged)
    pd.DataFrame([metrics]).to_csv(os.path.join(perf_dir, 'global_metrics.csv'), index=False)
    logger.info(
        f"Global metrics over {n_matched} pooled outer-test predictions: "
        f"AUC={metrics['AUC']:.3f}, AP={metrics['AP']:.3f}, MCC={metrics['MCC']:.3f}"
    )

    # --- Per-fold / per-group metrics ---
    # Each iteration held out one fold (a group in cv_mode='group'); label it
    # from fold_group.txt when available so the breakdown is by serotype name.
    if 'iteration' in merged.columns:
        fold_rows = []
        for it, sub in merged.groupby('iteration'):
            label = _read_fold_label(output_dir, int(it))
            row = {'iteration': int(it), 'fold': label}
            row.update(_metrics_for(sub))
            fold_rows.append(row)
        if fold_rows:
            per_fold = pd.DataFrame(fold_rows).sort_values('iteration')
            per_fold.to_csv(os.path.join(perf_dir, 'per_fold_metrics.csv'), index=False)
            logger.info(
                f"Saved per_fold_metrics.csv for {len(per_fold)} held-out fold(s)."
            )

    if not both_classes:
        logger.warning(
            "Pooled outer-test labels contain a single class; ROC/PR curves "
            "are not meaningful and were skipped."
        )
        return metrics

    # ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_score)
    plt.figure()
    plt.plot(fpr, tpr, label=f"ROC (AUC = {metrics['AUC']:.3f})")
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Outer-test ROC (pooled folds)')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(perf_dir, 'roc_curve.png'), dpi=150, bbox_inches='tight')
    plt.close()

    # Precision-recall curve
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    plt.figure()
    plt.step(recall, precision, where='post', label=f"PR (AP = {metrics['AP']:.3f})")
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Outer-test Precision-Recall (pooled folds)')
    plt.legend(loc='lower left')
    plt.savefig(os.path.join(perf_dir, 'pr_curve.png'), dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved ROC and PR curves to {perf_dir}")
    return metrics


def run_nested_cv_workflow(
    input_strain_dir,
    input_phage_dir,
    interaction_matrix,
    output_dir,
    cv_mode='kfold',
    n_folds=10,
    cv_rounds=1,
    group_metadata=None,
    group_column=None,
    group_ungrouped='drop',
    folds_file=None,
    strain_column='strain',
    suffix='faa',
    min_seq_id=0.4,
    coverage=0.8,
    sensitivity=7.5,
    threads=4,
    num_runs_fs=25,
    num_runs_modeling=50,
    use_dynamic_weights=False,
    weights_method='log10',
    use_clustering=False,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=2,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    filter_by_cluster_presence=False,
    min_cluster_presence=2,
    max_ram=40,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    duplicate_all=False,
    clear_tmp=False,
    shared_clustering=False,
):
    """Run nested cross-validation (deterministic k-fold or leave-one-group-out).

    Folds run sequentially in this process. Two split modes:

    - cv_mode='kfold' (default): deterministic (repeated) k-fold. The total
      number of folds is ``n_folds * cv_rounds`` (each round reshuffles the fold
      assignment with a new seed).
    - cv_mode='group': leave-one-group-out. One fold per unique group value in
      ``group_metadata[group_column]``; the held-out group is the validation set
      and all other grouped strains are the modeling set. ``n_folds``/
      ``cv_rounds`` are ignored. Strains with no group label are handled by
      ``group_ungrouped``: 'drop' (default) excludes them from every fold;
      'train' adds them to every fold's modeling set (never held out).
    - cv_mode='predefined': run the exact folds listed in ``folds_file`` (a
      long-format CSV of fold_label/strain/role). Folds may overlap (a strain may
      be held out in several folds) -- useful for leave-one-serotype-out with
      multi-call/cross-reactive antigens. The package applies no grouping logic;
      the caller controls all fold membership (including symmetric hold-out).

    A failing fold is logged and skipped; the run continues and aggregates
    whatever folds completed. Each iteration_N/ records its held-out fold/group
    name in ``fold_group.txt``.

    After aggregation, global performance metrics and pooled outer-test PR/ROC
    curves are written to ``output_dir/performance/`` (global_metrics.csv,
    roc_curve.png, pr_curve.png). The pooled predictions are merged back against
    the interaction matrix to recover true labels.

    Args:
        input_strain_dir (str): Directory of strain FASTA files.
        input_phage_dir (str): Directory of phage FASTA files.
        interaction_matrix (str): Path to the interaction/phenotype matrix CSV.
        output_dir (str): Directory for all results (per-fold under iteration_N/).
        cv_mode (str): 'kfold' (default), 'group' (leave-one-group-out), or
            'predefined' (explicit folds from folds_file; folds may overlap).
        n_folds (int): Number of folds per round, kfold mode (default: 10).
        cv_rounds (int): Number of repeated k-fold rounds, kfold mode (default: 1).
        group_metadata (str): Path to a strain->group CSV (required for cv_mode='group').
        group_column (str): Group/serotype column in group_metadata (required for cv_mode='group').
        folds_file (str): Path to a long-format CSV (fold_label, <strain_column>, role)
            of explicit folds (required for cv_mode='predefined').
        group_ungrouped (str): How to handle strains with no group label in
            cv_mode='group': 'drop' (default, exclude from CV) or 'train' (add to
            every fold's modeling set; never held out).
        strain_column (str): Strain identifier column (in the matrix and group metadata; default: 'strain').
        suffix (str): FASTA file suffix for strain files (default: 'faa').
        min_seq_id (float): Minimum sequence identity for MMseqs2 (default: 0.4).
        coverage (float): Minimum coverage for MMseqs2 (default: 0.8).
        sensitivity (float): MMseqs2 sensitivity (default: 7.5).
        threads (int): Threads per fold (default: 4).
        clear_tmp (bool): Delete per-fold tmp dirs on success (default: False).
        shared_clustering (bool): Cluster once over all genomes and reuse across
            folds; per-fold feature collapse stays training-only (default: False).
        (remaining args are passed through to run_protein_family_workflow /
        assign_predict_workflow; see those functions for details.)

    Returns:
        int: Number of folds successfully aggregated.
    """
    if cv_mode not in ('kfold', 'group', 'predefined'):
        raise ValueError(
            f"cv_mode must be 'kfold', 'group', or 'predefined', got '{cv_mode}'."
        )

    logger.info("=== GenoPHI nested cross-validation ===")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"CV mode: {cv_mode}")

    os.makedirs(output_dir, exist_ok=True)

    full_strain_list = _get_full_strain_list(
        interaction_matrix, input_strain_dir, strain_column, suffix=suffix
    )
    if not full_strain_list:
        raise RuntimeError(
            "No valid strains found (intersection of interaction matrix and "
            "strain directory is empty). Cannot run cross-validation."
        )

    # Build the ordered split plan: a list of (fold_label, modeling, validation).
    if cv_mode == 'group':
        if not group_metadata or not group_column:
            raise ValueError(
                "cv_mode='group' requires both group_metadata and group_column."
            )
        group_map, ungrouped = _load_group_map(
            group_metadata, full_strain_list, strain_column, group_column
        )
        raw_splits = split_strains_by_group(
            full_strain_list, group_map,
            ungrouped=ungrouped, ungrouped_mode=group_ungrouped,
        )
        split_plan = [
            (label, modeling, validation)
            for (label, modeling, validation) in raw_splits
        ]
        if group_ungrouped == 'drop' and ungrouped:
            logger.info(
                f"{len(ungrouped)} ungrouped strain(s) dropped from CV "
                f"(group_ungrouped='drop'). Use group_ungrouped='train' to keep "
                f"them in every fold's training set."
            )
        logger.info(
            f"Leave-one-group-out: {len(split_plan)} fold(s), one per group."
        )
    elif cv_mode == 'predefined':
        if not folds_file:
            raise ValueError("cv_mode='predefined' requires folds_file.")
        split_plan = load_predefined_folds(
            folds_file, full_strain_list, strain_column=strain_column
        )
        logger.info(
            f"Predefined folds: {len(split_plan)} fold(s) from {folds_file} "
            f"(folds may overlap)."
        )
    else:
        if len(full_strain_list) < n_folds:
            raise ValueError(
                f"Cannot create {n_folds} folds from {len(full_strain_list)} strains."
            )
        total_iterations = n_folds * cv_rounds
        split_plan = []
        for iteration in range(1, total_iterations + 1):
            modeling, validation = split_strains_kfold(
                full_strain_list, iteration=iteration, n_folds=n_folds
            )
            split_plan.append((f"fold_{iteration}", modeling, validation))
        logger.info(
            f"K-fold: {n_folds} x {cv_rounds} round(s) = {total_iterations} total fits"
        )

    total_iterations = len(split_plan)

    # Shared clustering: run MMseqs2 once over all genomes; folds reuse it.
    shared_dir = None
    if shared_clustering:
        shared_dir = os.path.join(output_dir, 'shared_clustering')
        logger.info("Shared clustering enabled: clustering once over all genomes.")
        run_shared_clustering(
            input_strain_dir=input_strain_dir,
            input_phage_dir=input_phage_dir,
            output_dir=shared_dir,
            min_seq_id=min_seq_id,
            coverage=coverage,
            sensitivity=sensitivity,
            suffix=suffix,
            threads=threads,
            strain_column=strain_column,
        )

    succeeded, failed = [], []
    for iteration, (fold_label, modeling_strains, validation_strains) in enumerate(split_plan, start=1):
        try:
            # Record which group/fold this iteration holds out (self-documenting,
            # resume-safe).
            iteration_dir = os.path.join(output_dir, f'iteration_{iteration}')
            os.makedirs(iteration_dir, exist_ok=True)
            with open(os.path.join(iteration_dir, 'fold_group.txt'), 'w') as fh:
                fh.write(f"{fold_label}\n")

            if shared_clustering:
                _run_fold_shared_wrapper(
                    iteration=iteration,
                    modeling_strains=modeling_strains,
                    validation_strains=validation_strains,
                    shared_dir=shared_dir,
                    input_strain_dir=input_strain_dir,
                    interaction_matrix=interaction_matrix,
                    output_dir=output_dir,
                    min_seq_id=min_seq_id,
                    coverage=coverage,
                    sensitivity=sensitivity,
                    threads=threads,
                    num_runs_fs=num_runs_fs,
                    num_runs_modeling=num_runs_modeling,
                    use_dynamic_weights=use_dynamic_weights,
                    weights_method=weights_method,
                    use_clustering=use_clustering,
                    cluster_method=cluster_method,
                    n_clusters=n_clusters,
                    min_cluster_size=min_cluster_size,
                    min_samples=min_samples,
                    cluster_selection_epsilon=cluster_selection_epsilon,
                    check_feature_presence=check_feature_presence,
                    filter_by_cluster_presence=filter_by_cluster_presence,
                    min_cluster_presence=min_cluster_presence,
                    max_ram=max_ram,
                    use_feature_clustering=use_feature_clustering,
                    feature_cluster_method=feature_cluster_method,
                    feature_n_clusters=feature_n_clusters,
                    feature_min_cluster_presence=feature_min_cluster_presence,
                    duplicate_all=True,
                )
                succeeded.append(iteration)
                if clear_tmp:
                    _clear_fold_tmp(output_dir, iteration)
                continue

            _run_single_fold(
                iteration=iteration,
                modeling_strains=modeling_strains,
                validation_strains=validation_strains,
                input_strain_dir=input_strain_dir,
                input_phage_dir=input_phage_dir,
                interaction_matrix=interaction_matrix,
                output_dir=output_dir,
                min_seq_id=min_seq_id,
                coverage=coverage,
                sensitivity=sensitivity,
                threads=threads,
                num_runs_fs=num_runs_fs,
                num_runs_modeling=num_runs_modeling,
                use_dynamic_weights=use_dynamic_weights,
                weights_method=weights_method,
                use_clustering=use_clustering,
                cluster_method=cluster_method,
                n_clusters=n_clusters,
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                cluster_selection_epsilon=cluster_selection_epsilon,
                check_feature_presence=check_feature_presence,
                filter_by_cluster_presence=filter_by_cluster_presence,
                min_cluster_presence=min_cluster_presence,
                max_ram=max_ram,
                use_feature_clustering=use_feature_clustering,
                feature_cluster_method=feature_cluster_method,
                feature_n_clusters=feature_n_clusters,
                feature_min_cluster_presence=feature_min_cluster_presence,
                duplicate_all=duplicate_all,
            )
            succeeded.append(iteration)
            if clear_tmp:
                _clear_fold_tmp(output_dir, iteration)
        except Exception as e:  # noqa: BLE001 - keep going across folds
            failed.append(iteration)
            logger.error(
                f"Iteration {iteration} ({fold_label}) failed: {e}", exc_info=True
            )

    logger.info(
        f"Folds complete: {len(succeeded)} succeeded, {len(failed)} failed."
    )
    if failed:
        logger.warning(f"Failed folds: {failed}")

    aggregated = _aggregate_results(output_dir, total_iterations)
    logger.info(f"Aggregated {aggregated} fold(s) into final predictions.")

    if aggregated > 0:
        try:
            _compute_global_metrics(
                output_dir=output_dir,
                interaction_matrix=interaction_matrix,
                strain_column=strain_column,
                phage_column='phage',
                phenotype_column='interaction',
            )
        except Exception as e:  # noqa: BLE001 - metrics are non-fatal
            logger.error(f"Failed to compute global metrics: {e}", exc_info=True)

    return aggregated


def _clear_fold_tmp(output_dir, iteration):
    """Remove the per-fold tmp directories after a successful fold.

    Only call this AFTER the fold's validation prediction has completed. The
    fold's MMseqs2 database (iteration_N/tmp/strain/mmseqs_db) is required to
    assign features to the held-out validation strains during prediction, so it
    must not be removed before _run_single_fold returns. Each fold builds its
    own mmseqs_db and nothing reuses it across folds, so clearing here is safe.
    """
    import shutil

    iteration_output_dir = os.path.join(output_dir, f'iteration_{iteration}')
    for tmp_path in (
        os.path.join(iteration_output_dir, 'tmp'),
        os.path.join(iteration_output_dir, 'model_validation', 'tmp'),
    ):
        if os.path.isdir(tmp_path):
            shutil.rmtree(tmp_path, ignore_errors=True)
            logger.info(f"Cleared tmp directory: {tmp_path}")
