"""
End-to-end test for the nested k-fold cross-validation workflow.

Exercises the full per-fold pipeline (clustering -> feature selection ->
modeling -> held-out prediction), aggregation across folds, and the pooled
outer-test performance metrics / PR-ROC curve generation.
"""

import pytest
import pandas as pd
from pathlib import Path

from genophi.workflows.nested_cv_workflow import run_nested_cv_workflow


@pytest.mark.e2e
@pytest.mark.slow
@pytest.mark.requires_mmseqs2
def test_nested_cv_kfold_end_to_end(medium_test_dataset, temp_output_dir):
    """Run a small repeated-free 5-fold CV and validate all expected outputs.

    Uses reduced run counts so the test stays within E2E time budgets while
    still exercising every stage on real data.
    """
    strain_dir = medium_test_dataset['strain_dir']
    phage_dir = medium_test_dataset['phage_dir']
    interaction_matrix = medium_test_dataset['interaction_matrix']
    output_dir = temp_output_dir / 'nested_cv_output'

    n_folds = 5
    cv_rounds = 1

    aggregated = run_nested_cv_workflow(
        input_strain_dir=str(strain_dir),
        input_phage_dir=str(phage_dir),
        interaction_matrix=str(interaction_matrix),
        output_dir=str(output_dir),
        n_folds=n_folds,
        cv_rounds=cv_rounds,
        num_runs_fs=5,
        num_runs_modeling=5,
        min_seq_id=0.4,
        coverage=0.8,
        sensitivity=7.5,
        threads=2,
    )

    total_iterations = n_folds * cv_rounds

    # At least one fold must complete for the run to be meaningful.
    assert aggregated > 0, "No folds produced predictions"

    # Per-fold output directories and held-out predictions exist for the
    # folds that succeeded.
    completed_folds = 0
    for i in range(1, total_iterations + 1):
        iteration_dir = output_dir / f'iteration_{i}'
        median_pred = (
            iteration_dir / 'model_validation' / 'predict_results'
            / 'strain_median_predictions.csv'
        )
        if median_pred.exists():
            completed_folds += 1
            # Each fold's split files should also be present.
            assert (iteration_dir / 'modeling_strains.csv').exists()
            assert (iteration_dir / 'validation_strains.csv').exists()
    assert completed_folds == aggregated

    # Aggregated predictions.
    final_predictions = output_dir / 'final_predictions.csv'
    assert final_predictions.exists(), "final_predictions.csv not written"
    final_df = pd.read_csv(final_predictions)
    assert len(final_df) > 0
    assert 'iteration' in final_df.columns
    assert {'strain', 'phage'}.issubset(final_df.columns)

    summary = output_dir / 'prediction_summary.csv'
    assert summary.exists(), "prediction_summary.csv not written"

    # Global performance metrics and pooled outer-test curves.
    perf_dir = output_dir / 'performance'
    metrics_csv = perf_dir / 'global_metrics.csv'
    assert metrics_csv.exists(), "global_metrics.csv not written"

    metrics = pd.read_csv(metrics_csv)
    assert len(metrics) == 1
    for col in ('AUC', 'AP', 'Accuracy', 'Precision', 'Recall', 'F1', 'MCC', 'n_predictions'):
        assert col in metrics.columns, f"Missing metric column: {col}"
    assert metrics['n_predictions'].iloc[0] > 0

    # Curves are produced when both classes are present in the pooled labels.
    # On real interaction data this is expected; assert when applicable.
    assert (perf_dir / 'roc_curve.png').exists(), "roc_curve.png not written"
    assert (perf_dir / 'pr_curve.png').exists(), "pr_curve.png not written"
