"""
Smoke tests for GenoPHI package.

Quick tests to verify package installation and basic functionality.
These tests should run in seconds and catch major installation issues.
"""

import pytest
import subprocess
import sys


@pytest.mark.smoke
def test_package_import():
    """Test that the main package can be imported."""
    import genophi
    assert genophi is not None


@pytest.mark.smoke
def test_workflow_imports():
    """Test that workflow modules can be imported."""
    from genophi.workflows import protein_family_workflow
    from genophi.workflows import kmer_full_workflow
    from genophi.workflows import full_workflow
    from genophi.workflows import nested_cv_workflow

    assert hasattr(protein_family_workflow, 'run_protein_family_workflow')
    assert hasattr(kmer_full_workflow, 'run_kmer_workflow')
    assert hasattr(full_workflow, 'run_full_workflow')
    assert hasattr(nested_cv_workflow, 'run_nested_cv_workflow')


@pytest.mark.smoke
def test_core_module_imports():
    """Test that core processing modules can be imported."""
    from genophi import mmseqs2_clustering
    from genophi import feature_selection
    from genophi import select_feature_modeling
    from genophi import feature_annotations
    from genophi import utils

    assert hasattr(mmseqs2_clustering, 'run_clustering_workflow')
    assert hasattr(feature_selection, 'run_feature_selection_iterations')
    assert hasattr(select_feature_modeling, 'model_testing_select_MCC')
    assert hasattr(feature_annotations, 'get_predictive_features')
    assert hasattr(utils, 'validate_phenotype_task_type')


@pytest.mark.smoke
def test_cli_module_import():
    """Test that CLI module can be imported."""
    from genophi import cli
    assert hasattr(cli, 'main')


@pytest.mark.smoke
def test_cli_executable():
    """Test that genophi CLI is available as executable."""
    result = subprocess.run(
        [sys.executable, '-m', 'genophi.cli', '--help'],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    assert 'genophi' in result.stdout.lower() or 'usage' in result.stdout.lower()


@pytest.mark.smoke
def test_required_dependencies():
    """Test that all required dependencies are importable."""
    dependencies = [
        'pandas',
        'numpy',
        'sklearn',
        'Bio',  # biopython
        'catboost',
        'matplotlib',
        'seaborn',
        'tqdm',
        'joblib',
        'plotnine',
        'shap',
        'psutil',
        'scipy',
        'hdbscan',
    ]

    missing = []
    for dep in dependencies:
        try:
            __import__(dep)
        except ImportError:
            missing.append(dep)

    assert not missing, f"Missing required dependencies: {missing}"


@pytest.mark.smoke
def test_validation_functions():
    """Test that validation utilities work with simple data."""
    from genophi.utils import validate_phenotype_task_type
    import pandas as pd

    # Test classification validation (should pass)
    y_classification = pd.Series([0, 1, 0, 1, 1, 0])
    validate_phenotype_task_type(y_classification, 'classification')

    # Test regression validation (should pass)
    y_regression = pd.Series([1.5, 2.3, 4.7, 3.2, 5.1])
    validate_phenotype_task_type(y_regression, 'regression')

    # Test that invalid task_type raises error
    with pytest.raises(ValueError, match="Invalid task_type"):
        validate_phenotype_task_type(y_classification, 'invalid')


@pytest.mark.smoke
def test_kfold_split_properties():
    """Test deterministic k-fold strain splitting: determinism, disjointness, coverage."""
    from genophi.workflows.nested_cv_workflow import split_strains_kfold

    strains = [f"strain_{i}" for i in range(23)]  # not divisible by n_folds
    n_folds = 5

    # Determinism: same iteration -> identical split regardless of input order.
    m1, v1 = split_strains_kfold(strains, iteration=1, n_folds=n_folds)
    m2, v2 = split_strains_kfold(list(reversed(strains)), iteration=1, n_folds=n_folds)
    assert sorted(m1) == sorted(m2)
    assert sorted(v1) == sorted(v2)

    # For one full round: every strain is validated exactly once, folds are
    # disjoint, and the union covers all strains.
    all_validation = []
    for iteration in range(1, n_folds + 1):
        modeling, validation = split_strains_kfold(strains, iteration=iteration, n_folds=n_folds)
        # modeling and validation partition the full set with no overlap.
        assert set(modeling).isdisjoint(set(validation))
        assert set(modeling) | set(validation) == set(strains)
        all_validation.extend(validation)

    assert sorted(all_validation) == sorted(strains)  # each strain held out once

    # A second round (iterations n_folds+1..2*n_folds) reshuffles: validation
    # assignment should differ from round 0 for at least one fold.
    _, v_round0_fold0 = split_strains_kfold(strains, iteration=1, n_folds=n_folds)
    _, v_round1_fold0 = split_strains_kfold(strains, iteration=n_folds + 1, n_folds=n_folds)
    assert sorted(v_round0_fold0) != sorted(v_round1_fold0)


@pytest.mark.smoke
def test_global_metrics_and_curves(tmp_path):
    """Test pooled outer-test metrics + PR/ROC curve generation from synthetic data."""
    from genophi.workflows.nested_cv_workflow import _compute_global_metrics
    import pandas as pd

    output_dir = tmp_path / "cv_run"
    output_dir.mkdir()

    # Pooled per-fold predictions: strain, phage, Confidence, Final_Prediction.
    # Scores are well-separated by the true label so AUC is high and finite.
    pd.DataFrame({
        'strain': ['s1', 's2', 's3', 's4', 's5', 's6'],
        'phage': ['p1'] * 6,
        'Confidence': [0.9, 0.8, 0.7, 0.2, 0.1, 0.3],
        'Final_Prediction': [1, 1, 1, 0, 0, 0],
    }).to_csv(output_dir / 'final_predictions.csv', index=False)

    # Ground-truth interaction matrix with the 'interaction' label column.
    interaction_matrix = tmp_path / "matrix.csv"
    pd.DataFrame({
        'strain': ['s1', 's2', 's3', 's4', 's5', 's6'],
        'phage': ['p1'] * 6,
        'interaction': [1, 1, 1, 0, 0, 0],
    }).to_csv(interaction_matrix, index=False)

    metrics = _compute_global_metrics(str(output_dir), str(interaction_matrix))

    assert metrics is not None
    assert metrics['n_predictions'] == 6
    # Perfectly separable labels -> AUC == 1.0, MCC == 1.0.
    assert metrics['AUC'] == 1.0
    assert metrics['MCC'] == 1.0

    perf_dir = output_dir / 'performance'
    assert (perf_dir / 'global_metrics.csv').exists()
    assert (perf_dir / 'roc_curve.png').exists()
    assert (perf_dir / 'pr_curve.png').exists()


@pytest.mark.smoke
def test_feature_annotation_functions():
    """Test that feature annotation functions work with simple data."""
    from genophi.feature_annotations import get_predictive_features
    import pandas as pd
    import tempfile

    # Create simple test feature table
    df = pd.DataFrame({
        'strain': ['s1', 's2', 's3'],
        'interaction': [0, 1, 0],
        'sc_0': [1, 0, 1],
        'sc_1': [0, 1, 0],
        'pc_0': [1, 1, 0],
    })

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        temp_path = f.name
        df.to_csv(temp_path, index=False)

    try:
        # Test with strain features only (default)
        strain_features = get_predictive_features(
            temp_path,
            sample_column='strain',
            phenotype_column='interaction',
            feature_type='strain'
        )

        assert 'sc_0' in strain_features
        assert 'sc_1' in strain_features
        assert 'pc_0' not in strain_features  # pc features filtered out when feature_type='strain'

        # Test with all features
        all_features = get_predictive_features(
            temp_path,
            sample_column='strain',
            phenotype_column='interaction',
            feature_type='all'
        )

        assert 'sc_0' in all_features
        assert 'sc_1' in all_features
        assert 'pc_0' in all_features
        assert 'strain' not in all_features
        assert 'interaction' not in all_features
    finally:
        # Clean up temp file
        import os
        if os.path.exists(temp_path):
            os.unlink(temp_path)
