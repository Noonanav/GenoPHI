"""
Unit tests for multi-output task detection and validation (Phase 1).

Covers:
  - genophi.utils.detect_task_type: the (target_mode x strategy) resolution.
  - genophi.utils.validate_phenotype_task_type: DataFrame (per-column) support
    and unchanged Series behavior.
  - genophi.cli.resolve_phenotype_columns: scalar-vs-list target resolution.
"""

import pandas as pd
import pytest

from genophi.utils import (
    detect_task_type,
    validate_phenotype_task_type,
    BINARY_CLASSIFICATION,
    MULTICLASS_CLASSIFICATION,
    MULTILABEL_CLASSIFICATION,
    SINGLE_REGRESSION,
    MULTITARGET_REGRESSION,
)


# --------------------------------------------------------------------------
# detect_task_type: auto-detection
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_auto_binary_classification():
    y = pd.Series([0, 1, 0, 1, 1])
    assert detect_task_type(y, 'classification') == (BINARY_CLASSIFICATION, 'joint')


@pytest.mark.unit
def test_auto_multiclass_classification():
    y = pd.Series([0, 1, 2, 3, 2, 1])
    task, _ = detect_task_type(y, 'classification')
    assert task == MULTICLASS_CLASSIFICATION


@pytest.mark.unit
def test_auto_single_regression():
    y = pd.Series([1.5, 2.3, 4.7, 3.2])
    assert detect_task_type(y, 'regression') == (SINGLE_REGRESSION, 'joint')


@pytest.mark.unit
def test_auto_multilabel_classification():
    y = pd.DataFrame({'res_A': [0, 1, 1], 'res_B': [1, 0, 1]})
    task, _ = detect_task_type(y, 'classification')
    assert task == MULTILABEL_CLASSIFICATION


@pytest.mark.unit
def test_auto_multitarget_regression():
    y = pd.DataFrame({'growth': [1.1, 2.2], 'yield': [3.3, 4.4]})
    task, _ = detect_task_type(y, 'regression')
    assert task == MULTITARGET_REGRESSION


@pytest.mark.unit
def test_single_column_dataframe_treated_as_single_target():
    # A 1-column DataFrame is single-target, not multilabel.
    y = pd.DataFrame({'label': [0, 1, 0, 1]})
    task, _ = detect_task_type(y, 'classification')
    assert task == BINARY_CLASSIFICATION


# --------------------------------------------------------------------------
# detect_task_type: strategy echo
# --------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize('strategy', ['joint', 'independent'])
def test_strategy_echoed(strategy):
    y = pd.DataFrame({'a': [0, 1], 'b': [1, 0]})
    _, resolved = detect_task_type(y, 'classification', strategy=strategy)
    assert resolved == strategy


# --------------------------------------------------------------------------
# detect_task_type: explicit target_mode overrides
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_override_multilabel_requires_multiple_columns():
    y = pd.Series([0, 1, 0])
    with pytest.raises(ValueError, match='requires multiple target columns'):
        detect_task_type(y, 'classification', target_mode='multilabel')


@pytest.mark.unit
def test_override_binary_rejects_multiple_columns():
    y = pd.DataFrame({'a': [0, 1], 'b': [1, 0]})
    with pytest.raises(ValueError, match='expects a single target'):
        detect_task_type(y, 'classification', target_mode='binary')


@pytest.mark.unit
def test_override_task_type_mismatch():
    y = pd.Series([0, 1])
    with pytest.raises(ValueError, match="implies task_type"):
        detect_task_type(y, 'classification', target_mode='multitarget')


@pytest.mark.unit
def test_override_multitarget_ok():
    y = pd.DataFrame({'a': [1.1, 2.2], 'b': [3.3, 4.4]})
    task, _ = detect_task_type(y, 'regression', target_mode='multitarget')
    assert task == MULTITARGET_REGRESSION


# --------------------------------------------------------------------------
# detect_task_type: invalid args
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_invalid_task_type():
    with pytest.raises(ValueError, match='Invalid task_type'):
        detect_task_type(pd.Series([0, 1]), 'clustering')


@pytest.mark.unit
def test_invalid_target_mode():
    with pytest.raises(ValueError, match='Invalid target_mode'):
        detect_task_type(pd.Series([0, 1]), 'classification', target_mode='nonsense')


@pytest.mark.unit
def test_invalid_strategy():
    with pytest.raises(ValueError, match='Invalid strategy'):
        detect_task_type(pd.Series([0, 1]), 'classification', strategy='ensemble')


# --------------------------------------------------------------------------
# validate_phenotype_task_type: DataFrame + Series
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_validate_dataframe_per_column_ok():
    y = pd.DataFrame({'a': [0, 1, 1], 'b': [1, 0, 1]})
    # Should not raise.
    validate_phenotype_task_type(y, 'classification')


@pytest.mark.unit
def test_validate_dataframe_rejects_continuous_for_classification():
    # One column is clearly continuous -> classification should reject it.
    y = pd.DataFrame({
        'ok': [0, 1, 1, 0],
        'bad': [0.13, 0.27, 0.55, 0.91] * 1,
    })
    # Pad 'bad' with enough distinct decimals to trip the continuous check.
    y = pd.DataFrame({
        'ok': list(range(2)) * 6,
        'bad': [i / 13.0 for i in range(12)],
    })
    with pytest.raises(ValueError):
        validate_phenotype_task_type(y, 'classification')


@pytest.mark.unit
def test_validate_series_still_works():
    validate_phenotype_task_type(pd.Series([0, 1, 0, 1]), 'classification')
    validate_phenotype_task_type(pd.Series([1.2, 3.4, 5.6]), 'regression')


@pytest.mark.unit
def test_validate_series_rejects_continuous_for_classification():
    y = pd.Series([i / 13.0 for i in range(20)])
    with pytest.raises(ValueError):
        validate_phenotype_task_type(y, 'classification', 'score')


# --------------------------------------------------------------------------
# cli.resolve_phenotype_columns
# --------------------------------------------------------------------------

@pytest.mark.unit
def test_resolve_single_value_stays_scalar():
    from genophi.cli import create_parser, resolve_phenotype_columns
    ns = create_parser().parse_args(['train', '-i', 'd', '-o', 'o'])
    assert resolve_phenotype_columns(ns) == 'interaction'


@pytest.mark.unit
def test_resolve_comma_separated_becomes_list():
    from genophi.cli import create_parser, resolve_phenotype_columns
    ns = create_parser().parse_args(
        ['select-and-train', '-i', 'x.csv', '-o', 'out',
         '--phenotype_column', 'a,b,c']
    )
    assert resolve_phenotype_columns(ns) == ['a', 'b', 'c']


@pytest.mark.unit
def test_resolve_single_comma_value_collapses_to_scalar():
    from genophi.cli import create_parser, resolve_phenotype_columns
    ns = create_parser().parse_args(
        ['select-and-train', '-i', 'x.csv', '-o', 'out',
         '--phenotype_column', 'only_one']
    )
    assert resolve_phenotype_columns(ns) == 'only_one'


@pytest.mark.unit
def test_resolve_columns_file(tmp_path):
    from genophi.cli import create_parser, resolve_phenotype_columns
    f = tmp_path / 'targets.txt'
    f.write_text('res_A\nres_B\nres_C\n')
    ns = create_parser().parse_args(
        ['select-and-train', '-i', 'x.csv', '-o', 'out',
         '--phenotype_columns_file', str(f)]
    )
    assert resolve_phenotype_columns(ns) == ['res_A', 'res_B', 'res_C']


@pytest.mark.unit
def test_resolve_columns_file_missing_raises(tmp_path):
    from genophi.cli import create_parser, resolve_phenotype_columns
    ns = create_parser().parse_args(
        ['select-and-train', '-i', 'x.csv', '-o', 'out',
         '--phenotype_columns_file', str(tmp_path / 'nope.txt')]
    )
    with pytest.raises(FileNotFoundError):
        resolve_phenotype_columns(ns)


@pytest.mark.unit
def test_strategy_and_target_mode_defaults():
    from genophi.cli import create_parser
    ns = create_parser().parse_args(['train', '-i', 'd', '-o', 'o'])
    assert ns.strategy == 'joint'
    assert ns.target_mode == 'auto'
