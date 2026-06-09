"""
Utility functions for GenoPHI package.

This module provides general utility and validation functions to ensure data
integrity and compatibility across different modeling tasks.
"""

import os
import logging
import pandas as pd
import numpy as np


def validate_phenotype_task_type(y, task_type, phenotype_column='interaction'):
    """
    Validates that the phenotype data type matches the specified task type.

    This function checks whether the phenotype data is appropriate for the
    specified task (classification or regression) and raises an error if there's
    a mismatch. This prevents users from accidentally running classification
    on continuous data or vice versa.

    Accepts either a single-column ``pd.Series`` (single-target) or a
    ``pd.DataFrame`` (multi-target). For a DataFrame, each column is validated
    independently with the same Series logic.

    Args:
        y (pd.Series or pd.DataFrame): The phenotype/target variable(s) to validate.
        task_type (str): The task type - must be 'classification' or 'regression'.
        phenotype_column (str): Name of the phenotype column for error messages.
            Default is 'interaction'. Ignored for DataFrame input (column names
            are used instead).

    Raises:
        ValueError: If task_type is invalid or if the phenotype data type doesn't
            match the task type.

    Examples:
        >>> import pandas as pd
        >>> # Continuous data - should work with regression
        >>> y = pd.Series([1.5, 2.3, 4.7, 3.2])
        >>> validate_phenotype_task_type(y, 'regression', 'score')

        >>> # Binary data - should work with classification
        >>> y = pd.Series([0, 1, 0, 1, 1])
        >>> validate_phenotype_task_type(y, 'classification', 'label')

        >>> # Multi-target: validate each column independently
        >>> y = pd.DataFrame({'res_A': [0, 1, 1], 'res_B': [1, 0, 1]})
        >>> validate_phenotype_task_type(y, 'classification')

        >>> # Continuous data with classification - should raise error
        >>> y = pd.Series([1.5, 2.3, 4.7, 3.2])
        >>> validate_phenotype_task_type(y, 'classification', 'score')
        Traceback (most recent call last):
            ...
        ValueError: Task type is set to 'classification', but the phenotype column...
    """
    # Validate task_type parameter (do this once, before per-column dispatch)
    if task_type not in ['classification', 'regression']:
        raise ValueError(
            f"Invalid task_type '{task_type}'. Must be 'classification' or 'regression'."
        )

    # Multi-target: validate each column independently with the Series logic.
    if isinstance(y, pd.DataFrame):
        for column in y.columns:
            _validate_phenotype_series(y[column], task_type, str(column))
        logging.info(
            f"Validation passed: {len(y.columns)} phenotype columns suitable for "
            f"{task_type} ({list(y.columns)})"
        )
        return

    _validate_phenotype_series(y, task_type, phenotype_column)


def _validate_phenotype_series(y, task_type, phenotype_column='interaction'):
    """
    Validate a single phenotype Series against a task type.

    Contains the original single-column validation logic, factored out so that
    both Series and per-column DataFrame validation share identical behavior.
    Assumes ``task_type`` has already been validated by the caller.
    """
    # Get basic statistics about the data
    unique_values = y.nunique()

    # Detect data type characteristics
    is_float_type = pd.api.types.is_float_dtype(y)
    is_integer_type = pd.api.types.is_integer_dtype(y)
    is_numeric = pd.api.types.is_numeric_dtype(y)

    if task_type == 'classification':
        # Classification should have discrete categorical values
        # Check 1: If data is float type with many unique values, likely continuous
        if is_float_type and unique_values > 10:
            raise ValueError(
                f"Task type is set to 'classification', but the phenotype column '{phenotype_column}' "
                f"contains continuous data with {unique_values} unique values.\n"
                f"Data type: {y.dtype}\n"
                f"Sample values: {y.head().tolist()}\n\n"
                f"For continuous/numeric phenotype data, please use:\n"
                f"  --task_type regression\n\n"
                f"If this is truly categorical data, please encode it as integers (0, 1, 2, etc.)."
            )

        # Check 2: If float values have decimal places, likely continuous
        if is_float_type:
            # Sample first 100 non-null values to check for decimals
            sample_values = y.dropna().head(100)
            if len(sample_values) > 0:
                has_decimals = any(val % 1 != 0 for val in sample_values if pd.notna(val))
                if has_decimals:
                    raise ValueError(
                        f"Task type is set to 'classification', but the phenotype column '{phenotype_column}' "
                        f"contains float values with decimal places.\n"
                        f"Sample values: {y.dropna().head(5).tolist()}\n\n"
                        f"Classification requires discrete categories (e.g., 0, 1, 2).\n"
                        f"For continuous/numeric phenotype data, please use:\n"
                        f"  --task_type regression"
                    )

        # Check 3: Warn if there are many unique values even for integer data
        if is_integer_type and unique_values > 20:
            logging.warning(
                f"Phenotype column '{phenotype_column}' has {unique_values} unique integer values "
                f"for classification. This is unusual for categorical data. "
                f"Please verify this is the correct task type."
            )

        logging.info(
            f"Validation passed: Phenotype data is suitable for classification "
            f"(unique classes: {unique_values}, dtype: {y.dtype})"
        )

    elif task_type == 'regression':
        # For regression, we expect continuous data or many unique values
        # Provide a warning if data looks categorical/binary
        if is_numeric and unique_values <= 2:
            unique_vals_list = sorted(y.unique())
            logging.warning(
                f"Task type is set to 'regression', but the phenotype column '{phenotype_column}' "
                f"appears to contain binary/categorical data with only {unique_values} unique value(s): {unique_vals_list}.\n"
                f"If this is a classification problem (predicting categories), consider using:\n"
                f"  --task_type classification"
            )

        # Warn if data looks like multi-class categorical
        if is_integer_type and 3 <= unique_values <= 10:
            logging.warning(
                f"Task type is set to 'regression', but the phenotype column '{phenotype_column}' "
                f"contains {unique_values} unique integer values, which may indicate categorical data. "
                f"If predicting categories, consider using --task_type classification."
            )

        logging.info(
            f"Validation passed: Phenotype data is suitable for regression "
            f"(unique values: {unique_values}, dtype: {y.dtype})"
        )


# Resolved task identifiers returned by detect_task_type().
BINARY_CLASSIFICATION = 'binary_classification'
MULTICLASS_CLASSIFICATION = 'multiclass_classification'
MULTILABEL_CLASSIFICATION = 'multilabel_classification'
SINGLE_REGRESSION = 'single_regression'
MULTITARGET_REGRESSION = 'multitarget_regression'

VALID_TARGET_MODES = (
    'auto', 'binary', 'multiclass', 'multilabel', 'single', 'multitarget'
)
VALID_STRATEGIES = ('joint', 'independent')


def detect_task_type(y, task_type, target_mode='auto', strategy='joint'):
    """
    Resolve the specific modeling task from the target data and user overrides.

    This determines *what* the targets are (binary/multiclass/multilabel/
    single/multitarget regression). It is orthogonal to ``strategy``, which
    determines *how* multiple targets are modeled ('joint' single multi-output
    model vs. 'independent' one model per target). ``strategy`` is echoed back
    unchanged for single-target tasks (where it is irrelevant).

    Args:
        y (pd.Series or pd.DataFrame): Target variable(s). A Series is treated
            as single-target; a DataFrame as multi-target.
        task_type (str): 'classification' or 'regression'.
        target_mode (str): One of 'auto', 'binary', 'multiclass', 'multilabel',
            'single', 'multitarget'. 'auto' (default) infers from the data.
        strategy (str): 'joint' or 'independent'. Echoed back; only meaningful
            for multi-target tasks.

    Returns:
        tuple[str, str]: ``(specific_task, strategy)`` where ``specific_task`` is
        one of the module-level ``*_CLASSIFICATION`` / ``*_REGRESSION`` constants.

    Raises:
        ValueError: If task_type, target_mode, or strategy is invalid, or if an
            explicit target_mode is incompatible with the shape of ``y``.
    """
    if task_type not in ('classification', 'regression'):
        raise ValueError(
            f"Invalid task_type '{task_type}'. Must be 'classification' or 'regression'."
        )
    if target_mode not in VALID_TARGET_MODES:
        raise ValueError(
            f"Invalid target_mode '{target_mode}'. Must be one of {VALID_TARGET_MODES}."
        )
    if strategy not in VALID_STRATEGIES:
        raise ValueError(
            f"Invalid strategy '{strategy}'. Must be one of {VALID_STRATEGIES}."
        )

    is_multi = isinstance(y, pd.DataFrame) and y.shape[1] > 1

    # --- Explicit target_mode override (skip auto-detection) ---
    if target_mode != 'auto':
        explicit = {
            'binary': (BINARY_CLASSIFICATION, 'classification', False),
            'multiclass': (MULTICLASS_CLASSIFICATION, 'classification', False),
            'multilabel': (MULTILABEL_CLASSIFICATION, 'classification', True),
            'single': (SINGLE_REGRESSION, 'regression', False),
            'multitarget': (MULTITARGET_REGRESSION, 'regression', True),
        }
        specific_task, expected_task, expects_multi = explicit[target_mode]
        if task_type != expected_task:
            raise ValueError(
                f"target_mode '{target_mode}' implies task_type '{expected_task}', "
                f"but task_type='{task_type}' was given."
            )
        if expects_multi and not is_multi:
            raise ValueError(
                f"target_mode '{target_mode}' requires multiple target columns, "
                f"but a single target was provided."
            )
        if not expects_multi and is_multi:
            raise ValueError(
                f"target_mode '{target_mode}' expects a single target, but "
                f"multiple target columns were provided."
            )
        return specific_task, strategy

    # --- Auto-detection ---
    if is_multi:
        # Multiple columns: multilabel (classification) or multitarget (regression).
        specific_task = (
            MULTILABEL_CLASSIFICATION if task_type == 'classification'
            else MULTITARGET_REGRESSION
        )
        return specific_task, strategy

    # Single target: collapse a 1-column DataFrame to a Series for inspection.
    y_series = y.iloc[:, 0] if isinstance(y, pd.DataFrame) else y

    if task_type == 'regression':
        return SINGLE_REGRESSION, strategy

    # Single-column classification: binary vs. multiclass by unique count.
    n_unique = y_series.nunique()
    specific_task = (
        BINARY_CLASSIFICATION if n_unique <= 2 else MULTICLASS_CLASSIFICATION
    )
    return specific_task, strategy


def validate_file(path, name):
    """
    Validate that a file exists and is accessible.

    Args:
        path (str): Path to the file to validate.
        name (str): Descriptive name of the file for error messages.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the path is not a file.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found: {path}")
    if not os.path.isfile(path):
        raise ValueError(f"{name} is not a file: {path}")


def validate_directory(path, name, create=False):
    """
    Validate or create a directory.

    Args:
        path (str): Path to the directory to validate.
        name (str): Descriptive name of the directory for error messages.
        create (bool): If True, create the directory if it doesn't exist.
            Default is False.

    Raises:
        FileNotFoundError: If the directory does not exist and create is False.
        ValueError: If the path exists but is not a directory.
    """
    if create:
        os.makedirs(path, exist_ok=True)
    elif not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found: {path}")
    elif not os.path.isdir(path):
        raise ValueError(f"{name} is not a directory: {path}")
