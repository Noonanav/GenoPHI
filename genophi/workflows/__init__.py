"""
GenoPHI Workflow Modules

High-level workflow functions for running complete analysis pipelines.
"""

# Importing workflows from their respective modules
from .feature_table_workflow import run_full_feature_workflow
from .feature_selection_workflow import run_feature_selection_workflow
from .modeling_workflow import run_modeling_workflow
from .full_workflow import run_full_workflow
from .assign_features_workflow import run_assign_features_workflow
from .prediction_workflow import run_prediction_workflow
from .assign_predict_workflow import assign_predict_workflow
from .select_and_model_workflow import run_modeling_workflow_from_feature_table
from .feature_annotations_workflow import run_predictive_proteins_workflow
from .kmer_assign_features_workflow import run_kmer_assign_features_workflow
from .kmer_assign_predict_workflow import kmer_assign_predict_workflow
from .kmer_full_workflow import run_kmer_workflow

__all__ = [
    'run_full_feature_workflow',
    'run_feature_selection_workflow',
    'run_modeling_workflow',
    'run_full_workflow',
    'run_assign_features_workflow',
    'run_prediction_workflow',
    'assign_predict_workflow',
    'run_modeling_workflow_from_feature_table',
    'run_predictive_proteins_workflow',
    'run_kmer_assign_features_workflow',
    'kmer_assign_predict_workflow',
    'run_kmer_workflow',
]