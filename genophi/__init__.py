"""
GenoPHI: Genome-to-Phenotype Phage-Host Interaction Prediction

A Python package for whole-genome genotype-to-phenotype modeling with a focus
on phage-host interaction prediction.
"""

__version__ = "0.1.0"
__author__ = "Avery Noonan"
__license__ = "MIT"
__url__ = "https://github.com/Noonanav/genophi"

# Core clustering and feature extraction functions
from .mmseqs2_clustering import (
    run_clustering_workflow,
    run_feature_assignment,
    merge_feature_tables,
    load_strains,
    create_mmseqs_database,
    create_contig_to_genome_dict,
    run_mmseqs_cluster,
    assign_sequences_to_clusters,
    generate_presence_absence_matrix,
    compare_cluster_and_search_results,
    filter_presence_absence,
    get_genome_assignments_tables,
    feature_selection_optimized,
    feature_assignment
)

# Feature selection functions
from .feature_selection import (
    perform_rfe,
    grid_search,
    load_and_prepare_data,
    filter_data,
    run_feature_selection_iterations,
    generate_feature_tables
)

# Modeling functions
from .select_feature_modeling import run_experiments

# Annotation functions
from .feature_annotations import (
    get_predictive_proteins,
    merge_annotation_table,
    parse_and_filter_aa_sequences
)

__all__ = [
    # Metadata
    '__version__',
    '__author__',
    '__license__',
    # Clustering and feature extraction
    'run_clustering_workflow',
    'run_feature_assignment',
    'merge_feature_tables',
    'load_strains',
    'create_mmseqs_database',
    'create_contig_to_genome_dict',
    'run_mmseqs_cluster',
    'assign_sequences_to_clusters',
    'generate_presence_absence_matrix',
    'compare_cluster_and_search_results',
    'filter_presence_absence',
    'get_genome_assignments_tables',
    'feature_selection_optimized',
    'feature_assignment',
    # Feature selection
    'perform_rfe',
    'grid_search',
    'load_and_prepare_data',
    'filter_data',
    'run_feature_selection_iterations',
    'generate_feature_tables',
    # Modeling
    'run_experiments',
    # Annotations
    'get_predictive_proteins',
    'merge_annotation_table',
    'parse_and_filter_aa_sequences'
]