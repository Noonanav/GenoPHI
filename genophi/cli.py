#!/usr/bin/env python
"""
GenoPHI Command-Line Interface

Unified entry point for all GenoPHI workflows.
"""

import sys
import os
import argparse
import logging
from pathlib import Path


def setup_logging(verbose=False):
    """Configure logging for the CLI."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_file(path, name):
    """Validate that a file exists."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found: {path}")
    if not os.path.isfile(path):
        raise ValueError(f"{name} is not a file: {path}")


def validate_directory(path, name, create=False):
    """Validate or create a directory."""
    if create:
        os.makedirs(path, exist_ok=True)
    elif not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found: {path}")
    elif not os.path.isdir(path):
        raise ValueError(f"{name} is not a directory: {path}")


def add_common_args(parser):
    """Add common arguments shared across commands."""
    parser.add_argument('--threads', type=int, default=4, 
                       help='Number of threads to use (default: 4)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')


def create_parser():
    """Create the argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog='genophi',
        description='GenoPHI: Genome-to-Phenotype Interaction Prediction',
        epilog='Use genophi <command> --help for detailed information on each command.'
    )
    
    parser.add_argument('--version', action='version', version='genophi 0.1.0')
    subparsers = parser.add_subparsers(dest='command', required=True,
                                       help='Available commands')
    
    # ==================== PROTEIN FAMILY WORKFLOWS ====================
    
    # CLUSTER - Protein family clustering and feature generation
    p = subparsers.add_parser(
        'cluster',
        help='Generate protein family clusters and feature tables',
        description='Run MMseqs2 clustering and generate feature tables for strains and optionally phages.'
    )
    p.add_argument('--input_strain', '-is', required=True,
                   help='Directory containing strain FASTA files')
    p.add_argument('--phenotype_matrix', '-pm', required=True,
                   help='Path to phenotype matrix CSV file')
    p.add_argument('--output', '-o', required=True,
                   help='Output directory for results')
    p.add_argument('--input_phage', '-ip',
                   help='Directory containing phage FASTA files (optional)')
    p.add_argument('--tmp', default='tmp',
                   help='Temporary directory for intermediate files (default: tmp)')
    p.add_argument('--min_seq_id', type=float, default=0.4,
                   help='Minimum sequence identity for clustering (default: 0.4)')
    p.add_argument('--coverage', type=float, default=0.8,
                   help='Minimum coverage for clustering (default: 0.8)')
    p.add_argument('--sensitivity', type=float, default=7.5,
                   help='MMseqs2 sensitivity parameter (default: 7.5)')
    p.add_argument('--suffix', default='faa',
                   help='File suffix for FASTA files (default: faa)')
    p.add_argument('--strain_list',
                   help='Path to file with list of strains to include')
    p.add_argument('--strain_column', default='strain',
                   help='Column name for strain identifiers (default: strain)')
    p.add_argument('--phage_list',
                   help='Path to file with list of phages to include')
    p.add_argument('--phage_column', default='phage',
                   help='Column name for phage identifiers (default: phage)')
    p.add_argument('--source_strain', default='strain',
                   help='Prefix for strain features (default: strain)')
    p.add_argument('--source_phage', default='phage',
                   help='Prefix for phage features (default: phage)')
    p.add_argument('--compare', action='store_true',
                   help='Compare original clusters with assigned clusters')
    p.add_argument('--max_ram', type=float, default=16,
                   help='Maximum RAM usage in GB (default: 16)')
    p.add_argument('--use_feature_clustering', action='store_true',
                   help='Enable pre-processing cluster-based feature filtering')
    p.add_argument('--feature_cluster_method', default='hierarchical',
                   choices=['hierarchical'],
                   help='Pre-processing clustering method (default: hierarchical)')
    p.add_argument('--feature_n_clusters', type=int, default=20,
                   help='Number of clusters for pre-processing (default: 20)')
    p.add_argument('--feature_min_cluster_presence', type=int, default=2,
                   help='Min clusters a feature must appear in (default: 2)')
    add_common_args(p)
    
    # SELECT-FEATURES - Feature selection
    p = subparsers.add_parser(
        'select-features',
        help='Perform feature selection on feature table',
        description='Run iterative feature selection to identify important features.'
    )
    p.add_argument('--input', '-i', required=True,
                   help='Path to input feature table CSV')
    p.add_argument('--output', '-o', required=True,
                   help='Output directory for results')
    p.add_argument('--task_type', default='classification', choices=['classification', 'regression'],
                   help='Type of modeling task')
    p.add_argument('--method', default='rfe',
                   choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                   help='Feature selection method (default: rfe)')
    p.add_argument('--num_features', default='none',
                   help='Number of features to select or "none" for automatic (default: none)')
    p.add_argument('--num_runs', type=int, default=50,
                   help='Number of feature selection iterations (default: 50)')
    p.add_argument('--filter_type', default='none', choices=['none', 'strain', 'phage'],
                   help='Filter type for input data (default: none)')
    p.add_argument('--phenotype_column',
                   help='Column name for phenotype/target variable')
    p.add_argument('--sample_column', default='strain',
                   help='Column name for sample identifiers (default: strain)')
    p.add_argument('--phage_column', default='phage',
                   help='Column name for phage identifiers (default: phage)')
    p.add_argument('--binary_data', action='store_true',
                   help='Feature values are binary (0/1)')
    p.add_argument('--max_features', default='none',
                   help='Maximum number of features to include (default: none)')
    p.add_argument('--use_dynamic_weights', action='store_true',
                   help='Use dynamic feature weights')
    p.add_argument('--weights_method', default='inverse_frequency',
                   choices=['log10', 'inverse_frequency', 'balanced'],
                   help='Method for calculating weights (default: inverse_frequency)')
    p.add_argument('--no-clustering', dest='use_clustering', action='store_false', 
               default=True, help='Disable clustering (default: enabled)')
    p.add_argument('--cluster_method', default='hierarchical',
                   choices=['hdbscan', 'hierarchical'],
                   help='Clustering method (default: hierarchical)')
    p.add_argument('--n_clusters', type=int, default=20,
                   help='Number of clusters for hierarchical clustering (default: 20)')
    p.add_argument('--min_cluster_size', type=int, default=5,
                   help='Minimum cluster size for HDBSCAN (default: 5)')
    p.add_argument('--min_samples', type=int,
                   help='Minimum samples for HDBSCAN (default: same as min_cluster_size)')
    p.add_argument('--cluster_selection_epsilon', type=float, default=0.0,
                   help='Epsilon for HDBSCAN clustering (default: 0.0)')
    p.add_argument('--check_feature_presence', action='store_true',
                   help='Check feature presence during selection')
    p.add_argument('--filter_by_cluster_presence', action='store_true',
                   help='Filter features by cluster presence')
    p.add_argument('--min_cluster_presence', type=int, default=2,
                   help='Minimum cluster presence for filtering (default: 2)')
    p.add_argument('--max_ram', type=float, default=16,
                   help='Maximum RAM usage in GB (default: 16)')
    add_common_args(p)
    
    # TRAIN - Model training
    p = subparsers.add_parser(
        'train',
        help='Train predictive models',
        description='Train machine learning models on selected feature tables.'
    )
    p.add_argument('--input', '-i', required=True,
                   help='Directory containing feature table(s)')
    p.add_argument('--output', '-o', required=True,
                   help='Output directory for model results')
    p.add_argument('--num_runs', type=int, default=50,
                   help='Number of training runs per feature table (default: 50)')
    p.add_argument('--task_type', default='classification',
                   choices=['classification', 'regression'],
                   help='Type of modeling task (default: classification)')
    p.add_argument('--set_filter', default='strain', choices=['none', 'strain', 'phage'],
                   help='Dataset filter (default: strain)')
    p.add_argument('--sample_column', default='strain',
                   help='Column name for sample identifiers (default: strain)')
    p.add_argument('--phenotype_column', default='interaction',
                   help='Column name for phenotype (default: interaction)')
    p.add_argument('--phage_column', default='phage',
                   help='Column name for phage identifiers (default: phage)')
    p.add_argument('--binary_data', action='store_true',
                   help='Use binary data for SHAP plots')
    p.add_argument('--use_dynamic_weights', action='store_true',
                   help='Use dynamic feature weights')
    p.add_argument('--weights_method', default='inverse_frequency',
                   choices=['log10', 'inverse_frequency', 'balanced'],
                   help='Method for calculating weights (default: inverse_frequency)')
    p.add_argument('--no-clustering', dest='use_clustering', action='store_false', 
                   default=True, help='Disable clustering (default: enabled)')
    p.add_argument('--cluster_method', default='hierarchical',
                   choices=['hdbscan', 'hierarchical'],
                   help='Clustering method (default: hierarchical)')
    p.add_argument('--n_clusters', type=int, default=20,
                   help='Number of clusters (default: 20)')
    p.add_argument('--min_cluster_size', type=int, default=5,
                   help='Minimum cluster size (default: 5)')
    p.add_argument('--min_samples', type=int,
                   help='Minimum samples for HDBSCAN')
    p.add_argument('--cluster_selection_epsilon', type=float, default=0.0,
                   help='Epsilon for HDBSCAN (default: 0.0)')
    add_common_args(p)
    
    # PREDICT - Make predictions
    p = subparsers.add_parser(
        'predict',
        help='Make predictions using trained models',
        description='Generate predictions for new genome combinations.'
    )
    p.add_argument('--input_dir', required=True,
                   help='Directory with strain-specific feature tables')
    p.add_argument('--model_dir', required=True,
                   help='Directory containing trained models')
    p.add_argument('--output_dir', required=True,
                   help='Output directory for predictions')
    p.add_argument('--phage_feature_table',
                   help='Path to phage feature table (optional for single-strain mode)')
    p.add_argument('--feature_table',
                   help='Path to combined feature table for filtering')
    p.add_argument('--strain_source', default='strain',
                   help='Prefix for strain features (default: strain)')
    p.add_argument('--phage_source', default='phage',
                   help='Prefix for phage features (default: phage)')
    add_common_args(p)
    
    # ASSIGN-FEATURES - Assign features to new genomes
    p = subparsers.add_parser(
        'assign-features',
        help='Assign protein family features to new genomes',
        description='Map new genome sequences to existing protein family clusters.'
    )
    p.add_argument('--input_dir', required=True,
                   help='Directory containing genome FASTA files')
    p.add_argument('--mmseqs_db', required=True,
                   help='Path to existing MMseqs2 database')
    p.add_argument('--clusters_tsv', required=True,
                   help='Path to clusters TSV file')
    p.add_argument('--feature_map', required=True,
                   help='Path to feature mapping CSV (selected_features.csv)')
    p.add_argument('--output_dir', required=True,
                   help='Output directory for results')
    p.add_argument('--tmp_dir', required=True,
                   help='Temporary directory for intermediate files')
    p.add_argument('--genome_type', default='strain', choices=['strain', 'phage'],
                   help='Type of genomes to process (default: strain)')
    p.add_argument('--genome_list',
                   help='Path to file with list of genomes to process')
    p.add_argument('--sensitivity', type=float, default=7.5,
                   help='MMseqs2 sensitivity (default: 7.5)')
    p.add_argument('--coverage', type=float, default=0.8,
                   help='Minimum coverage for assignment (default: 0.8)')
    p.add_argument('--min_seq_id', type=float, default=0.4,
                   help='Minimum sequence identity (default: 0.4)')
    p.add_argument('--suffix', default='faa',
                   help='Suffix for FASTA files (default: faa)')
    p.add_argument('--duplicate_all', action='store_true',
                   help='Process all genomes even if duplicates found')
    add_common_args(p)
    
    # ANNOTATE - Annotate predictive features
    p = subparsers.add_parser(
        'annotate',
        help='Extract and annotate predictive protein features',
        description='Generate detailed annotations for predictive features from models.'
    )
    p.add_argument('--feature_file_path', required=True,
                   help='Path to file with predictive features')
    p.add_argument('--feature2cluster_path', required=True,
                   help='Path to feature-to-cluster mapping CSV')
    p.add_argument('--cluster2protein_path', required=True,
                   help='Path to cluster-to-protein mapping file')
    p.add_argument('--fasta_dir_or_file', required=True,
                   help='Path to FASTA file or directory with FASTA files')
    p.add_argument('--modeling_dir', required=True,
                   help='Directory containing modeling runs')
    p.add_argument('--output_dir', default='.',
                   help='Output directory (default: current directory)')
    p.add_argument('--output_fasta', default='predictive_AA_seqs.faa',
                   help='Name of output FASTA file (default: predictive_AA_seqs.faa)')
    p.add_argument('--protein_id_col', default='protein_ID',
                   help='Column name for protein IDs (default: protein_ID)')
    p.add_argument('--annotation_table_path',
                   help='Path to optional annotation table CSV')
    p.add_argument('--feature_assignments_path',
                   help='Path to feature assignments CSV')
    p.add_argument('--strain_column', default='strain',
                   help='Column for strain information (default: strain)')
    p.add_argument('--feature_type', default='strain', choices=['strain', 'phage'],
                   help='Type of features to extract (default: strain)')
    p.add_argument('--phenotype_column', default='interaction',
                   help='Column name for phenotype (default: interaction)')
    add_common_args(p)
    
    # SELECT-AND-TRAIN - Combined feature selection and training
    p = subparsers.add_parser(
        'select-and-train',
        help='Run feature selection and model training sequentially',
        description='Automated pipeline: feature selection → model training → extract predictive proteins.'
    )
    p.add_argument('--full_feature_table', '-i', required=True,
                   help='Path to full feature table')
    p.add_argument('--output', '-o', required=True,
                   help='Output directory for all results')
    p.add_argument('--num_features', default='none',
                   help='Number of features to select (default: none)')
    p.add_argument('--filter_type', default='strain', choices=['none', 'strain', 'phage'],
                   help='Filter type (default: strain)')
    p.add_argument('--num_runs_fs', type=int, default=25,
                   help='Number of feature selection runs (default: 25)')
    p.add_argument('--num_runs_modeling', type=int, default=50,
                   help='Number of modeling runs (default: 50)')
    p.add_argument('--sample_column', default='strain',
                   help='Sample identifier column (default: strain)')
    p.add_argument('--phage_column', default='phage',
                   help='Phage identifier column (default: phage)')
    p.add_argument('--phenotype_column', default='interaction',
                   help='Phenotype column (default: interaction)')
    p.add_argument('--method', default='rfe',
                   choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                   help='Feature selection method (default: rfe)')
    p.add_argument('--task_type', default='classification',
                   choices=['classification', 'regression'],
                   help='Modeling task type (default: classification)')
    p.add_argument('--max_features', default='none',
                   help='Maximum features to include (default: none)')
    p.add_argument('--max_ram', type=float, default=16,
                   help='Maximum RAM in GB (default: 16)')
    p.add_argument('--binary_data', action='store_true',
                   help='Convert to binary data')
    p.add_argument('--use_dynamic_weights', action='store_true',
                   help='Use dynamic weights')
    p.add_argument('--weights_method', default='inverse_frequency',
                   choices=['log10', 'inverse_frequency', 'balanced'],
                   help='Weight calculation method (default: inverse_frequency)')
    p.add_argument('--no-clustering', dest='use_clustering', action='store_false', 
                   default=True, help='Disable clustering (default: enabled)')
    p.add_argument('--cluster_method', default='hierarchical',
                   choices=['hdbscan', 'hierarchical'],
                   help='Clustering method (default: hierarchical)')
    p.add_argument('--n_clusters', type=int, default=20,
                   help='Number of clusters (default: 20)')
    p.add_argument('--min_cluster_size', type=int, default=5,
                   help='Minimum cluster size (default: 5)')
    p.add_argument('--min_samples', type=int,
                   help='Minimum samples for HDBSCAN')
    p.add_argument('--cluster_selection_epsilon', type=float, default=0.0,
                   help='HDBSCAN epsilon (default: 0.0)')
    p.add_argument('--check_feature_presence', action='store_true',
                   help='Check feature presence')
    p.add_argument('--filter_by_cluster_presence', action='store_true',
                   help='Filter by cluster presence')
    p.add_argument('--min_cluster_presence', type=int, default=2,
                   help='Min cluster presence (default: 2)')
    # Predictive proteins arguments
    p.add_argument('--run_predictive_proteins', action='store_true',
                   help='Extract predictive proteins after modeling')
    p.add_argument('--feature2cluster_path',
                   help='Path to strain feature-to-cluster mapping')
    p.add_argument('--cluster2protein_path',
                   help='Path to strain cluster-to-protein mapping')
    p.add_argument('--fasta_dir_or_file',
                   help='Path to strain FASTA file(s)')
    p.add_argument('--phage_feature2cluster_path',
                   help='Path to phage feature-to-cluster mapping')
    p.add_argument('--phage_cluster2protein_path',
                   help='Path to phage cluster-to-protein mapping')
    p.add_argument('--phage_fasta_dir_or_file',
                   help='Path to phage FASTA file(s)')
    p.add_argument('--annotation_table_path',
                   help='Path to annotation table CSV')
    p.add_argument('--protein_id_col', default='protein_ID',
                   help='Protein ID column name (default: protein_ID)')
    add_common_args(p)

    # PROTEIN-FAMILY-WORKFLOW - Complete protein family workflow
    p = subparsers.add_parser(
        'protein-family-workflow',
        help='Run complete protein family analysis workflow (clustering, feature selection, modeling)',
        description='Complete pipeline: clustering → feature selection → training → extract predictive proteins.'
    )
    # Required arguments
    p.add_argument('--input_path_strain', '-is', required=True,
                help='Path to the input directory or file for strain clustering')
    p.add_argument('--phenotype_matrix', '-pm', required=True,
                help='Path to phenotype matrix CSV file')
    p.add_argument('--output_dir', '-o', required=True,
                help='Output directory for results')

    # Optional inputs
    p.add_argument('--input_path_phage', '-ip',
                help='Path to the input directory or file for phage clustering (optional)')
    p.add_argument('--clustering_dir',
                help='Reuse existing clustering results from this directory')
    p.add_argument('--tmp_dir', default='tmp',
                help='Temporary directory for intermediate files (default: tmp)')

    # Clustering parameters
    p.add_argument('--min_seq_id', type=float, default=0.4,
                help='Minimum sequence identity for clustering (default: 0.4)')
    p.add_argument('--coverage', type=float, default=0.8,
                help='Minimum coverage for clustering (default: 0.8)')
    p.add_argument('--sensitivity', type=float, default=7.5,
                help='MMseqs2 sensitivity parameter (default: 7.5)')
    p.add_argument('--suffix', default='faa',
                help='File suffix for FASTA files (default: faa)')

    # Sample selection
    p.add_argument('--strain_list', default='none',
                help='Path to file with list of strains to include (default: none)')
    p.add_argument('--phage_list', default='none',
                help='Path to file with list of phages to include (default: none)')
    p.add_argument('--strain_column', default='strain',
                help='Column name for strain identifiers (default: strain)')
    p.add_argument('--phage_column', default='phage',
                help='Column name for phage identifiers (default: phage)')
    p.add_argument('--sample_column', default='strain',
                help='Column name for sample identifiers in phenotype matrix (default: strain)')
    p.add_argument('--phenotype_column', default='interaction',
                help='Column name for phenotype values (default: interaction)')

    # Feature parameters
    p.add_argument('--source_strain', default='strain',
                help='Prefix for strain features (default: strain)')
    p.add_argument('--source_phage', default='phage',
                help='Prefix for phage features (default: phage)')
    p.add_argument('--compare', action='store_true',
                help='Compare original clusters with assigned clusters')

    # Feature selection parameters
    p.add_argument('--num_features', default='none',
                help='Number of features to select (default: none)')
    p.add_argument('--filter_type', default='none',
                help='Type of feature filtering to apply (default: none)')
    p.add_argument('--method', default='rfe',
                choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                help='Feature selection method (default: rfe)')
    p.add_argument('--num_runs_fs', type=int, default=25,
                help='Number of runs for feature selection (default: 25)')

    # Modeling parameters
    p.add_argument('--num_runs_modeling', type=int, default=50,
                help='Number of runs for model training (default: 50)')
    p.add_argument('--task_type', default='classification',
                choices=['classification', 'regression'],
                help='Type of prediction task (default: classification)')
    p.add_argument('--max_features', default='none',
                help='Maximum number of features for modeling (default: none)')

    # Clustering-based filtering (for samples)
    p.add_argument('--use_clustering', action='store_true',
                help='Enable sample clustering for feature filtering')
    p.add_argument('--cluster_method', default='hdbscan',
                choices=['hdbscan', 'hierarchical'],
                help='Clustering method for samples (default: hdbscan)')
    p.add_argument('--n_clusters', type=int, default=20,
                help='Number of clusters for hierarchical (default: 20)')
    p.add_argument('--min_cluster_size', type=int, default=5,
                help='Minimum cluster size for HDBSCAN (default: 5)')
    p.add_argument('--min_samples', type=int,
                help='Min samples for core points (HDBSCAN)')
    p.add_argument('--cluster_selection_epsilon', type=float, default=0.0,
                help='Cluster selection epsilon (HDBSCAN) (default: 0.0)')
    p.add_argument('--check_feature_presence', action='store_true',
                help='Check if features are present in clusters')
    p.add_argument('--filter_by_cluster_presence', action='store_true',
                help='Filter features by cluster presence')
    p.add_argument('--min_cluster_presence', type=int, default=2,
                help='Minimum number of clusters a feature must appear in (default: 2)')

    # Feature clustering (pre-processing)
    p.add_argument('--use_feature_clustering', action='store_true',
                help='Enable pre-processing cluster-based feature filtering')
    p.add_argument('--feature_cluster_method', default='hierarchical',
                choices=['hierarchical'],
                help='Pre-processing clustering method (default: hierarchical)')
    p.add_argument('--feature_n_clusters', type=int, default=20,
                help='Number of clusters for pre-processing (default: 20)')
    p.add_argument('--feature_min_cluster_presence', type=int, default=2,
                help='Min clusters a feature must appear in pre-processing (default: 2)')

    # Annotation
    p.add_argument('--annotation_table_path',
                help='Path to annotation table for protein families')
    p.add_argument('--protein_id_col', default='protein_ID',
                help='Column name for protein IDs in annotation table (default: protein_ID)')

    # Advanced options
    p.add_argument('--use_dynamic_weights', action='store_true',
                help='Use dynamic class weights in modeling')
    p.add_argument('--weights_method', default='log10',
                choices=['log10', 'inverse_frequency', 'balanced'],
                help='Method for calculating class weights (default: log10)')
    p.add_argument('--use_shap', action='store_true',
                help='Calculate SHAP values for model interpretation')
    p.add_argument('--bootstrapping', action='store_true',
                help='Use bootstrapping in clustering')
    p.add_argument('--clear_tmp', action='store_true',
                help='Clear temporary files after completion')

    # System parameters
    p.add_argument('--max_ram', type=int, default=8,
                help='Maximum RAM usage in GB (default: 8)')
    add_common_args(p)
    
    # FULL-WORKFLOW - Complete protein family workflow
    p = subparsers.add_parser(
        'full-workflow',
        help='Run complete protein family analysis workflow',
        description='Complete pipeline: clustering → feature selection → training → prediction.'
    )
    p.add_argument('--input_strain', '-is', required=True,
                   help='Directory with strain genomes')
    p.add_argument('--input_phage', '-ip', required=True,
                   help='Directory with phage genomes')
    p.add_argument('--phenotype_matrix', '-pm', required=True,
                   help='Path to phenotype matrix')
    p.add_argument('--output', '-o', required=True,
                   help='Output directory')
    p.add_argument('--clustering_dir',
                   help='Existing clustering output directory (optional)')
    p.add_argument('--min_seq_id', type=float, default=0.4,
                   help='Minimum sequence identity (default: 0.4)')
    p.add_argument('--coverage', type=float, default=0.8,
                   help='Minimum coverage (default: 0.8)')
    p.add_argument('--sensitivity', type=float, default=7.5,
                   help='MMseqs2 sensitivity (default: 7.5)')
    p.add_argument('--suffix', default='faa',
                   help='FASTA file suffix (default: faa)')
    p.add_argument('--strain_list', default='none',
                   help='Strain list file or "none" (default: none)')
    p.add_argument('--phage_list', default='none',
                   help='Phage list file or "none" (default: none)')
    p.add_argument('--strain_column', default='strain',
                   help='Strain column name (default: strain)')
    p.add_argument('--phage_column', default='phage',
                   help='Phage column name (default: phage)')
    p.add_argument('--source_strain', default='strain',
                   help='Strain feature prefix (default: strain)')
    p.add_argument('--source_phage', default='phage',
                   help='Phage feature prefix (default: phage)')
    p.add_argument('--compare', action='store_true',
                   help='Compare clusters')
    p.add_argument('--num_features', default='none',
                   help='Number of features (default: none)')
    p.add_argument('--filter_type', default='strain',
                   help='Filter type (default: strain)')
    p.add_argument('--num_runs_fs', type=int, default=25,
                   help='Feature selection runs (default: 25)')
    p.add_argument('--num_runs_modeling', type=int, default=50,
                   help='Modeling runs (default: 50)')
    p.add_argument('--sample_column', default='strain',
                   help='Sample column (default: strain)')
    p.add_argument('--phenotype_column', default='interaction',
                   help='Phenotype column (default: interaction)')
    p.add_argument('--method', default='rfe',
                   help='Feature selection method (default: rfe)')
    p.add_argument('--annotation_table_path',
                   help='Annotation table path')
    p.add_argument('--protein_id_col', default='protein_ID',
                   help='Protein ID column (default: protein_ID)')
    p.add_argument('--task_type', default='classification',
                   choices=['classification', 'regression'],
                   help='Task type (default: classification)')
    p.add_argument('--max_features', default='none',
                   help='Max features (default: none)')
    p.add_argument('--max_ram', type=float, default=16,
                   help='Max RAM in GB (default: 16)')
    p.add_argument('--use_dynamic_weights', action='store_true')
    p.add_argument('--weights_method', default='inverse_frequency')
    p.add_argument('--no-clustering', dest='use_clustering', action='store_false', 
                   default=True, help='Disable clustering (default: enabled)')
    p.add_argument('--cluster_method', default='hierarchical')
    p.add_argument('--n_clusters', type=int, default=20)
    p.add_argument('--min_cluster_size', type=int, default=5)
    p.add_argument('--min_samples', type=int)
    p.add_argument('--cluster_selection_epsilon', type=float, default=0.0)
    p.add_argument('--check_feature_presence', action='store_true')
    p.add_argument('--filter_by_cluster_presence', action='store_true')
    p.add_argument('--min_cluster_presence', type=int, default=2)
    p.add_argument('--use_shap', action='store_true')
    p.add_argument('--clear_tmp', action='store_true')
    p.add_argument('--ignore_families', action='store_true')
    p.add_argument('--modeling', action='store_true')
    p.add_argument('--use_feature_clustering', action='store_true')
    p.add_argument('--feature_cluster_method', default='hierarchical')
    p.add_argument('--feature_n_clusters', type=int, default=20)
    p.add_argument('--feature_min_cluster_presence', type=int, default=2)
    # K-mer parameters for combined workflow
    p.add_argument('--k', type=int, default=5)
    p.add_argument('--k_range', action='store_true')
    p.add_argument('--remove_suffix', action='store_true')
    p.add_argument('--one_gene', action='store_true')
    add_common_args(p)
    
    # ==================== K-MER WORKFLOWS ====================
    
    # KMER-WORKFLOW - Complete k-mer analysis
    p = subparsers.add_parser(
        'kmer-workflow',
        help='Run complete k-mer based analysis workflow',
        description='Complete k-mer pipeline: feature generation → selection → training.'
    )
    p.add_argument('--input_strain_dir', '-is', required=True,
                   help='Directory with strain FASTA files')
    p.add_argument('--input_phage_dir', '-ip',
                   help='Directory with phage FASTA files (optional)')
    p.add_argument('--phenotype_matrix', '-pm', required=True,
                   help='Path to phenotype matrix')
    p.add_argument('--output', '-o', required=True,
                   help='Output directory')
    p.add_argument('--k', type=int, default=5,
                   help='K-mer size (default: 5)')
    p.add_argument('--k_range', action='store_true',
                   help='Use k-mer range instead of single k')
    p.add_argument('--suffix', default='faa',
                   help='FASTA file suffix (default: faa)')
    p.add_argument('--one_gene', action='store_true',
                   help='Use one gene per protein family')
    p.add_argument('--remove_suffix', action='store_true',
                   help='Remove suffix from strain names')
    p.add_argument('--strain_list', default='none',
                   help='Strain list or "none" (default: none)')
    p.add_argument('--phage_list', default='none',
                   help='Phage list or "none" (default: none)')
    p.add_argument('--strain_column', default='strain',
                   help='Strain column name (default: strain)')
    p.add_argument('--phage_column', default='phage',
                   help='Phage column name (default: phage)')
    p.add_argument('--sample_column', default='strain',
                   help='Sample column (default: strain)')
    p.add_argument('--phenotype_column', default='interaction',
                   help='Phenotype column (default: interaction)')
    p.add_argument('--num_features', default='none',
                   help='Number of features (default: none)')
    p.add_argument('--filter_type', default='strain',
                   help='Filter type (default: strain)')
    p.add_argument('--num_runs_fs', type=int, default=25,
                   help='Feature selection runs (default: 25)')
    p.add_argument('--num_runs_modeling', type=int, default=50,
                   help='Modeling runs (default: 50)')
    p.add_argument('--method', default='rfe',
                   help='Feature selection method (default: rfe)')
    p.add_argument('--task_type', default='classification',
                   choices=['classification', 'regression'],
                   help='Task type (default: classification)')
    p.add_argument('--max_features', default='none',
                   help='Max features (default: none)')
    p.add_argument('--max_ram', type=float, default=16,
                   help='Max RAM in GB (default: 16)')
    p.add_argument('--use_shap', action='store_true',
                   help='Use SHAP analysis')
    p.add_argument('--no-clustering', dest='use_clustering', action='store_false', 
                   default=True, help='Disable clustering (default: enabled)')
    p.add_argument('--cluster_method', default='hierarchical')
    p.add_argument('--n_clusters', type=int, default=20)
    p.add_argument('--min_cluster_size', type=int, default=5)
    p.add_argument('--min_samples', type=int)
    p.add_argument('--cluster_selection_epsilon', type=float, default=0.0)
    p.add_argument('--use_dynamic_weights', action='store_true')
    p.add_argument('--weights_method', default='inverse_frequency')
    p.add_argument('--check_feature_presence', action='store_true')
    p.add_argument('--filter_by_cluster_presence', action='store_true')
    p.add_argument('--min_cluster_presence', type=int, default=2)
    add_common_args(p)
    
    # KMER-ASSIGN-FEATURES - Assign k-mer features to new genomes
    p = subparsers.add_parser(
        'kmer-assign-features',
        help='Assign k-mer features to new genomes',
        description='Map new genomes to existing k-mer feature space.'
    )
    p.add_argument('--input_dir', required=True,
                   help='Directory with genome FASTA files')
    p.add_argument('--mmseqs_db', required=True,
                   help='Path to MMseqs2 database')
    p.add_argument('--clusters_tsv', required=True,
                   help='Path to clusters TSV')
    p.add_argument('--feature_map', required=True,
                   help='Path to feature map CSV')
    p.add_argument('--filtered_kmers', required=True,
                   help='Path to filtered k-mers CSV')
    p.add_argument('--aa_sequence_file', required=True,
                   help='Path to reference FASTA file')
    p.add_argument('--output_dir', required=True,
                   help='Output directory')
    p.add_argument('--tmp_dir', required=True,
                   help='Temporary directory')
    p.add_argument('--genome_type', default='strain', choices=['strain', 'phage'],
                   help='Genome type (default: strain)')
    p.add_argument('--genome_list',
                   help='Path to genome list file')
    p.add_argument('--sensitivity', type=float, default=7.5,
                   help='MMseqs2 sensitivity (default: 7.5)')
    p.add_argument('--coverage', type=float, default=0.8,
                   help='Minimum coverage (default: 0.8)')
    p.add_argument('--min_seq_id', type=float, default=0.4,
                   help='Minimum sequence identity (default: 0.4)')
    p.add_argument('--suffix', default='faa',
                   help='FASTA suffix (default: faa)')
    p.add_argument('--threshold', type=float, default=0.001,
                   help='K-mer matching threshold (default: 0.001)')
    p.add_argument('--reuse_existing', action='store_true',
                   help='Reuse existing outputs if available')
    add_common_args(p)
    
    # KMER-ASSIGN-PREDICT - Assign and predict with k-mers
    p = subparsers.add_parser(
        'kmer-assign-predict',
        help='Assign k-mer features and make predictions',
        description='Combined workflow: assign k-mer features → predict interactions.'
    )
    p.add_argument('--input_dir', required=True,
                   help='Directory with genome FASTA files')
    p.add_argument('--mmseqs_db', required=True,
                   help='Path to MMseqs2 database')
    p.add_argument('--clusters_tsv', required=True,
                   help='Path to clusters TSV')
    p.add_argument('--feature_map', required=True,
                   help='Path to feature map CSV')
    p.add_argument('--filtered_kmers', required=True,
                   help='Path to filtered k-mers CSV')
    p.add_argument('--aa_sequence_file', required=True,
                   help='Path to reference FASTA')
    p.add_argument('--tmp_dir', required=True,
                   help='Temporary directory')
    p.add_argument('--output_dir', required=True,
                   help='Output directory')
    p.add_argument('--model_dir', required=True,
                   help='Directory with trained models')
    p.add_argument('--feature_table',
                   help='Path to combined feature table')
    p.add_argument('--phage_feature_table',
                   help='Path to phage feature table')
    p.add_argument('--genome_type', default='strain', choices=['strain', 'phage'],
                   help='Genome type (default: strain)')
    p.add_argument('--genome_list',
                   help='Genome list file')
    p.add_argument('--sensitivity', type=float, default=7.5,
                   help='MMseqs2 sensitivity (default: 7.5)')
    p.add_argument('--coverage', type=float, default=0.8,
                   help='Minimum coverage (default: 0.8)')
    p.add_argument('--min_seq_id', type=float, default=0.4,
                   help='Min sequence identity (default: 0.4)')
    p.add_argument('--suffix', default='faa',
                   help='FASTA suffix (default: faa)')
    p.add_argument('--threshold', type=float, default=0.001,
                   help='K-mer threshold (default: 0.001)')
    p.add_argument('--reuse_existing', action='store_true',
                   help='Reuse existing outputs')
    add_common_args(p)
    
    # ASSIGN-PREDICT - Assign protein families and predict
    p = subparsers.add_parser(
        'assign-predict',
        help='Assign protein family features and make predictions',
        description='Combined workflow: assign features → predict interactions.'
    )
    p.add_argument('--input_dir', required=True,
                   help='Directory with genome FASTA files')
    p.add_argument('--mmseqs_db', required=True,
                   help='Path to MMseqs2 database')
    p.add_argument('--clusters_tsv', required=True,
                   help='Path to clusters TSV')
    p.add_argument('--feature_map', required=True,
                   help='Path to feature map')
    p.add_argument('--tmp_dir', required=True,
                   help='Temporary directory')
    p.add_argument('--output_dir', required=True,
                   help='Output directory')
    p.add_argument('--model_dir', required=True,
                   help='Model directory')
    p.add_argument('--feature_table',
                   help='Feature selection table for filtering')
    p.add_argument('--strain_feature_table',
                   help='Existing strain features (for phage prediction)')
    p.add_argument('--phage_feature_table',
                   help='Existing phage features (for strain prediction)')
    p.add_argument('--genome_type', default='strain', choices=['strain', 'phage'],
                   help='Genome type (default: strain)')
    p.add_argument('--genome_list',
                   help='Genome list file')
    p.add_argument('--sensitivity', type=float, default=7.5,
                   help='MMseqs2 sensitivity (default: 7.5)')
    p.add_argument('--coverage', type=float, default=0.8,
                   help='Minimum coverage (default: 0.8)')
    p.add_argument('--min_seq_id', type=float, default=0.4,
                   help='Min sequence identity (default: 0.4)')
    p.add_argument('--suffix', default='faa',
                   help='FASTA suffix (default: faa)')
    p.add_argument('--duplicate_all', action='store_true',
                   help='Process all genomes even with duplicates')
    add_common_args(p)
    
    return parser


def run_cluster(args):
    """Run protein family clustering workflow."""
    from genophi.workflows.feature_table_workflow import run_full_feature_workflow
    
    validate_directory(args.input_strain, "Strain input directory")
    validate_file(args.phenotype_matrix, "Phenotype matrix")
    if args.input_phage:
        validate_directory(args.input_phage, "Phage input directory")
    validate_directory(args.output, "Output directory", create=True)
    
    run_full_feature_workflow(
        input_path_strain=args.input_strain,
        input_path_phage=args.input_phage,
        phenotype_matrix=args.phenotype_matrix,
        output_dir=args.output,
        tmp_dir=args.tmp,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        strain_column=args.strain_column,
        phage_list=args.phage_list,
        phage_column=args.phage_column,
        compare=args.compare,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        max_ram=args.max_ram,
        use_feature_clustering=args.use_feature_clustering,
        feature_cluster_method=args.feature_cluster_method,
        feature_n_clusters=args.feature_n_clusters,
        feature_min_cluster_presence=args.feature_min_cluster_presence
    )


def run_select_features(args):
    """Run feature selection workflow."""
    from genophi.workflows.feature_selection_workflow import run_feature_selection_workflow
    
    validate_file(args.input, "Input feature table")
    validate_directory(args.output, "Output directory", create=True)
    
    run_feature_selection_workflow(
        input_path=args.input,
        base_output_dir=args.output,
        threads=args.threads,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs=args.num_runs,
        method=args.method,
        task_type=args.task_type,
        phenotype_column=args.phenotype_column,
        sample_column=args.sample_column,
        phage_column=args.phage_column,
        binary_data=args.binary_data,
        max_features=args.max_features,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        filter_by_cluster_presence=args.filter_by_cluster_presence,
        min_cluster_presence=args.min_cluster_presence,
        max_ram=args.max_ram
    )


def run_train(args):
    """Run model training workflow."""
    from genophi.workflows.modeling_workflow import run_modeling_workflow
    
    validate_directory(args.input, "Input directory")
    validate_directory(args.output, "Output directory", create=True)
    
    run_modeling_workflow(
        input_dir=args.input,
        base_output_dir=args.output,
        threads=args.threads,
        num_runs=args.num_runs,
        set_filter=args.set_filter,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        phage_column=args.phage_column,
        task_type=args.task_type,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        binary_data=args.binary_data
    )


def run_predict(args):
    """Run prediction workflow."""
    from genophi.workflows.prediction_workflow import run_prediction_workflow
    
    validate_directory(args.input_dir, "Input directory")
    validate_directory(args.model_dir, "Model directory")
    validate_directory(args.output_dir, "Output directory", create=True)
    
    run_prediction_workflow(
        input_dir=args.input_dir,
        model_dir=args.model_dir,
        output_dir=args.output_dir,
        phage_feature_table_path=args.phage_feature_table,
        feature_table=args.feature_table,
        strain_source=args.strain_source,
        phage_source=args.phage_source,
        threads=args.threads
    )


def run_assign_features(args):
    """Run feature assignment workflow."""
    from genophi.workflows.assign_features_workflow import run_assign_features_workflow
    
    validate_directory(args.input_dir, "Input directory")
    validate_file(args.mmseqs_db, "MMseqs2 database")
    validate_file(args.clusters_tsv, "Clusters TSV")
    validate_file(args.feature_map, "Feature map")
    validate_directory(args.output_dir, "Output directory", create=True)
    validate_directory(args.tmp_dir, "Temporary directory", create=True)
    
    run_assign_features_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        feature_map=args.feature_map,
        clusters_tsv=args.clusters_tsv,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        duplicate_all=args.duplicate_all
    )


def run_annotate(args):
    """Run annotation workflow."""
    from genophi.workflows.feature_annotations_workflow import run_predictive_proteins_workflow
    
    validate_file(args.feature_file_path, "Feature file")
    validate_file(args.feature2cluster_path, "Feature to cluster mapping")
    validate_file(args.cluster2protein_path, "Cluster to protein mapping")
    validate_directory(args.modeling_dir, "Modeling directory")
    validate_directory(args.output_dir, "Output directory", create=True)
    
    run_predictive_proteins_workflow(
        feature_file_path=args.feature_file_path,
        feature2cluster_path=args.feature2cluster_path,
        cluster2protein_path=args.cluster2protein_path,
        fasta_dir_or_file=args.fasta_dir_or_file,
        modeling_dir=args.modeling_dir,
        output_dir=args.output_dir,
        output_fasta=args.output_fasta,
        protein_id_col=args.protein_id_col,
        annotation_table_path=args.annotation_table_path,
        feature_assignments_path=args.feature_assignments_path,
        strain_column=args.strain_column,
        feature_type=args.feature_type,
        phenotype_column=args.phenotype_column
    )


def run_select_and_train(args):
    """Run combined feature selection and training workflow."""
    from genophi.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table
    
    validate_file(args.full_feature_table, "Full feature table")
    validate_directory(args.output, "Output directory", create=True)
    
    run_modeling_workflow_from_feature_table(
        full_feature_table=args.full_feature_table,
        output_dir=args.output,
        threads=args.threads,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phage_column=args.phage_column,
        phenotype_column=args.phenotype_column,
        method=args.method,
        annotation_table_path=args.annotation_table_path,
        protein_id_col=args.protein_id_col,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram,
        binary_data=args.binary_data,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        filter_by_cluster_presence=args.filter_by_cluster_presence,
        min_cluster_presence=args.min_cluster_presence,
        run_predictive_proteins=args.run_predictive_proteins,
        feature2cluster_path=args.feature2cluster_path,
        cluster2protein_path=args.cluster2protein_path,
        fasta_dir_or_file=args.fasta_dir_or_file,
        phage_feature2cluster_path=args.phage_feature2cluster_path,
        phage_cluster2protein_path=args.phage_cluster2protein_path,
        phage_fasta_dir_or_file=args.phage_fasta_dir_or_file
    )


def run_full_workflow(args):
    """Run complete protein family workflow."""
    from genophi.workflows.full_workflow import run_full_workflow
    
    validate_directory(args.input_strain, "Strain input directory")
    validate_directory(args.input_phage, "Phage input directory")
    validate_file(args.phenotype_matrix, "Phenotype matrix")
    validate_directory(args.output, "Output directory", create=True)
    
    run_full_workflow(
        input_strain=args.input_strain,
        input_phage=args.input_phage,
        phenotype_matrix=args.phenotype_matrix,
        output=args.output,
        clustering_dir=args.clustering_dir,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        compare=args.compare,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        method=args.method,
        annotation_table_path=args.annotation_table_path,
        protein_id_col=args.protein_id_col,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram,
        threads=args.threads,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        filter_by_cluster_presence=args.filter_by_cluster_presence,
        min_cluster_presence=args.min_cluster_presence,
        use_shap=args.use_shap,
        clear_tmp=args.clear_tmp,
        k=args.k,
        k_range=args.k_range,
        remove_suffix=args.remove_suffix,
        one_gene=args.one_gene,
        ignore_families=args.ignore_families,
        modeling=args.modeling,
        use_feature_clustering=args.use_feature_clustering,
        feature_cluster_method=args.feature_cluster_method,
        feature_n_clusters=args.feature_n_clusters,
        feature_min_cluster_presence=args.feature_min_cluster_presence
    )


def run_kmer_workflow(args):
    """Run complete k-mer workflow."""
    from genophi.workflows.kmer_full_workflow import run_kmer_workflow
    
    validate_directory(args.input_strain_dir, "Strain input directory")
    if args.input_phage_dir:
        validate_directory(args.input_phage_dir, "Phage input directory")
    validate_file(args.phenotype_matrix, "Phenotype matrix")
    validate_directory(args.output, "Output directory", create=True)
    
    run_kmer_workflow(
        input_strain_dir=args.input_strain_dir,
        input_phage_dir=args.input_phage_dir,
        phenotype_matrix=args.phenotype_matrix,
        output_dir=args.output,
        k=args.k,
        k_range=args.k_range,
        suffix=args.suffix,
        one_gene=args.one_gene,
        remove_suffix=args.remove_suffix,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        method=args.method,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram,
        threads=args.threads,
        use_shap=args.use_shap,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        check_feature_presence=args.check_feature_presence,
        filter_by_cluster_presence=args.filter_by_cluster_presence,
        min_cluster_presence=args.min_cluster_presence
    )


def run_kmer_assign_features(args):
    """Run k-mer feature assignment workflow."""
    from genophi.workflows.kmer_assign_features_workflow import run_kmer_assign_features_workflow
    
    validate_directory(args.input_dir, "Input directory")
    validate_file(args.mmseqs_db, "MMseqs2 database")
    validate_file(args.clusters_tsv, "Clusters TSV")
    validate_file(args.feature_map, "Feature map")
    validate_file(args.filtered_kmers, "Filtered k-mers")
    validate_file(args.aa_sequence_file, "AA sequence file")
    validate_directory(args.output_dir, "Output directory", create=True)
    validate_directory(args.tmp_dir, "Temporary directory", create=True)
    
    run_kmer_assign_features_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        feature_map=args.feature_map,
        filtered_kmers=args.filtered_kmers,
        aa_sequence_file=args.aa_sequence_file,
        clusters_tsv=args.clusters_tsv,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        threshold=args.threshold,
        reuse_existing=args.reuse_existing
    )


def run_kmer_assign_predict(args):
    """Run k-mer assignment and prediction workflow."""
    from genophi.workflows.kmer_assign_predict_workflow import kmer_assign_predict_workflow
    
    validate_directory(args.input_dir, "Input directory")
    validate_file(args.mmseqs_db, "MMseqs2 database")
    validate_file(args.clusters_tsv, "Clusters TSV")
    validate_file(args.feature_map, "Feature map")
    validate_file(args.filtered_kmers, "Filtered k-mers")
    validate_file(args.aa_sequence_file, "AA sequence file")
    validate_directory(args.model_dir, "Model directory")
    validate_directory(args.output_dir, "Output directory", create=True)
    validate_directory(args.tmp_dir, "Temporary directory", create=True)
    
    kmer_assign_predict_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        feature_map=args.feature_map,
        filtered_kmers=args.filtered_kmers,
        aa_sequence_file=args.aa_sequence_file,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        feature_table=args.feature_table,
        model_dir=args.model_dir,
        phage_feature_table_path=args.phage_feature_table,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        threshold=args.threshold,
        reuse_existing=args.reuse_existing
    )


def run_assign_predict(args):
    """Run protein family assignment and prediction workflow."""
    from genophi.workflows.assign_predict_workflow import assign_predict_workflow
    
    validate_directory(args.input_dir, "Input directory")
    validate_file(args.mmseqs_db, "MMseqs2 database")
    validate_file(args.clusters_tsv, "Clusters TSV")
    validate_file(args.feature_map, "Feature map")
    validate_directory(args.model_dir, "Model directory")
    validate_directory(args.output_dir, "Output directory", create=True)
    validate_directory(args.tmp_dir, "Temporary directory", create=True)
    
    assign_predict_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        feature_map=args.feature_map,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        model_dir=args.model_dir,
        feature_table=args.feature_table,
        strain_feature_table_path=args.strain_feature_table,  # Fixed: map to correct param name
        phage_feature_table_path=args.phage_feature_table,    # Fixed: map to correct param name
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        duplicate_all=args.duplicate_all
    )

def protein_family_workflow_command(args):
    """Execute protein family workflow command."""
    from genophi.workflows.protein_family_workflow import setup_logging
    from genophi.workflows.protein_family_workflow import run_protein_family_workflow
    setup_logging(args.output_dir)
    
    run_protein_family_workflow(
        input_path_strain=args.input_path_strain,
        input_path_phage=args.input_path_phage,
        phenotype_matrix=args.phenotype_matrix,
        output_dir=args.output_dir,
        clustering_dir=args.clustering_dir,
        tmp_dir=args.tmp_dir,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        compare=args.compare,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        method=args.method,
        annotation_table_path=args.annotation_table_path,
        protein_id_col=args.protein_id_col,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        filter_by_cluster_presence=args.filter_by_cluster_presence,
        min_cluster_presence=args.min_cluster_presence,
        use_shap=args.use_shap,
        bootstrapping=args.bootstrapping,
        clear_tmp=args.clear_tmp,
        use_feature_clustering=args.use_feature_clustering,
        feature_cluster_method=args.feature_cluster_method,
        feature_n_clusters=args.feature_n_clusters,
        feature_min_cluster_presence=args.feature_min_cluster_presence
    )

def main():
    """Main entry point for the CLI."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging if verbose flag is present
    if hasattr(args, 'verbose'):
        setup_logging(args.verbose)
    else:
        setup_logging(False)
    
    try:
        # Route to appropriate workflow based on command
        if args.command == 'cluster':
            run_cluster(args)
        elif args.command == 'select-features':
            run_select_features(args)
        elif args.command == 'train':
            run_train(args)
        elif args.command == 'predict':
            run_predict(args)
        elif args.command == 'assign-features':
            run_assign_features(args)
        elif args.command == 'annotate':
            run_annotate(args)
        elif args.command == 'select-and-train':
            run_select_and_train(args)
        elif args.command == 'protein-family-workflow':
            protein_family_workflow_command(args)
        elif args.command == 'full-workflow':
            run_full_workflow(args)
        elif args.command == 'kmer-workflow':
            run_kmer_workflow(args)
        elif args.command == 'kmer-assign-features':
            run_kmer_assign_features(args)
        elif args.command == 'kmer-assign-predict':
            run_kmer_assign_predict(args)
        elif args.command == 'assign-predict':
            run_assign_predict(args)
        else:
            parser.print_help()
            sys.exit(1)
            
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        sys.exit(1)
    except ValueError as e:
        logging.error(f"Invalid value: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("\nOperation cancelled by user")
        sys.exit(130)
    except Exception as e:
        logging.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()