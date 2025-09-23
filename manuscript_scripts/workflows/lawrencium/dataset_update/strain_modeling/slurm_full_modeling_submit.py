#!/usr/bin/env python3
"""
Modified SLURM workflow submission for full_workflow.py
Breaks the workflow into 5 sequential SLURM jobs with proper dependencies.
Stage 5 now includes both k-mer generation AND modeling.
Includes capability to resume from specific stages.
"""

import os
import sys
import argparse
import subprocess
import time

def check_stage_completion(args, stage):
    """Check if a stage has been completed based on expected output files."""
    output_dir = args.output
    
    completion_markers = {
        1: f"{output_dir}/strain/features/feature_table.csv",
        2: f"{output_dir}/feature_selection/filtered_feature_tables",
        3: f"{output_dir}/modeling_results/model_performance/model_performance_metrics.csv", 
        4: f"{output_dir}/modeling_results/model_performance/predictive_proteins/strain_predictive_AA_seqs.faa",
        5: f"{output_dir}/kmer_modeling/modeling/model_performance/model_performance_metrics.csv"
    }
    
    marker = completion_markers.get(stage)
    if marker and os.path.exists(marker):
        print(f"✅ Stage {stage} appears complete (found: {marker})")
        return True
    return False

def submit_job(script_path):
    """Submit a SLURM job and return job ID."""
    try:
        result = subprocess.run(['sbatch', '--parsable', script_path], 
                              capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error submitting {script_path}: {e}")
        print(f"Error output: {e.stderr}")
        return None

def create_stage1_clustering(args, run_dir):
    """Stage 1: MMSeqs2 Clustering"""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=stage1_clustering
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=6:00:00
#SBATCH --output=logs/stage1_%j.out
#SBATCH --error=logs/stage1_%j.err

echo "=== Stage 1: Clustering ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import sys
sys.path.append('{os.path.dirname(os.path.abspath(__file__))}')
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment
import os

# Strain clustering
if not os.path.exists('{args.output}/strain/features/feature_table.csv'):
    print('Running strain clustering...')
    run_clustering_workflow(
        '{args.input_strain}', '{args.output}/strain', '{args.output}/tmp/strain',
        {args.min_seq_id}, {args.coverage}, {args.sensitivity}, '{args.suffix}', 
        {args.threads}, '{args.strain_list}', '{args.strain_column}', {args.compare}
    )
    
    run_feature_assignment(
        '{args.output}/strain/presence_absence_matrix.csv',
        '{args.output}/strain/features',
        source='{args.source_strain}',
        select='{args.strain_list}',
        select_column='{args.strain_column}',
        max_ram={args.max_ram}
    )

# Phage clustering if provided
if '{args.input_phage}' and '{args.input_phage}' != 'None':
    if not os.path.exists('{args.output}/phage/features/feature_table.csv'):
        print('Running phage clustering...')
        run_clustering_workflow(
            '{args.input_phage}', '{args.output}/phage', '{args.output}/tmp/phage',
            {args.min_seq_id}, {args.coverage}, {args.sensitivity}, '{args.suffix}',
            {args.threads}, '{args.phage_list}', '{args.phage_column}', {args.compare}
        )
        
        run_feature_assignment(
            '{args.output}/phage/presence_absence_matrix.csv',
            '{args.output}/phage/features',
            source='{args.source_phage}',
            select='{args.phage_list}',
            select_column='{args.phage_column}',
            max_ram={args.max_ram}
        )

print('Stage 1 completed successfully')
"

echo "Stage 1 completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "stage1_clustering.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_stage2_feature_selection(args, run_dir, dependency):
    """Stage 2: Feature Selection"""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=stage2_feature_sel
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=6:00:00
#SBATCH --dependency=afterok:{dependency}
#SBATCH --output=logs/stage2_%j.out
#SBATCH --error=logs/stage2_%j.err

echo "=== Stage 2: Feature Selection ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import sys
sys.path.append('{os.path.dirname(os.path.abspath(__file__))}')
from phage_modeling.mmseqs2_clustering import merge_feature_tables
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
import pandas as pd
import os

# Merge feature tables with phenotype
if '{args.input_phage}' and '{args.input_phage}' != 'None':
    merged_dir = '{args.output}/merged'
    os.makedirs(merged_dir, exist_ok=True)
    feature_input = merge_feature_tables(
        strain_features='{args.output}/strain/features/feature_table.csv',
        phenotype_matrix='{args.phenotype_matrix}',
        output_dir=merged_dir,
        sample_column='{args.sample_column}',
        phage_features='{args.output}/phage/features/feature_table.csv',
        remove_suffix=False,
        use_feature_clustering={args.use_feature_clustering},
        feature_cluster_method='{args.feature_cluster_method}',
        feature_n_clusters={args.feature_n_clusters},
        feature_min_cluster_presence={args.feature_min_cluster_presence}
    )
else:
    feature_input = merge_feature_tables(
        strain_features='{args.output}/strain/features/feature_table.csv',
        phenotype_matrix='{args.phenotype_matrix}',
        output_dir='{args.output}',
        sample_column='{args.sample_column}',
        remove_suffix=False,
        use_feature_clustering={args.use_feature_clustering},
        feature_cluster_method='{args.feature_cluster_method}',
        feature_n_clusters={args.feature_n_clusters},
        feature_min_cluster_presence={args.feature_min_cluster_presence}
    )

# Determine num_features
df = pd.read_csv(feature_input)
num_rows = len(df)
if '{args.num_features}' == 'none':
    if num_rows < 500:
        num_features = 50
    elif num_rows < 2000:
        num_features = 100
    else:
        num_features = int(num_rows/20)
else:
    num_features = int('{args.num_features}')

# Run feature selection
run_feature_selection_iterations(
    input_path=feature_input,
    base_output_dir='{args.output}/feature_selection',
    threads={args.threads},
    num_features=num_features,
    filter_type='{args.filter_type}',
    num_runs={args.num_runs_fs},
    method='{args.method}',
    sample_column='{args.sample_column}',
    phenotype_column='{args.phenotype_column}',
    phage_column='{args.phage_column}',
    task_type='{args.task_type}',
    use_dynamic_weights={args.use_dynamic_weights},
    weights_method='{args.weights_method}',
    use_clustering={args.use_clustering},
    cluster_method='{args.cluster_method}',
    n_clusters={args.n_clusters},
    min_cluster_size={args.min_cluster_size},
    min_samples={args.min_samples},
    cluster_selection_epsilon={args.cluster_selection_epsilon},
    check_feature_presence={args.check_feature_presence},
    filter_by_cluster_presence={args.filter_by_cluster_presence},
    min_cluster_presence={args.min_cluster_presence},
    max_ram={args.max_ram}
)

# Generate feature tables
max_features = num_features if '{args.max_features}' == 'none' else int('{args.max_features}')
generate_feature_tables(
    model_testing_dir='{args.output}/feature_selection',
    full_feature_table_file=feature_input,
    filter_table_dir='{args.output}/feature_selection/filtered_feature_tables',
    phenotype_column='{args.phenotype_column}',
    sample_column='{args.sample_column}',
    cut_offs=[3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50],
    binary_data=True,
    max_features=max_features,
    filter_type='{args.filter_type}'
)

print('Stage 2 completed successfully')
"

echo "Stage 2 completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "stage2_feature_selection.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_stage3_modeling(args, run_dir, dependency):
    """Stage 3: Modeling"""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=stage3_modeling
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --dependency=afterok:{dependency}
#SBATCH --output=logs/stage3_%j.out
#SBATCH --error=logs/stage3_%j.err

echo "=== Stage 3: Modeling ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import sys
sys.path.append('{os.path.dirname(os.path.abspath(__file__))}')
from phage_modeling.select_feature_modeling import run_experiments

run_experiments(
    input_dir='{args.output}/feature_selection/filtered_feature_tables',
    base_output_dir='{args.output}/modeling_results',
    threads={args.threads},
    num_runs={args.num_runs_modeling},
    set_filter='{args.filter_type}',
    sample_column='{args.sample_column}',
    phenotype_column='{args.phenotype_column}',
    phage_column='{args.phage_column}',
    use_dynamic_weights={args.use_dynamic_weights},
    weights_method='{args.weights_method}',
    task_type='{args.task_type}',
    binary_data=True,
    max_ram={args.max_ram},
    use_clustering={args.use_clustering},
    cluster_method='{args.cluster_method}',
    n_clusters={args.n_clusters},
    min_cluster_size={args.min_cluster_size},
    min_samples={args.min_samples},
    cluster_selection_epsilon={args.cluster_selection_epsilon},
    use_shap={args.use_shap}
)

print('Stage 3 completed successfully')
"

echo "Stage 3 completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "stage3_modeling.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_stage4_predictive_proteins(args, run_dir, dependency):
    """Stage 4: Predictive Proteins"""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=stage4_pred_proteins
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH --dependency=afterok:{dependency}
#SBATCH --output=logs/stage4_%j.out
#SBATCH --error=logs/stage4_%j.err

echo "=== Stage 4: Predictive Proteins ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import sys
sys.path.append('{os.path.dirname(os.path.abspath(__file__))}')
from phage_modeling.workflows.feature_annotations_workflow import run_predictive_proteins_workflow
import pandas as pd
import os

# Get top cutoff
metrics_file = '{args.output}/modeling_results/model_performance/model_performance_metrics.csv'
performance_df = pd.read_csv(metrics_file)
top_cutoff = performance_df.iloc[0]['cut_off'].split('_')[-1]

# Set up paths
feature_file_path = f'{args.output}/feature_selection/filtered_feature_tables/select_feature_table_cutoff_{{top_cutoff}}.csv'
modeling_dir = f'{args.output}/modeling_results/cutoff_{{top_cutoff}}'
output_dir = '{args.output}/modeling_results/model_performance/predictive_proteins'

# Strain predictive proteins
modified_AA_path = '{args.output}/strain/modified_AAs/strain'
fasta_input = modified_AA_path if os.path.exists(modified_AA_path) else '{args.input_strain}'

run_predictive_proteins_workflow(
    feature_file_path=feature_file_path,
    feature2cluster_path='{args.output}/strain/features/selected_features.csv',
    cluster2protein_path='{args.output}/strain/clusters.tsv',
    fasta_dir_or_file=fasta_input,
    modeling_dir=modeling_dir,
    output_dir=output_dir,
    output_fasta='predictive_AA_seqs.faa',
    protein_id_col='{args.protein_id_col}',
    annotation_table_path={repr(args.annotation_table_path)},
    feature_assignments_path='{args.output}/strain/features/feature_assignments.csv',
    strain_column='strain',
    feature_type='strain'
)

# Phage predictive proteins if phage data provided
if '{args.input_phage}' and '{args.input_phage}' != 'None':
    modified_AA_path_phage = '{args.output}/phage/modified_AAs/phage'
    fasta_input_phage = modified_AA_path_phage if os.path.exists(modified_AA_path_phage) else '{args.input_phage}'
    
    run_predictive_proteins_workflow(
        feature_file_path=feature_file_path,
        feature2cluster_path='{args.output}/phage/features/selected_features.csv',
        cluster2protein_path='{args.output}/phage/clusters.tsv',
        fasta_dir_or_file=fasta_input_phage,
        modeling_dir=modeling_dir,
        output_dir=output_dir,
        output_fasta='predictive_AA_seqs.faa',
        protein_id_col='{args.protein_id_col}',
        feature_assignments_path='{args.output}/phage/features/feature_assignments.csv',
        strain_column='phage',
        feature_type='phage'
    )

print('Stage 4 completed successfully')
"

echo "Stage 4 completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "stage4_predictive_proteins.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_stage5_kmer_generation_and_modeling(args, run_dir, dependency):
    """Stage 5: K-mer Generation AND Modeling (Combined)"""
    script_content = f"""#!/bin/bash
#SBATCH --job-name=stage5_kmer_complete
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=10:00:00
#SBATCH --dependency=afterok:{dependency}
#SBATCH --output=logs/stage5_%j.out
#SBATCH --error=logs/stage5_%j.err

echo "=== Stage 5: K-mer Generation & Modeling ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import sys
sys.path.append('{os.path.dirname(os.path.abspath(__file__))}')
from phage_modeling.workflows.kmer_table_workflow import run_kmer_table_workflow
import pandas as pd
import os

# Get paths from stage 4 results
metrics_file = '{args.output}/modeling_results/model_performance/model_performance_metrics.csv'
performance_df = pd.read_csv(metrics_file)
top_cutoff = performance_df.iloc[0]['cut_off'].split('_')[-1]
feature_file_path = f'{args.output}/feature_selection/filtered_feature_tables/select_feature_table_cutoff_{{top_cutoff}}.csv'

strain_fasta = '{args.output}/modeling_results/model_performance/predictive_proteins/strain_predictive_AA_seqs.faa'
protein_csv = '{args.output}/modeling_results/model_performance/predictive_proteins/strain_predictive_feature_overview.csv'
phage_fasta = '{args.output}/modeling_results/model_performance/predictive_proteins/phage_predictive_AA_seqs.faa'
protein_csv_phage = '{args.output}/modeling_results/model_performance/predictive_proteins/phage_predictive_feature_overview.csv'

# Validate required files exist
for file_path in [strain_fasta, protein_csv]:
    if not os.path.exists(file_path):
        print(f'Error: Required file {{file_path}} not found')
        sys.exit(1)

print('Running k-mer table workflow with modeling=True...')
run_kmer_table_workflow(
    strain_fasta=strain_fasta,
    protein_csv=protein_csv,
    k={args.k},
    id_col='strain',
    one_gene={args.one_gene},
    output_dir='{args.output}/kmer_modeling',
    k_range={args.k_range},
    phenotype_matrix='{args.phenotype_matrix}',
    phage_fasta=phage_fasta if os.path.exists(phage_fasta) else None,
    protein_csv_phage=protein_csv_phage if os.path.exists(protein_csv_phage) else None,
    remove_suffix={args.remove_suffix},
    sample_column='{args.sample_column}',
    phenotype_column='{args.phenotype_column}',
    modeling=True,  # CHANGED: Now True instead of False
    filter_type='{args.filter_type}',
    num_features='{args.num_features}',
    num_runs_fs={args.num_runs_fs},
    num_runs_modeling={args.num_runs_modeling},
    method='{args.method}',
    strain_list=feature_file_path,
    phage_list=feature_file_path,
    threads={args.threads},
    task_type='{args.task_type}',
    max_features='{args.max_features}',
    ignore_families={args.ignore_families},
    max_ram={args.max_ram},
    use_shap={args.use_shap},
    use_dynamic_weights={args.use_dynamic_weights},
    weights_method='{args.weights_method}',
    use_clustering={args.use_clustering},
    cluster_method='{args.cluster_method}',
    n_clusters={args.n_clusters},
    min_cluster_size={args.min_cluster_size},
    min_samples={args.min_samples},
    cluster_selection_epsilon={args.cluster_selection_epsilon},
    check_feature_presence={args.check_feature_presence},
    filter_by_cluster_presence={args.filter_by_cluster_presence},
    min_cluster_presence={args.min_cluster_presence}
)

print('Stage 5 (K-mer generation & modeling) completed successfully')
"

echo "Stage 5 completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "stage5_kmer_complete.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Submit full_workflow.py as 5 sequential SLURM jobs")
    
    # Copy ALL arguments from full_workflow.py
    parser.add_argument('--input_strain', required=True, help='Input strain FASTA path.')
    parser.add_argument('--input_phage', help='Input phage FASTA path.')
    parser.add_argument('--phenotype_matrix', required=True, help='Phenotype matrix file path.')
    parser.add_argument('--output', required=True, help='Output directory.')
    parser.add_argument('--clustering_dir', help='Path to an existing strain clustering directory.')
    parser.add_argument('--min_seq_id', type=float, default=0.4, help='Minimum sequence identity for clustering (default: 0.4).')
    parser.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering (default: 0.8).')
    parser.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering (default: 7.5).')
    parser.add_argument('--suffix', default='faa', help='Suffix for input FASTA files (default: faa).')
    parser.add_argument('--strain_list', default='none', help='List of strains for filtering (default: none).')
    parser.add_argument('--phage_list', default='none', help='List of phages for filtering (default: none).')
    parser.add_argument('--strain_column', default='strain', help='Column name for strain data (default: strain).')
    parser.add_argument('--phage_column', default='phage', help='Column name for phage data (default: phage).')
    parser.add_argument('--source_strain', default='strain', help='Source prefix for strain (default: strain).')
    parser.add_argument('--source_phage', default='phage', help='Source prefix for phage (default: phage).')
    parser.add_argument('--compare', action='store_true', help='Compare clustering results.')
    parser.add_argument('--num_features', default='none', help='Number of features for selection (default: none).')
    parser.add_argument('--filter_type', default='none', choices=['none', 'strain', 'phage'], help='Filter type for feature selection. (default: none)')
    parser.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection runs (default: 10).')
    parser.add_argument('--num_runs_modeling', type=int, default=20, help='Number of modeling runs (default: 20).')
    parser.add_argument('--sample_column', default='strain', help='Sample column name (default: strain).')
    parser.add_argument('--phenotype_column', default='interaction', help='Phenotype column name (default: interaction).')
    parser.add_argument('--method', default='rfe', choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'], help='Feature selection method (default: rfe)')
    parser.add_argument('--annotation_table_path', help='Path to annotation table.')
    parser.add_argument('--protein_id_col', default='protein_ID', help='Protein ID column name (default: protein_ID).')
    parser.add_argument('--task_type', default='classification', choices=['classification', 'regression'], help='Task type for modeling (default: classification)')
    parser.add_argument('--max_features', default='none', help='Maximum number of features for modeling (default: none).')
    parser.add_argument('--max_ram', type=float, default=8, help='Maximum RAM usage in GB (default: 8).')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')
    parser.add_argument('--use_dynamic_weights', action='store_true', help='Use dynamic weights for feature selection and modeling.')
    parser.add_argument('--weights_method', default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help='Method to calculate class weights (default: log10)')
    parser.add_argument('--use_clustering', action='store_true', help='Use clustering for feature selection.')
    parser.add_argument('--cluster_method', default='hdbscan', choices=['hdbscan', 'hierarchical'], help='Clustering method for feature selection (default: hdbscan)')
    parser.add_argument('--n_clusters', type=int, default=20, help='Number of clusters for hierarchical clustering (default: 20)')
    parser.add_argument('--min_cluster_size', type=int, default=5, help='Minimum cluster size for HDBSCAN clustering (default: 5)')
    parser.add_argument('--min_samples', type=int, default=None, help='Min samples parameter for HDBSCAN')
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help='Cluster selection epsilon for HDBSCAN (default: 0.0)')
    parser.add_argument('--check_feature_presence', action='store_true', help='Check for presence of features during train-test split.')
    parser.add_argument('--filter_by_cluster_presence', action='store_true', help='Filter features by cluster/group presence instead of train/test presence.')
    parser.add_argument('--min_cluster_presence', type=int, default=2, help='Minimum number of clusters/groups a feature must be present in (default: 2).')
    
    # Feature clustering parameters (pre-processing)
    parser.add_argument('--use_feature_clustering', action='store_true', help='Enable pre-processing cluster-based feature filtering')
    parser.add_argument('--feature_cluster_method', default='hierarchical', choices=['hierarchical'], help='Pre-processing clustering method')
    parser.add_argument('--feature_n_clusters', type=int, default=20, help='Number of clusters for pre-processing feature clustering')
    parser.add_argument('--feature_min_cluster_presence', type=int, default=2, help='Min clusters a feature must appear in during pre-processing')
    
    parser.add_argument('--use_shap', action='store_true', help='Use SHAP values for analysis (default: False).')
    parser.add_argument('--clear_tmp', action='store_true', help='Clear temporary files after workflow.')
    parser.add_argument('--k', type=int, default=5, help='K-mer length (default: 5).')
    parser.add_argument('--k_range', action='store_true', help='Use range of k-mer lengths.')
    parser.add_argument('--remove_suffix', action='store_true', help='Remove suffix from genome names.')
    parser.add_argument('--one_gene', action='store_true', help='Include features with one gene.')
    parser.add_argument('--ignore_families', action='store_true', help='Ignore protein families.')
    parser.add_argument('--modeling', action='store_true', help='Run modeling workflow.')
    
    # SLURM-specific arguments
    parser.add_argument('--account', default='ac_mak', help='SLURM account (default: ac_mak).')
    parser.add_argument('--partition', default='lr7', help='SLURM partition (default: lr7).')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS (default: lr_normal).')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment name (default: phage_modeling).')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts but do not submit jobs.')
    
    # NEW: Resume capability
    parser.add_argument('--start_from_stage', type=int, default=1, choices=[1, 2, 3, 4, 5], 
                       help='Stage to start from (default: 1). Use this to resume from a specific stage.')
    
    args = parser.parse_args()
    
    # Create timestamped run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"slurm_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== SLURM Workflow Submission (5 Stages) ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output}")
    print(f"Starting from stage: {args.start_from_stage}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print()
    
    # Check which stages are already complete
    print("Checking for completed stages...")
    completed_stages = []
    for stage in range(1, 6):
        if check_stage_completion(args, stage):
            completed_stages.append(stage)
    
    if completed_stages:
        suggested_start = max(completed_stages) + 1
        if suggested_start <= 5:
            print(f"💡 Suggestion: You could start from stage {suggested_start} using --start_from_stage {suggested_start}")
            if args.start_from_stage == 1 and suggested_start > 1:
                print(f"⚠️  Warning: Starting from stage 1 will resubmit completed stages (wastes compute time)")
                print(f"   Consider using: --start_from_stage {suggested_start}")
    print()
    
    # Create all scripts
    print("Creating SLURM job scripts...")
    scripts = {}
    scripts[1] = create_stage1_clustering(args, run_dir)
    scripts[2] = create_stage2_feature_selection(args, run_dir, "PLACEHOLDER")
    scripts[3] = create_stage3_modeling(args, run_dir, "PLACEHOLDER")
    scripts[4] = create_stage4_predictive_proteins(args, run_dir, "PLACEHOLDER")
    scripts[5] = create_stage5_kmer_generation_and_modeling(args, run_dir, "PLACEHOLDER")  # CHANGED: Now includes modeling
    
    if args.dry_run:
        print("Dry run - scripts created but not submitted")
        print("Scripts:")
        for i, script in scripts.items():
            print(f"  Stage {i}: {script}")
        return
    
    # Change to run directory first
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    # Submit jobs with dependencies, starting from specified stage
    job_ids = {}
    
    print("Submitting jobs...")
    
    stage_names = {
        1: "Clustering", 
        2: "Feature Selection", 
        3: "Modeling", 
        4: "Predictive Proteins", 
        5: "K-mer Generation & Modeling"  # CHANGED: Updated name
    }
    
    # Submit starting stage
    start_stage = args.start_from_stage
    script_filename = f"stage{start_stage}_" + {
        1: "clustering.sh",
        2: "feature_selection.sh", 
        3: "modeling.sh",
        4: "predictive_proteins.sh",
        5: "kmer_complete.sh"  # CHANGED: New filename
    }[start_stage]
    
    job_ids[start_stage] = submit_job(script_filename)
    print(f"Stage {start_stage} ({stage_names[start_stage]}): {job_ids[start_stage]}")
    
    if job_ids[start_stage]:
        # Submit remaining stages with dependencies
        stage_filenames = {
            2: "stage2_feature_selection.sh",
            3: "stage3_modeling.sh", 
            4: "stage4_predictive_proteins.sh",
            5: "stage5_kmer_complete.sh"  # CHANGED: New filename
        }
        
        for stage in range(start_stage + 1, 6):  # CHANGED: Now goes to 6 instead of 7
            # Update dependency in script
            script_filename = stage_filenames[stage]
            with open(script_filename, 'r') as f:
                content = f.read()
            content = content.replace("PLACEHOLDER", job_ids[stage-1])
            with open(script_filename, 'w') as f:
                f.write(content)
            
            job_ids[stage] = submit_job(script_filename)
            print(f"Stage {stage} ({stage_names[stage]}): {job_ids[stage]}")
    
    # Change back to original directory
    os.chdir(original_dir)
    
    print(f"\n=== Job Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Stages submitted: {start_stage}-5")
    print("Monitor with: squeue -u $USER")
    print("View logs: tail -f logs/stage*_*.out")
    print("\n🎯 Key change: Stage 5 now includes BOTH k-mer generation AND modeling")
    print("Expected total runtime: 6-10 hours (5 sequential stages)")

if __name__ == "__main__":
    main()