#!/usr/bin/env python3
"""
SLURM workflow submission for phi_bootstrap_workflow_update.py
Creates a job array where each job handles one iteration of the bootstrap validation.
Iterations run in parallel, but steps within each iteration are sequential.
Takes the same arguments as phi_bootstrap_workflow_update.py
"""

import os
import sys
import argparse
import subprocess
import time

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

def create_bootstrap_job_array(args, run_dir):
    """Create SLURM job array script for bootstrap iterations"""
    
    # Get the absolute path to the original script directory for imports
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=bootstrap_validation
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_per_job}G
#SBATCH --time={args.time_limit}
#SBATCH --array=1-{args.n_iterations}
#SBATCH --output=logs/bootstrap_%A_%a.out
#SBATCH --error=logs/bootstrap_%A_%a.err

echo "=== Bootstrap Validation - Iteration $SLURM_ARRAY_TASK_ID ==="
echo "Job: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

# Run the bootstrap iteration logic directly

# Create a custom single-iteration script
python3 -c "
import sys
import os
sys.path.insert(0, '{script_dir}')

# Import the required functions
import pandas as pd
import shutil
import random
from phage_modeling.workflows.protein_family_workflow import run_protein_family_workflow
from phage_modeling.workflows.assign_predict_workflow import assign_predict_workflow

def get_full_strain_list(interaction_matrix, input_strain_dir, strain_column):
    print(f'Reading interaction matrix: {{interaction_matrix}}')
    interaction_df = pd.read_csv(interaction_matrix)
    print(f'Interaction matrix shape: {{interaction_df.shape}}')
    print(f'Interaction matrix columns: {{list(interaction_df.columns)}}')
    
    if strain_column not in interaction_df.columns:
        print(f'ERROR: Column \\'{{strain_column}}\\' not found in interaction matrix')
        print(f'Available columns: {{list(interaction_df.columns)}}')
        return []
    
    strains_in_matrix = interaction_df[strain_column].unique()
    strains_in_matrix = [str(s) for s in strains_in_matrix]
    print(f'Found {{len(strains_in_matrix)}} unique strains in interaction matrix')
    print(f'First 10 strains from matrix: {{list(strains_in_matrix[:10])}}')
    
    print(f'Reading strain directory: {{input_strain_dir}}')
    strain_files = [f for f in os.listdir(input_strain_dir) if f.endswith('.faa')]
    strains_in_dir = ['.'.join(f.split('.')[:-1]) for f in strain_files]
    print(f'Found {{len(strains_in_dir)}} strain files in directory')
    print(f'First 10 strains from directory: {{strains_in_dir[:10]}}')
    
    full_strain_list = list(set(strains_in_matrix).intersection(set(strains_in_dir)))
    print(f'Intersection: {{len(full_strain_list)}} strains found in both matrix and directory')
    
    if len(full_strain_list) == 0:
        print('ERROR: No strains found in both interaction matrix and input directory!')
        print('This might be due to naming mismatches between the files and matrix.')
    
    return full_strain_list

def split_strains(full_strain_list, iteration, validation_percentage=0.1):
    if len(full_strain_list) == 0:
        print('ERROR: Cannot split empty strain list!')
        return [], []
        
    print(f'Splitting {{len(full_strain_list)}} strains for iteration {{iteration}}')
    random.seed(iteration)
    strain_list_copy = full_strain_list.copy()  # Don't modify original list
    random.shuffle(strain_list_copy)
    
    split_index = int(len(strain_list_copy) * (1 - validation_percentage))
    modeling_strains = strain_list_copy[:split_index]
    validation_strains = strain_list_copy[split_index:]
    
    print(f'Created {{len(modeling_strains)}} modeling strains and {{len(validation_strains)}} validation strains')
    return modeling_strains, validation_strains

def select_best_cutoff(output_dir):
    metrics_file = os.path.join(output_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    best_cutoff = metrics_df['cut_off'].values[0]
    return best_cutoff

# Get iteration number from SLURM array task ID
iteration = int(os.environ['SLURM_ARRAY_TASK_ID'])
print(f'Processing iteration {{iteration}}')

# Parameters from command line arguments
output_dir = '{args.output_dir}'
iteration_output_dir = os.path.join(output_dir, f'iteration_{{iteration}}')
median_predictions_file = os.path.join(iteration_output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')

# Skip if already completed
if os.path.exists(median_predictions_file):
    print(f'Iteration {{iteration}} is already complete, skipping.')
    sys.exit(0)

print(f'Starting iteration {{iteration}}...')
os.makedirs(iteration_output_dir, exist_ok=True)
modeling_tmp_dir = os.path.join(iteration_output_dir, 'tmp')

# Get full strain list
full_strain_list = get_full_strain_list('{args.interaction_matrix}', '{args.input_strain_dir}', '{args.strain_column}')

if len(full_strain_list) == 0:
    print('FATAL ERROR: No valid strains found. Cannot proceed with iteration.')
    sys.exit(1)

# Verify phage directory exists (required for this workflow)
if not os.path.exists('{args.input_phage_dir}'):
    print('FATAL ERROR: Phage directory not found: {args.input_phage_dir}')
    sys.exit(1)

# Handle clustering directory and strain splits
clustering_dir = {repr(args.clustering_dir) if args.clustering_dir else 'None'}
if clustering_dir:
    modeling_strains_old = os.path.join(clustering_dir, f'iteration_{{iteration}}', 'modeling_strains.csv')
    validation_strains_old = os.path.join(clustering_dir, f'iteration_{{iteration}}', 'validation_strains.csv')
    
    modeling_strains_new = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    validation_strains_new = os.path.join(iteration_output_dir, 'validation_strains.csv')
    
    if not os.path.exists(modeling_strains_new):
        os.symlink(modeling_strains_old, modeling_strains_new)
    if not os.path.exists(validation_strains_new):
        os.symlink(validation_strains_old, validation_strains_new)
        
    modeling_strains = pd.read_csv(modeling_strains_new)['strain'].tolist()
    validation_strains = pd.read_csv(validation_strains_new)['strain'].tolist()
else:
    modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')
    
    if not os.path.exists(modeling_strains_path) or not os.path.exists(validation_strains_path):
        modeling_strains, validation_strains = split_strains(full_strain_list, iteration=iteration)
        
        if len(modeling_strains) == 0 or len(validation_strains) == 0:
            print(f'ERROR: Empty strain lists created for iteration {{iteration}}')
            print(f'Modeling strains: {{len(modeling_strains)}}, Validation strains: {{len(validation_strains)}}')
            sys.exit(1)
        
        # Save strain lists
        print(f'Saving {{len(modeling_strains)}} modeling strains to {{modeling_strains_path}}')
        pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
        
        print(f'Saving {{len(validation_strains)}} validation strains to {{validation_strains_path}}')
        pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
        
        # Verify files were written correctly
        test_modeling = pd.read_csv(modeling_strains_path)
        test_validation = pd.read_csv(validation_strains_path)
        print(f'Verification - Modeling CSV has {{len(test_modeling)}} rows, Validation CSV has {{len(test_validation)}} rows')
        
    else:
        print('Strain lists already exist. Loading...')
        modeling_strains_df = pd.read_csv(modeling_strains_path)
        validation_strains_df = pd.read_csv(validation_strains_path)
        
        if len(modeling_strains_df) == 0 or len(validation_strains_df) == 0:
            print(f'ERROR: Existing strain CSV files are empty!')
            print(f'Modeling CSV: {{len(modeling_strains_df)}} rows, Validation CSV: {{len(validation_strains_df)}} rows')
            # Regenerate the files
            modeling_strains, validation_strains = split_strains(full_strain_list, iteration=iteration)
            pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
            pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
        else:
            modeling_strains = modeling_strains_df['strain'].tolist()
            validation_strains = validation_strains_df['strain'].tolist()
            print(f'Loaded {{len(modeling_strains)}} modeling strains and {{len(validation_strains)}} validation strains')

# Step 1: Run protein_family_workflow with modeling strains
metrics_file = os.path.join(iteration_output_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
iteration_clustering_dir = None
if clustering_dir:
    iteration_clustering_dir = os.path.join(clustering_dir, f'iteration_{{iteration}}')

if not os.path.exists(metrics_file):
    print('Running protein family workflow...')
    run_protein_family_workflow(
        input_path_strain='{args.input_strain_dir}',
        input_path_phage='{args.input_phage_dir}',
        clustering_dir=iteration_clustering_dir,
        phenotype_matrix='{args.interaction_matrix}',
        output_dir=iteration_output_dir,
        tmp_dir=modeling_tmp_dir,
        min_seq_id=0.4,
        coverage=0.8,
        sensitivity=7.5,
        threads={args.threads},
        phenotype_column='interaction',
        phage_column='phage',
        strain_list=os.path.join(iteration_output_dir, 'modeling_strains.csv'),
        phage_list='{args.interaction_matrix}',
        filter_type='strain',
        num_runs_fs={args.num_runs_fs},
        num_runs_modeling={args.num_runs_modeling},
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
        bootstrapping={args.bootstrapping},
        max_ram={args.max_ram}
    )
else:
    print('Modeling results already exist. Skipping protein_family_workflow.')

# Step 2: Get the cutoff with the highest MCC
best_cutoff = select_best_cutoff(iteration_output_dir)
model_dir = os.path.join(iteration_output_dir, f'modeling_results', f'{{best_cutoff}}')

# Step 3: Predict interactions for validation strains
validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
os.makedirs(validation_output_dir, exist_ok=True)
validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')

# Check for modified AA directory
modified_aa_dir = os.path.join(iteration_output_dir, 'strain', 'modified_AAs', 'strain')
input_strain_dir = modified_aa_dir if os.path.exists(modified_aa_dir) else '{args.input_strain_dir}'

select_feature_table = os.path.join(iteration_output_dir, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_{{best_cutoff}}.csv')

print('Running prediction workflow...')
assign_predict_workflow(
    input_dir=input_strain_dir,
    genome_list=os.path.join(iteration_output_dir, 'validation_strains.csv'),
    mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
    clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
    feature_map=os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv'),
    tmp_dir=validation_tmp_dir,
    suffix='faa',
    model_dir=model_dir,
    feature_table=select_feature_table,
    phage_feature_table_path=os.path.join(iteration_output_dir, 'phage', 'features', 'feature_table.csv'),
    output_dir=validation_output_dir,
    threads={args.threads},
    genome_type='strain',
    sensitivity=7.5,
    coverage=0.8,
    min_seq_id=0.4,
    duplicate_all={args.duplicate_all}
)

print(f'Iteration {{iteration}} completed successfully')
"

echo "Iteration $SLURM_ARRAY_TASK_ID completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "bootstrap_job_array.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_final_aggregation_job(args, run_dir, dependency):
    """Create job to aggregate all iteration results"""
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=bootstrap_aggregate
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --dependency=afterok:{dependency}
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err

echo "=== Bootstrap Results Aggregation ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import pandas as pd
import os

# Aggregate final predictions from all iterations
output_dir = '{args.output_dir}'
final_predictions = pd.DataFrame()

print('Aggregating results from {args.n_iterations} iterations...')
completed_iterations = 0

for i in range(1, {args.n_iterations} + 1):
    iteration_output_dir = os.path.join(output_dir, f'iteration_{{i}}')
    median_predictions_file = os.path.join(iteration_output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')
    
    if os.path.exists(median_predictions_file):
        median_predictions = pd.read_csv(median_predictions_file)
        median_predictions['iteration'] = i
        final_predictions = pd.concat([final_predictions, median_predictions], ignore_index=True)
        completed_iterations += 1
        print(f'Added results from iteration {{i}}')
    else:
        print(f'Warning: Results missing for iteration {{i}}')

# Save final concatenated predictions
if len(final_predictions) > 0:
    final_predictions.to_csv(os.path.join(output_dir, 'final_predictions.csv'), index=False)
    print(f'Final predictions saved with {{len(final_predictions)}} total predictions from {{completed_iterations}} iterations')
    
    # Generate summary statistics
    summary_stats = final_predictions.groupby(['strain', 'phage']).agg({{
        'prediction': ['mean', 'std', 'count']
    }}).round(4)
    summary_stats.columns = ['mean_prediction', 'std_prediction', 'n_iterations']
    summary_stats = summary_stats.reset_index()
    summary_stats.to_csv(os.path.join(output_dir, 'prediction_summary.csv'), index=False)
    print(f'Prediction summary saved with {{len(summary_stats)}} unique strain-phage pairs')
else:
    print('ERROR: No results found to aggregate!')
    exit(1)
"

echo "Aggregation completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "aggregate_results.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Submit bootstrap validation workflow as SLURM job array")
    
    # Copy ALL arguments from phi_bootstrap_workflow_update.py
    parser.add_argument('--input_strain_dir', type=str, required=True, help="Directory containing strain FASTA files.")
    parser.add_argument('--input_phage_dir', type=str, required=True, help="Directory containing phage FASTA files.")
    parser.add_argument('--interaction_matrix', type=str, required=True, help="Path to the interaction matrix.")
    parser.add_argument('--clustering_dir', type=str, help="Directory containing strain clustering results.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--n_iterations', type=int, default=10, help="Number of iterations.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads to use.")
    parser.add_argument('--strain_column', type=str, default='strain', help="Column in the interaction matrix containing strain names.")
    parser.add_argument('--num_runs_fs', type=int, default=25, help="Number of runs for feature selection.")
    parser.add_argument('--num_runs_modeling', type=int, default=50, help="Number of runs for modeling.")
    parser.add_argument('--use_dynamic_weights', action='store_true', help="Use dynamic weights for feature selection.")
    parser.add_argument('--weights_method', type=str, default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help="Method for calculating dynamic weights.")
    parser.add_argument('--use_clustering', action='store_true', help="Use clustering results for feature selection.")
    parser.add_argument('--cluster_method', type=str, default='hdbscan', choices=['hdbscan', 'hierarchical'], help="Clustering method to use.")
    parser.add_argument('--n_clusters', type=int, default=20, help="Number of clusters for clustering feature selection.")
    parser.add_argument('--min_cluster_size', type=int, default=2, help="Minimum cluster size for clustering feature selection.")
    parser.add_argument('--min_samples', type=int, help="Minimum number of samples for clustering feature selection.")
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help="Epsilon value for clustering feature selection.")
    parser.add_argument('--check_feature_presence', action='store_true', help="Check for feature presence for train-test splits.")
    parser.add_argument('--filter_by_cluster_presence', action='store_true', help="Filter features by cluster presence across train/test splits.")
    parser.add_argument('--min_cluster_presence', type=float, default=2, help="Minimum fraction of clusters that must contain the feature.")
    parser.add_argument('--duplicate_all', action='store_true', help="Duplicate all genomes in the feature table for predictions.")
    parser.add_argument('--max_ram', type=float, default=40, help="Maximum RAM usage in GB.")
    parser.add_argument('--bootstrapping', action='store_true', help="Enable bootstrapping for feature selection and modeling.")
    
    # SLURM-specific arguments
    parser.add_argument('--account', default='ac_mak', help='SLURM account (default: ac_mak).')
    parser.add_argument('--partition', default='lr7', help='SLURM partition (default: lr7).')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS (default: lr_normal).')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment name (default: phage_modeling).')
    parser.add_argument('--mem_per_job', type=int, default=60, help='Memory per job in GB (default: 60).')
    parser.add_argument('--time_limit', default='12:00:00', help='Time limit per job (default: 12:00:00).')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts but do not submit jobs.')
    
    args = parser.parse_args()
    
    # Validate critical inputs exist
    if not os.path.exists(args.input_strain_dir):
        print(f"Error: Input strain directory not found: {args.input_strain_dir}")
        return 1
    if not os.path.exists(args.input_phage_dir):
        print(f"Error: Input phage directory not found: {args.input_phage_dir}")
        return 1
    if not os.path.exists(args.interaction_matrix):
        print(f"Error: Interaction matrix not found: {args.interaction_matrix}")
        return 1
    if args.clustering_dir and not os.path.exists(args.clustering_dir):
        print(f"Error: Clustering directory not found: {args.clustering_dir}")
        return 1
    
    # Check if output directory parent exists
    output_parent = os.path.dirname(args.output_dir)
    if output_parent and not os.path.exists(output_parent):
        print(f"Error: Output parent directory not found: {output_parent}")
        return 1
    
    # Validate interaction matrix has expected column
    try:
        import pandas as pd
        df = pd.read_csv(args.interaction_matrix)
        if args.strain_column not in df.columns:
            print(f"Error: Column '{args.strain_column}' not found in interaction matrix")
            print(f"Available columns: {list(df.columns)}")
            return 1
        print(f"Interaction matrix validation: ✓ Found {len(df)} rows with column '{args.strain_column}'")
    except Exception as e:
        print(f"Error reading interaction matrix: {e}")
        return 1
    
    # Check strain files exist
    strain_files = [f for f in os.listdir(args.input_strain_dir) if f.endswith('.faa')]
    if len(strain_files) == 0:
        print(f"Error: No .faa files found in strain directory: {args.input_strain_dir}")
        return 1
    print(f"Strain directory validation: ✓ Found {len(strain_files)} .faa files")
    
    # Check phage files exist  
    phage_files = [f for f in os.listdir(args.input_phage_dir) if f.endswith('.faa')]
    if len(phage_files) == 0:
        print(f"Error: No .faa files found in phage directory: {args.input_phage_dir}")
        return 1
    print(f"Phage directory validation: ✓ Found {len(phage_files)} .faa files")
    
    # Create timestamped run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"bootstrap_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== Bootstrap Validation SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of iterations: {args.n_iterations}")
    print(f"Threads per job: {args.threads}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    
    # Resource validation
    if args.mem_per_job < 32:
        print(f"WARNING: {args.mem_per_job}GB memory may be insufficient for complex datasets")
    if args.n_iterations > 50:
        print(f"WARNING: {args.n_iterations} iterations will create many parallel jobs")
    
    print()
    
    # Create scripts
    print("Creating SLURM job scripts...")
    bootstrap_script = create_bootstrap_job_array(args, run_dir)
    aggregate_script = create_final_aggregation_job(args, run_dir, "PLACEHOLDER")
    
    if args.dry_run:
        print("Dry run - scripts created but not submitted")
        print("Scripts:")
        print(f"  Bootstrap array: {bootstrap_script}")
        print(f"  Aggregation: {aggregate_script}")
        return
    
    # Change to run directory
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    # Submit bootstrap job array
    print("Submitting bootstrap job array...")
    bootstrap_job_id = submit_job("bootstrap_job_array.sh")
    print(f"Bootstrap job array: {bootstrap_job_id}")
    
    if bootstrap_job_id:
        # Update aggregation script dependency and submit
        with open("aggregate_results.sh", 'r') as f:
            content = f.read()
        content = content.replace("PLACEHOLDER", bootstrap_job_id)
        with open("aggregate_results.sh", 'w') as f:
            f.write(content)
        
        aggregate_job_id = submit_job("aggregate_results.sh")
        print(f"Aggregation job: {aggregate_job_id}")
    
    # Change back to original directory
    os.chdir(original_dir)
    
    print(f"\n=== Job Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Bootstrap array: {args.n_iterations} parallel jobs")
    print(f"Expected total runtime: {args.time_limit} per iteration (parallel)")
    print("\nMonitor with:")
    print("  squeue -u $USER")
    print("  squeue -u $USER --name=bootstrap_validation  # Just bootstrap jobs")
    print("  tail -f logs/bootstrap_*.out")
    print(f"\nResults will be in:")
    print(f"  Individual iterations: {args.output_dir}/iteration_N/")
    print(f"  Final aggregated: {args.output_dir}/final_predictions.csv")
    print(f"  Summary stats: {args.output_dir}/prediction_summary.csv")

if __name__ == "__main__":
    main()