#!/usr/bin/env python3
"""
SLURM workflow submission for phi_kmer_bootstrap_workflow_update.py
Creates a job array where each job handles one iteration of the k-mer bootstrap validation.
Iterations run in parallel, but steps within each iteration are sequential.
Takes the same arguments as phi_kmer_bootstrap_workflow_update.py
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

def create_kmer_bootstrap_job_array(args, run_dir):
    """Create SLURM job array script for k-mer bootstrap iterations"""
    
    # Get the absolute path to the original script directory for imports
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=kmer_bootstrap_validation
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_per_job}G
#SBATCH --time={args.time_limit}
#SBATCH --array=1-{args.n_iterations}
#SBATCH --output=logs/kmer_bootstrap_%A_%a.out
#SBATCH --error=logs/kmer_bootstrap_%A_%a.err

echo "=== K-mer Bootstrap Validation - Iteration $SLURM_ARRAY_TASK_ID ==="
echo "Job: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

# Run the k-mer bootstrap iteration logic directly

# Create a custom single-iteration script
python3 -c "
import sys
import os
sys.path.insert(0, '{script_dir}')

# Import the required functions
import pandas as pd
import shutil
import random
import logging
from Bio import SeqIO
from collections import defaultdict
from phage_modeling.workflows.full_workflow import run_full_workflow
from phage_modeling.workflows.kmer_assign_predict_workflow import kmer_assign_predict_workflow
from phage_modeling.workflows.kmer_analysis_workflow import kmer_analysis_workflow
from phage_modeling.workflows.kmer_table_workflow import is_fasta_empty
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

def detect_and_modify_duplicates(input_dir, output_dir, suffix='faa', strains_to_process=None, duplicate_all=False):
    '''Detect and resolve duplicate protein IDs by prefixing them with strain names.'''
    logging.info('Detecting and resolving duplicate protein IDs...')
    duplicate_found = False
    protein_id_tracker = defaultdict(set)
    modified_dir = os.path.join(output_dir, 'modified_AAs')
    
    # Check if modified directory already exists and contains files
    if os.path.exists(modified_dir) and os.listdir(modified_dir):
        logging.info(f'Modified AA directory already exists at {{modified_dir}}. Reusing it.')
        return modified_dir
    
    # First pass: detect duplicates by checking all files
    logging.info('Checking for duplicate protein IDs...')
    for file_name in os.listdir(input_dir):
        if file_name.endswith(suffix):
            strain_name = file_name.replace(f'.{{suffix}}', '')
            if strains_to_process and strain_name not in strains_to_process:
                continue
            
            file_path = os.path.join(input_dir, file_name)
            for record in SeqIO.parse(file_path, 'fasta'):
                if record.id in protein_id_tracker:
                    duplicate_found = True
                    logging.warning(f'Duplicate protein ID found: {{record.id}} in strain {{strain_name}}')
                    protein_id_tracker[record.id].add(strain_name)
                else:
                    protein_id_tracker[record.id].add(strain_name)
            
            if duplicate_found and not duplicate_all:
                break

    # If duplicates were found, modify the files
    if duplicate_found:
        os.makedirs(modified_dir, exist_ok=True)
        logging.info(f'Duplicate IDs detected. Creating modified files...')
        
        # Determine which files to process
        files_to_process = []
        for file_name in os.listdir(input_dir):
            if file_name.endswith(suffix):
                strain_name = file_name.replace(f'.{{suffix}}', '')
                if duplicate_all:  # Process all files when duplicate_all=True
                    if not strains_to_process or strain_name in strains_to_process:
                        files_to_process.append((file_name, strain_name))
                else:  # Otherwise only process files with duplicates
                    if not strains_to_process or strain_name in strains_to_process:
                        # Check if this strain contains any duplicated IDs
                        strain_has_duplicates = any(strain_name in strains for protein_id, strains in 
                                                  protein_id_tracker.items() if len(strains) > 1)
                        if strain_has_duplicates:
                            files_to_process.append((file_name, strain_name))
        
        # Process the determined files
        for file_name, strain_name in files_to_process:
            file_path = os.path.join(input_dir, file_name)
            modified_file_path = os.path.join(modified_dir, file_name)

            with open(modified_file_path, 'w') as modified_file:
                for record in SeqIO.parse(file_path, 'fasta'):
                    record.id = f'{{strain_name}}::{{record.id}}'
                    record.description = ''
                    SeqIO.write(record, modified_file, 'fasta')
        
        if files_to_process:
            logging.info(f'Resolved duplicate IDs and saved {{len(files_to_process)}} modified files.')
            return modified_dir
        else:
            logging.info('No files needed modification after duplicate check.')
            return input_dir
    else:
        logging.info('No duplicate protein IDs found. Using original files.')
        return input_dir

def get_full_strain_list(interaction_matrix, input_strain_dir, strain_column):
    logging.info(f'Reading interaction matrix: {{interaction_matrix}}')
    interaction_df = pd.read_csv(interaction_matrix)
    logging.info(f'Interaction matrix shape: {{interaction_df.shape}}')
    logging.info(f'Interaction matrix columns: {{list(interaction_df.columns)}}')
    
    if strain_column not in interaction_df.columns:
        logging.error(f'ERROR: Column \\'{{strain_column}}\\' not found in interaction matrix')
        logging.error(f'Available columns: {{list(interaction_df.columns)}}')
        return []
    
    strains_in_matrix = interaction_df[strain_column].unique()
    strains_in_matrix = [str(s) for s in strains_in_matrix]
    logging.info(f'Found {{len(strains_in_matrix)}} unique strains in interaction matrix')
    
    logging.info(f'Reading strain directory: {{input_strain_dir}}')
    strain_files = [f for f in os.listdir(input_strain_dir) if f.endswith('.faa')]
    strains_in_dir = ['.'.join(f.split('.')[:-1]) for f in strain_files]
    logging.info(f'Found {{len(strains_in_dir)}} strain files in directory')
    
    full_strain_list = list(set(strains_in_matrix).intersection(set(strains_in_dir)))
    logging.info(f'Intersection: {{len(full_strain_list)}} strains found in both matrix and directory')
    
    if len(full_strain_list) == 0:
        logging.error('ERROR: No strains found in both interaction matrix and input directory!')
        logging.error('This might be due to naming mismatches between the files and matrix.')
    
    return full_strain_list

def split_strains(full_strain_list, iteration, validation_percentage=0.1):
    if len(full_strain_list) == 0:
        logging.error('ERROR: Cannot split empty strain list!')
        return [], []
        
    logging.info(f'Splitting {{len(full_strain_list)}} strains for iteration {{iteration}}')
    random.seed(iteration)
    strain_list_copy = full_strain_list.copy()
    random.shuffle(strain_list_copy)
    
    split_index = int(len(strain_list_copy) * (1 - validation_percentage))
    modeling_strains = strain_list_copy[:split_index]
    validation_strains = strain_list_copy[split_index:]
    
    logging.info(f'Created {{len(modeling_strains)}} modeling strains and {{len(validation_strains)}} validation strains')
    return modeling_strains, validation_strains

def select_best_cutoff(output_dir):
    metrics_file = os.path.join(output_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    best_cutoff = metrics_df['cut_off'].values[0]
    return best_cutoff

# Get iteration number from SLURM array task ID
iteration = int(os.environ['SLURM_ARRAY_TASK_ID'])
logging.info(f'Processing iteration {{iteration}}')

# Configure logging for this iteration
output_dir = '{args.output_dir}'
iteration_output_dir = os.path.join(output_dir, f'iteration_{{iteration}}')
os.makedirs(iteration_output_dir, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(iteration_output_dir, f'kmer_bootstrap_iteration_{{iteration}}.log'))
    ]
)

median_predictions_file = os.path.join(iteration_output_dir, 'kmer_modeling', 'model_validation', 'predict_results', 'strain_median_predictions.csv')

# Skip if already completed
if os.path.exists(median_predictions_file):
    logging.info(f'Iteration {{iteration}} is already complete, skipping.')
    sys.exit(0)

logging.info(f'Starting iteration {{iteration}}...')
modeling_tmp_dir = os.path.join(iteration_output_dir, 'tmp')

# Get full strain list
full_strain_list = get_full_strain_list('{args.interaction_matrix}', '{args.input_strain_dir}', '{args.strain_column}')

if len(full_strain_list) == 0:
    logging.error('FATAL ERROR: No valid strains found. Cannot proceed with iteration.')
    sys.exit(1)

# Verify phage directory exists (required for this workflow)
if not os.path.exists('{args.input_phage_dir}'):
    logging.error('FATAL ERROR: Phage directory not found: {args.input_phage_dir}')
    sys.exit(1)

# Handle strain splits
modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')

if not os.path.exists(modeling_strains_path) or not os.path.exists(validation_strains_path):
    modeling_strains, validation_strains = split_strains(full_strain_list, iteration=iteration)
    
    if len(modeling_strains) == 0 or len(validation_strains) == 0:
        logging.error(f'ERROR: Empty strain lists created for iteration {{iteration}}')
        logging.error(f'Modeling strains: {{len(modeling_strains)}}, Validation strains: {{len(validation_strains)}}')
        sys.exit(1)
    
    # Save strain lists
    logging.info(f'Saving {{len(modeling_strains)}} modeling strains to {{modeling_strains_path}}')
    pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
    
    logging.info(f'Saving {{len(validation_strains)}} validation strains to {{validation_strains_path}}')
    pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
    
    # Verify files were written correctly
    test_modeling = pd.read_csv(modeling_strains_path)
    test_validation = pd.read_csv(validation_strains_path)
    logging.info(f'Verification - Modeling CSV has {{len(test_modeling)}} rows, Validation CSV has {{len(test_validation)}} rows')
    
else:
    logging.info('Strain lists already exist. Loading...')
    modeling_strains_df = pd.read_csv(modeling_strains_path)
    validation_strains_df = pd.read_csv(validation_strains_path)
    
    if len(modeling_strains_df) == 0 or len(validation_strains_df) == 0:
        logging.error(f'ERROR: Existing strain CSV files are empty!')
        logging.error(f'Modeling CSV: {{len(modeling_strains_df)}} rows, Validation CSV: {{len(validation_strains_df)}} rows')
        # Regenerate the files
        modeling_strains, validation_strains = split_strains(full_strain_list, iteration=iteration)
        pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
        pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
    else:
        modeling_strains = modeling_strains_df['strain'].tolist()
        validation_strains = validation_strains_df['strain'].tolist()
        logging.info(f'Loaded {{len(modeling_strains)}} modeling strains and {{len(validation_strains)}} validation strains')

# Detect and modify duplicates
modified_input_dir = detect_and_modify_duplicates(
    input_dir='{args.input_strain_dir}',
    output_dir=os.path.join(iteration_output_dir, 'strain'),
    suffix='faa',
    strains_to_process=modeling_strains + validation_strains,
    duplicate_all={args.duplicate_all}
)

if modified_input_dir != '{args.input_strain_dir}':
    logging.info(f'Using modified AA directory: {{modified_input_dir}}')

# Step 1: Run full workflow (includes k-mer modeling) with modeling strains
metrics_file = os.path.join(iteration_output_dir, 'kmer_modeling/modeling/modeling_results/model_performance/model_performance_metrics.csv')

if not os.path.exists(metrics_file):
    logging.info('Running full k-mer workflow...')
    run_full_workflow(
        input_strain=modified_input_dir,
        input_phage='{args.input_phage_dir}',
        phenotype_matrix='{args.interaction_matrix}',
        output=iteration_output_dir,
        clustering_dir=None,
        min_seq_id=0.4,
        coverage=0.8,
        sensitivity=7.5,
        suffix='faa',
        strain_list=modeling_strains_path,
        phage_list='{args.interaction_matrix}',
        strain_column='{args.strain_column}',
        phage_column='phage',
        source_strain='strain',
        source_phage='phage',
        compare=False,
        num_features='none',
        filter_type='strain',
        num_runs_fs={args.num_runs_fs},
        num_runs_modeling={args.num_runs_modeling},
        sample_column='{args.strain_column}',
        phenotype_column='interaction',
        method='rfe',
        annotation_table_path=None,
        protein_id_col='protein_ID',
        task_type='classification',
        max_features='none',
        max_ram={args.max_ram},
        threads={args.threads},
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
        use_shap=False,
        clear_tmp=False,
        k=5,
        k_range=False,
        remove_suffix=False,
        one_gene=False,
        ignore_families=False,
        modeling=True,
        use_feature_clustering={args.use_feature_clustering},
        feature_cluster_method='{args.feature_cluster_method}',
        feature_n_clusters={args.feature_n_clusters},
        feature_min_cluster_presence={args.feature_min_cluster_presence}
    )
else:
    logging.info('K-mer modeling results already exist. Skipping full workflow.')

# Step 2: Generate filtered_kmers.csv using kmer_analysis_workflow
best_cutoff = select_best_cutoff(os.path.join(iteration_output_dir, 'kmer_modeling', 'modeling'))
model_dir = os.path.join(iteration_output_dir, 'kmer_modeling', 'modeling', f'modeling_results', f'{{best_cutoff}}')
cutoff_num = best_cutoff.split('_')[1]

feature_file_path = os.path.join(iteration_output_dir, 'kmer_modeling', 'modeling', 'feature_selection', 'filtered_feature_tables', f'select_feature_table_cutoff_{{cutoff_num}}.csv')
feature2cluster_path = os.path.join(iteration_output_dir, 'kmer_modeling', 'feature_tables', 'selected_features.csv')
protein_families_file = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'predictive_proteins', 'strain_predictive_feature_overview.csv')
aa_sequence_file = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'predictive_proteins', 'strain_predictive_AA_seqs.faa')
kmer_analysis_dir = os.path.join(iteration_output_dir, 'kmer_modeling', 'kmer_analysis')

filtered_kmers_path = os.path.join(kmer_analysis_dir, 'strain', 'filtered_kmers.csv')

if is_fasta_empty(aa_sequence_file):
    logging.warning('No predictive sequences found. Skipping k-mer analysis.')
elif not os.path.exists(filtered_kmers_path):
    logging.info('Running k-mer analysis workflow...')
    kmer_analysis_workflow(
        aa_sequence_file=aa_sequence_file,
        feature_file_path=feature_file_path,
        feature2cluster_path=feature2cluster_path,
        protein_families_file=protein_families_file,
        output_dir=kmer_analysis_dir,
        feature_type='strain'
    )
else:
    logging.info('Filtered kmers already exist. Skipping kmer analysis workflow.')

# Step 3: Predict interactions for validation strains using kmer-based workflow
validation_output_dir = os.path.join(iteration_output_dir, 'kmer_modeling', 'model_validation')
validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')
os.makedirs(validation_output_dir, exist_ok=True)

# Create symbolic links to the necessary directories
protein_tmp_dir = os.path.abspath(os.path.join(iteration_output_dir, 'model_validation', 'tmp'))
validation_tmp_dir = os.path.abspath(os.path.join(validation_output_dir, 'tmp'))

# Check if the protein tmp directory exists
if os.path.exists(protein_tmp_dir):
    # If validation tmp directory already exists, don't try to create it
    if os.path.exists(validation_tmp_dir):
        logging.info(f'Validation tmp directory already exists at {{validation_tmp_dir}}, not creating symlink.')
    else:
        # Make sure the parent directory exists
        os.makedirs(os.path.dirname(validation_tmp_dir), exist_ok=True)
        
        # Check if parent directory exists after creation
        if os.path.exists(os.path.dirname(validation_tmp_dir)):
            try:
                os.symlink(protein_tmp_dir, validation_tmp_dir, target_is_directory=True)
                logging.info(f'Created symlink from {{protein_tmp_dir}} to {{validation_tmp_dir}}')
            except FileExistsError:
                logging.info(f'Validation tmp directory already exists at {{validation_tmp_dir}}, not creating symlink.')
            except Exception as e:
                logging.error(f'Error creating symlink: {{e}}')
                # If symlink creation fails, create a new directory
                os.makedirs(validation_tmp_dir, exist_ok=True)
                logging.info(f'Created new directory {{validation_tmp_dir}} instead of symlink')
        else:
            logging.error(f'Failed to create parent directory for {{validation_tmp_dir}}')

if is_fasta_empty(aa_sequence_file):
    logging.warning('No predictive strain sequences found. Creating empty strain feature table and proceeding directly to prediction.')
    
    # Create empty strain feature table
    validation_strains_df = pd.read_csv(validation_strains_path)
    strain_list = validation_strains_df['strain'].tolist()
    
    # Create empty strain feature table with just the strain column
    empty_feature_table = pd.DataFrame({{'strain': strain_list}})
    
    # Save the empty feature table in the expected location
    assign_results_dir = os.path.join(validation_output_dir, 'assign_results')
    os.makedirs(assign_results_dir, exist_ok=True)
    empty_feature_table.to_csv(os.path.join(assign_results_dir, 'strain_combined_feature_table.csv'), index=False)
    
    # Run just the prediction part
    predict_output_dir = os.path.join(validation_output_dir, 'predict_results')
    os.makedirs(predict_output_dir, exist_ok=True)
    
    run_prediction_workflow(
        input_dir=assign_results_dir,
        phage_feature_table_path=os.path.join(iteration_output_dir, 'kmer_modeling', 'feature_tables', 'phage_final_feature_table.csv'),
        model_dir=model_dir,
        output_dir=predict_output_dir
    )
else:
    # Proceed with the kmer assignment and prediction workflow
    logging.info('Running kmer assignment and prediction workflow...')
    kmer_assign_predict_workflow(
        input_dir=modified_input_dir,
        genome_list=validation_strains_path,
        mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
        clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
        feature_map=os.path.join(iteration_output_dir, 'kmer_modeling', 'feature_tables', 'selected_features.csv'),
        filtered_kmers=filtered_kmers_path,
        aa_sequence_file=aa_sequence_file,
        tmp_dir=validation_tmp_dir,
        output_dir=validation_output_dir,
        model_dir=model_dir,
        phage_feature_table_path=os.path.join(iteration_output_dir, 'kmer_modeling', 'feature_tables', 'phage_final_feature_table.csv'),
        threads={args.threads},
        genome_type='strain',
        sensitivity=7.5,
        coverage=0.8,
        min_seq_id=0.4,
        threshold=0.01
    )

logging.info(f'Iteration {{iteration}} completed successfully')
"

echo "Iteration $SLURM_ARRAY_TASK_ID completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "kmer_bootstrap_job_array.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_final_aggregation_job(args, run_dir, dependency):
    """Create job to aggregate all iteration results"""
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=kmer_bootstrap_aggregate
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

echo "=== K-mer Bootstrap Results Aggregation ==="
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
    median_predictions_file = os.path.join(iteration_output_dir, 'kmer_modeling', 'model_validation', 'predict_results', 'strain_median_predictions.csv')
    
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
    final_predictions.to_csv(os.path.join(output_dir, 'kmer_final_predictions.csv'), index=False)
    print(f'K-mer final predictions saved with {{len(final_predictions)}} total predictions from {{completed_iterations}} iterations')
    
    # Generate summary statistics
    summary_stats = final_predictions.groupby(['strain', 'phage']).agg({{
        'prediction': ['mean', 'std', 'count']
    }}).round(4)
    summary_stats.columns = ['mean_prediction', 'std_prediction', 'n_iterations']
    summary_stats = summary_stats.reset_index()
    summary_stats.to_csv(os.path.join(output_dir, 'kmer_prediction_summary.csv'), index=False)
    print(f'K-mer prediction summary saved with {{len(summary_stats)}} unique strain-phage pairs')
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
    parser = argparse.ArgumentParser(description="Submit k-mer bootstrap validation workflow as SLURM job array")
    
    # Copy ALL arguments from phi_kmer_bootstrap_workflow_update.py
    parser.add_argument('--input_strain_dir', type=str, required=True, help="Directory containing strain FASTA files.")
    parser.add_argument('--input_phage_dir', type=str, required=True, help="Directory containing phage FASTA files.")
    parser.add_argument('--interaction_matrix', type=str, required=True, help="Path to the interaction matrix.")
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
    parser.add_argument('--min_cluster_presence', type=int, default=2, help="Minimum number of clusters a feature must appear in.")
    parser.add_argument('--duplicate_all', action='store_true', help="Duplicate all genomes in the feature table for predictions.")
    parser.add_argument('--max_ram', type=float, default=40, help="Maximum RAM usage in GB.")
    
    # Pre-processing feature clustering parameters
    parser.add_argument('--use_feature_clustering', action='store_true', help="Enable pre-processing cluster-based feature filtering.")
    parser.add_argument('--feature_cluster_method', type=str, default='hierarchical', choices=['hierarchical'], help="Pre-processing clustering method ('hierarchical' only for now).")
    parser.add_argument('--feature_n_clusters', type=int, default=20, help="Number of clusters for pre-processing feature clustering.")
    parser.add_argument('--feature_min_cluster_presence', type=int, default=2, help="Min clusters a feature must appear in during pre-processing.")
    
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
    run_dir = f"kmer_bootstrap_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== K-mer Bootstrap Validation SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of iterations: {args.n_iterations}")
    print(f"Threads per job: {args.threads}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    print(f"Duplicate all: {args.duplicate_all}")
    print(f"Filter by cluster presence: {args.filter_by_cluster_presence}")
    print(f"Use feature clustering: {args.use_feature_clustering}")
    if args.use_feature_clustering:
        print(f"Feature cluster method: {args.feature_cluster_method}")
        print(f"Feature n_clusters: {args.feature_n_clusters}")
        print(f"Feature min presence: {args.feature_min_cluster_presence}")
    
    # Resource validation
    if args.mem_per_job < 32:
        print(f"WARNING: {args.mem_per_job}GB memory may be insufficient for complex datasets")
    if args.n_iterations > 50:
        print(f"WARNING: {args.n_iterations} iterations will create many parallel jobs")
    
    print()
    
    # Create scripts
    print("Creating SLURM job scripts...")
    bootstrap_script = create_kmer_bootstrap_job_array(args, run_dir)
    aggregate_script = create_final_aggregation_job(args, run_dir, "PLACEHOLDER")
    
    if args.dry_run:
        print("Dry run - scripts created but not submitted")
        print("Scripts:")
        print(f"  K-mer bootstrap array: {bootstrap_script}")
        print(f"  Aggregation: {aggregate_script}")
        return
    
    # Change to run directory
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    # Submit bootstrap job array
    print("Submitting k-mer bootstrap job array...")
    bootstrap_job_id = submit_job("kmer_bootstrap_job_array.sh")
    print(f"K-mer bootstrap job array: {bootstrap_job_id}")
    
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
    print("  squeue -u $USER --name=kmer_bootstrap_validation  # Just k-mer bootstrap jobs")
    print("  tail -f logs/kmer_bootstrap_*.out")
    print(f"\nResults will be in:")
    print(f"  Individual iterations: {args.output_dir}/iteration_N/")
    print(f"  K-mer modeling: {args.output_dir}/iteration_N/kmer_modeling/")
    print(f"  Final aggregated: {args.output_dir}/kmer_final_predictions.csv")
    print(f"  Summary stats: {args.output_dir}/kmer_prediction_summary.csv")

if __name__ == "__main__":
    main()