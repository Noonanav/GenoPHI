#!/usr/bin/env python3
"""
SLURM workflow submission for k-mer-based bootstrap cross-validation.
Creates a job array where each job handles one iteration of the bootstrap validation.
Iterations run in parallel, but steps within each iteration are sequential.
Uses separate strain and phage k-mer processing with proper merging.
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
    """Create SLURM job array script for k-mer bootstrap iterations"""
    
    # Get the absolute path to the original script directory for imports
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=kmer_bootstrap
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
from phage_modeling.workflows.kmer_table_workflow import run_kmer_table_workflow
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow
from Bio import SeqIO
import multiprocessing
from functools import partial
import gc
from collections import defaultdict

def generate_protein_mapping_csv(fasta_dir, output_csv, genome_col_name, file_extension='.faa'):
    '''Generate protein-to-genome mapping CSV from FASTA directory'''
    data = []
    
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(file_extension):
            genome_name = fasta_file.replace(file_extension, '')
            fasta_path = os.path.join(fasta_dir, fasta_file)
            
            for record in SeqIO.parse(fasta_path, 'fasta'):
                data.append({{
                    'protein_ID': record.id,
                    genome_col_name: genome_name
                }})
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f'Generated {{output_csv}} with {{len(df)}} protein mappings')
    return output_csv

def detect_duplicate_ids(input_path, suffix='faa', strains_to_process=None, input_type='directory'):
    '''Detects duplicate protein IDs across relevant .faa files in the input directory.'''
    print('Detecting duplicate protein IDs...')
    duplicate_found = False
    protein_id_tracker = defaultdict(set)

    if input_type == 'directory':
        file_list = os.listdir(input_path)
    elif input_type == 'file':
        file_list = [input_path]
    else:
        print(f'Invalid input type: {{input_type}}')
        return duplicate_found

    # Identify and track protein IDs in specified strains
    for file_name in file_list:
        if file_name.endswith(suffix):
            strain_name = file_name.replace(f'.{{suffix}}', '')
            if strains_to_process and strain_name not in strains_to_process:
                continue  # Skip strains not in the strain list
            
            file_path = os.path.join(input_path, file_name)
            for record in SeqIO.parse(file_path, 'fasta'):
                if record.id in protein_id_tracker:
                    duplicate_found = True
                    print(f'Duplicate protein ID found: {{record.id}} in strain {{strain_name}}')
                    break
                protein_id_tracker[record.id].add(strain_name)
    
    return duplicate_found

def modify_duplicate_ids(input_path, output_dir, suffix='faa', strains_to_process=None, strain_column='strain'):
    '''Detects duplicate protein IDs and modifies all protein IDs in relevant .faa files 
    by prefixing them with genome names to ensure uniqueness.'''
    print(f'Duplicate protein IDs found; modifying protein IDs and saving to {{output_dir}}/modified_AAs/{{strain_column}}')
    modified_file_dir = os.path.join(output_dir, 'modified_AAs', strain_column)
    os.makedirs(modified_file_dir, exist_ok=True)

    for file_name in os.listdir(input_path):
        if file_name.endswith(suffix):
            strain_name = file_name.replace(f'.{{suffix}}', '')
            if strains_to_process and strain_name not in strains_to_process:
                continue
            
            file_path = os.path.join(input_path, file_name)
            modified_file_path = os.path.join(modified_file_dir, file_name)
            
            with open(modified_file_path, 'w') as modified_file:
                for record in SeqIO.parse(file_path, 'fasta'):
                    # Update ID with <genome_id>::<protein_ID> format
                    record.id = f'{{strain_name}}::{{record.id}}'
                    record.description = ''  # Clear description to avoid duplication
                    SeqIO.write(record, modified_file, 'fasta')

            print(f'Modified protein IDs in file: {{modified_file_path}}')
    
    return modified_file_dir

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

def get_full_phage_list(interaction_matrix, input_phage_dir):
    print(f'Reading interaction matrix for phages: {{interaction_matrix}}')
    interaction_df = pd.read_csv(interaction_matrix)
    
    # Assume phage column is named 'phage' - adjust if different
    phage_column = 'phage'
    if phage_column not in interaction_df.columns:
        print(f'WARNING: Column \\'{{phage_column}}\\' not found in interaction matrix')
        # Try to find a phage-related column
        potential_cols = [col for col in interaction_df.columns if 'phage' in col.lower()]
        if potential_cols:
            phage_column = potential_cols[0]
            print(f'Using column: {{phage_column}}')
        else:
            print('No phage column found, using all phage files from directory')
            phage_files = [f for f in os.listdir(input_phage_dir) if f.endswith('.faa')]
            return ['.'.join(f.split('.')[:-1]) for f in phage_files]
    
    phages_in_matrix = interaction_df[phage_column].unique()
    phages_in_matrix = [str(s) for s in phages_in_matrix]
    print(f'Found {{len(phages_in_matrix)}} unique phages in interaction matrix')
    
    print(f'Reading phage directory: {{input_phage_dir}}')
    phage_files = [f for f in os.listdir(input_phage_dir) if f.endswith('.faa')]
    phages_in_dir = ['.'.join(f.split('.')[:-1]) for f in phage_files]
    print(f'Found {{len(phages_in_dir)}} phage files in directory')
    
    full_phage_list = list(set(phages_in_matrix).intersection(set(phages_in_dir)))
    print(f'Intersection: {{len(full_phage_list)}} phages found in both matrix and directory')
    
    return full_phage_list

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
    metrics_file = os.path.join(output_dir, 'modeling', 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    best_cutoff = metrics_df['cut_off'].values[0]
    return best_cutoff

def create_modeling_strain_fasta(strain_dir, strain_list, output_path):
    '''Create FASTA file with only modeling strains'''
    print(f'Creating modeling strain FASTA from {{len(strain_list)}} strains')
    
    with open(output_path, 'w') as outfile:
        for strain in strain_list:
            strain_file = os.path.join(strain_dir, f'{{strain}}.faa')
            if os.path.exists(strain_file):
                for record in SeqIO.parse(strain_file, 'fasta'):
                    SeqIO.write(record, outfile, 'fasta')
    
    print(f'Modeling strain FASTA created at {{output_path}}')

def create_modeling_phage_fasta(phage_dir, phage_list, output_path):
    '''Create FASTA file with modeling phages'''
    print(f'Creating modeling phage FASTA from {{len(phage_list)}} phages')
    
    with open(output_path, 'w') as outfile:
        for phage in phage_list:
            phage_file = os.path.join(phage_dir, f'{{phage}}.faa')
            if os.path.exists(phage_file):
                for record in SeqIO.parse(phage_file, 'fasta'):
                    SeqIO.write(record, outfile, 'fasta')
    
    print(f'Modeling phage FASTA created at {{output_path}}')

def create_modeling_strain_csv(strain_csv, strain_list, output_path):
    '''Create filtered strain protein CSV with only modeling strains'''
    print('Creating filtered strain protein CSV')
    
    strain_df = pd.read_csv(strain_csv)
    filtered_df = strain_df[strain_df['strain'].isin(strain_list)]
    filtered_df.to_csv(output_path, index=False)
    
    print(f'Filtered strain CSV created at {{output_path}} with {{len(filtered_df)}} proteins')

def create_modeling_phage_csv(phage_csv, phage_list, output_path):
    '''Create filtered phage protein CSV with modeling phages'''
    print('Creating filtered phage protein CSV')
    
    phage_df = pd.read_csv(phage_csv)
    filtered_df = phage_df[phage_df['phage'].isin(phage_list)]
    filtered_df.to_csv(output_path, index=False)
    
    print(f'Filtered phage CSV created at {{output_path}} with {{len(filtered_df)}} proteins')

def load_strain_sequences(input_dir, strain_list, suffix='faa'):
    '''Load amino acid sequences from strain FASTA files.'''
    strain_sequences = {{}}
    for strain_name in strain_list:
        fasta_file = f'{{strain_name}}.{{suffix}}'
        file_path = os.path.join(input_dir, fasta_file)
        
        if os.path.exists(file_path):
            try:
                sequences = [str(record.seq) for record in SeqIO.parse(file_path, 'fasta')]
                strain_sequences[strain_name] = sequences
                print(f'Loaded {{len(sequences)}} sequences for strain {{strain_name}}')
            except Exception as e:
                print(f'Error loading {{fasta_file}}: {{e}}')
        else:
            print(f'Warning: FASTA file not found for strain {{strain_name}}: {{file_path}}')
    
    print(f'Loaded {{len(strain_sequences)}} strains')
    return strain_sequences

def process_strain_features(strain_name, sequences, feature_kmers):
    '''Process a single strain for k-mer feature assignment.'''
    feature_presence = {{}}
    for feature, kmers in feature_kmers.items():
        feature_present = any(any(kmer in seq for seq in sequences) for kmer in kmers)
        feature_presence[feature] = 1 if feature_present else 0
    return strain_name, feature_presence

def assign_kmer_features_to_strains(strain_sequences, feature_map_file, predictive_feature_table, threads=4):
    '''Assign k-mer features to strains, limited to strain features only.'''
    print(f'Loading feature map from {{feature_map_file}}')
    feature_map = pd.read_csv(feature_map_file)
    
    print(f'Loading predictive features from {{predictive_feature_table}}')
    predictive_features_df = pd.read_csv(predictive_feature_table)
    
    # Only assign strain features (sc_*) to strains
    strain_features = [col for col in predictive_features_df.columns if col.startswith('sc_')]
    predictive_features = strain_features
    
    # Filter feature_map to only predictive features
    feature_map = feature_map[feature_map['Feature'].isin(predictive_features)]
    feature_to_kmers = feature_map.groupby('Feature')['Cluster_Label'].apply(list).to_dict()
    
    print(f'Processing {{len(feature_to_kmers)}} predictive k-mer features')
    
    # Process strains
    if threads > 1 and len(strain_sequences) > 1:
        print(f'Processing {{len(strain_sequences)}} strains with {{threads}} threads')
        process_func = partial(process_strain_features, feature_kmers=feature_to_kmers)
        with multiprocessing.Pool(processes=threads) as pool:
            results = pool.starmap(process_func, strain_sequences.items())
    else:
        print(f'Processing {{len(strain_sequences)}} strains sequentially')
        results = [process_strain_features(strain, seqs, feature_to_kmers) 
                  for strain, seqs in strain_sequences.items()]
    
    # Create feature presence DataFrame
    feature_presence_dict = {{strain: features for strain, features in results}}
    strain_feature_df = pd.DataFrame.from_dict(feature_presence_dict, orient='index')
    strain_feature_df.index.name = 'strain'
    
    # Ensure all predictive features are present, filling missing with 0
    for feature in predictive_features:
        if feature not in strain_feature_df.columns:
            strain_feature_df[feature] = 0
    
    # Reset index to make 'strain' a column and reorder
    strain_feature_df = strain_feature_df.reset_index()
    strain_feature_df = strain_feature_df[['strain'] + predictive_features]
    
    print(f'K-mer feature assignment completed for {{len(strain_feature_df)}} strains')
    return strain_feature_df

# Get iteration number from SLURM array task ID
iteration = int(os.environ['SLURM_ARRAY_TASK_ID'])
print(f'Processing iteration {{iteration}}')

# Parameters from command line arguments
output_dir = '{args.output_dir}'
iteration_output_dir = os.path.join(output_dir, f'iteration_{{iteration}}')
final_predictions_file = os.path.join(iteration_output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')

# Skip if already completed
if os.path.exists(final_predictions_file):
    print(f'Iteration {{iteration}} is already complete, skipping.')
    sys.exit(0)

print(f'Starting iteration {{iteration}}...')
os.makedirs(iteration_output_dir, exist_ok=True)
modeling_tmp_dir = os.path.join(iteration_output_dir, 'tmp')

# Get full strain and phage lists
full_strain_list = get_full_strain_list('{args.interaction_matrix}', '{args.input_strain_dir}', '{args.strain_column}')
full_phage_list = get_full_phage_list('{args.interaction_matrix}', '{args.input_phage_dir}')

if len(full_strain_list) == 0:
    print('FATAL ERROR: No valid strains found. Cannot proceed with iteration.')
    sys.exit(1)

if len(full_phage_list) == 0:
    print('FATAL ERROR: No valid phages found. Cannot proceed with iteration.')
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
        
    else:
        print('Strain lists already exist. Loading...')
        modeling_strains_df = pd.read_csv(modeling_strains_path)
        validation_strains_df = pd.read_csv(validation_strains_path)
        
        if len(modeling_strains_df) == 0 or len(validation_strains_df) == 0:
            print(f'ERROR: Existing strain CSV files are empty!')
            modeling_strains, validation_strains = split_strains(full_strain_list, iteration=iteration)
            pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
            pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
        else:
            modeling_strains = modeling_strains_df['strain'].tolist()
            validation_strains = validation_strains_df['strain'].tolist()

# For phages, use all phages for modeling (can be refined later)
modeling_phages = full_phage_list.copy()
validation_phages = full_phage_list.copy()

 # Check for duplicate protein IDs and modify if necessary
strain_dir_to_use = '{args.input_strain_dir}'
phage_dir_to_use = '{args.input_phage_dir}'

# Check strain duplicates
strain_duplicate_found = detect_duplicate_ids('{args.input_strain_dir}', 'faa', None, 'directory')
if strain_duplicate_found:
    print('Duplicate protein IDs found in strain directory; modifying all protein IDs.')
    strain_dir_to_use = modify_duplicate_ids('{args.input_strain_dir}', iteration_output_dir, 'faa', None, 'strain')

# Check phage duplicates  
phage_duplicate_found = detect_duplicate_ids('{args.input_phage_dir}', 'faa', None, 'directory')
if phage_duplicate_found:
    print('Duplicate protein IDs found in phage directory; modifying all protein IDs.')
    phage_dir_to_use = modify_duplicate_ids('{args.input_phage_dir}', iteration_output_dir, 'faa', None, 'phage')

# Step 1: Run k-mer table workflow with modeling strains and phages
merged_feature_file = os.path.join(iteration_output_dir, 'full_feature_table.csv')

modeling_performance_path = os.path.join(iteration_output_dir, 'modeling', 'modeling_results', 'model_performance', 'model_performance_metrics.csv')

if not os.path.exists(modeling_performance_path):
    print('Running k-mer table workflow...')
    
    # Create filtered FASTA and protein CSV files for modeling genomes
    modeling_strain_fasta = os.path.join(iteration_output_dir, 'modeling_strains.faa')
    modeling_phage_fasta = os.path.join(iteration_output_dir, 'modeling_phages.faa')
    modeling_strain_csv = os.path.join(iteration_output_dir, 'modeling_strain_proteins.csv')
    modeling_phage_csv = os.path.join(iteration_output_dir, 'modeling_phage_proteins.csv')

    create_modeling_strain_fasta(strain_dir_to_use, modeling_strains, modeling_strain_fasta)
    create_modeling_phage_fasta(phage_dir_to_use, modeling_phages, modeling_phage_fasta)
    
    # Generate protein mapping CSV files using potentially modified directories
    feature_tables_dir = os.path.join(iteration_output_dir, 'feature_tables')
    os.makedirs(feature_tables_dir, exist_ok=True)
    
    strain_protein_csv_path = generate_protein_mapping_csv(strain_dir_to_use, 
        os.path.join(feature_tables_dir, 'strain_protein_features.csv'), 'strain')
    phage_protein_csv_path = generate_protein_mapping_csv(phage_dir_to_use, 
        os.path.join(feature_tables_dir, 'phage_protein_features.csv'), 'phage')
    
    create_modeling_strain_csv(strain_protein_csv_path, modeling_strains, modeling_strain_csv)
    create_modeling_phage_csv(phage_protein_csv_path, modeling_phages, modeling_phage_csv)
    
    # Run k-mer workflow with separate strain and phage processing
    returned_path = run_kmer_table_workflow(
        strain_fasta=modeling_strain_fasta,
        protein_csv=modeling_strain_csv,
        k={args.k},
        id_col='strain',
        one_gene={args.one_gene},
        output_dir=iteration_output_dir,
        k_range={args.k_range},
        phenotype_matrix='{args.interaction_matrix}',
        phage_fasta=modeling_phage_fasta,
        protein_csv_phage=modeling_phage_csv,
        remove_suffix=False,
        sample_column='{args.strain_column}',
        phenotype_column='interaction',
        modeling=True,
        filter_type='none',
        num_features=100,
        num_runs_fs={args.num_runs_fs},
        num_runs_modeling={args.num_runs_modeling},
        method='rfe',
        strain_list=None,  # Already filtered in FASTA files
        phage_list=None,   # Already filtered in FASTA files
        threads={args.threads},
        task_type='classification',
        max_features='none',
        ignore_families={args.ignore_families},
        max_ram={args.max_ram},
        use_shap=False,
        use_clustering={args.use_clustering},
        cluster_method='{args.cluster_method}',
        n_clusters={args.n_clusters},
        min_cluster_size={args.min_cluster_size},
        min_samples={args.min_samples},
        cluster_selection_epsilon={args.cluster_selection_epsilon},
        use_dynamic_weights={args.use_dynamic_weights},
        weights_method='{args.weights_method}',
        check_feature_presence={args.check_feature_presence},
        filter_by_cluster_presence={args.filter_by_cluster_presence},
        min_cluster_presence={args.min_cluster_presence},
        use_feature_clustering={args.use_feature_clustering},
        feature_cluster_method='{args.feature_cluster_method}',
        feature_n_clusters={args.feature_n_clusters},
        feature_min_cluster_presence={args.feature_min_cluster_presence}
    )
    
    # The workflow returns the path to the merged feature table
    if returned_path and os.path.exists(returned_path):
        merged_feature_file = returned_path
        print(f'K-mer workflow completed. Merged feature table: {{merged_feature_file}}')
    else:
        print('K-mer workflow failed to produce merged feature table.')
        sys.exit(1)
else:
    print('K-mer modeling results already exist. Skipping k-mer workflow.')

# Step 2: Get the cutoff with the highest MCC
best_cutoff = select_best_cutoff(iteration_output_dir)
model_dir = os.path.join(iteration_output_dir, 'modeling', 'modeling_results', str(best_cutoff))

# Step 3: Predict interactions for validation strains and phages
validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
predict_output_dir = os.path.join(validation_output_dir, 'predict_results')

# Check if validation is already completed
if os.path.exists(final_predictions_file):
    print('Validation predictions already exist. Skipping validation workflow.')
else:
    print('Running validation workflow...')
    os.makedirs(validation_output_dir, exist_ok=True)
    os.makedirs(predict_output_dir, exist_ok=True)

    # Get required paths from the k-mer workflow results
    feature_map = os.path.join(iteration_output_dir, 'feature_tables', 'selected_features.csv')
    select_feature_table = os.path.join(iteration_output_dir, 'modeling', 'feature_selection', 
                                       'filtered_feature_tables', f'select_feature_table_{{best_cutoff}}.csv')
    phage_feature_table = os.path.join(iteration_output_dir, 'feature_tables', 'phage_final_feature_table.csv')

    # Check for required files
    required_files = [
        (feature_map, 'Feature map (selected_features.csv)'),
        (select_feature_table, 'Selected feature table'),
        (phage_feature_table, 'Phage feature table')
    ]

    missing_files = []
    for file_path, description in required_files:
        if not os.path.exists(file_path):
            missing_files.append(f'{{description}}: {{file_path}}')

    if missing_files:
        print('Missing required files for k-mer assignment:')
        for missing in missing_files:
            print(f'  - {{missing}}')
        
        # Try alternative paths
        alt_phage_table = os.path.join(iteration_output_dir, 'phage', 'features', 'feature_table.csv')
        if not os.path.exists(phage_feature_table) and os.path.exists(alt_phage_table):
            phage_feature_table = alt_phage_table
            print(f'Using alternative phage feature table: {{phage_feature_table}}')
        
        # Check again after alternatives
        still_missing = []
        for file_path, description in required_files:
            if not os.path.exists(file_path):
                still_missing.append(f'{{description}}: {{file_path}}')
        
        if still_missing:
            print('Still missing critical files. Cannot proceed with prediction.')
            for missing in still_missing:
                print(f'  - {{missing}}')
            sys.exit(1)

    # Load validation strain sequences
    strain_sequences = load_strain_sequences(strain_dir_to_use, validation_strains, 'faa')

    if not strain_sequences:
        print('ERROR: No validation strain sequences loaded')
        sys.exit(1)

    # Assign k-mer features to validation strains
    print('Assigning k-mer features to validation strains...')
    validation_strain_features = assign_kmer_features_to_strains(
        strain_sequences=strain_sequences,
        feature_map_file=feature_map,
        predictive_feature_table=select_feature_table,
        threads={args.threads}
    )

    # Save validation strain feature table
    validation_strain_feature_path = os.path.join(validation_output_dir, 'validation_feature_table.csv')
    validation_strain_features.to_csv(validation_strain_feature_path, index=False)
    print(f'Validation strain feature table saved to {{validation_strain_feature_path}}')

    # Run prediction workflow
    print('Running prediction workflow...')
    run_prediction_workflow(
        input_dir=validation_output_dir,
        phage_feature_table_path=phage_feature_table,
        model_dir=model_dir,
        output_dir=predict_output_dir,
        strain_source='strain',
        phage_source='phage',
        threads={args.threads}
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

print('Aggregating k-mer bootstrap results from {args.n_iterations} iterations...')
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
    final_predictions.to_csv(os.path.join(output_dir, 'final_kmer_predictions.csv'), index=False)
    print(f'Final k-mer predictions saved with {{len(final_predictions)}} total predictions from {{completed_iterations}} iterations')
    
    # Generate summary statistics
    summary_stats = final_predictions.groupby(['strain', 'phage']).agg({{
        'prediction': ['mean', 'std', 'count']
    }}).round(4)
    summary_stats.columns = ['mean_prediction', 'std_prediction', 'n_iterations']
    summary_stats = summary_stats.reset_index()
    summary_stats.to_csv(os.path.join(output_dir, 'kmer_prediction_summary.csv'), index=False)
    print(f'K-mer prediction summary saved with {{len(summary_stats)}} unique strain-phage pairs')
    
    # Generate performance summary by iteration
    iteration_performance = final_predictions.groupby('iteration').agg({{
        'prediction': ['count', 'mean', 'std']
    }}).round(4)
    iteration_performance.columns = ['n_predictions', 'mean_prediction', 'std_prediction']
    iteration_performance = iteration_performance.reset_index()
    iteration_performance.to_csv(os.path.join(output_dir, 'iteration_performance_summary.csv'), index=False)
    print(f'Iteration performance summary saved')
else:
    print('ERROR: No results found to aggregate!')
    exit(1)
"

echo "K-mer bootstrap aggregation completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "aggregate_results.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Submit k-mer bootstrap validation workflow as SLURM job array")
    
    # Required input arguments
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
    parser.add_argument('--method', type=str, default='rfe', choices=['rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'], help="Feature selection method.")
    
    # K-mer specific parameters
    parser.add_argument('--k', type=int, default=5, help="K-mer length.")
    parser.add_argument('--k_range', action='store_true', help="Generate k-mers from 3 to k.")
    parser.add_argument('--one_gene', action='store_true', help="Include features with only one gene.")
    parser.add_argument('--ignore_families', action='store_true', help="Ignore protein families when defining k-mer features.")
    
    # Feature selection parameters
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
    
    # Feature clustering parameters
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
    run_dir = f"kmer_bootstrap_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== K-mer Bootstrap Validation SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of iterations: {args.n_iterations}")
    print(f"Threads per job: {args.threads}")
    print(f"K-mer length: {args.k}")
    print(f"K-mer range: {args.k_range}")
    print(f"Ignore families: {args.ignore_families}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    print(f"Use feature clustering: {args.use_feature_clustering}")
    if args.use_feature_clustering:
        print(f"Feature cluster method: {args.feature_cluster_method}")
        print(f"Feature n_clusters: {args.feature_n_clusters}")
        print(f"Feature min presence: {args.feature_min_cluster_presence}")
    
    # Resource validation
    if args.mem_per_job < 32:
        print(f"WARNING: {args.mem_per_job}GB memory may be insufficient for k-mer generation on complex datasets")
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
    print("Submitting k-mer bootstrap job array...")
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
    
    print(f"\n=== K-mer Bootstrap Job Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Bootstrap array: {args.n_iterations} parallel jobs")
    print(f"Expected total runtime: {args.time_limit} per iteration (parallel)")
    print("\nMonitor with:")
    print("  squeue -u $USER")
    print("  squeue -u $USER --name=kmer_bootstrap  # Just k-mer bootstrap jobs")
    print("  tail -f logs/bootstrap_*.out")
    print(f"\nResults will be in:")
    print(f"  Individual iterations: {args.output_dir}/iteration_N/")
    print(f"  Final aggregated: {args.output_dir}/final_kmer_predictions.csv")
    print(f"  Summary stats: {args.output_dir}/kmer_prediction_summary.csv")

if __name__ == "__main__":
    main()