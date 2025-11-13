#!/usr/bin/env python3
"""
SLURM workflow submission for hybrid bootstrap cross-validation.
Uses PROTEIN FAMILIES for strains and K-MER FEATURES for phages.
Creates a job array where each job handles one iteration of the bootstrap validation.
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
    """Create SLURM job array script for hybrid bootstrap iterations"""
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=hybrid_bootstrap
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

echo "=== Hybrid Bootstrap Validation - Iteration $SLURM_ARRAY_TASK_ID ==="
echo "Job: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

# Run the hybrid bootstrap iteration logic
python3 -c "
import sys
import os
sys.path.insert(0, '{script_dir}')

import pandas as pd
import shutil
import random
from Bio import SeqIO
from collections import defaultdict
from phage_modeling.workflows.kmer_table_workflow import construct_feature_table
from phage_modeling.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table
from phage_modeling.workflows.assign_predict_workflow import assign_predict_workflow
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables

def get_full_strain_list(interaction_matrix, input_strain_dir, strain_column):
    print(f'Reading interaction matrix: {{interaction_matrix}}')
    interaction_df = pd.read_csv(interaction_matrix)
    
    if strain_column not in interaction_df.columns:
        print(f'ERROR: Column \\'{{strain_column}}\\' not found')
        return []
    
    strains_in_matrix = [str(s) for s in interaction_df[strain_column].unique()]
    strain_files = [f for f in os.listdir(input_strain_dir) if f.endswith('.faa')]
    strains_in_dir = ['.'.join(f.split('.')[:-1]) for f in strain_files]
    
    full_strain_list = list(set(strains_in_matrix).intersection(set(strains_in_dir)))
    print(f'Found {{len(full_strain_list)}} strains in both matrix and directory')
    return full_strain_list

def get_full_phage_list(interaction_matrix, input_phage_dir):
    print(f'Reading interaction matrix for phages')
    interaction_df = pd.read_csv(interaction_matrix)
    
    phage_column = 'phage'
    if phage_column not in interaction_df.columns:
        potential_cols = [col for col in interaction_df.columns if 'phage' in col.lower()]
        if potential_cols:
            phage_column = potential_cols[0]
    
    phages_in_matrix = [str(s) for s in interaction_df[phage_column].unique()]
    phage_files = [f for f in os.listdir(input_phage_dir) if f.endswith('.faa')]
    phages_in_dir = ['.'.join(f.split('.')[:-1]) for f in phage_files]
    
    full_phage_list = list(set(phages_in_matrix).intersection(set(phages_in_dir)))
    print(f'Found {{len(full_phage_list)}} phages')
    return full_phage_list

def generate_protein_mapping_csv(fasta_dir, output_csv, genome_col_name, file_extension='.faa'):
    '''Generate protein-to-genome mapping CSV from FASTA directory'''
    data = []
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(file_extension):
            genome_name = fasta_file.replace(file_extension, '')
            fasta_path = os.path.join(fasta_dir, fasta_file)
            for record in SeqIO.parse(fasta_path, 'fasta'):
                data.append({{'protein_ID': record.id, genome_col_name: genome_name}})
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f'Generated {{output_csv}} with {{len(df)}} protein mappings')
    return output_csv

def detect_duplicate_ids(input_path, suffix='faa', genomes_to_process=None, input_type='directory'):
    '''Detects duplicate protein IDs across .faa files in the input directory.'''
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

    for file_name in file_list:
        if file_name.endswith(suffix):
            genome_name = file_name.replace(f'.{{suffix}}', '')
            if genomes_to_process and genome_name not in genomes_to_process:
                continue
            
            file_path = os.path.join(input_path, file_name)
            for record in SeqIO.parse(file_path, 'fasta'):
                if record.id in protein_id_tracker:
                    duplicate_found = True
                    print(f'Duplicate protein ID found: {{record.id}} in genome {{genome_name}}')
                    break
                protein_id_tracker[record.id].add(genome_name)
    
    return duplicate_found

def modify_duplicate_ids(input_path, output_dir, suffix='faa', genomes_to_process=None, genome_column='genome'):
    '''Modifies all protein IDs by prefixing them with genome names to ensure uniqueness.'''
    print(f'Duplicate protein IDs found; modifying protein IDs and saving to {{output_dir}}/modified_AAs/{{genome_column}}')
    modified_file_dir = os.path.join(output_dir, 'modified_AAs', genome_column)
    os.makedirs(modified_file_dir, exist_ok=True)

    for file_name in os.listdir(input_path):
        if file_name.endswith(suffix):
            genome_name = file_name.replace(f'.{{suffix}}', '')
            if genomes_to_process and genome_name not in genomes_to_process:
                continue
            
            file_path = os.path.join(input_path, file_name)
            modified_file_path = os.path.join(modified_file_dir, file_name)
            
            with open(modified_file_path, 'w') as modified_file:
                for record in SeqIO.parse(file_path, 'fasta'):
                    # Update ID with <genome_id>::<protein_ID> format
                    record.id = f'{{genome_name}}::{{record.id}}'
                    record.description = ''  # Clear description to avoid duplication
                    SeqIO.write(record, modified_file, 'fasta')

            print(f'Modified protein IDs in file: {{modified_file_path}}')
    
    return modified_file_dir

def split_strains(full_strain_list, iteration, validation_percentage=0.1):
    if len(full_strain_list) == 0:
        return [], []
    
    random.seed(iteration)
    strain_list_copy = full_strain_list.copy()
    random.shuffle(strain_list_copy)
    
    split_index = int(len(strain_list_copy) * (1 - validation_percentage))
    modeling_strains = strain_list_copy[:split_index]
    validation_strains = strain_list_copy[split_index:]
    
    print(f'Split: {{len(modeling_strains)}} modeling, {{len(validation_strains)}} validation')
    return modeling_strains, validation_strains

def select_best_cutoff(output_dir):
    metrics_file = os.path.join(output_dir, 'modeling', 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    return metrics_df['cut_off'].values[0]

def feature_selection_optimized(presence_absence, source, genome_column_name, output_dir=None, prefix=None):
    '''Optimizes feature selection by identifying perfect co-occurrence of features using hashing.'''
    print('Optimizing feature selection using hashing...')
    
    presence_absence.set_index(genome_column_name, inplace=True)
    presence_absence = presence_absence.applymap(lambda x: 1 if x > 0 else 0)
    
    column_hashes = presence_absence.apply(lambda col: hash(tuple(col)), axis=0)
    
    hash_to_columns = {{}}
    for col, col_hash in column_hashes.items():
        hash_to_columns.setdefault(col_hash, []).append(col)
    
    unique_clusters = list(hash_to_columns.values())
    
    data = [
        (f'{{source[0]}}c_{{idx}}', cluster)
        for idx, cluster_group in enumerate(unique_clusters)
        for cluster in cluster_group
    ]
    selected_features = pd.DataFrame(data, columns=['Feature', 'Cluster_Label'])
    
    if output_dir:
        selected_features_output = os.path.join(output_dir, f'{{prefix}}_selected_features.csv' if prefix else 'selected_features.csv')
        selected_features.to_csv(selected_features_output, index=False)
        print(f'Selected features saved to {{selected_features_output}}')
    
    return selected_features

# Get iteration number
iteration = int(os.environ['SLURM_ARRAY_TASK_ID'])
print(f'Processing iteration {{iteration}}')

# Parameters
output_dir = '{args.output_dir}'
iteration_output_dir = os.path.join(output_dir, f'iteration_{{iteration}}')
final_predictions_file = os.path.join(iteration_output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')

# Skip if already completed
if os.path.exists(final_predictions_file):
    print(f'Iteration {{iteration}} already complete, skipping')
    sys.exit(0)

print(f'Starting iteration {{iteration}}...')
os.makedirs(iteration_output_dir, exist_ok=True)

# Get full lists
full_strain_list = get_full_strain_list('{args.interaction_matrix}', '{args.input_strain_dir}', '{args.strain_column}')
full_phage_list = get_full_phage_list('{args.interaction_matrix}', '{args.input_phage_dir}')

if len(full_strain_list) == 0 or len(full_phage_list) == 0:
    print('FATAL ERROR: No valid strains or phages found')
    sys.exit(1)

# Handle strain splits
clustering_dir = {repr(args.clustering_dir) if args.clustering_dir else 'None'}
if clustering_dir:
    modeling_strains_path = os.path.join(clustering_dir, f'iteration_{{iteration}}', 'modeling_strains.csv')
    validation_strains_path = os.path.join(clustering_dir, f'iteration_{{iteration}}', 'validation_strains.csv')
    
    new_modeling_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    new_validation_path = os.path.join(iteration_output_dir, 'validation_strains.csv')
    
    if not os.path.exists(new_modeling_path):
        os.symlink(modeling_strains_path, new_modeling_path)
    if not os.path.exists(new_validation_path):
        os.symlink(validation_strains_path, new_validation_path)
    
    modeling_strains = pd.read_csv(new_modeling_path)['strain'].tolist()
    validation_strains = pd.read_csv(new_validation_path)['strain'].tolist()
else:
    modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
    validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')
    
    if not os.path.exists(modeling_strains_path):
        modeling_strains, validation_strains = split_strains(full_strain_list, iteration)
        pd.DataFrame(modeling_strains, columns=['strain']).to_csv(modeling_strains_path, index=False)
        pd.DataFrame(validation_strains, columns=['strain']).to_csv(validation_strains_path, index=False)
    else:
        modeling_strains = pd.read_csv(modeling_strains_path)['strain'].tolist()
        validation_strains = pd.read_csv(validation_strains_path)['strain'].tolist()

# Check for duplicate protein IDs and modify if necessary
print('\\n=== Checking for duplicate protein IDs ===')
strain_dir_to_use = '{args.input_strain_dir}'
phage_dir_to_use = '{args.input_phage_dir}'

# Check and handle strain duplicates
strain_duplicate_found = detect_duplicate_ids('{args.input_strain_dir}', 'faa', None, 'directory')
if strain_duplicate_found:
    print('Duplicate protein IDs found in strain directory; modifying all protein IDs.')
    strain_dir_to_use = modify_duplicate_ids('{args.input_strain_dir}', iteration_output_dir, 'faa', None, 'strain')
else:
    print('No duplicate protein IDs found in strain directory.')

# Check and handle phage duplicates  
phage_duplicate_found = detect_duplicate_ids('{args.input_phage_dir}', 'faa', None, 'directory')
if phage_duplicate_found:
    print('Duplicate protein IDs found in phage directory; modifying all protein IDs.')
    phage_dir_to_use = modify_duplicate_ids('{args.input_phage_dir}', iteration_output_dir, 'faa', None, 'phage')
else:
    print('No duplicate protein IDs found in phage directory.')

# STEP 1: Generate STRAIN protein family features for modeling strains
print('\\n=== STEP 1: Generating strain protein family features ===')
strain_feature_table = os.path.join(iteration_output_dir, 'strain', 'features', 'feature_table.csv')

if not os.path.exists(strain_feature_table):
    print('Running protein family clustering for strains...')
    
    # Import the lower-level functions
    
    strain_output_dir = os.path.join(iteration_output_dir, 'strain')
    strain_tmp_dir = os.path.join(iteration_output_dir, 'tmp', 'strain')
    
    # Run clustering
    run_clustering_workflow(
        strain_dir_to_use,
        strain_output_dir,
        strain_tmp_dir,
        min_seq_id={args.min_seq_id},
        coverage={args.coverage},
        sensitivity={args.sensitivity},
        suffix='faa',
        threads={args.threads},
        strain_list=modeling_strains_path,
        strain_column='strain',
        compare=False,
        bootstrapping=False,
        clear_tmp=False
    )
    
    # Run feature assignment
    run_feature_assignment(
        input_file=os.path.join(strain_output_dir, 'presence_absence_matrix.csv'),
        output_dir=os.path.join(strain_output_dir, 'features'),
        source='strain',
        select=modeling_strains_path,
        select_column='strain',
        max_ram={args.max_ram},
        threads={args.threads}
    )
    
    print(f'Strain protein family features saved to {{strain_feature_table}}')
else:
    print('Strain features already exist')

# STEP 2: Generate PHAGE k-mer features for ALL phages
print('\\n=== STEP 2: Generating phage k-mer features ===')
phage_raw_feature_table = os.path.join(iteration_output_dir, 'phage', 'features', 'phage_feature_table.csv')
phage_selected_features_path = os.path.join(iteration_output_dir, 'phage', 'features', 'phage_selected_features.csv')
phage_genome_assignments_path = os.path.join(iteration_output_dir, 'phage', 'features', 'phage_genome_assignments.csv')
phage_final_feature_table = os.path.join(iteration_output_dir, 'phage', 'features', 'phage_final_feature_table.csv')

if not os.path.exists(phage_final_feature_table):
    print('Generating k-mer features for all phages...')
    
    # Generate phage protein CSV
    phage_features_dir = os.path.join(iteration_output_dir, 'phage', 'features')
    os.makedirs(phage_features_dir, exist_ok=True)
    
    phage_protein_csv = os.path.join(phage_features_dir, 'phage_protein_mapping.csv')
    if not os.path.exists(phage_protein_csv):
        generate_protein_mapping_csv(phage_dir_to_use, phage_protein_csv, 'phage')
    
    # Create single FASTA with all phages
    phage_combined_fasta = os.path.join(iteration_output_dir, 'phage', 'all_phages.faa')
    if not os.path.exists(phage_combined_fasta):
        with open(phage_combined_fasta, 'w') as outfile:
            for phage in full_phage_list:
                phage_file = os.path.join(phage_dir_to_use, f'{{phage}}.faa')
                if os.path.exists(phage_file):
                    for record in SeqIO.parse(phage_file, 'fasta'):
                        SeqIO.write(record, outfile, 'fasta')
    
    # Checkpoint 1: Generate raw k-mer feature table
    if not os.path.exists(phage_raw_feature_table):
        print(f'Generating k-mers with k={args.k}')
        construct_feature_table(
            fasta_file=phage_combined_fasta,
            protein_csv=phage_protein_csv,
            k={args.k},
            id_col='phage',
            one_gene={args.one_gene},
            output_dir=phage_features_dir,
            output_name='phage',
            k_range=False,
            ignore_families=True,
            genome_list=full_phage_list
        )
    else:
        print(f'Reusing existing raw phage feature table: {{phage_raw_feature_table}}')
    
    # Load raw phage feature table
    phage_feature_table_df = pd.read_csv(phage_raw_feature_table)
    
    # Checkpoint 2: Generate genome assignments and feature selection
    if not os.path.exists(phage_selected_features_path) or not os.path.exists(phage_genome_assignments_path):
        print('Generating phage genome assignments and feature selection...')
        
        # Get genome assignments
        phage_feature_table_df.rename(columns={{'Genome': 'phage'}}, inplace=True)
        phage_genome_assignments = phage_feature_table_df.melt(
            id_vars='phage', 
            var_name='Cluster_Label', 
            value_name='Presence'
        )
        phage_genome_assignments = phage_genome_assignments[phage_genome_assignments['Presence'] == 1]
        phage_genome_assignments = phage_genome_assignments.drop(columns=['Presence'])
        phage_genome_assignments.to_csv(phage_genome_assignments_path, index=False)
        print(f'Phage genome assignments saved to {{phage_genome_assignments_path}}')
        
        # Optimize feature selection
        phage_selected_features = feature_selection_optimized(
            phage_feature_table_df.copy(),
            'phage',
            'phage',
            phage_features_dir,
            prefix='phage'
        )
    else:
        print('Reusing existing phage assignments and selected features')
        phage_genome_assignments = pd.read_csv(phage_genome_assignments_path)
        phage_selected_features = pd.read_csv(phage_selected_features_path)
    
    del phage_feature_table_df
    import gc
    gc.collect()
    
    # Checkpoint 3: Create final phage feature table with pc_* features
    print('Creating final phage feature table with pc_* features...')
    
    # Merge assignments with selected features to get pc_* mapping
    phage_assignment_df = phage_genome_assignments.merge(
        phage_selected_features, 
        on='Cluster_Label', 
        how='inner'
    )
    phage_assignment_df = phage_assignment_df.drop(columns=['Cluster_Label']).drop_duplicates()
    
    # Create final feature table in wide format
    phage_final_df = phage_assignment_df.pivot_table(
        index='phage', 
        columns='Feature', 
        aggfunc=lambda x: 1, 
        fill_value=0
    ).reset_index()
    
    # Ensure all phages are included (even if they have no features)
    missing_phages = set(full_phage_list) - set(phage_final_df['phage'])
    if missing_phages:
        print(f'Adding {{len(missing_phages)}} phages with zero values')
        missing_df = pd.DataFrame({{'phage': list(missing_phages)}})
        for feature in phage_final_df.columns:
            if feature != 'phage':
                missing_df[feature] = 0
        phage_final_df = pd.concat([phage_final_df, missing_df], ignore_index=True)
    
    # Save final feature table
    phage_final_df.to_csv(phage_final_feature_table, index=False)
    print(f'Phage k-mer features with pc_* aggregation saved to {{phage_final_feature_table}}')
    
    del phage_assignment_df, phage_final_df, phage_genome_assignments, phage_selected_features
    gc.collect()
else:
    print('Phage k-mer features already exist')

# Set the phage feature table path for later use
phage_feature_table = phage_final_feature_table

# STEP 3: Merge strain and phage features with interaction matrix
print('\\n=== STEP 3: Merging feature tables ===')
merged_dir = os.path.join(iteration_output_dir, 'merged')
merged_feature_table = os.path.join(merged_dir, 'full_feature_table.csv')

if not os.path.exists(merged_feature_table):
    print('Merging strain protein families with phage k-mers...')
    os.makedirs(merged_dir, exist_ok=True)
    
    merge_feature_tables(
        strain_features=strain_feature_table,
        phenotype_matrix='{args.interaction_matrix}',
        output_dir=merged_dir,
        sample_column='{args.strain_column}',
        phage_features=phage_feature_table,
        remove_suffix=False,
        use_feature_clustering={args.use_feature_clustering},
        feature_cluster_method='{args.feature_cluster_method}',
        feature_n_clusters={args.feature_n_clusters},
        feature_min_cluster_presence={args.feature_min_cluster_presence}
    )
    print(f'Merged feature table saved to {{merged_feature_table}}')
else:
    print('Merged feature table already exists')

# STEP 4: Feature selection and modeling
print('\\n=== STEP 4: Feature selection and modeling ===')
modeling_dir = os.path.join(iteration_output_dir, 'modeling')
metrics_file = os.path.join(modeling_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')

if not os.path.exists(metrics_file):
    print('Running feature selection and modeling...')
    
    # Get paths for predictive proteins
    strain_feature2cluster = os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv')
    strain_cluster2protein = os.path.join(iteration_output_dir, 'strain', 'clusters.tsv')
    
    run_modeling_workflow_from_feature_table(
        full_feature_table=merged_feature_table,
        output_dir=modeling_dir,
        threads={args.threads},
        num_features=100,
        filter_type='strain',
        num_runs_fs={args.num_runs_fs},
        num_runs_modeling={args.num_runs_modeling},
        sample_column='{args.strain_column}',
        phage_column='phage',
        phenotype_column='interaction',
        method='rfe',
        annotation_table_path=None,
        protein_id_col='protein_ID',
        feature2cluster_path=strain_feature2cluster,
        cluster2protein_path=strain_cluster2protein,
        fasta_dir_or_file='{args.input_strain_dir}',
        run_predictive_proteins=False,
        phage_feature2cluster_path=None,
        phage_cluster2protein_path=None,
        phage_fasta_dir_or_file=None,
        task_type='classification',
        binary_data=True,
        max_features='none',
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
        max_ram={args.max_ram},
        use_shap=False
    )
    print('Modeling completed')
else:
    print('Modeling results already exist')

# STEP 5: Validation - Assign features to validation strains and predict
print('\\n=== STEP 5: Validation predictions ===')
validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')

if not os.path.exists(final_predictions_file):
    print('Running validation workflow...')
    
    # Get best cutoff
    best_cutoff = select_best_cutoff(iteration_output_dir)
    model_dir = os.path.join(modeling_dir, 'modeling_results', str(best_cutoff))
    
    # Get feature selection table
    select_feature_table = os.path.join(modeling_dir, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_{{best_cutoff}}.csv')
    
    # Check for modified AA directory (from protein family workflow)
    modified_aa_dir = os.path.join(iteration_output_dir, 'strain', 'modified_AAs', 'strain')
    if os.path.exists(modified_aa_dir):
        validation_strain_dir = modified_aa_dir
        print(f'Using modified strain directory for validation: {{modified_aa_dir}}')
    else:
        validation_strain_dir = strain_dir_to_use  # Use potentially modified directory from duplicate check
        print(f'Using strain directory for validation: {{validation_strain_dir}}')
    
    # Run assign_predict_workflow
    validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')
    
    assign_predict_workflow(
        input_dir=validation_strain_dir,
        genome_list=validation_strains_path,
        mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
        clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
        feature_map=os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv'),
        tmp_dir=validation_tmp_dir,
        suffix='faa',
        model_dir=model_dir,
        feature_table=select_feature_table,
        strain_feature_table_path=None,
        phage_feature_table_path=phage_feature_table,  # Use same phage features from training
        output_dir=validation_output_dir,
        threads={args.threads},
        genome_type='strain',
        sensitivity={args.sensitivity},
        coverage={args.coverage},
        min_seq_id={args.min_seq_id},
        duplicate_all={args.duplicate_all}
    )
    print('Validation completed')
else:
    print('Validation predictions already exist')

print(f'\\nIteration {{iteration}} completed successfully')
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
#SBATCH --job-name=hybrid_aggregate
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

echo "=== Hybrid Bootstrap Results Aggregation ==="
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

output_dir = '{args.output_dir}'
final_predictions = pd.DataFrame()

print('Aggregating hybrid bootstrap results from {args.n_iterations} iterations...')
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

# Save final predictions
if len(final_predictions) > 0:
    final_predictions.to_csv(os.path.join(output_dir, 'final_hybrid_predictions.csv'), index=False)
    print(f'Final hybrid predictions saved with {{len(final_predictions)}} total predictions from {{completed_iterations}} iterations')
    
    # Generate summary statistics
    summary_stats = final_predictions.groupby(['strain', 'phage']).agg({{
        'Confidence': ['mean', 'std', 'count']
    }}).round(4)
    summary_stats.columns = ['mean_confidence', 'std_confidence', 'n_iterations']
    summary_stats = summary_stats.reset_index()
    summary_stats.to_csv(os.path.join(output_dir, 'hybrid_prediction_summary.csv'), index=False)
    print(f'Hybrid prediction summary saved with {{len(summary_stats)}} unique strain-phage pairs')
else:
    print('ERROR: No results found to aggregate!')
    exit(1)
"

echo "Hybrid bootstrap aggregation completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "aggregate_results.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Submit hybrid bootstrap validation workflow")
    
    # Required inputs
    parser.add_argument('--input_strain_dir', type=str, required=True)
    parser.add_argument('--input_phage_dir', type=str, required=True)
    parser.add_argument('--interaction_matrix', type=str, required=True)
    parser.add_argument('--clustering_dir', type=str, help="Pre-computed clustering directory")
    parser.add_argument('--output_dir', type=str, required=True)
    
    # Core parameters
    parser.add_argument('--n_iterations', type=int, default=10)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--strain_column', type=str, default='strain')
    parser.add_argument('--max_ram', type=float, default=40)
    
    # Feature selection and modeling
    parser.add_argument('--num_runs_fs', type=int, default=25)
    parser.add_argument('--num_runs_modeling', type=int, default=50)
    parser.add_argument('--use_dynamic_weights', action='store_true')
    parser.add_argument('--weights_method', type=str, default='log10')
    parser.add_argument('--use_clustering', action='store_true')
    parser.add_argument('--cluster_method', type=str, default='hdbscan')
    parser.add_argument('--n_clusters', type=int, default=20)
    parser.add_argument('--min_cluster_size', type=int, default=2)
    parser.add_argument('--min_samples', type=int)
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0)
    parser.add_argument('--check_feature_presence', action='store_true')
    parser.add_argument('--filter_by_cluster_presence', action='store_true')
    parser.add_argument('--min_cluster_presence', type=float, default=2)
    parser.add_argument('--duplicate_all', action='store_true')
    parser.add_argument('--bootstrapping', action='store_true')
    
    # Protein family parameters (strains)
    parser.add_argument('--min_seq_id', type=float, default=0.4)
    parser.add_argument('--coverage', type=float, default=0.8)
    parser.add_argument('--sensitivity', type=float, default=7.5)
    
    # K-mer parameters (phages)
    parser.add_argument('--k', type=int, default=5)
    parser.add_argument('--k_range', action='store_true')
    parser.add_argument('--one_gene', action='store_true')
    parser.add_argument('--ignore_families', action='store_true')
    
    # Feature clustering
    parser.add_argument('--use_feature_clustering', action='store_true')
    parser.add_argument('--feature_cluster_method', type=str, default='hierarchical')
    parser.add_argument('--feature_n_clusters', type=int, default=20)
    parser.add_argument('--feature_min_cluster_presence', type=int, default=2)
    
    # SLURM configuration
    parser.add_argument('--account', default='ac_mak')
    parser.add_argument('--partition', default='lr7')
    parser.add_argument('--qos', default='lr_normal')
    parser.add_argument('--environment', default='phage_modeling')
    parser.add_argument('--mem_per_job', type=int, default=60)
    parser.add_argument('--time_limit', default='12:00:00')
    parser.add_argument('--dry_run', action='store_true')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input_strain_dir):
        print(f"Error: Strain directory not found: {args.input_strain_dir}")
        return 1
    if not os.path.exists(args.input_phage_dir):
        print(f"Error: Phage directory not found: {args.input_phage_dir}")
        return 1
    if not os.path.exists(args.interaction_matrix):
        print(f"Error: Interaction matrix not found: {args.interaction_matrix}")
        return 1
    
    # Validate interaction matrix
    try:
        import pandas as pd
        df = pd.read_csv(args.interaction_matrix)
        if args.strain_column not in df.columns:
            print(f"Error: Column '{args.strain_column}' not found in interaction matrix")
            return 1
        print(f"Interaction matrix validation: ✓ Found {len(df)} rows")
    except Exception as e:
        print(f"Error reading interaction matrix: {e}")
        return 1
    
    # Create run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"hybrid_bootstrap_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== Hybrid Bootstrap Validation SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"HYBRID APPROACH:")
    print(f"  Strain features: Protein families (sc_*)")
    print(f"  Phage features: K-mers (pc_*), k={args.k}")
    print(f"Iterations: {args.n_iterations}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print(f"Memory: {args.mem_per_job}GB, Time: {args.time_limit}")
    print()
    
    # Create scripts
    print("Creating SLURM job scripts...")
    bootstrap_script = create_bootstrap_job_array(args, run_dir)
    aggregate_script = create_final_aggregation_job(args, run_dir, "PLACEHOLDER")
    
    if args.dry_run:
        print("Dry run - scripts created but not submitted")
        return 0
    
    # Submit jobs
    original_dir = os.getcwd()
    os.chdir(run_dir)
    
    print("Submitting bootstrap job array...")
    bootstrap_job_id = submit_job("bootstrap_job_array.sh")
    print(f"Bootstrap job array: {bootstrap_job_id}")
    
    if bootstrap_job_id:
        # Update and submit aggregation
        with open("aggregate_results.sh", 'r') as f:
            content = f.read()
        content = content.replace("PLACEHOLDER", bootstrap_job_id)
        with open("aggregate_results.sh", 'w') as f:
            f.write(content)
        
        aggregate_job_id = submit_job("aggregate_results.sh")
        print(f"Aggregation job: {aggregate_job_id}")
    
    os.chdir(original_dir)
    
    print(f"\n=== Job Submission Summary ===")
    print(f"Bootstrap array: {args.n_iterations} parallel jobs")
    print(f"Results: {args.output_dir}/final_hybrid_predictions.csv")
    print("\nMonitor: squeue -u $USER")

if __name__ == "__main__":
    main()
