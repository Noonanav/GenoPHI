#!/usr/bin/env python3
"""
K-mer bootstrap workflow with ALL protein families support.
Designed to work with existing protein family workflow results.
Generates all_AA_seqs.faa from clusters.tsv and runs k-mer modeling.
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
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Determine output directory name based on mode
    if args.use_all_families_for_kmers:
        kmer_output_dirname = "kmer_modeling_all_families"
        workflow_mode_comment = "# MODE: Using ALL protein families for k-mer modeling"
    else:
        kmer_output_dirname = "kmer_modeling"
        workflow_mode_comment = "# MODE: Using PREDICTIVE protein families only"
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=kmer_bootstrap_all
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
echo "{workflow_mode_comment}"
echo "Job: $SLURM_JOB_ID, Task: $SLURM_ARRAY_TASK_ID, Node: $SLURMD_NODENAME, Started: $(date)"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

python3 -c "
import sys
import os
sys.path.insert(0, '{script_dir}')

import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import workflow components
from phage_modeling.workflows.kmer_table_workflow import run_kmer_table_workflow, is_fasta_empty
from phage_modeling.workflows.kmer_assign_predict_workflow import kmer_assign_predict_workflow
from phage_modeling.workflows.kmer_analysis_workflow import kmer_analysis_workflow
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

def generate_all_aa_seqs_fasta(clusters_tsv_path, source_fasta_dir, output_fasta_path, genome_type='strain'):
    '''
    Generate all_AA_seqs.faa from clusters.tsv by extracting all sequences.
    
    Args:
        clusters_tsv_path: Path to clusters.tsv (protein_family, protein_ID)
        source_fasta_dir: Directory containing source FASTA files (modified_AAs or original)
        output_fasta_path: Where to save the combined all_AA_seqs.faa
        genome_type: 'strain' or 'phage'
    '''
    logging.info(f'Generating all_AA_seqs.faa for {{genome_type}}...')
    
    # Read clusters.tsv to get all protein IDs
    clusters_df = pd.read_csv(clusters_tsv_path, sep='\\t', names=['protein_family', 'protein_ID'])
    all_protein_ids = set(clusters_df['protein_ID'].unique())
    n_protein_ids = len(all_protein_ids)
    logging.info(f'Found {{n_protein_ids}} unique protein IDs in clusters.tsv')
    
    # Read all sequences from source FASTA files
    all_sequences = {{}}
    
    if not os.path.exists(source_fasta_dir):
        logging.error(f'Source FASTA directory not found: {{source_fasta_dir}}')
        return False
    
    fasta_files = [f for f in os.listdir(source_fasta_dir) if f.endswith('.faa')]
    n_fasta_files = len(fasta_files)
    logging.info(f'Reading sequences from {{n_fasta_files}} FASTA files...')
    
    for fasta_file in fasta_files:
        fasta_path = os.path.join(source_fasta_dir, fasta_file)
        for record in SeqIO.parse(fasta_path, 'fasta'):
            all_sequences[record.id] = str(record.seq)
    
    n_sequences = len(all_sequences)
    logging.info(f'Loaded {{n_sequences}} total sequences from FASTA files')
    
    # Match protein IDs from clusters.tsv with sequences
    matched_sequences = []
    unmatched_ids = []
    
    for protein_id in all_protein_ids:
        if protein_id in all_sequences:
            matched_sequences.append(SeqRecord(
                Seq(all_sequences[protein_id]),
                id=protein_id,
                description=''
            ))
        else:
            unmatched_ids.append(protein_id)
    
    if unmatched_ids:
        n_unmatched = len(unmatched_ids)
        logging.warning(f'Could not find sequences for {{n_unmatched}} protein IDs')
        if n_unmatched <= 10:
            logging.warning(f'Unmatched IDs: {{unmatched_ids}}')
    
    # Write all sequences to output file
    os.makedirs(os.path.dirname(output_fasta_path), exist_ok=True)
    SeqIO.write(matched_sequences, output_fasta_path, 'fasta')
    n_matched = len(matched_sequences)
    logging.info(f'Wrote {{n_matched}} sequences to {{output_fasta_path}}')
    
    return n_matched > 0

def create_protein_overview_from_clusters(clusters_tsv_path, output_csv_path):
    '''
    Create protein overview CSV from clusters.tsv in the format expected by kmer workflow.
    
    Format: Feature, protein_family, protein_ID, strain/phage
    '''
    logging.info(f'Creating protein overview from {{clusters_tsv_path}}...')
    
    clusters_df = pd.read_csv(clusters_tsv_path, sep='\\t', names=['protein_family', 'protein_ID'])
    
    # Create Feature column (all families are features in this mode)
    clusters_df['Feature'] = clusters_df['protein_family']
    
    # Extract genome name from protein_ID (format: strain::protein or just protein)
    if '::' in clusters_df['protein_ID'].iloc[0]:
        clusters_df['strain'] = clusters_df['protein_ID'].str.split('::').str[0]
    else:
        # If no :: delimiter, try to infer from protein_ID
        # This is a fallback - may need adjustment based on your data
        clusters_df['strain'] = 'unknown'
    
    # Reorder columns to match expected format
    output_df = clusters_df[['Feature', 'protein_family', 'protein_ID']].copy()
    
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    output_df.to_csv(output_csv_path, index=False)
    n_proteins = len(output_df)
    n_families = output_df['protein_family'].nunique()
    logging.info(f'Created protein overview with {{n_proteins}} proteins from {{n_families}} families')
    
    return output_csv_path

def select_best_cutoff(modeling_dir):
    '''Select best cutoff from model performance.'''
    performance_file = os.path.join(modeling_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
    
    if not os.path.exists(performance_file):
        logging.error(f'Performance metrics not found: {{performance_file}}')
        return None
    
    df = pd.read_csv(performance_file)
    best_cutoff = df.iloc[0]['cut_off']
    logging.info(f'Selected best cutoff: {{best_cutoff}}')
    return best_cutoff

# Configure logging
iteration = int(os.environ['SLURM_ARRAY_TASK_ID'])
iteration_output_dir = os.path.join('{args.output_dir}', f'iteration_{{iteration}}')

log_file = os.path.join(iteration_output_dir, f'kmer_all_families_iteration_{{iteration}}.log')
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

logging.info(f'Starting k-mer modeling for iteration {{iteration}}')
logging.info('Mode: {'ALL FAMILIES' if {args.use_all_families_for_kmers} else 'PREDICTIVE ONLY'}')

# Check if iteration directory exists
if not os.path.exists(iteration_output_dir):
    logging.error(f'Iteration directory not found: {{iteration_output_dir}}')
    logging.error('Please run the protein family workflow first!')
    sys.exit(1)

# Check for required existing files
modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')
metrics_file = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')

if not os.path.exists(modeling_strains_path):
    logging.error(f'Modeling strains file not found: {{modeling_strains_path}}')
    sys.exit(1)
if not os.path.exists(validation_strains_path):
    logging.error(f'Validation strains file not found: {{validation_strains_path}}')
    sys.exit(1)
if not os.path.exists(metrics_file):
    logging.error(f'Metrics file not found: {{metrics_file}}')
    logging.error('Protein family modeling must be complete first!')
    sys.exit(1)

logging.info('✓ Found existing protein family workflow results')
logging.info('✓ Found existing train/test splits')

# Load strain lists
modeling_strains_df = pd.read_csv(modeling_strains_path)
validation_strains_df = pd.read_csv(validation_strains_path)
modeling_strains = modeling_strains_df['strain'].tolist()
validation_strains = validation_strains_df['strain'].tolist()
n_modeling = len(modeling_strains)
n_validation = len(validation_strains)
logging.info(f'Using {{n_modeling}} modeling strains, {{n_validation}} validation strains')

# Determine source directory for sequences
# Try multiple possible locations for modified_AAs
strain_source_dir = None
possible_strain_dirs = [
    os.path.join(iteration_output_dir, 'strain', 'modified_AAs', 'strain'),  # Nested structure
    os.path.join(iteration_output_dir, 'strain', 'modified_AAs'),            # Direct structure
    '{args.input_strain_dir}'                                                 # Original input
]

for dir_path in possible_strain_dirs:
    if os.path.exists(dir_path):
        fasta_files = [f for f in os.listdir(dir_path) if f.endswith('.faa')]
        if fasta_files:
            strain_source_dir = dir_path
            n_files = len(fasta_files)
            logging.info(f'Using strain directory: {{strain_source_dir}} ({{n_files}} files)')
            break

if not strain_source_dir:
    logging.error('Could not find strain FASTA files in any expected location')
    sys.exit(1)

# Try multiple possible locations for phage modified_AAs
phage_source_dir = None
possible_phage_dirs = [
    os.path.join(iteration_output_dir, 'phage', 'modified_AAs', 'phage'),   # Nested structure
    os.path.join(iteration_output_dir, 'phage', 'modified_AAs'),            # Direct structure
    '{args.input_phage_dir}'                                                 # Original input
]

for dir_path in possible_phage_dirs:
    if os.path.exists(dir_path):
        fasta_files = [f for f in os.listdir(dir_path) if f.endswith('.faa')]
        if fasta_files:
            phage_source_dir = dir_path
            n_files = len(fasta_files)
            logging.info(f'Using phage directory: {{phage_source_dir}} ({{n_files}} files)')
            break

if not phage_source_dir:
    logging.error('Could not find phage FASTA files in any expected location')
    sys.exit(1)

# ==================================================================================
# PREPARE INPUT FILES BASED ON MODE
# ==================================================================================

if {args.use_all_families_for_kmers}:
    logging.info('=== Preparing files for ALL FAMILIES mode ===')
    
    # Create directory for generated files
    all_families_dir = os.path.join(iteration_output_dir, 'all_families_data')
    os.makedirs(all_families_dir, exist_ok=True)
    
    # STRAIN: Generate all_AA_seqs.faa and protein overview
    strain_clusters_tsv = os.path.join(iteration_output_dir, 'strain', 'clusters.tsv')
    strain_all_aa_fasta = os.path.join(all_families_dir, 'strain_all_AA_seqs.faa')
    strain_protein_csv = os.path.join(all_families_dir, 'strain_all_families_overview.csv')
    
    if not os.path.exists(strain_all_aa_fasta):
        success = generate_all_aa_seqs_fasta(
            clusters_tsv_path=strain_clusters_tsv,
            source_fasta_dir=strain_source_dir,
            output_fasta_path=strain_all_aa_fasta,
            genome_type='strain'
        )
        if not success:
            logging.error('Failed to generate strain all_AA_seqs.faa')
            sys.exit(1)
    else:
        logging.info(f'Using existing strain all_AA_seqs.faa: {{strain_all_aa_fasta}}')
    
    if not os.path.exists(strain_protein_csv):
        create_protein_overview_from_clusters(strain_clusters_tsv, strain_protein_csv)
    else:
        logging.info(f'Using existing strain overview: {{strain_protein_csv}}')
    
    # PHAGE: Generate all_AA_seqs.faa and protein overview
    phage_clusters_tsv = os.path.join(iteration_output_dir, 'phage', 'clusters.tsv')
    phage_all_aa_fasta = os.path.join(all_families_dir, 'phage_all_AA_seqs.faa')
    phage_protein_csv = os.path.join(all_families_dir, 'phage_all_families_overview.csv')
    
    if not os.path.exists(phage_all_aa_fasta):
        success = generate_all_aa_seqs_fasta(
            clusters_tsv_path=phage_clusters_tsv,
            source_fasta_dir=phage_source_dir,
            output_fasta_path=phage_all_aa_fasta,
            genome_type='phage'
        )
        if not success:
            logging.error('Failed to generate phage all_AA_seqs.faa')
            sys.exit(1)
    else:
        logging.info(f'Using existing phage all_AA_seqs.faa: {{phage_all_aa_fasta}}')
    
    if not os.path.exists(phage_protein_csv):
        create_protein_overview_from_clusters(phage_clusters_tsv, phage_protein_csv)
    else:
        logging.info(f'Using existing phage overview: {{phage_protein_csv}}')
    
    strain_fasta = strain_all_aa_fasta
    protein_csv = strain_protein_csv
    phage_fasta = phage_all_aa_fasta
    protein_csv_phage = phage_protein_csv
    protein_families_file = strain_protein_csv
    
else:
    logging.info('=== Using PREDICTIVE families only (standard mode) ===')
    
    strain_fasta = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'predictive_proteins', 'strain_predictive_AA_seqs.faa')
    protein_csv = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'predictive_proteins', 'strain_predictive_feature_overview.csv')
    phage_fasta = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'predictive_proteins', 'phage_predictive_AA_seqs.faa')
    protein_csv_phage = os.path.join(iteration_output_dir, 'modeling_results', 'model_performance', 'predictive_proteins', 'phage_predictive_feature_overview.csv')
    protein_families_file = protein_csv
    
    # Validate predictive files exist
    for filepath in [strain_fasta, protein_csv, phage_fasta, protein_csv_phage]:
        if not os.path.exists(filepath):
            logging.error(f'Required predictive file not found: {{filepath}}')
            sys.exit(1)

# Validate prepared files exist
logging.info('Validating prepared input files...')
for desc, filepath in [('strain_fasta', strain_fasta), ('protein_csv', protein_csv), 
                       ('phage_fasta', phage_fasta), ('protein_csv_phage', protein_csv_phage)]:
    if not os.path.exists(filepath):
        logging.error(f'Required file {{desc}} not found: {{filepath}}')
        sys.exit(1)
    logging.info(f'✓ {{desc}}: {{filepath}}')

# Get best cutoff from protein family modeling
performance_df = pd.read_csv(metrics_file)
top_cutoff = performance_df.iloc[0]['cut_off'].split('_')[-1]
feature_file_path = os.path.join(iteration_output_dir, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_cutoff_{{top_cutoff}}.csv')

if not os.path.exists(feature_file_path):
    logging.error(f'Feature file not found: {{feature_file_path}}')
    sys.exit(1)

# ==================================================================================
# RUN K-MER TABLE WORKFLOW
# ==================================================================================
kmer_output_dir = os.path.join(iteration_output_dir, '{kmer_output_dirname}')
kmer_metrics_file = os.path.join(kmer_output_dir, 'modeling', 'modeling_results', 'model_performance', 'model_performance_metrics.csv')

if not os.path.exists(kmer_metrics_file):
    logging.info(f'Running k-mer table workflow in {{kmer_output_dir}}...')
    run_kmer_table_workflow(
        strain_fasta=strain_fasta,
        protein_csv=protein_csv,
        k=5,
        id_col='strain',
        one_gene=False,
        output_dir=kmer_output_dir,
        k_range=False,
        phenotype_matrix='{args.interaction_matrix}',
        phage_fasta=phage_fasta,
        protein_csv_phage=protein_csv_phage,
        remove_suffix=False,
        sample_column='{args.strain_column}',
        phenotype_column='interaction',
        modeling=True,
        filter_type='strain',
        num_features='none',
        num_runs_fs={args.num_runs_fs},
        num_runs_modeling={args.num_runs_modeling},
        method='rfe',
        strain_list=feature_file_path,
        phage_list=feature_file_path,
        threads={args.threads},
        task_type='classification',
        max_features='none',
        ignore_families=False,
        max_ram={args.max_ram},
        use_shap=False,
        use_dynamic_weights={args.use_dynamic_weights},
        weights_method='{args.weights_method}',
        use_clustering={args.use_clustering},
        cluster_method='{args.cluster_method}',
        n_clusters={args.n_clusters},
        min_cluster_size={args.min_cluster_size},
        min_samples={args.min_samples if args.min_samples else 'None'},
        cluster_selection_epsilon={args.cluster_selection_epsilon},
        check_feature_presence={args.check_feature_presence},
        filter_by_cluster_presence={args.filter_by_cluster_presence},
        min_cluster_presence={args.min_cluster_presence},
        use_feature_clustering={args.use_feature_clustering},
        feature_cluster_method='{args.feature_cluster_method}',
        feature_n_clusters={args.feature_n_clusters},
        feature_min_cluster_presence={args.feature_min_cluster_presence}
    )
else:
    logging.info(f'K-mer workflow already complete: {{kmer_metrics_file}}')

# ==================================================================================
# RUN K-MER ANALYSIS WORKFLOW
# ==================================================================================
best_cutoff = select_best_cutoff(os.path.join(kmer_output_dir, 'modeling'))
if not best_cutoff:
    logging.error('Could not determine best cutoff')
    sys.exit(1)

model_dir = os.path.join(kmer_output_dir, 'modeling', 'modeling_results', f'{{best_cutoff}}')
cutoff_num = best_cutoff.split('_')[1]

kmer_feature_file_path = os.path.join(kmer_output_dir, 'modeling', 'feature_selection', 'filtered_feature_tables', f'select_feature_table_cutoff_{{cutoff_num}}.csv')
feature2cluster_path = os.path.join(kmer_output_dir, 'feature_tables', 'selected_features.csv')

kmer_analysis_dir = os.path.join(kmer_output_dir, 'kmer_analysis')
filtered_kmers_path = os.path.join(kmer_analysis_dir, 'strain', 'filtered_kmers.csv')

aa_sequence_file = strain_fasta

if is_fasta_empty(aa_sequence_file):
    logging.warning('No sequences found. Skipping k-mer analysis.')
elif not os.path.exists(filtered_kmers_path):
    logging.info('Running k-mer analysis workflow...')
    kmer_analysis_workflow(
        aa_sequence_file=aa_sequence_file,
        feature_file_path=kmer_feature_file_path,
        feature2cluster_path=feature2cluster_path,
        protein_families_file=protein_families_file,
        output_dir=kmer_analysis_dir,
        feature_type='strain',
        quick_run={str(args.use_all_families_for_kmers).lower()}
    )
else:
    logging.info(f'K-mer analysis already complete: {{filtered_kmers_path}}')

# ==================================================================================
# RUN PREDICTION ON VALIDATION STRAINS
# ==================================================================================
validation_output_dir = os.path.join(kmer_output_dir, 'model_validation')
validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')
os.makedirs(validation_output_dir, exist_ok=True)

# Create symlink to tmp directory
protein_tmp_dir = os.path.abspath(os.path.join(iteration_output_dir, 'tmp'))
validation_tmp_dir_abs = os.path.abspath(validation_tmp_dir)

if os.path.exists(protein_tmp_dir) and not os.path.exists(validation_tmp_dir_abs):
    try:
        os.makedirs(os.path.dirname(validation_tmp_dir_abs), exist_ok=True)
        os.symlink(protein_tmp_dir, validation_tmp_dir_abs, target_is_directory=True)
        logging.info(f'Created symlink: {{validation_tmp_dir_abs}} -> {{protein_tmp_dir}}')
    except Exception as e:
        logging.warning(f'Could not create symlink: {{e}}, creating directory instead')
        os.makedirs(validation_tmp_dir_abs, exist_ok=True)
else:
    os.makedirs(validation_tmp_dir_abs, exist_ok=True)

median_predictions_file = os.path.join(validation_output_dir, 'predict_results', 'strain_median_predictions.csv')

if os.path.exists(median_predictions_file):
    logging.info(f'Predictions already complete: {{median_predictions_file}}')
elif is_fasta_empty(aa_sequence_file):
    logging.warning('No sequences. Creating empty feature table for prediction.')
    
    validation_strains_list = validation_strains_df['strain'].tolist()
    empty_feature_table = pd.DataFrame({{'strain': validation_strains_list}})
    
    assign_results_dir = os.path.join(validation_output_dir, 'assign_results')
    os.makedirs(assign_results_dir, exist_ok=True)
    empty_feature_table.to_csv(os.path.join(assign_results_dir, 'strain_combined_feature_table.csv'), index=False)
    
    predict_output_dir = os.path.join(validation_output_dir, 'predict_results')
    os.makedirs(predict_output_dir, exist_ok=True)
    
    run_prediction_workflow(
        input_dir=assign_results_dir,
        phage_feature_table_path=os.path.join(kmer_output_dir, 'feature_tables', 'phage_final_feature_table.csv'),
        model_dir=model_dir,
        output_dir=predict_output_dir
    )
else:
    logging.info('Running k-mer assignment and prediction workflow...')
    kmer_assign_predict_workflow(
        input_dir=strain_source_dir,
        genome_list=validation_strains_path,
        mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
        clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
        feature_map=os.path.join(kmer_output_dir, 'feature_tables', 'selected_features.csv'),
        filtered_kmers=filtered_kmers_path,
        aa_sequence_file=aa_sequence_file,
        tmp_dir=validation_tmp_dir_abs,
        output_dir=validation_output_dir,
        model_dir=model_dir,
        phage_feature_table_path=os.path.join(kmer_output_dir, 'feature_tables', 'phage_final_feature_table.csv'),
        genome_type='strain',
        threads={args.threads},
        reuse_existing=True
    )

logging.info(f'Iteration {{iteration}} completed successfully!')
logging.info(f'Results: {{median_predictions_file}}')
"

echo "Iteration $SLURM_ARRAY_TASK_ID completed at $(date)"
"""
    
    script_path = os.path.join(run_dir, "kmer_bootstrap_job_array.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    return script_path

def create_final_aggregation_job(args, run_dir, dependency_job_id):
    """Create job to aggregate results from all iterations"""
    
    # Determine output directory name based on mode
    if args.use_all_families_for_kmers:
        kmer_output_dirname = "kmer_modeling_all_families"
        output_prefix = "all_families_"
    else:
        kmer_output_dirname = "kmer_modeling"
        output_prefix = "kmer_"
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=kmer_aggregate
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:{dependency_job_id}
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err

echo "=== Aggregating K-mer Bootstrap Results ==="
echo "Started: $(date)"

module load anaconda3
conda activate {args.environment}

python3 << 'EOF'
import os
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

output_dir = '{args.output_dir}'
n_iterations = {args.n_iterations}
kmer_dir = '{kmer_output_dirname}'

all_predictions = []
for i in range(1, n_iterations + 1):
    pred_file = os.path.join(output_dir, f'iteration_{{i}}', kmer_dir, 'model_validation', 
                            'predict_results', 'strain_median_predictions.csv')
    
    if os.path.exists(pred_file):
        df = pd.read_csv(pred_file)
        df['iteration'] = i
        all_predictions.append(df)
        n_preds = len(df)
        logging.info(f'Loaded iteration {{i}}: {{n_preds}} predictions')
    else:
        logging.warning(f'Missing iteration {{i}}: {{pred_file}}')

if not all_predictions:
    logging.error('No prediction files found!')
    exit(1)

combined_df = pd.concat(all_predictions, ignore_index=True)
n_combined = len(combined_df)
n_iters = len(all_predictions)
logging.info(f'Combined {{n_combined}} total predictions from {{n_iters}} iterations')

combined_file = os.path.join(output_dir, '{output_prefix}final_predictions.csv')
combined_df.to_csv(combined_file, index=False)
logging.info(f'Saved combined predictions: {{combined_file}}')

summary_stats = combined_df.groupby(['strain', 'phage']).agg({{
    'pred_interaction': ['mean', 'std', 'min', 'max', 'count']
}}).reset_index()

summary_stats.columns = ['strain', 'phage', 'mean_prediction', 'std_prediction', 
                         'min_prediction', 'max_prediction', 'n_predictions']

summary_file = os.path.join(output_dir, '{output_prefix}prediction_summary.csv')
summary_stats.to_csv(summary_file, index=False)
logging.info(f'Saved summary statistics: {{summary_file}}')

print(f"\\n=== Aggregation Complete ===")
n_pairs = len(summary_stats)
avg_pred_per_pair = summary_stats['n_predictions'].mean()
print(f"Strain-phage pairs: {{n_pairs}}")
print(f"Avg predictions per pair: {{avg_pred_per_pair:.1f}}")
print(f"\\nOutput files:")
print(f"  All predictions: {{combined_file}}")
print(f"  Summary stats: {{summary_file}}")

EOF

echo "Aggregation completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "aggregate_results.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    return script_path

def main():
    parser = argparse.ArgumentParser(
        description="K-mer bootstrap validation with ALL families support (works with existing results)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script runs k-mer modeling on iterations where protein family workflow is already complete.

Example usage:
  # Run with ALL families mode:
  python3 kmer_bootstrap_all_families_FINAL.py \\
    --input_strain_dir /path/to/strains \\
    --input_phage_dir /path/to/phages \\
    --interaction_matrix /path/to/matrix.csv \\
    --output_dir /path/to/existing/output \\
    --use_all_families_for_kmers \\
    --mem_per_job 120 \\
    --time_limit "24:00:00"

  # Run with predictive families only:
  python3 kmer_bootstrap_all_families_FINAL.py \\
    --input_strain_dir /path/to/strains \\
    --input_phage_dir /path/to/phages \\
    --interaction_matrix /path/to/matrix.csv \\
    --output_dir /path/to/existing/output
        """
    )
    
    # Required
    parser.add_argument('--input_strain_dir', required=True, help='Original strain FASTA directory')
    parser.add_argument('--input_phage_dir', required=True, help='Original phage FASTA directory')
    parser.add_argument('--interaction_matrix', required=True, help='Interaction matrix CSV')
    parser.add_argument('--output_dir', required=True, help='Output directory (should already contain iteration_N/ dirs)')
    
    # Bootstrap
    parser.add_argument('--n_iterations', type=int, default=20, help='Number of iterations (default: 20)')
    parser.add_argument('--strain_column', default='strain', help='Strain column name (default: strain)')
    parser.add_argument('--threads', type=int, default=16, help='Threads per job (default: 16)')
    
    # Feature selection/modeling
    parser.add_argument('--num_runs_fs', type=int, default=25, help='Feature selection runs (default: 25)')
    parser.add_argument('--num_runs_modeling', type=int, default=50, help='Modeling runs (default: 50)')
    parser.add_argument('--use_dynamic_weights', action='store_true', help='Use dynamic weights')
    parser.add_argument('--weights_method', default='log10', choices=['log10', 'inverse_frequency', 'balanced'])
    parser.add_argument('--use_clustering', action='store_true', help='Use clustering')
    parser.add_argument('--cluster_method', default='hdbscan', choices=['hdbscan', 'hierarchical'])
    parser.add_argument('--n_clusters', type=int, default=20, help='Number of clusters (default: 20)')
    parser.add_argument('--min_cluster_size', type=int, default=2, help='Min cluster size (default: 2)')
    parser.add_argument('--min_samples', type=int, help='Min samples for clustering')
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0)
    parser.add_argument('--check_feature_presence', action='store_true')
    parser.add_argument('--filter_by_cluster_presence', action='store_true')
    parser.add_argument('--min_cluster_presence', type=int, default=2)
    parser.add_argument('--duplicate_all', action='store_true')
    parser.add_argument('--max_ram', type=float, default=40, help='Max RAM in GB (default: 40)')
    
    # Pre-processing clustering
    parser.add_argument('--use_feature_clustering', action='store_true')
    parser.add_argument('--feature_cluster_method', default='hierarchical')
    parser.add_argument('--feature_n_clusters', type=int, default=20)
    parser.add_argument('--feature_min_cluster_presence', type=int, default=2)
    
    # 🆕 KEY PARAMETER
    parser.add_argument('--use_all_families_for_kmers', action='store_true',
                       help='🆕 Use ALL protein families for k-mer modeling (creates kmer_modeling_all_families/ dir)')
    
    # SLURM
    parser.add_argument('--account', default='ac_mak', help='SLURM account')
    parser.add_argument('--partition', default='lr7', help='SLURM partition')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment')
    parser.add_argument('--mem_per_job', type=int, default=60, help='Memory per job in GB (default: 60)')
    parser.add_argument('--time_limit', default='12:00:00', help='Time limit (default: 12:00:00)')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts without submitting')
    
    args = parser.parse_args()
    
    # Validation
    if not os.path.exists(args.output_dir):
        print(f"Error: Output directory not found: {args.output_dir}")
        print("This script expects an existing output directory with completed protein family workflows.")
        return 1
    
    # Check for iteration directories
    iteration_dirs = [d for d in os.listdir(args.output_dir) if d.startswith('iteration_')]
    if not iteration_dirs:
        print(f"Error: No iteration directories found in {args.output_dir}")
        print("Expected to find directories like iteration_1/, iteration_2/, etc.")
        return 1
    
    print(f"✓ Found {len(iteration_dirs)} iteration directories")
    
    # Check a sample iteration for required files
    sample_iter = os.path.join(args.output_dir, 'iteration_1')
    required_files = {
        'modeling_strains.csv': os.path.join(sample_iter, 'modeling_strains.csv'),
        'validation_strains.csv': os.path.join(sample_iter, 'validation_strains.csv'),
        'metrics_file': os.path.join(sample_iter, 'modeling_results', 'model_performance', 'model_performance_metrics.csv'),
        'strain_clusters': os.path.join(sample_iter, 'strain', 'clusters.tsv'),
        'phage_clusters': os.path.join(sample_iter, 'phage', 'clusters.tsv')
    }
    
    missing_files = []
    for desc, filepath in required_files.items():
        if not os.path.exists(filepath):
            missing_files.append(f"{desc}: {filepath}")
    
    if missing_files:
        print(f"\\nError: Required files missing in {sample_iter}:")
        for mf in missing_files:
            print(f"  ✗ {mf}")
        print("\\nPlease ensure protein family workflow has completed successfully.")
        return 1
    
    print("✓ Protein family workflow results validated")
    
    if args.use_all_families_for_kmers:
        print("\\n⚠️  ALL FAMILIES MODE ENABLED")
        print(f"   Will create: {args.output_dir}/iteration_N/kmer_modeling_all_families/")
        print("   Processing ALL protein families for k-mer modeling")
        print(f"   Output files: {args.output_dir}/all_families_final_predictions.csv")
        if args.mem_per_job < 120:
            print(f"   ⚠️  Recommend: --mem_per_job 120 (currently {args.mem_per_job}GB)")
        if args.time_limit < "24:00:00":
            print(f"   ⚠️  Recommend: --time_limit 24:00:00 (currently {args.time_limit})")
    else:
        print("\\n📋 PREDICTIVE FAMILIES MODE")
        print(f"   Will create: {args.output_dir}/iteration_N/kmer_modeling/")
        print(f"   Output files: {args.output_dir}/kmer_final_predictions.csv")
    
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"kmer_bootstrap_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"\\n=== K-mer Bootstrap Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Iterations: {args.n_iterations}")
    print(f"Mode: {'ALL FAMILIES' if args.use_all_families_for_kmers else 'PREDICTIVE ONLY'}")
    print(f"Memory: {args.mem_per_job}GB, Time: {args.time_limit}")
    
    print("\\nCreating SLURM scripts...")
    bootstrap_script = create_kmer_bootstrap_job_array(args, run_dir)
    aggregate_script = create_final_aggregation_job(args, run_dir, "PLACEHOLDER")
    
    if args.dry_run:
        print("\\n✓ Dry run - scripts created but not submitted")
        print(f"Scripts in: {run_dir}/")
        return 0
    
    original_dir = os.getcwd()
    os.chdir(run_dir)
    
    print("\\nSubmitting jobs...")
    bootstrap_job_id = submit_job("kmer_bootstrap_job_array.sh")
    print(f"Bootstrap array: {bootstrap_job_id}")
    
    if bootstrap_job_id:
        with open("aggregate_results.sh", 'r') as f:
            content = f.read()
        content = content.replace("PLACEHOLDER", bootstrap_job_id)
        with open("aggregate_results.sh", 'w') as f:
            f.write(content)
        
        aggregate_job_id = submit_job("aggregate_results.sh")
        print(f"Aggregation: {aggregate_job_id}")
    
    os.chdir(original_dir)
    
    output_prefix = "all_families_" if args.use_all_families_for_kmers else "kmer_"
    
    print(f"\\n=== Submission Complete ===")
    print(f"Monitor: squeue -u $USER")
    print(f"Logs: {run_dir}/logs/")
    print(f"Results will be in:")
    print(f"  {args.output_dir}/{output_prefix}final_predictions.csv")
    print(f"  {args.output_dir}/{output_prefix}prediction_summary.csv")

if __name__ == "__main__":
    main()