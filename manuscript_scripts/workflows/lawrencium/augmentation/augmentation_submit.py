#!/usr/bin/env python3
"""
SLURM workflow submission for augmentation testing on existing CV results.
Creates a job array where each job tests one augmentation fraction on one CV iteration.
Jobs run in parallel, testing all fraction-iteration combinations.
"""

import os
import sys
import argparse
import subprocess
import time
import pandas as pd

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

def detect_iterations(cv_dir):
    """Auto-detect iteration directories from CV output"""
    iterations = []
    if not os.path.exists(cv_dir):
        return []
    
    for d in os.listdir(cv_dir):
        if d.startswith('iteration_') and os.path.isdir(os.path.join(cv_dir, d)):
            try:
                iter_num = int(d.split('_')[1])
                iterations.append(iter_num)
            except:
                pass
    return sorted(iterations)

def validate_iteration_structure(cv_dir, iterations, feature_table_path):
    """Validate that required files exist in iteration directories"""
    required_files = [
        feature_table_path,
        'validation_strains.csv',
        'strain/clusters.tsv', 
        'strain/features/selected_features.csv',
        'phage/features/feature_table.csv',
        'tmp/strain/mmseqs_db'
    ]
    
    missing_files = []
    valid_iterations = []
    
    for iteration in iterations:
        iter_dir = os.path.join(cv_dir, f'iteration_{iteration}')
        iter_valid = True
        
        # Check required files
        for req_file in required_files:
            file_path = os.path.join(iter_dir, req_file)
            if not os.path.exists(file_path):
                missing_files.append(f"iteration_{iteration}/{req_file}")
                iter_valid = False
        
        if iter_valid:
            valid_iterations.append(iteration)
    
    return valid_iterations, missing_files

def create_fractions_file(fractions, run_dir):
    """Create file with augmentation fractions for array job"""
    fractions_file = os.path.join(run_dir, 'augmentation_fractions.txt')
    with open(fractions_file, 'w') as f:
        for frac in fractions:
            f.write(f"{frac}\n")
    return fractions_file

def create_augmentation_job_array(args, run_dir, valid_iterations):
    """Create SLURM job array script for augmentation testing"""
    
    # Get the absolute path to the script directory for imports
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Calculate total tasks
    n_fractions = len(args.fractions)
    n_iterations = len(valid_iterations)
    total_tasks = n_fractions * n_iterations
    
    # Create fractions file
    fractions_file = create_fractions_file(args.fractions, run_dir)
    
    # Create iterations file
    iterations_file = os.path.join(run_dir, 'valid_iterations.txt')
    with open(iterations_file, 'w') as f:
        for iteration in valid_iterations:
            f.write(f"{iteration}\n")
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=aug_test
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_per_job}G
#SBATCH --time={args.time_limit}
#SBATCH --array=1-{total_tasks}
#SBATCH --output=logs/aug_%A_%a.out
#SBATCH --error=logs/aug_%A_%a.err

echo "=== Augmentation Test - Task $SLURM_ARRAY_TASK_ID of {total_tasks} ==="
echo "Job: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID, Node: $SLURMD_NODENAME, Started: $(date)"

# Calculate iteration and fraction indices
n_fractions={n_fractions}
n_iterations={n_iterations}
iteration_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) / $n_fractions + 1 ))
fraction_idx=$(( ($SLURM_ARRAY_TASK_ID - 1) % $n_fractions + 1 ))

# Get actual iteration number and fraction value  
iteration=$(sed -n "${iteration_idx}p" /{run_dir_abs}/valid_iterations.txt)
fraction=$(sed -n "${fraction_idx}p" /{run_dir_abs}/augmentation_fractions.txt)

echo "Processing iteration $iteration with augmentation fraction $fraction"

module load anaconda3
conda activate {args.environment} 2>&1 || {{
    echo "Direct activation failed, trying with conda init..."
    conda init bash >/dev/null 2>&1
    source ~/.bashrc >/dev/null 2>&1
    conda activate {args.environment}
}}

# Run the augmentation test directly
python3 -c "
import sys
import os
sys.path.insert(0, '{script_dir}')

import pandas as pd
from phage_modeling.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table
from phage_modeling.workflows.assign_predict_workflow import assign_predict_workflow

# Get parameters from environment
iteration = int('$iteration')
fraction = float('$fraction')
cv_dir = '{args.cv_dir}'

print(f'Running augmentation test for iteration {{iteration}}, fraction {{fraction}}')

# Set up paths - write directly into CV iteration directory
iteration_dir = os.path.join(cv_dir, f'iteration_{{iteration}}')
feature_table = os.path.join(iteration_dir, '{args.feature_table_path}')

# Create output directory with fraction percentage
frac_pct = int(fraction * 100)
output_dir = os.path.join(iteration_dir, 'augmentation', f'frac_{{frac_pct:02d}}pct')
os.makedirs(output_dir, exist_ok=True)

# Check if already completed
final_predictions_file = os.path.join(output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')
if os.path.exists(final_predictions_file):
    print(f'Task already completed for iteration {{iteration}}, fraction {{fraction}}. Skipping.')
    sys.exit(0)

# Validate feature table exists
if not os.path.exists(feature_table):
    print(f'ERROR: Feature table not found: {{feature_table}}')
    sys.exit(1)

print(f'Using feature table: {{feature_table}}')
print(f'Output directory: {{output_dir}}')

# Step 1: Run modeling workflow with augmentation
print('Running modeling workflow with augmentation...')
try:
    run_modeling_workflow_from_feature_table(
        full_feature_table=feature_table,
        output_dir=output_dir,
        threads={args.threads},
        num_features={args.num_features},
        filter_type='strain',
        num_runs_fs={args.num_runs_fs},
        num_runs_modeling={args.num_runs_modeling},
        sample_column='strain',
        phage_column='phage',
        phenotype_column='interaction',
        method='rfe',
        task_type='classification',
        binary_data=True,
        max_features='none',
        max_ram={args.max_ram},
        use_dynamic_weights={args.use_dynamic_weights},
        weights_method='{args.weights_method}',
        use_clustering={args.use_clustering},
        cluster_method='{args.cluster_method}',
        n_clusters={args.n_clusters},
        min_cluster_size={args.min_cluster_size},
        min_samples={args.min_samples if args.min_samples is not None else None},
        cluster_selection_epsilon={args.cluster_selection_epsilon},
        check_feature_presence={args.check_feature_presence},
        filter_by_cluster_presence={args.filter_by_cluster_presence},
        min_cluster_presence={args.min_cluster_presence},
        use_augmentation=True,
        augmentation_strain_fraction=fraction,
        augmentation_phage_fraction=fraction,
        augmentation_fold_increase={args.fold_increase},
        run_predictive_proteins=False  # Skip predictive proteins for speed
    )
except Exception as e:
    print(f'ERROR in modeling workflow: {{e}}')
    sys.exit(1)

# Step 2: Select best cutoff
print('Selecting best cutoff from modeling results...')
metrics_file = os.path.join(output_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
if not os.path.exists(metrics_file):
    print(f'ERROR: Model performance metrics not found: {{metrics_file}}')
    sys.exit(1)

metrics_df = pd.read_csv(metrics_file)
if len(metrics_df) == 0:
    print('ERROR: No model performance metrics found')
    sys.exit(1)

# Sort by MCC (descending) then by cutoff (descending for tie-breaking)
metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
best_cutoff = metrics_df['cut_off'].values[0]
best_mcc = metrics_df['MCC'].values[0]
print(f'Selected best cutoff: {{best_cutoff}} (MCC: {{best_mcc:.4f}})')

# Step 3: Prepare validation prediction
print('Setting up validation prediction...')
validation_output_dir = os.path.join(output_dir, 'model_validation')
os.makedirs(validation_output_dir, exist_ok=True)
validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')

# Get required paths from original iteration
validation_strains = os.path.join(iteration_dir, 'validation_strains.csv')
mmseqs_db = os.path.join(iteration_dir, 'tmp/strain/mmseqs_db')
clusters_tsv = os.path.join(iteration_dir, 'strain/clusters.tsv')
feature_map = os.path.join(iteration_dir, 'strain/features/selected_features.csv')
phage_feature_table = os.path.join(iteration_dir, 'phage/features/feature_table.csv')

# Check for modified AA directory
modified_aa_dir = os.path.join(iteration_dir, 'strain/modified_AAs/strain')
input_strain_dir = modified_aa_dir if os.path.exists(modified_aa_dir) else '{args.input_strain_dir}'

# Construct model directory
model_dir = os.path.join(output_dir, f'modeling_results/{{best_cutoff}}')
if not os.path.exists(model_dir):
    print(f'ERROR: Model directory not found: {{model_dir}}')
    sys.exit(1)

# Construct feature table path
select_feature_table = os.path.join(output_dir, f'feature_selection/filtered_feature_tables/select_feature_table_{{best_cutoff}}.csv')
if not os.path.exists(select_feature_table):
    print(f'ERROR: Selected feature table not found: {{select_feature_table}}')
    sys.exit(1)

print(f'Model directory: {{model_dir}}')
print(f'Feature table: {{select_feature_table}}')
print(f'Validation strains: {{validation_strains}}')

# Step 4: Run prediction workflow
print('Running prediction workflow...')
try:
    assign_predict_workflow(
        input_dir=input_strain_dir,
        genome_list=validation_strains,
        mmseqs_db=mmseqs_db,
        clusters_tsv=clusters_tsv,
        feature_map=feature_map,
        tmp_dir=validation_tmp_dir,
        suffix='faa',
        model_dir=model_dir,
        feature_table=select_feature_table,
        phage_feature_table_path=phage_feature_table,
        output_dir=validation_output_dir,
        threads={args.threads},
        genome_type='strain',
        sensitivity=7.5,
        coverage=0.8,
        min_seq_id=0.4,
        duplicate_all={args.duplicate_all}
    )
except Exception as e:
    print(f'ERROR in prediction workflow: {{e}}')
    sys.exit(1)

print(f'Augmentation test completed successfully for iteration {{iteration}}, fraction {{fraction}}')

# Save task summary
task_summary = {{
    'iteration': iteration,
    'augmentation_fraction': fraction,
    'best_cutoff': best_cutoff,
    'best_mcc': best_mcc,
    'completed': True
}}

summary_file = os.path.join(output_dir, 'task_summary.csv')
pd.DataFrame([task_summary]).to_csv(summary_file, index=False)
print(f'Task summary saved to: {{summary_file}}')
"

echo "Task $SLURM_ARRAY_TASK_ID completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "augmentation_job_array.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_aggregation_job(args, run_dir, dependency, valid_iterations):
    """Create job to aggregate results from all augmentation tests"""
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=aug_aggregate
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

echo "=== Augmentation Results Aggregation ==="
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
import numpy as np
import os
from sklearn.metrics import balanced_accuracy_score, f1_score, precision_score, recall_score, matthews_corrcoef

# Parameters
cv_dir = '{args.cv_dir}'
fractions = {args.fractions}
iterations = {valid_iterations}
n_fractions = len(fractions)
n_iterations = len(iterations)

print(f'Aggregating results from {{n_iterations}} iterations and {{n_fractions}} fractions...')

# Collect all results
all_results = []
completed_tasks = 0
missing_tasks = []

for iteration in iterations:
    for fraction in fractions:
        frac_pct = int(fraction * 100)
        
        # Paths for this iteration-fraction combination - now in CV dir
        result_dir = os.path.join(cv_dir, f'iteration_{{iteration}}', 'augmentation', f'frac_{{frac_pct:02d}}pct')
        
        # Check for validation predictions
        predictions_file = os.path.join(result_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')
        
        # Check for model metrics
        metrics_file = os.path.join(result_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
        
        # Check for task summary
        task_summary_file = os.path.join(result_dir, 'task_summary.csv')
        
        if os.path.exists(predictions_file) and os.path.exists(metrics_file):
            try:
                # Load predictions
                predictions = pd.read_csv(predictions_file)
                
                # Load model metrics (get best performing model)
                metrics = pd.read_csv(metrics_file)
                best_metrics = metrics.sort_values(['MCC', 'cut_off'], ascending=[False, False]).iloc[0]
                
                # Calculate validation performance
                y_true = predictions['interaction'] if 'interaction' in predictions.columns else predictions[predictions.columns[-1]]
                
                # Handle different prediction column names
                if 'Final_Prediction' in predictions.columns:
                    y_pred = predictions['Final_Prediction']
                elif 'Prediction' in predictions.columns:
                    y_pred = predictions['Prediction'] 
                else:
                    # Use median prediction based on confidence
                    y_pred = (predictions['Confidence'] > 0.5).astype(int) if 'Confidence' in predictions.columns else predictions.iloc[:, 1]
                
                if 'Confidence' in predictions.columns:
                    y_conf = predictions['Confidence']
                elif 'Final_Prediction' in predictions.columns:
                    y_conf = predictions['Final_Prediction'] 
                else:
                    y_conf = y_pred
                
                # Calculate metrics
                val_performance = {{
                    'iteration': iteration,
                    'augmentation_fraction': fraction,
                    'best_cutoff': best_metrics['cut_off'],
                    'train_accuracy': best_metrics.get('Accuracy', np.nan),
                    'train_f1': best_metrics.get('F1', np.nan),
                    'train_mcc': best_metrics.get('MCC', np.nan),
                    'val_balanced_accuracy': balanced_accuracy_score(y_true, y_pred),
                    'val_f1': f1_score(y_true, y_pred, zero_division=0),
                    'val_precision': precision_score(y_true, y_pred, zero_division=0),
                    'val_recall': recall_score(y_true, y_pred, zero_division=0),
                    'val_mcc': matthews_corrcoef(y_true, y_pred),
                    'n_val_samples': len(predictions),
                    'mean_confidence': y_conf.mean(),
                    'std_confidence': y_conf.std()
                }}
                
                all_results.append(val_performance)
                completed_tasks += 1
                print(f'✓ Added results from iteration {{iteration}}, fraction {{fraction:.3f}}')
                
            except Exception as e:
                print(f'✗ Error processing iteration {{iteration}}, fraction {{fraction:.3f}}: {{e}}')
                missing_tasks.append(f'iteration_{{iteration}}_frac_{{fraction:.3f}}')
        else:
            print(f'✗ Missing results for iteration {{iteration}}, fraction {{fraction:.3f}}')
            missing_tasks.append(f'iteration_{{iteration}}_frac_{{fraction:.3f}}')

print(f'\\nCompleted: {{completed_tasks}}/{{n_iterations * n_fractions}} tasks')
if missing_tasks:
    print(f'Missing: {{len(missing_tasks)}} tasks')
    print('Missing tasks:', missing_tasks[:10], '...' if len(missing_tasks) > 10 else '')

# Save results if we have any
if all_results:
    # Create main results DataFrame
    final_results = pd.DataFrame(all_results)
    final_results = final_results.sort_values(['augmentation_fraction', 'iteration'])
    
    # Save detailed results
    detailed_output = os.path.join(cv_dir, 'augmentation_test_detailed_results.csv')
    final_results.to_csv(detailed_output, index=False)
    print(f'\\nDetailed results saved to: {{detailed_output}}')
    
    # Create summary statistics by fraction
    summary_stats = final_results.groupby('augmentation_fraction').agg({{
        'val_balanced_accuracy': ['mean', 'std', 'count'],
        'val_f1': ['mean', 'std'],
        'val_mcc': ['mean', 'std'],
        'val_precision': ['mean', 'std'],
        'val_recall': ['mean', 'std'],
        'train_mcc': ['mean', 'std'],
        'n_val_samples': ['mean']
    }}).round(4)
    
    # Flatten column names
    summary_stats.columns = ['_'.join(col).strip() for col in summary_stats.columns.values]
    summary_stats = summary_stats.reset_index()
    
    # Save summary
    summary_output = os.path.join(cv_dir, 'augmentation_test_summary.csv')
    summary_stats.to_csv(summary_output, index=False)
    print(f'Summary statistics saved to: {{summary_output}}')
    
    # Print summary to console
    print('\\n' + '='*80)
    print('AUGMENTATION TEST SUMMARY')
    print('='*80)
    print('Summary by augmentation fraction:')
    print(summary_stats[['augmentation_fraction', 'val_balanced_accuracy_mean', 'val_balanced_accuracy_std', 
                       'val_f1_mean', 'val_f1_std', 'val_mcc_mean', 'val_mcc_std']].to_string(index=False))
    
    # Find best performing fraction
    best_fraction = summary_stats.loc[summary_stats['val_mcc_mean'].idxmax(), 'augmentation_fraction']
    best_mcc = summary_stats.loc[summary_stats['val_mcc_mean'].idxmax(), 'val_mcc_mean']
    print(f'\\nBest performing fraction: {{best_fraction:.3f}} (Mean MCC: {{best_mcc:.4f}})')
    
    print(f'\\nTotal results processed: {{len(final_results)}} out of {{n_iterations * n_fractions}} expected')
    
else:
    print('\\nERROR: No complete results found to aggregate!')
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
    parser = argparse.ArgumentParser(description="Submit augmentation testing workflow as SLURM job array")
    
    # Required paths
    parser.add_argument('--cv_dir', type=str, required=True, help="Directory containing CV iteration results")
    parser.add_argument('--input_strain_dir', type=str, required=True, help="Directory containing strain FASTA files")
    parser.add_argument('--input_phage_dir', type=str, required=True, help="Directory containing phage FASTA files")
    
    # Augmentation parameters
    parser.add_argument('--fractions', nargs='+', type=float, required=True, help="Augmentation fractions to test")
    parser.add_argument('--fold_increase', type=int, default=3, help="Fold increase for augmentation (default: 3)")
    
    # CV parameters
    parser.add_argument('--n_iterations', type=int, help="Number of iterations (auto-detected if not specified)")
    parser.add_argument('--feature_table_path', default='merged/full_feature_table.csv', 
                       help="Path to feature table within iteration directory (default: merged/full_feature_table.csv)")
    
    # Modeling parameters
    parser.add_argument('--threads', type=int, default=16, help="Number of threads (default: 16)")
    parser.add_argument('--num_runs_fs', type=int, default=25, help="Number of feature selection runs (default: 25)")
    parser.add_argument('--num_runs_modeling', type=int, default=50, help="Number of modeling runs (default: 50)")
    parser.add_argument('--max_ram', type=int, default=60, help="Maximum RAM in GB (default: 60)")
    parser.add_argument('--num_features', type=int, default=100, help="Number of features to select (default: 100)")
    
    # Feature selection parameters
    parser.add_argument('--use_dynamic_weights', action='store_true', help="Use dynamic weights for feature selection")
    parser.add_argument('--weights_method', type=str, default='log10', choices=['log10', 'inverse_frequency', 'balanced'], 
                       help="Method for calculating dynamic weights (default: log10)")
    parser.add_argument('--use_clustering', action='store_true', help="Use clustering for feature selection")
    parser.add_argument('--cluster_method', type=str, default='hdbscan', choices=['hdbscan', 'hierarchical'], 
                       help="Clustering method (default: hdbscan)")
    parser.add_argument('--n_clusters', type=int, default=20, help="Number of clusters for hierarchical clustering (default: 20)")
    parser.add_argument('--min_cluster_size', type=int, default=2, help="Minimum cluster size (default: 2)")
    parser.add_argument('--min_samples', type=int, help="Minimum number of samples for clustering")
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help="Epsilon value for clustering (default: 0.0)")
    parser.add_argument('--check_feature_presence', action='store_true', help="Check feature presence in train-test splits")
    parser.add_argument('--filter_by_cluster_presence', action='store_true', help="Filter features by cluster presence")
    parser.add_argument('--min_cluster_presence', type=int, default=2, help="Minimum cluster presence (default: 2)")
    parser.add_argument('--bootstrapping', action='store_true', help="Enable bootstrapping")
    parser.add_argument('--duplicate_all', action='store_true', help="Duplicate all genomes in predictions")
    
    # SLURM configuration
    parser.add_argument('--account', default='ac_mak', help='SLURM account (default: ac_mak)')
    parser.add_argument('--partition', default='lr7', help='SLURM partition (default: lr7)')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS (default: lr_normal)')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment (default: phage_modeling)')
    parser.add_argument('--mem_per_job', type=int, default=60, help='Memory per job in GB (default: 60)')
    parser.add_argument('--time_limit', default='12:00:00', help='Time limit per job (default: 12:00:00)')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts but do not submit jobs')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.cv_dir):
        print(f"Error: CV directory not found: {args.cv_dir}")
        return 1
    if not os.path.exists(args.input_strain_dir):
        print(f"Error: Strain directory not found: {args.input_strain_dir}")
        return 1
    if not os.path.exists(args.input_phage_dir):
        print(f"Error: Phage directory not found: {args.input_phage_dir}")
        return 1
    
    # Detect iterations
    if args.n_iterations:
        iterations = list(range(1, args.n_iterations + 1))
    else:
        iterations = detect_iterations(args.cv_dir)
        if not iterations:
            print(f"Error: No iterations found in {args.cv_dir}")
            return 1
        args.n_iterations = len(iterations)
    
    print(f"Found {len(iterations)} iterations: {iterations}")
    
    # Validate iteration structure
    valid_iterations, missing_files = validate_iteration_structure(args.cv_dir, iterations, args.feature_table_path)
    
    if not valid_iterations:
        print("Error: No valid iterations found!")
        print("Missing files:")
        for missing in missing_files[:20]:  # Show first 20
            print(f"  {missing}")
        return 1
    
    if len(valid_iterations) < len(iterations):
        print(f"Warning: Only {len(valid_iterations)}/{len(iterations)} iterations are valid")
        print(f"Valid iterations: {valid_iterations}")
    
    # Create timestamped run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"augmentation_test_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"=== Augmentation Testing SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"CV directory: {args.cv_dir}")
    print(f"Valid iterations: {len(valid_iterations)}")
    print(f"Augmentation fractions: {args.fractions}")
    print(f"Total tasks: {len(valid_iterations) * len(args.fractions)}")
    print(f"Threads per job: {args.threads}")
    print(f"Account: {args.account}, Environment: {args.environment}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    print()
    
    # Create scripts
    print("Creating SLURM job scripts...")
    array_script = create_augmentation_job_array(args, run_dir, valid_iterations)
    aggregate_script = create_aggregation_job(args, run_dir, "PLACEHOLDER", valid_iterations)
    
    if args.dry_run:
        print("Dry run - scripts created but not submitted")
        print("Scripts:")
        print(f"  Array job: {array_script}")
        print(f"  Aggregation: {aggregate_script}")
        return 0
    
    # Change to run directory
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    # Submit array job
    print("Submitting augmentation test array job...")
    array_job_id = submit_job("augmentation_job_array.sh")
    print(f"Array job: {array_job_id}")
    
    if array_job_id:
        # Update aggregation script dependency and submit
        with open("aggregate_results.sh", 'r') as f:
            content = f.read()
        content = content.replace("PLACEHOLDER", array_job_id)
        with open("aggregate_results.sh", 'w') as f:
            f.write(content)
        
        aggregate_job_id = submit_job("aggregate_results.sh")
        print(f"Aggregation job: {aggregate_job_id}")
    
    # Change back to original directory
    os.chdir(original_dir)
    
    print(f"\n=== Job Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Array job: {len(valid_iterations) * len(args.fractions)} parallel tasks")
    print(f"Expected total runtime: {args.time_limit} per task (parallel)")
    print("\nMonitor with:")
    print("  squeue -u $USER")
    print("  squeue -u $USER --name=aug_test  # Just augmentation jobs")
    print("  tail -f logs/aug_*.out")
    print(f"\nResults will be in:")
    print(f"  Individual tests: {args.cv_dir}/iteration_N/augmentation/frac_*pct/")
    print(f"  Summary: {args.cv_dir}/augmentation_test_summary.csv")
    print(f"  Detailed: {args.cv_dir}/augmentation_test_detailed_results.csv")

if __name__ == "__main__":
    main()