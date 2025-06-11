#!/usr/bin/env python3
"""
Run only the prediction step for phage-based bootstrap cross-validation iterations.
This script runs the validation prediction step for iterations that have completed modeling
but failed or haven't completed the prediction step.
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

def find_incomplete_iterations(output_dir, n_iterations):
    """Find iterations that have completed modeling but not prediction."""
    incomplete_iterations = []
    
    for i in range(1, n_iterations + 1):
        iteration_dir = os.path.join(output_dir, f'iteration_{i}')
        
        # Check if modeling is complete
        metrics_file = os.path.join(iteration_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
        modeling_complete = os.path.exists(metrics_file)
        
        # Check if prediction is complete
        prediction_file = os.path.join(iteration_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')
        prediction_complete = os.path.exists(prediction_file)
        
        if modeling_complete and not prediction_complete:
            incomplete_iterations.append(i)
            print(f"Iteration {i}: Modeling ✓, Prediction ✗ - Will run prediction")
        elif modeling_complete and prediction_complete:
            print(f"Iteration {i}: Modeling ✓, Prediction ✓ - Complete")
        else:
            print(f"Iteration {i}: Modeling ✗ - Cannot run prediction without completed modeling")
    
    return incomplete_iterations

def create_prediction_job_array(args, run_dir, incomplete_iterations):
    """Create SLURM job array script for prediction-only iterations"""
    
    if not incomplete_iterations:
        print("No incomplete iterations found. Nothing to run.")
        return None
    
    # Convert list to comma-separated string for SLURM array
    array_indices = ','.join(map(str, incomplete_iterations))
    
    # Get the absolute path to the original script directory for imports
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=prediction_only
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_per_job}G
#SBATCH --time={args.time_limit}
#SBATCH --array={array_indices}
#SBATCH --output=logs/prediction_only_%A_%a.out
#SBATCH --error=logs/prediction_only_%A_%a.err

echo "=== Prediction-Only Step - Iteration $SLURM_ARRAY_TASK_ID ==="
echo "Job: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID, Node: $SLURMD_NODENAME, Started: $(date)"

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
import shutil

def select_best_cutoff(output_dir):
    '''Select the cutoff with highest MCC from modeling results'''
    metrics_file = os.path.join(output_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
    metrics_df = pd.read_csv(metrics_file)
    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    best_cutoff = metrics_df['cut_off'].values[0]
    return best_cutoff

# Get iteration number from SLURM array task ID
iteration = int(os.environ['SLURM_ARRAY_TASK_ID'])
print(f'Running prediction for iteration {{iteration}}')

# Parameters
output_dir = '{args.output_dir}'
iteration_output_dir = os.path.join(output_dir, f'iteration_{{iteration}}')
median_predictions_file = os.path.join(iteration_output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')

# Skip if already completed
if os.path.exists(median_predictions_file):
    print(f'Iteration {{iteration}} prediction is already complete, skipping.')
    sys.exit(0)

# Verify modeling is complete
metrics_file = os.path.join(iteration_output_dir, 'modeling_results/model_performance/model_performance_metrics.csv')
if not os.path.exists(metrics_file):
    print(f'ERROR: Modeling not complete for iteration {{iteration}}. Cannot run prediction.')
    sys.exit(1)

print(f'Starting prediction for iteration {{iteration}}...')

# Get the cutoff with the highest MCC
best_cutoff = select_best_cutoff(iteration_output_dir)
model_dir = os.path.join(iteration_output_dir, f'modeling_results', f'{{best_cutoff}}')
print(f'Using best cutoff: {{best_cutoff}}')

# Import workflow functions
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

# Verify required files exist
validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
assign_output_dir = os.path.join(validation_output_dir, 'assign_results')
validation_phage_features_path = os.path.join(assign_output_dir, 'phage_combined_feature_table.csv')
strain_feature_table_path = os.path.join(iteration_output_dir, 'strain', 'features', 'feature_table.csv')
select_feature_table = os.path.join(iteration_output_dir, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_{{best_cutoff}}.csv')

# Safety checks - ensure required files exist
missing_files = []
if not os.path.exists(validation_phage_features_path):
    missing_files.append(f'Phage features: {{validation_phage_features_path}}')
if not os.path.exists(strain_feature_table_path):
    missing_files.append(f'Strain features: {{strain_feature_table_path}}')
if not os.path.exists(select_feature_table):
    missing_files.append(f'Selected features: {{select_feature_table}}')

if missing_files:
    print('ERROR: Required files not found for prediction:')
    for missing in missing_files:
        print(f'  {{missing}}')
    print('Complete the full bootstrap workflow first.')
    sys.exit(1)

print('✓ All required files found. Running prediction...')

# Run prediction
strain_feature_table_path = os.path.join(iteration_output_dir, 'strain', 'features', 'feature_table.csv')

# Create properly named strain features for prediction workflow
predict_input_dir = os.path.join(validation_output_dir, 'prediction_input')
os.makedirs(predict_input_dir, exist_ok=True)
strain_features_for_prediction = os.path.join(predict_input_dir, 'strain_feature_table.csv')

# Only copy if not already exists (safety check)
if not os.path.exists(strain_features_for_prediction):
    shutil.copy2(strain_feature_table_path, strain_features_for_prediction)
    print(f'Copied strain features for prediction workflow')

predict_output_dir = os.path.join(validation_output_dir, 'predict_results')
run_prediction_workflow(
    input_dir=predict_input_dir,
    phage_feature_table_path=validation_phage_features_path,
    model_dir=model_dir,
    output_dir=predict_output_dir,
    feature_table=select_feature_table,
    strain_source='strain',
    phage_source='phage', 
    threads={args.threads}
)

print(f'Prediction for iteration {{iteration}} completed successfully')
"

echo "Iteration $SLURM_ARRAY_TASK_ID prediction completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, "prediction_only_job_array.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def main():
    parser = argparse.ArgumentParser(description="Run prediction-only step for incomplete bootstrap validation iterations")
    
    # Required arguments
    parser.add_argument('--output_dir', type=str, required=True, help="Directory with bootstrap validation results.")
    parser.add_argument('--n_iterations', type=int, required=True, help="Total number of iterations to check.")
    
    # SLURM-specific arguments  
    parser.add_argument('--account', default='ac_mak', help='SLURM account (default: ac_mak).')
    parser.add_argument('--partition', default='lr7', help='SLURM partition (default: lr7).')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS (default: lr_normal).')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment name (default: phage_modeling).')
    parser.add_argument('--mem_per_job', type=int, default=32, help='Memory per job in GB (default: 32).')
    parser.add_argument('--time_limit', default='2:00:00', help='Time limit per job (default: 2:00:00).')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads to use (default: 8).')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts but do not submit jobs.')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.output_dir):
        print(f"Error: Output directory not found: {args.output_dir}")
        return 1
    
    print(f"=== Prediction-Only Job Submission ===")
    print(f"Output directory: {args.output_dir}")
    print(f"Checking {args.n_iterations} iterations...")
    print()
    
    # Find iterations that need prediction
    incomplete_iterations = find_incomplete_iterations(args.output_dir, args.n_iterations)
    
    if not incomplete_iterations:
        print("All iterations are complete or don't have required modeling results.")
        print("Nothing to submit.")
        return 0
    
    print(f"\nFound {len(incomplete_iterations)} iterations needing prediction:")
    print(f"Iterations: {incomplete_iterations}")
    print()
    
    # Create timestamped run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"prediction_only_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"Run directory: {run_dir}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    print()
    
    # Create script
    script_path = create_prediction_job_array(args, run_dir, incomplete_iterations)
    
    if script_path is None:
        return 0
    
    if args.dry_run:
        print("Dry run - script created but not submitted")
        print(f"Script: {script_path}")
        return 0
    
    # Change to run directory and submit
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    print("Submitting prediction-only job array...")
    job_id = submit_job("prediction_only_job_array.sh")
    print(f"Job array submitted: {job_id}")
    
    # Change back to original directory
    os.chdir(original_dir)
    
    print(f"\n=== Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Iterations to process: {incomplete_iterations}")
    print(f"Jobs submitted: {len(incomplete_iterations)}")
    print(f"Expected runtime: {args.time_limit} per job (parallel)")
    print("\nMonitor with:")
    print("  squeue -u $USER")
    print("  squeue -u $USER --name=prediction_only")
    print("  tail -f logs/prediction_only_*.out")
    print(f"\nResults will be saved to existing iteration directories")
    print(f"Final predictions can be aggregated using the main bootstrap aggregation script")

if __name__ == "__main__":
    main()