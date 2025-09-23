#!/usr/bin/env python3
"""
SLURM workflow submission for strain prediction using existing trained models.
Creates batched jobs to assign features to new strains and predict interactions.
"""

import os
import sys
import argparse
import subprocess
import time
import pandas as pd
import math

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

def validate_training_results(training_dir):
    """Validate that training results directory has expected structure."""
    required_paths = {
        'mmseqs_db': os.path.join(training_dir, 'tmp', 'strain', 'mmseqs_db'),
        'clusters_tsv': os.path.join(training_dir, 'strain', 'clusters.tsv'),
        'feature_map': os.path.join(training_dir, 'strain', 'features', 'selected_features.csv'),
        'phage_features': os.path.join(training_dir, 'phage', 'features', 'feature_table.csv'),
        'model_performance': os.path.join(training_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv'),
        'modeling_results': os.path.join(training_dir, 'modeling_results')
    }
    
    missing = []
    for name, path in required_paths.items():
        if not os.path.exists(path):
            missing.append(f"{name}: {path}")
    
    if missing:
        print("❌ Missing required files from training results:")
        for item in missing:
            print(f"   {item}")
        return False, {}
    
    print("✅ Training results validation passed")
    return True, required_paths

def get_best_cutoff(model_performance_file):
    """Get the best performing cutoff from model performance metrics."""
    try:
        df = pd.read_csv(model_performance_file)
        if 'MCC' in df.columns:
            # Sort by MCC (higher is better) then by cutoff
            df_sorted = df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
            best_cutoff = df_sorted.iloc[0]['cut_off']
            best_mcc = df_sorted.iloc[0]['MCC']
            print(f"✅ Best cutoff: {best_cutoff} (MCC: {best_mcc:.4f})")
            return best_cutoff
        else:
            print("❌ No MCC column found in model performance file")
            return None
    except Exception as e:
        print(f"❌ Error reading model performance file: {e}")
        return None

def get_strain_list(strain_dir):
    """Get list of strain names from .faa files in directory."""
    try:
        strain_files = [f for f in os.listdir(strain_dir) if f.endswith('.faa')]
        strains = ['.'.join(f.split('.')[:-1]) for f in strain_files]
        print(f"✅ Found {len(strains)} strain files in {strain_dir}")
        return strains
    except Exception as e:
        print(f"❌ Error reading strain directory: {e}")
        return []

def create_strain_batches(strains, strains_per_batch, max_parallel_jobs):
    """Split strains into batches for parallel processing."""
    total_strains = len(strains)
    ideal_batches = math.ceil(total_strains / strains_per_batch)
    actual_batches = min(ideal_batches, max_parallel_jobs)
    
    # Recalculate strains per batch if we're limited by max_parallel_jobs
    if actual_batches < ideal_batches:
        strains_per_batch = math.ceil(total_strains / actual_batches)
        print(f"⚠️  Adjusted to {strains_per_batch} strains per batch due to max_parallel_jobs limit")
    
    batches = []
    for i in range(0, total_strains, strains_per_batch):
        batch = strains[i:i + strains_per_batch]
        batches.append(batch)
    
    print(f"✅ Created {len(batches)} batches:")
    for i, batch in enumerate(batches):
        print(f"   Batch {i+1}: {len(batch)} strains")
    
    return batches

def create_batch_job_script(args, run_dir, batch_id, strain_batch, required_paths, cutoff):
    """Create SLURM job script for a single batch of strains."""
    
    # Get the absolute path to the original script directory for imports
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=strain_pred_batch_{batch_id}
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={args.threads}
#SBATCH --mem={args.mem_per_job}G
#SBATCH --time={args.time_limit}
#SBATCH --output=logs/batch_{batch_id}_%j.out
#SBATCH --error=logs/batch_{batch_id}_%j.err

echo "=== Strain Prediction Batch {batch_id} ==="
echo "Job: $SLURM_JOB_ID, Node: $SLURMD_NODENAME, Started: $(date)"
echo "Processing {len(strain_batch)} strains"

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

# Import required functions
from phage_modeling.workflows.assign_predict_workflow import assign_predict_workflow
import pandas as pd

# Create strain list file for this batch
batch_output_dir = '{args.output_dir}/batch_{batch_id}'
os.makedirs(batch_output_dir, exist_ok=True)

strain_list = {strain_batch}
strain_list_file = os.path.join(batch_output_dir, 'strain_list.csv')
strain_df = pd.DataFrame(strain_list, columns=['strain'])
strain_df.to_csv(strain_list_file, index=False)
print(f'Created strain list file with {{len(strain_list)}} strains: {{strain_list_file}}')

# Set up paths
model_dir = '{required_paths["modeling_results"]}/{cutoff}'
feature_table = None

# Check for feature selection table
feature_selection_dir = '{args.training_results_dir}/feature_selection/filtered_feature_tables'
cutoff_num = '{cutoff}'.split('_')[-1]
feature_table_path = os.path.join(feature_selection_dir, f'select_feature_table_cutoff_{{cutoff_num}}.csv')
if os.path.exists(feature_table_path):
    feature_table = feature_table_path
    print(f'Using feature selection table: {{feature_table}}')
else:
    print(f'No feature selection table found at {{feature_table_path}}, using all features')

print(f'Model directory: {{model_dir}}')
print(f'MMSeqs2 database: {required_paths["mmseqs_db"]}')
print(f'Clusters TSV: {required_paths["clusters_tsv"]}')
print(f'Feature map: {required_paths["feature_map"]}')
print(f'Phage features: {required_paths["phage_features"]}')

# Run assign_predict workflow for this batch
tmp_dir = os.path.join(batch_output_dir, 'tmp')
try:
    assign_predict_workflow(
        input_dir='{args.new_strain_dir}',
        mmseqs_db='{required_paths["mmseqs_db"]}',
        clusters_tsv='{required_paths["clusters_tsv"]}',
        feature_map='{required_paths["feature_map"]}',
        tmp_dir=tmp_dir,
        output_dir=batch_output_dir,
        model_dir=model_dir,
        feature_table=feature_table,
        strain_feature_table_path=None,
        phage_feature_table_path='{required_paths["phage_features"]}',
        genome_type='strain',
        genome_list=strain_list_file,
        sensitivity={args.sensitivity},
        coverage={args.coverage},
        min_seq_id={args.min_seq_id},
        threads={args.threads},
        suffix='faa',
        duplicate_all={args.duplicate_all}
    )
    print(f'Batch {batch_id} completed successfully')
except Exception as e:
    print(f'Error in batch {batch_id}: {{e}}')
    sys.exit(1)
"

echo "Batch {batch_id} completed: $(date)"
"""
    
    script_path = os.path.join(run_dir, f"batch_{batch_id}.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    return script_path

def create_aggregation_job(args, run_dir, dependencies, num_batches):
    """Create job to aggregate results from all batches."""
    
    dependency_str = ":".join(dependencies) if dependencies else ""
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=strain_pred_aggregate
#SBATCH --account={args.account}
#SBATCH --partition={args.partition}
#SBATCH --qos={args.qos}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --dependency=afterok:{dependency_str}
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err

echo "=== Strain Prediction Results Aggregation ==="
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
import glob

# Aggregate median predictions from all batches
output_dir = '{args.output_dir}'
final_predictions = pd.DataFrame()

print('Aggregating results from {num_batches} batches...')
completed_batches = 0

for batch_id in range(1, {num_batches} + 1):
    batch_dir = os.path.join(output_dir, f'batch_{{batch_id}}')
    pred_file = os.path.join(batch_dir, 'predict_results', 'strain_median_predictions.csv')
    
    if os.path.exists(pred_file):
        batch_predictions = pd.read_csv(pred_file)
        batch_predictions['batch_id'] = batch_id
        final_predictions = pd.concat([final_predictions, batch_predictions], ignore_index=True)
        completed_batches += 1
        print(f'Added results from batch {{batch_id}}: {{len(batch_predictions)}} predictions')
    else:
        print(f'Warning: Results missing for batch {{batch_id}} at {{pred_file}}')

# Save final aggregated predictions
if len(final_predictions) > 0:
    # Remove batch_id column for final output
    final_predictions_clean = final_predictions.drop(columns=['batch_id'])
    final_predictions_clean.to_csv(os.path.join(output_dir, 'final_strain_predictions.csv'), index=False)
    
    print(f'\\n✅ Final predictions saved with {{len(final_predictions_clean)}} total predictions from {{completed_batches}} batches')
    
    # Generate summary statistics
    if 'strain' in final_predictions_clean.columns and 'phage' in final_predictions_clean.columns:
        n_strains = final_predictions_clean['strain'].nunique()
        n_phages = final_predictions_clean['phage'].nunique()
        n_interactions = len(final_predictions_clean)
        avg_confidence = final_predictions_clean['Confidence'].mean()
        
        print(f'\\n📊 Summary Statistics:')
        print(f'   Unique strains: {{n_strains}}')
        print(f'   Unique phages: {{n_phages}}')
        print(f'   Total predictions: {{n_interactions}}')
        print(f'   Average confidence: {{avg_confidence:.4f}}')
        
        # Save summary
        summary_stats = {{
            'n_strains': n_strains,
            'n_phages': n_phages,
            'n_predictions': n_interactions,
            'avg_confidence': avg_confidence,
            'completed_batches': completed_batches,
            'total_batches': {num_batches}
        }}
        summary_df = pd.DataFrame([summary_stats])
        summary_df.to_csv(os.path.join(output_dir, 'prediction_summary.csv'), index=False)
        print(f'   Summary saved to prediction_summary.csv')
        
        # Save strain-level summary (how many phages each strain can infect)
        strain_summary = final_predictions_clean.groupby('strain').agg({{
            'Confidence': ['count', 'mean', 'std'],
            'Final_Prediction': 'sum'
        }}).round(4)
        strain_summary.columns = ['total_phages_tested', 'avg_confidence', 'std_confidence', 'predicted_infections']
        strain_summary = strain_summary.reset_index()
        strain_summary.to_csv(os.path.join(output_dir, 'strain_summary.csv'), index=False)
        print(f'   Strain-level summary saved to strain_summary.csv')
    
else:
    print('❌ ERROR: No results found to aggregate!')
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
    parser = argparse.ArgumentParser(description="Submit strain prediction workflow as batched SLURM jobs")
    
    # Required arguments
    parser.add_argument('--new_strain_dir', type=str, required=True, help="Directory containing new strain .faa files")
    parser.add_argument('--training_results_dir', type=str, required=True, help="Directory containing training results")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save prediction results")
    
    # SLURM configuration
    parser.add_argument('--account', default='ac_mak', help='SLURM account (default: ac_mak)')
    parser.add_argument('--partition', default='lr7', help='SLURM partition (default: lr7)')
    parser.add_argument('--qos', default='lr_normal', help='SLURM QOS (default: lr_normal)')
    parser.add_argument('--environment', default='phage_modeling', help='Conda environment (default: phage_modeling)')
    parser.add_argument('--mem_per_job', type=int, default=60, help='Memory per job in GB (default: 60)')
    parser.add_argument('--time_limit', default='4:00:00', help='Time limit per job (default: 4:00:00)')
    
    # Batching parameters
    parser.add_argument('--strains_per_batch', type=int, default=300, help='Number of strains per batch (default: 300)')
    parser.add_argument('--max_parallel_jobs', type=int, default=20, help='Maximum parallel jobs (default: 20)')
    
    # Workflow parameters
    parser.add_argument('--threads', type=int, default=16, help='Number of threads per job (default: 16)')
    parser.add_argument('--sensitivity', type=float, default=7.5, help='MMseqs2 sensitivity (default: 7.5)')
    parser.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage (default: 0.8)')
    parser.add_argument('--min_seq_id', type=float, default=0.4, help='Minimum sequence identity (default: 0.4)')
    parser.add_argument('--use_best_cutoff', action='store_true', help='Use best cutoff from training results')
    parser.add_argument('--specific_cutoff', type=str, help='Use specific cutoff (e.g., cutoff_15)')
    parser.add_argument('--duplicate_all', action='store_true', help='Duplicate all genomes for predictions')
    parser.add_argument('--dry_run', action='store_true', help='Create scripts but do not submit jobs')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.new_strain_dir):
        print(f"❌ New strain directory not found: {args.new_strain_dir}")
        return 1
    
    if not os.path.exists(args.training_results_dir):
        print(f"❌ Training results directory not found: {args.training_results_dir}")
        return 1
    
    # Validate training results structure
    valid, required_paths = validate_training_results(args.training_results_dir)
    if not valid:
        return 1
    
    # Get strain list
    strains = get_strain_list(args.new_strain_dir)
    if not strains:
        print("❌ No strain files found")
        return 1
    
    # Determine cutoff to use
    if args.specific_cutoff:
        cutoff = args.specific_cutoff
        model_dir = os.path.join(required_paths['modeling_results'], cutoff)
        if not os.path.exists(model_dir):
            print(f"❌ Specified cutoff directory not found: {model_dir}")
            return 1
        print(f"✅ Using specified cutoff: {cutoff}")
    elif args.use_best_cutoff:
        cutoff = get_best_cutoff(required_paths['model_performance'])
        if not cutoff:
            return 1
        model_dir = os.path.join(required_paths['modeling_results'], cutoff)
        if not os.path.exists(model_dir):
            print(f"❌ Best cutoff directory not found: {model_dir}")
            return 1
    else:
        print("❌ Must specify either --use_best_cutoff or --specific_cutoff")
        return 1
    
    # Create batches
    batches = create_strain_batches(strains, args.strains_per_batch, args.max_parallel_jobs)
    if not batches:
        print("❌ Could not create strain batches")
        return 1
    
    # Create timestamped run directory
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    run_dir = f"strain_prediction_run_{timestamp}"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
    
    print(f"\n=== Strain Prediction SLURM Submission ===")
    print(f"Run directory: {run_dir}")
    print(f"New strains: {len(strains)} strains in {len(batches)} batches")
    print(f"Using cutoff: {cutoff}")
    print(f"Output directory: {args.output_dir}")
    print(f"Memory per job: {args.mem_per_job}GB, Time limit: {args.time_limit}")
    print(f"Threads per job: {args.threads}")
    print()
    
    # Estimate resources
    total_core_hours = len(batches) * float(args.time_limit.split(':')[0]) * args.threads
    estimated_cost = total_core_hours * 0.10  # Rough estimate at $0.10/core-hour
    
    print(f"📊 Resource Estimates:")
    print(f"   Parallel batch jobs: {len(batches)}")
    print(f"   Total core-hours: ~{total_core_hours:.0f}")
    print(f"   Estimated cost: ~${estimated_cost:.0f}")
    print()
    
    if args.dry_run:
        print("🧪 DRY RUN MODE - Scripts will be created but not submitted")
        print()
    
    # Create scripts
    print("Creating SLURM job scripts...")
    batch_scripts = []
    for i, batch in enumerate(batches, 1):
        script = create_batch_job_script(args, run_dir, i, batch, required_paths, cutoff)
        batch_scripts.append(script)
    
    aggregate_script = create_aggregation_job(args, run_dir, ["PLACEHOLDER"] * len(batches), len(batches))
    
    if args.dry_run:
        print("Dry run completed - scripts created but not submitted")
        print("Scripts:")
        for i, script in enumerate(batch_scripts, 1):
            print(f"  Batch {i}: {script}")
        print(f"  Aggregation: {aggregate_script}")
        return 0
    
    # Change to run directory
    original_dir = os.getcwd()
    run_dir_abs = os.path.abspath(run_dir)
    os.chdir(run_dir)
    
    # Submit batch jobs
    print("Submitting batch jobs...")
    batch_job_ids = []
    for i, script_path in enumerate(batch_scripts, 1):
        script_name = os.path.basename(script_path)
        job_id = submit_job(script_name)
        if job_id:
            batch_job_ids.append(job_id)
            print(f"  Batch {i}: {job_id}")
        else:
            print(f"  Batch {i}: Failed to submit")
    
    # Submit aggregation job if batch jobs were submitted
    if batch_job_ids:
        # Update aggregation script with dependencies
        with open("aggregate_results.sh", 'r') as f:
            content = f.read()
        dependency_str = ":".join(batch_job_ids)
        content = content.replace("PLACEHOLDER", dependency_str)
        with open("aggregate_results.sh", 'w') as f:
            f.write(content)
        
        aggregate_job_id = submit_job("aggregate_results.sh")
        print(f"  Aggregation: {aggregate_job_id}")
    
    # Change back to original directory
    os.chdir(original_dir)
    
    print(f"\n=== Job Submission Summary ===")
    print(f"Run directory: {run_dir_abs}")
    print(f"Batch jobs submitted: {len(batch_job_ids)}/{len(batches)}")
    print(f"Expected completion time: {args.time_limit}")
    print("\nMonitor with:")
    print("  squeue -u $USER")
    print("  tail -f logs/batch_*.out")
    print(f"\nResults will be in:")
    print(f"  Individual batches: {args.output_dir}/batch_N/")
    print(f"  Final predictions: {args.output_dir}/final_strain_predictions.csv")
    print(f"  Summary statistics: {args.output_dir}/prediction_summary.csv")
    print(f"  Strain-level summary: {args.output_dir}/strain_summary.csv")
    print(f"\n💡 This workflow:")
    print(f"  1. Assigns features to new strains using existing clustering")
    print(f"  2. Combines with existing phage features")
    print(f"  3. Predicts strain-phage interactions using trained models")
    print(f"  4. Aggregates results from all batches")

if __name__ == "__main__":
    main()