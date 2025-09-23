#!/usr/bin/env python3
"""
Configuration script for strain prediction workflow.
Predicts interactions for new strains using existing trained models and phage features.
Edit the paths and parameters below, then run: python3 strain_prediction_usage.py
"""

import subprocess
import sys

def main():
    # =============================================
    # YOUR PATHS - EDIT THESE
    # =============================================
    
    # New strain genomes to predict on
    new_strain_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/PSU_predictions/PSU_AAs"  # Directory with .faa files
    
    # Existing training results directory (contains strain/, phage/, modeling_results/, etc.)
    training_results_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/dataset_update/strain_modeling/ecoli_hierarchical_inv_freq_filter"
    
    # Output directory for predictions
    output_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/PSU_predictions/predictions"
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"                    # SLURM partition 
    qos = "lr_normal"                    # SLURM QOS
    environment = "phage_modeling"       # Conda environment name
    
    # =============================================
    # BATCHING PARAMETERS
    # =============================================
    strains_per_batch = 300             # Number of strains per batch job (adjust based on memory)
    max_parallel_jobs = 100               # Maximum number of parallel batch jobs
    
    # =============================================
    # PREDICTION PARAMETERS
    # =============================================
    
    # Resource parameters
    threads = "32"                       # Number of threads per job
    mem_per_job = "120"                   # Memory per job in GB
    time_limit = "6:00:00"               # Time limit per batch job
    
    # MMSeqs2 parameters for feature assignment
    sensitivity = "7.5"                  # Sensitivity for MMseqs2 search
    coverage = "0.8"                     # Minimum coverage for assignment
    min_seq_id = "0.4"                   # Minimum sequence identity for assignment
    
    # Workflow parameters
    use_best_cutoff = True               # Automatically use the best performing cutoff from training
    specific_cutoff = None               # Or specify a cutoff manually (e.g., "cutoff_15")
    duplicate_all = True                 # Duplicate all genomes in feature table for predictions
    
    # Debug options
    dry_run = False                      # Create scripts but don't submit jobs
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "strain_prediction_slurm.py",
        
        # Required arguments
        "--new_strain_dir", new_strain_dir,
        "--training_results_dir", training_results_dir,
        "--output_dir", output_dir,
        
        # SLURM configuration  
        "--account", account,
        "--partition", partition,
        "--qos", qos,
        "--environment", environment,
        "--mem_per_job", mem_per_job,
        "--time_limit", time_limit,
        
        # Batching parameters
        "--strains_per_batch", str(strains_per_batch),
        "--max_parallel_jobs", str(max_parallel_jobs),
        
        # Resource parameters
        "--threads", threads,
        
        # MMSeqs2 parameters
        "--sensitivity", sensitivity,
        "--coverage", coverage,
        "--min_seq_id", min_seq_id,
    ]
    
    # Add optional arguments
    if specific_cutoff:
        cmd.extend(["--specific_cutoff", specific_cutoff])
    
    # Add boolean flags
    if use_best_cutoff:
        cmd.append("--use_best_cutoff")
    if duplicate_all:
        cmd.append("--duplicate_all")
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # SUBMIT WORKFLOW
    # =============================================
    print("=" * 70)
    print("Strain Prediction Workflow Submission")
    print("=" * 70)
    print(f"New strain directory:    {new_strain_dir}")
    print(f"Training results:        {training_results_dir}")
    print(f"Output directory:        {output_dir}")
    print(f"Strains per batch:       {strains_per_batch}")
    print(f"Max parallel jobs:       {max_parallel_jobs}")
    print(f"Memory per job:          {mem_per_job}GB")
    print(f"Time limit per job:      {time_limit}")
    print(f"Use best cutoff:         {use_best_cutoff}")
    if specific_cutoff:
        print(f"Specific cutoff:         {specific_cutoff}")
    print()
    print(f"SLURM account:           {account}")
    print(f"Environment:             {environment}")
    print()
    
    if dry_run:
        print("🧪 DRY RUN MODE - Scripts will be created but not submitted")
        print()
    
    print("Submitting workflow with command:")
    print(" ".join(cmd))
    print()
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created in strain_prediction_run_* directory")
        else:
            print("\n✅ Strain prediction workflow submitted successfully!")
            print("\n📋 Monitor progress with:")
            print("  squeue -u $USER")
            print("  tail -f strain_prediction_run_*/logs/batch_*.out")
            print("\n⏱️  Expected workflow completion:")
            print(f"   Batch jobs will run in parallel (up to {max_parallel_jobs} at once)")
            print(f"   Each batch processes {strains_per_batch} strains")
            print(f"   Time per batch: ~{time_limit}")
            print(f"   Final results: {output_dir}/final_strain_predictions.csv")
            print("\n💡 What this does:")
            print("   1. Assigns features to new strains using existing clustering")
            print("   2. Combines assigned strain features with existing phage features")
            print("   3. Predicts strain-phage interactions using trained models")
            print("   4. Aggregates predictions from all batches")
            print("   5. Outputs median predictions for each strain-phage pair")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting workflow: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())