#!/usr/bin/env python3
"""
Usage script for running prediction-only step on incomplete bootstrap iterations.
Edit the paths below, then run: python3 run_prediction_only_usage.py
"""

import subprocess
import sys
import os

def main():
    # =============================================
    # YOUR PATHS - EDIT THESE
    # =============================================
    output_dir = "/global/scratch/users/anoonan/BRaVE/ecoli/ecoli_modeling_phage_cv_hierarchical"
    n_iterations = 20  # Total number of iterations to check
    
    # =============================================
    # SLURM CONFIGURATION
    # =============================================
    account = "ac_mak"
    partition = "lr7"
    qos = "lr_normal"
    environment = "phage_modeling"
    mem_per_job = "32"          # Reduced memory for prediction-only
    time_limit = "2:00:00"      # Shorter time for prediction-only
    threads = "8"               # Fewer threads needed
    
    # =============================================
    # OTHER OPTIONS
    # =============================================
    dry_run = False  # Set to True to just create scripts without submitting
    
    # =============================================
    # BUILD COMMAND
    # =============================================
    cmd = [
        "python3", "run_prediction_only.py",
        
        # Required arguments
        "--output_dir", output_dir,
        "--n_iterations", str(n_iterations),
        
        # SLURM configuration
        "--account", account,
        "--partition", partition,
        "--qos", qos,
        "--environment", environment,
        "--mem_per_job", mem_per_job,
        "--time_limit", time_limit,
        "--threads", threads,
    ]
    
    # Add boolean flags
    if dry_run:
        cmd.append("--dry_run")
    
    # =============================================
    # SUBMIT JOBS
    # =============================================
    print("=" * 60)
    print("Prediction-Only Job Submission")
    print("=" * 60)
    print(f"Output directory: {output_dir}")
    print(f"Checking {n_iterations} iterations for incomplete predictions")
    print()
    print(f"SLURM settings:")
    print(f"  Account:     {account}")
    print(f"  Partition:   {partition}")
    print(f"  Environment: {environment}")
    print(f"  Memory:      {mem_per_job}GB per job")
    print(f"  Time limit:  {time_limit} per job")
    print(f"  Threads:     {threads} per job")
    print()
    
    if dry_run:
        print("🧪 DRY RUN MODE - Scripts will be created but not submitted")
        print()
    
    print("Running command:")
    print(" ".join(cmd))
    print()
    
    # Validate paths
    if not os.path.exists(output_dir):
        print(f"❌ ERROR: Output directory not found: {output_dir}")
        return 1
    
    try:
        subprocess.run(cmd, check=True)
        
        if dry_run:
            print("\n✅ Dry run completed successfully!")
            print("Scripts created but not submitted")
        else:
            print("\n✅ Prediction-only jobs submitted successfully!")
            print("\n📋 Monitor progress with:")
            print("   squeue -u $USER")
            print("   squeue -u $USER --name=prediction_only")
            print("   tail -f prediction_only_run_*/logs/prediction_only_*.out")
            print("\n💡 What this does:")
            print("   - Checks all iterations for completed modeling but incomplete prediction")
            print("   - Submits lightweight jobs to run only the prediction step")
            print("   - Much faster than re-running entire iterations")
            print("   - Results integrate seamlessly with existing bootstrap results")
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error submitting jobs: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())