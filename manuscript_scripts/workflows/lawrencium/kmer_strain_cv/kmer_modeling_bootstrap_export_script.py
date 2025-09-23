#!/usr/bin/env python3
"""
Script to export k-mer bootstrap validation results from HPC cluster.
Extracts key results and performance metrics while avoiding intermediate files.
Handles k-mer-based modeling workflow from all iterations.
"""

import argparse
import shutil
import sys
import pandas as pd
from datetime import datetime
from pathlib import Path


def get_directory_size(path):
    """Calculate total size of directory in human-readable format."""
    try:
        total_size = sum(f.stat().st_size for f in path.rglob('*') if f.is_file())
        # Convert to human-readable format
        for unit in ['B', 'KB', 'MB', 'GB']:
            if total_size < 1024.0:
                return f"{total_size:.1f} {unit}"
            total_size /= 1024.0
        return f"{total_size:.1f} TB"
    except Exception:
        return "Unknown"


def copy_file_safe(src, dest, file_description=""):
    """Safely copy a file with error handling."""
    try:
        if src.exists():
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            return True
        else:
            print(f"  ⚠ Missing: {file_description or src.name}")
            return False
    except Exception as e:
        print(f"  ✗ Error copying {src.name}: {e}")
        return False


def copy_directory_safe(src, dest, dir_description=""):
    """Safely copy a directory with error handling."""
    try:
        if src.exists() and src.is_dir():
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, dest, dirs_exist_ok=True)
            return True
        else:
            print(f"  ⚠ Missing: {dir_description or src.name}")
            return False
    except Exception as e:
        print(f"  ✗ Error copying directory {src.name}: {e}")
        return False


def get_best_cutoff(model_perf_metrics_file):
    """
    Read model performance metrics and return the best cutoff.
    
    Args:
        model_perf_metrics_file (Path): Path to model_performance_metrics.csv
        
    Returns:
        str: Best cutoff identifier (e.g., "cutoff_5") or None if not found
    """
    try:
        if not model_perf_metrics_file.exists():
            return None
            
        df = pd.read_csv(model_perf_metrics_file)
        
        # Try to determine the best cutoff based on available metrics
        if 'MCC' in df.columns:
            # For classification, use MCC (higher is better)
            best_row = df.loc[df['MCC'].idxmax()]
            metric_name = 'MCC'
            metric_value = best_row['MCC']
        elif 'r2' in df.columns:
            # For regression, use R2 (higher is better)
            best_row = df.loc[df['r2'].idxmax()]
            metric_name = 'r2'
            metric_value = best_row['r2']
        else:
            print(f"  ⚠ Warning: No recognized performance metric found in {model_perf_metrics_file}")
            return None
        
        best_cutoff = best_row['cut_off']
        print(f"  ✓ Best cutoff: {best_cutoff} (best {metric_name}: {metric_value:.4f})")
        return best_cutoff
        
    except Exception as e:
        print(f"  ✗ Error reading performance metrics: {e}")
        return None


def process_kmer_iteration(iteration_dir, dest_iteration, minimal=False):
    """
    Process k-mer workflow files for a single iteration.
    
    Args:
        iteration_dir (Path): Source iteration directory
        dest_iteration (Path): Destination iteration directory
        minimal (bool): If True, skip large files and only copy best cutoff
        
    Returns:
        int: Number of files/directories copied
    """
    print(f"  Processing k-mer workflow...")
    files_copied = 0
    
    # Define source and destination paths for k-mer modeling
    model_perf_src = iteration_dir / "modeling" / "modeling_results" / "model_performance"
    feature_tables_src = iteration_dir / "modeling" / "feature_selection" / "filtered_feature_tables"
    modeling_results_src = iteration_dir / "modeling" / "modeling_results"
    predict_src = iteration_dir / "model_validation" / "predict_results"
    
    model_perf_dest = dest_iteration / "modeling" / "modeling_results" / "model_performance"
    feature_tables_dest = dest_iteration / "modeling" / "feature_selection" / "filtered_feature_tables"
    modeling_results_dest = dest_iteration / "modeling" / "modeling_results"
    predict_dest = dest_iteration / "model_validation" / "predict_results"
    
    # Copy model performance files
    performance_files = [
        "model_performance_metrics.csv",
        "pr_curve.png", 
        "roc_curve.png"
    ]
    
    for filename in performance_files:
        src_file = model_perf_src / filename
        dest_file = model_perf_dest / filename
        if copy_file_safe(src_file, dest_file, f"kmer/{filename}"):
            files_copied += 1
    
    # Copy filtered_feature_tables (best cutoff only)
    metrics_file = model_perf_src / "model_performance_metrics.csv"
    best_cutoff = get_best_cutoff(metrics_file)
    
    if best_cutoff:
        # Handle different cutoff formats
        if best_cutoff.startswith('cutoff_'):
            cutoff_num = best_cutoff.split('_')[-1]
            best_feature_file = f"select_feature_table_cutoff_{cutoff_num}.csv"
        else:
            best_feature_file = f"select_feature_table_{best_cutoff}.csv"
            
        src_file = feature_tables_src / best_feature_file
        dest_file = feature_tables_dest / best_feature_file
        
        if copy_file_safe(src_file, dest_file, f"kmer/best feature table ({best_feature_file})"):
            files_copied += 1
    else:
        print("    ⚠ Warning: Could not determine best cutoff, skipping feature tables")
    
    # Copy top-level modeling results files
    top_level_files = [
        "select_features_model_performance.csv",
        "select_features_model_predictions.csv"
    ]
    
    for filename in top_level_files:
        src_file = modeling_results_src / filename
        dest_file = modeling_results_dest / filename
        if copy_file_safe(src_file, dest_file, f"kmer/{filename}"):
            files_copied += 1
    
    # Copy best cutoff directory with importance files (only in standard mode)
    if not minimal and best_cutoff:
        if best_cutoff.startswith('cutoff_'):
            cutoff_num = best_cutoff.split('_')[-1]
            best_cutoff_src = modeling_results_src / f"cutoff_{cutoff_num}"
        else:
            best_cutoff_src = modeling_results_src / str(best_cutoff)
            
        best_cutoff_dest = modeling_results_dest / best_cutoff_src.name
        
        if best_cutoff_src.exists():
            # Find all run directories in the best cutoff
            run_dirs = sorted([d for d in best_cutoff_src.iterdir() 
                              if d.is_dir() and d.name.startswith('run_')])
            
            cutoff_files_copied = 0
            for run_dir in run_dirs:
                run_name = run_dir.name
                dest_run_dir = best_cutoff_dest / run_name
                
                # Copy importance files
                importance_files = ["shap_importances.csv", "feature_importances.csv"]
                
                for filename in importance_files:
                    src_file = run_dir / filename
                    dest_file = dest_run_dir / filename
                    if copy_file_safe(src_file, dest_file, f"kmer/{best_cutoff}/{run_name}/{filename}"):
                        cutoff_files_copied += 1
            
            if cutoff_files_copied > 0:
                files_copied += cutoff_files_copied
                print(f"    ✓ Copied {cutoff_files_copied} importance files from {best_cutoff}")
    
    # Copy k-mer feature_tables directory
    if minimal:
        # Only copy essential files from feature_tables
        essential_files = [
            "final_feature_table.csv",
            "selected_features.csv",
            "phage_final_feature_table.csv"
        ]
        
        feature_tables_src_dir = iteration_dir / "feature_tables"
        feature_tables_dest_dir = dest_iteration / "feature_tables"
        
        for filename in essential_files:
            src_file = feature_tables_src_dir / filename
            dest_file = feature_tables_dest_dir / filename
            if copy_file_safe(src_file, dest_file, f"kmer/feature_tables/{filename}"):
                files_copied += 1
    else:
        # Copy entire feature_tables directory
        feature_tables_src_dir = iteration_dir / "feature_tables"
        feature_tables_dest_dir = dest_iteration / "feature_tables"
        if copy_directory_safe(feature_tables_src_dir, feature_tables_dest_dir, "kmer/feature_tables directory"):
            files_copied += 1
    
    # Copy full_feature_table.csv (at iteration root level)
    full_feature_src = iteration_dir / "full_feature_table.csv"
    full_feature_dest = dest_iteration / "full_feature_table.csv"
    if copy_file_safe(full_feature_src, full_feature_dest, "kmer/full_feature_table.csv"):
        files_copied += 1
    
    # Copy model validation results
    validation_files = [
        "strain_median_predictions.csv",
        "strain_all_predictions.csv"
    ]
    
    for filename in validation_files:
        src_file = predict_src / filename
        dest_file = predict_dest / filename
        if copy_file_safe(src_file, dest_file, f"kmer/{filename}"):
            files_copied += 1
    
    # Copy top-level iteration files
    top_level_files = [
        "modeling_strains.csv", 
        "validation_strains.csv",
        "workflow_report.txt"
    ]
    
    for filename in top_level_files:
        src_file = iteration_dir / filename
        dest_file = dest_iteration / filename
        if copy_file_safe(src_file, dest_file, f"kmer/{filename}"):
            files_copied += 1
    
    return files_copied


def export_kmer_bootstrap_results(source_dir, dest_dir, minimal=False):
    """
    Export k-mer bootstrap validation results.
    
    Args:
        source_dir (Path): Source directory containing bootstrap_validation
        dest_dir (Path): Destination directory for export
        minimal (bool): If True, skip large files and only copy best cutoff
    """
    source_dir = Path(source_dir)
    dest_dir = Path(dest_dir)
    
    # Validate source directory exists
    if not source_dir.exists():
        print(f"✗ Error: Source directory '{source_dir}' does not exist")
        sys.exit(1)
    
    # Create destination directory
    dest_dir.mkdir(parents=True, exist_ok=True)
    
    mode_text = "MINIMAL" if minimal else "STANDARD"
    print(f"Exporting k-mer bootstrap validation results ({mode_text} mode)...")
    print(f"Source: {source_dir}")
    print(f"Destination: {dest_dir}")
    if minimal:
        print("📦 Minimal mode: Only copying essential k-mer feature files, skipping importance files")
    else:
        print("📦 Standard mode: Copying all k-mer results including best cutoffs with importance files")
    print("")
    
    # Copy top-level final predictions
    final_predictions_files = [
        "final_kmer_predictions.csv",
        "kmer_prediction_summary.csv",
        "iteration_performance_summary.csv"
    ]
    
    for filename in final_predictions_files:
        src_file = source_dir / filename
        dest_file = dest_dir / filename
        if copy_file_safe(src_file, dest_file, filename):
            print(f"✓ Copied {filename}")
    
    # Initialize counters
    total_iterations = 0
    successful_iterations = 0
    total_files_copied = 0
    
    # Find all iteration directories
    iteration_dirs = sorted([d for d in source_dir.iterdir() 
                           if d.is_dir() and d.name.startswith('iteration_')])
    
    if not iteration_dirs:
        print("⚠ Warning: No iteration directories found")
        return
    
    # Process each iteration
    for iteration_dir in iteration_dirs:
        iteration_name = iteration_dir.name
        total_iterations += 1
        
        print(f"Processing {iteration_name}...")
        
        # Create iteration directory structure in destination
        dest_iteration = dest_dir / iteration_name
        
        # Track files copied for this iteration
        iteration_files_copied = process_kmer_iteration(iteration_dir, dest_iteration, minimal)
        
        # Update counters
        if iteration_files_copied > 0:
            successful_iterations += 1
            total_files_copied += iteration_files_copied
            print(f"  ✓ Copied {iteration_files_copied} items from {iteration_name}")
        else:
            print(f"  ✗ No files found in {iteration_name}")
    
    # Print summary
    print("")
    print("Export complete!")
    print(f"Successfully processed: {successful_iterations}/{total_iterations} iterations")
    print(f"Total files/directories copied: {total_files_copied}")
    
    # Calculate and display total size
    total_size = get_directory_size(dest_dir)
    print(f"Total exported size: {total_size}")
    print(f"Results saved to: {dest_dir}")
    
    # Create manifest file
    create_manifest(dest_dir, source_dir, successful_iterations, total_iterations, minimal)


def create_manifest(dest_dir, source_dir, successful_iterations, total_iterations, minimal=False):
    """Create a manifest file documenting the export."""
    manifest_file = dest_dir / "export_manifest.txt"
    
    mode_text = "MINIMAL" if minimal else "STANDARD"
    
    manifest_content = f"""K-mer Bootstrap Validation Export Manifest ({mode_text} Mode)
=======================================================
Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Source Directory: {source_dir}
Destination Directory: {dest_dir}
Export Mode: {mode_text}
Total Iterations Processed: {successful_iterations}/{total_iterations}

Files Exported per Iteration:

K-MER MODELING WORKFLOW:
- modeling/modeling_results/model_performance/model_performance_metrics.csv
- modeling/modeling_results/model_performance/pr_curve.png
- modeling/modeling_results/model_performance/roc_curve.png
- modeling/feature_selection/filtered_feature_tables/select_feature_table_[BEST_CUTOFF].csv (best cutoff only)
- modeling/modeling_results/select_features_model_performance.csv
- modeling/modeling_results/select_features_model_predictions.csv"""

    if not minimal:
        manifest_content += """
- modeling/modeling_results/[BEST_CUTOFF]/run_*/[shap_importances.csv, feature_importances.csv]
- feature_tables/ (complete directory with all k-mer feature assignments)"""
    else:
        manifest_content += """
- feature_tables/final_feature_table.csv
- feature_tables/selected_features.csv  
- feature_tables/phage_final_feature_table.csv
  (Note: Complete feature_tables directory excluded in minimal mode)"""

    manifest_content += """
- model_validation/predict_results/strain_median_predictions.csv
- model_validation/predict_results/strain_all_predictions.csv
- full_feature_table.csv
- modeling_strains.csv
- validation_strains.csv
- workflow_report.txt

Top-level Files:
- final_kmer_predictions.csv
- kmer_prediction_summary.csv  
- iteration_performance_summary.csv

Usage:
This export contains the key results and performance metrics from
each k-mer bootstrap validation iteration.

STANDARD MODE (default):
- Includes all k-mer feature tables and assignments
- Copies best cutoff feature tables only (not all cutoffs)
- Includes importance files from best cutoffs
- Complete k-mer modeling results

MINIMAL MODE (--minimal flag):
- From feature_tables, only copies essential files (final_feature_table.csv, selected_features.csv, phage_final_feature_table.csv)
- Excludes importance files from cutoff directories
- Focuses on core k-mer modeling performance and predictions only
- Significantly reduces storage requirements
"""
    
    if minimal:
        manifest_content += """
MINIMAL MODE NOTES:
- feature_tables/ directory reduced to only essential CSV files
- Importance files from cutoff directories completely excluded
- This provides the core k-mer modeling results needed for analysis with minimal storage
"""
    
    try:
        with open(manifest_file, 'w') as f:
            f.write(manifest_content)
        print(f"Manifest created: {manifest_file}")
    except Exception as e:
        print(f"⚠ Warning: Could not create manifest file: {e}")


def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Export k-mer bootstrap validation results from HPC cluster",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python export_kmer_bootstrap_results.py -i ./kmer_bootstrap_validation -o ./export
  python export_kmer_bootstrap_results.py --input /scratch/user/kmer_results --output ~/results --minimal

  Standard mode (default): Copies all k-mer results, includes importance files
  Minimal mode (--minimal): Only copies essential k-mer files, skips importance files
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Input directory containing k-mer bootstrap validation results'
    )
    
    parser.add_argument(
        '-o', '--output', 
        type=str,
        required=True,
        help='Output directory for exported results'
    )
    
    parser.add_argument(
        '--minimal',
        action='store_true',
        help='Minimal export: skip importance files, only copy essential k-mer feature files'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be copied without actually copying files'
    )
    
    args = parser.parse_args()
    
    if args.dry_run:
        mode_text = "MINIMAL" if args.minimal else "STANDARD"
        print(f"DRY RUN MODE ({mode_text}) - No files will be copied")
        print("")
    
    # Convert to Path objects and export
    source_path = Path(args.input).resolve()
    dest_path = Path(args.output).resolve()
    
    if not args.dry_run:
        export_kmer_bootstrap_results(source_path, dest_path, minimal=args.minimal)
    else:
        print(f"Would export from: {source_path}")
        print(f"Would export to: {dest_path}")
        print(f"Minimal mode: {args.minimal}")
        print("Run without --dry-run to perform actual export")


if __name__ == "__main__":
    main()