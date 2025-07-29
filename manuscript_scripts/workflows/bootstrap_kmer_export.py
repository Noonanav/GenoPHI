#!/usr/bin/env python3
"""
Script to export bootstrap validation results from HPC cluster.
Extracts key results and performance metrics while avoiding intermediate files.
Handles both main workflow and kmer_modeling workflow from all iterations.
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


def process_main_workflow(iteration_dir, dest_iteration, minimal=False):
    """
    Process main workflow files for a single iteration.
    
    Args:
        iteration_dir (Path): Source iteration directory
        dest_iteration (Path): Destination iteration directory
        minimal (bool): If True, skip large files and only copy best cutoff
        
    Returns:
        int: Number of files/directories copied
    """
    print("  Processing main workflow...")
    files_copied = 0
    
    # Define source and destination paths
    model_perf_src = iteration_dir / "modeling_results" / "model_performance"
    predict_src = iteration_dir / "model_validation" / "predict_results"
    feature_tables_src = iteration_dir / "feature_selection" / "filtered_feature_tables"
    
    model_perf_dest = dest_iteration / "modeling_results" / "model_performance"
    predict_dest = dest_iteration / "model_validation" / "predict_results"
    feature_tables_dest = dest_iteration / "feature_selection" / "filtered_feature_tables"
    
    # Files to copy from model_performance
    performance_files = [
        "model_performance_metrics.csv",
        "pr_curve.png", 
        "roc_curve.png"
    ]
    
    # Copy model performance files
    for filename in performance_files:
        src_file = model_perf_src / filename
        dest_file = model_perf_dest / filename
        if copy_file_safe(src_file, dest_file, f"main/{filename}"):
            files_copied += 1
    
    # Copy predictive_proteins directory (only if not minimal)
    if not minimal:
        predictive_proteins_src = model_perf_src / "predictive_proteins"
        predictive_proteins_dest = model_perf_dest / "predictive_proteins"
        if copy_directory_safe(predictive_proteins_src, predictive_proteins_dest, "main/predictive_proteins directory"):
            files_copied += 1
    else:
        print("    📦 Skipping predictive_proteins directory (minimal mode)")
    
    # Copy filtered_feature_tables (best cutoff only - for both standard and ultra-minimal)
    metrics_file = model_perf_src / "model_performance_metrics.csv"
    best_cutoff = get_best_cutoff(metrics_file)
    
    if best_cutoff:
        # Extract cutoff number from "cutoff_X" format
        if best_cutoff.startswith('cutoff_'):
            cutoff_num = best_cutoff.split('_')[-1]
            best_feature_file = f"select_feature_table_cutoff_{cutoff_num}.csv"
            
            src_file = feature_tables_src / best_feature_file
            dest_file = feature_tables_dest / best_feature_file
            
            if copy_file_safe(src_file, dest_file, f"main/best feature table ({best_feature_file})"):
                files_copied += 1
        else:
            print(f"    ⚠ Warning: Unexpected cutoff format: {best_cutoff}")
    else:
        print("    ⚠ Warning: Could not determine best cutoff, skipping feature tables")
    
    # Copy strain median predictions
    strain_pred_src = predict_src / "strain_median_predictions.csv"
    strain_pred_dest = predict_dest / "strain_median_predictions.csv"
    if copy_file_safe(strain_pred_src, strain_pred_dest, "main/strain_median_predictions.csv"):
        files_copied += 1
    
    # Copy main workflow top-level modeling_results files
    main_modeling_files = [
        "select_features_model_performance.csv",
        "select_features_model_predictions.csv"
    ]
    
    modeling_results_src = iteration_dir / "modeling_results"
    modeling_results_dest = dest_iteration / "modeling_results"
    
    for filename in main_modeling_files:
        src_file = modeling_results_src / filename
        dest_file = modeling_results_dest / filename
        if copy_file_safe(src_file, dest_file, f"main/modeling_results/{filename}"):
            files_copied += 1
    
    # Copy the best cutoff directory with importance files (only in standard mode)
    if not minimal:
        metrics_file = modeling_results_src / "model_performance" / "model_performance_metrics.csv"
        best_cutoff = get_best_cutoff(metrics_file)
        
        if best_cutoff:
            # Extract cutoff number from "cutoff_X" format
            if best_cutoff.startswith('cutoff_'):
                cutoff_num = best_cutoff.split('_')[-1]
                best_cutoff_src = modeling_results_src / f"cutoff_{cutoff_num}"
                best_cutoff_dest = modeling_results_dest / f"cutoff_{cutoff_num}"
                
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
                            if copy_file_safe(src_file, dest_file, f"main/cutoff_{cutoff_num}/{run_name}/{filename}"):
                                cutoff_files_copied += 1
                    
                    if cutoff_files_copied > 0:
                        files_copied += cutoff_files_copied
                        print(f"    ✓ Copied {cutoff_files_copied} importance files from cutoff_{cutoff_num}")
    
    # Copy top-level iteration files
    top_level_files = [
        "combined_workflow_summary.txt",
        "modeling_strains.csv", 
        "validation_strains.csv",
        "workflow_report.txt",
        "workflow_report.csv",  # Both .txt and .csv versions
        "workflow_section_metrics.csv"
    ]
    
    for filename in top_level_files:
        src_file = iteration_dir / filename
        dest_file = dest_iteration / filename
        if copy_file_safe(src_file, dest_file, f"main/{filename}"):
            files_copied += 1
    
    # Copy phage and strain cluster and feature files (skip in minimal mode)
    if not minimal:
        for feature_type in ["phage", "strain"]:
            # Copy features directory
            features_src = iteration_dir / feature_type / "features"
            features_dest = dest_iteration / feature_type / "features"
            if copy_directory_safe(features_src, features_dest, f"main/{feature_type}/features directory"):
                files_copied += 1
            
            # Copy cluster and presence/absence files
            cluster_files = [
                "clusters.tsv",
                "presence_absence_matrix.csv",
                "assigned_clusters.tsv"
            ]
            
            for filename in cluster_files:
                src_file = iteration_dir / feature_type / filename
                dest_file = dest_iteration / feature_type / filename
                if copy_file_safe(src_file, dest_file, f"main/{feature_type}/{filename}"):
                    files_copied += 1
    else:
        print("    📦 Skipping phage/ and strain/ directories (minimal mode)")
    
    # Copy merged directory files
    merged_src = iteration_dir / "merged" / "full_feature_table.csv"
    merged_dest = dest_iteration / "merged" / "full_feature_table.csv"
    if copy_file_safe(merged_src, merged_dest, "main/merged/full_feature_table.csv"):
        files_copied += 1
    
    # Copy features_occurrence.csv from feature_selection
    features_occ_src = iteration_dir / "feature_selection" / "features_occurrence.csv"
    features_occ_dest = dest_iteration / "feature_selection" / "features_occurrence.csv"
    if copy_file_safe(features_occ_src, features_occ_dest, "main/feature_selection/features_occurrence.csv"):
        files_copied += 1
    
    return files_copied


def process_kmer_workflow(iteration_dir, dest_iteration, minimal=False):
    """
    Process kmer_modeling workflow files for a single iteration.
    
    Args:
        iteration_dir (Path): Source iteration directory
        dest_iteration (Path): Destination iteration directory
        minimal (bool): If True, skip large files and only copy best cutoff
        
    Returns:
        int: Number of files/directories copied
    """
    kmer_src = iteration_dir / "kmer_modeling"
    
    if not kmer_src.exists():
        return 0
    
    print("  Processing kmer_modeling workflow...")
    files_copied = 0
    
    # Define source and destination paths for kmer_modeling
    kmer_model_perf_src = kmer_src / "modeling" / "modeling_results" / "model_performance"
    kmer_feature_tables_src = kmer_src / "modeling" / "feature_selection" / "filtered_feature_tables"
    kmer_modeling_results_src = kmer_src / "modeling" / "modeling_results"
    
    kmer_dest = dest_iteration / "kmer_modeling"
    kmer_model_perf_dest = kmer_dest / "modeling" / "modeling_results" / "model_performance"
    kmer_feature_tables_dest = kmer_dest / "modeling" / "feature_selection" / "filtered_feature_tables"
    kmer_modeling_results_dest = kmer_dest / "modeling" / "modeling_results"
    
    # Copy kmer model performance files
    performance_files = [
        "model_performance_metrics.csv",
        "pr_curve.png", 
        "roc_curve.png"
    ]
    
    for filename in performance_files:
        src_file = kmer_model_perf_src / filename
        dest_file = kmer_model_perf_dest / filename
        if copy_file_safe(src_file, dest_file, f"kmer/{filename}"):
            files_copied += 1
    
    # Copy kmer filtered_feature_tables (best cutoff only - for both standard and ultra-minimal)
    metrics_file = kmer_model_perf_src / "model_performance_metrics.csv"
    best_cutoff = get_best_cutoff(metrics_file)
    
    if best_cutoff:
        # Extract cutoff number from "cutoff_X" format
        if best_cutoff.startswith('cutoff_'):
            cutoff_num = best_cutoff.split('_')[-1]
            best_feature_file = f"select_feature_table_cutoff_{cutoff_num}.csv"
            
            src_file = kmer_feature_tables_src / best_feature_file
            dest_file = kmer_feature_tables_dest / best_feature_file
            
            if copy_file_safe(src_file, dest_file, f"kmer/best feature table ({best_feature_file})"):
                files_copied += 1
        else:
            print(f"    ⚠ Warning: Unexpected kmer cutoff format: {best_cutoff}")
    else:
        print("    ⚠ Warning: Could not determine kmer best cutoff, skipping feature tables")
    
    # Copy kmer top-level modeling results files (keep in both standard and minimal modes)
    kmer_top_level_files = [
        "select_features_model_performance.csv",
        "select_features_model_predictions.csv"
    ]
    
    for filename in kmer_top_level_files:
        src_file = kmer_modeling_results_src / filename
        dest_file = kmer_modeling_results_dest / filename
        if copy_file_safe(src_file, dest_file, f"kmer/{filename}"):
            files_copied += 1
    
    # Copy kmer feature_tables directory (all files or just final_feature_table.csv)
    if minimal:
        # Only copy final_feature_table.csv from feature_tables
        kmer_final_feature_src = kmer_src / "feature_tables" / "final_feature_table.csv"
        kmer_final_feature_dest = kmer_dest / "feature_tables" / "final_feature_table.csv"
        if copy_file_safe(kmer_final_feature_src, kmer_final_feature_dest, "kmer/feature_tables/final_feature_table.csv"):
            files_copied += 1
    else:
        # Copy entire feature_tables directory
        kmer_feature_tables_src_dir = kmer_src / "feature_tables"
        kmer_feature_tables_dest_dir = kmer_dest / "feature_tables"
        if copy_directory_safe(kmer_feature_tables_src_dir, kmer_feature_tables_dest_dir, "kmer/feature_tables directory"):
            files_copied += 1
    
    # Copy kmer full_feature_table.csv
    kmer_full_feature_src = kmer_src / "full_feature_table.csv"
    kmer_full_feature_dest = kmer_dest / "full_feature_table.csv"
    if copy_file_safe(kmer_full_feature_src, kmer_full_feature_dest, "kmer/full_feature_table.csv"):
        files_copied += 1
    
    # Copy kmer workflow_report.txt
    kmer_report_src = kmer_src / "workflow_report.txt"
    kmer_report_dest = kmer_dest / "workflow_report.txt"
    if copy_file_safe(kmer_report_src, kmer_report_dest, "kmer/workflow_report.txt"):
        files_copied += 1
    
    # Copy kmer model validation results (keep in both standard and minimal modes)
    kmer_predict_src = kmer_src / "model_validation" / "predict_results"
    kmer_predict_dest = kmer_dest / "model_validation" / "predict_results"
    
    # Look for strain median predictions in kmer workflow
    kmer_strain_pred_src = kmer_predict_src / "strain_median_predictions.csv"
    kmer_strain_pred_dest = kmer_predict_dest / "strain_median_predictions.csv"
    if copy_file_safe(kmer_strain_pred_src, kmer_strain_pred_dest, "kmer/strain_median_predictions.csv"):
        files_copied += 1
    
    # Copy kmer best cutoff directory with importance files (only in standard mode)
    if not minimal:
        kmer_metrics_file = kmer_model_perf_src / "model_performance_metrics.csv"
        kmer_best_cutoff = get_best_cutoff(kmer_metrics_file)
        
        if kmer_best_cutoff:
            # Extract cutoff number from "cutoff_X" format
            if kmer_best_cutoff.startswith('cutoff_'):
                cutoff_num = kmer_best_cutoff.split('_')[-1]
                kmer_best_cutoff_src = kmer_modeling_results_src / f"cutoff_{cutoff_num}"
                kmer_best_cutoff_dest = kmer_modeling_results_dest / f"cutoff_{cutoff_num}"
                
                if kmer_best_cutoff_src.exists():
                    # Find all run directories in the best cutoff
                    run_dirs = sorted([d for d in kmer_best_cutoff_src.iterdir() 
                                      if d.is_dir() and d.name.startswith('run_')])
                    
                    kmer_cutoff_files_copied = 0
                    for run_dir in run_dirs:
                        run_name = run_dir.name
                        dest_run_dir = kmer_best_cutoff_dest / run_name
                        
                        # Copy importance files
                        importance_files = ["shap_importances.csv", "feature_importances.csv"]
                        
                        for filename in importance_files:
                            src_file = run_dir / filename
                            dest_file = dest_run_dir / filename
                            if copy_file_safe(src_file, dest_file, f"kmer/cutoff_{cutoff_num}/{run_name}/{filename}"):
                                kmer_cutoff_files_copied += 1
                    
                    if kmer_cutoff_files_copied > 0:
                        files_copied += kmer_cutoff_files_copied
                        print(f"    ✓ Copied {kmer_cutoff_files_copied} kmer importance files from cutoff_{cutoff_num}")
    
    return files_copied


def export_bootstrap_results(source_dir, dest_dir, minimal=False):
    """
    Export bootstrap validation results including both main and kmer_modeling workflows.
    
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
    
    mode_text = "ULTRA-MINIMAL" if minimal else "STANDARD"
    print(f"Exporting bootstrap validation results ({mode_text} mode)...")
    print(f"Source: {source_dir}")
    print(f"Destination: {dest_dir}")
    if minimal:
        print("📦 Ultra-minimal mode: Skipping phage/, strain/, kmer modeling_results, only copying kmer final_feature_table.csv")
    else:
        print("📦 Standard mode: Skipping only predictive_proteins, copying best cutoffs only")
    print("")
    
    # Copy top-level final predictions
    final_predictions = source_dir / "final_predictions.csv"
    if copy_file_safe(final_predictions, dest_dir / "final_predictions.csv", "final_predictions.csv"):
        print("✓ Copied final_predictions.csv")
    
    # Initialize counters
    total_iterations = 0
    successful_iterations = 0
    total_files_copied = 0
    kmer_iterations_found = 0
    
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
        iteration_files_copied = 0
        
        # Process main workflow
        main_files = process_main_workflow(iteration_dir, dest_iteration, minimal)
        iteration_files_copied += main_files
        
        # Process kmer_modeling workflow if it exists
        kmer_files = process_kmer_workflow(iteration_dir, dest_iteration, minimal)
        if kmer_files > 0:
            kmer_iterations_found += 1
        iteration_files_copied += kmer_files
        
        # Update counters
        if iteration_files_copied > 0:
            successful_iterations += 1
            total_files_copied += iteration_files_copied
            print(f"  ✓ Copied {iteration_files_copied} items from {iteration_name} (main: {main_files}, kmer: {kmer_files})")
        else:
            print(f"  ✗ No files found in {iteration_name}")
    
    # Print summary
    print("")
    print("Export complete!")
    print(f"Successfully processed: {successful_iterations}/{total_iterations} iterations")
    print(f"Iterations with kmer_modeling: {kmer_iterations_found}/{total_iterations}")
    print(f"Total files/directories copied: {total_files_copied}")
    
    # Calculate and display total size
    total_size = get_directory_size(dest_dir)
    print(f"Total exported size: {total_size}")
    print(f"Results saved to: {dest_dir}")
    
    # Create manifest file
    create_manifest(dest_dir, source_dir, successful_iterations, total_iterations, kmer_iterations_found, minimal)


def create_manifest(dest_dir, source_dir, successful_iterations, total_iterations, kmer_iterations_found, minimal=False):
    """Create a manifest file documenting the export."""
    manifest_file = dest_dir / "export_manifest.txt"
    
    mode_text = "ULTRA-MINIMAL" if minimal else "STANDARD"
    
    manifest_content = f"""Bootstrap Validation Export Manifest ({mode_text} Mode)
=====================================
Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Source Directory: {source_dir}
Destination Directory: {dest_dir}
Export Mode: {mode_text}
Total Iterations Processed: {successful_iterations}/{total_iterations}
Iterations with kmer_modeling: {kmer_iterations_found}/{total_iterations}

Files Exported per Iteration:

MAIN WORKFLOW:
- modeling_results/model_performance/model_performance_metrics.csv
- modeling_results/model_performance/pr_curve.png
- modeling_results/model_performance/roc_curve.png"""

    if not minimal:
        manifest_content += """
- modeling_results/model_performance/predictive_proteins/ (directory)"""
    else:
        manifest_content += """
  (Note: predictive_proteins directory excluded in ultra-minimal mode)"""

    manifest_content += """
- feature_selection/filtered_feature_tables/select_feature_table_cutoff_X.csv (best cutoff only)
- model_validation/predict_results/strain_median_predictions.csv
- modeling_results/select_features_model_performance.csv
- modeling_results/select_features_model_predictions.csv"""

    if not minimal:
        manifest_content += """
- modeling_results/cutoff_[BEST]/run_*/[shap_importances.csv, feature_importances.csv]"""

    manifest_content += """
- combined_workflow_summary.txt
- modeling_strains.csv
- validation_strains.csv
- workflow_report.txt
- workflow_report.csv
- workflow_section_metrics.csv"""

    if not minimal:
        manifest_content += """
- phage/features/ (directory)
- strain/features/ (directory)
- phage/clusters.tsv, presence_absence_matrix.csv, assigned_clusters.tsv
- strain/clusters.tsv, presence_absence_matrix.csv, assigned_clusters.tsv"""
    else:
        manifest_content += """
  (Note: phage/ and strain/ directories excluded in ultra-minimal mode)"""

    manifest_content += """
- merged/full_feature_table.csv
- feature_selection/features_occurrence.csv

KMER_MODELING WORKFLOW (if present):
- kmer_modeling/modeling/modeling_results/model_performance/model_performance_metrics.csv
- kmer_modeling/modeling/modeling_results/model_performance/pr_curve.png
- kmer_modeling/modeling/modeling_results/model_performance/roc_curve.png
- kmer_modeling/modeling/feature_selection/filtered_feature_tables/select_feature_table_cutoff_X.csv (best cutoff only)
- kmer_modeling/modeling/modeling_results/select_features_model_performance.csv
- kmer_modeling/modeling/modeling_results/select_features_model_predictions.csv"""

    if not minimal:
        manifest_content += """
- kmer_modeling/modeling/modeling_results/cutoff_[BEST]/run_*/[shap_importances.csv, feature_importances.csv]
- kmer_modeling/feature_tables/ (directory with feature assignments)"""
    else:
        manifest_content += """
- kmer_modeling/feature_tables/final_feature_table.csv (only this file)
  (Note: kmer cutoff directories excluded in ultra-minimal mode)"""

    manifest_content += """
- kmer_modeling/model_validation/predict_results/strain_median_predictions.csv (if present)
- kmer_modeling/full_feature_table.csv
- kmer_modeling/workflow_report.txt

Top-level Files:
- final_predictions.csv

Usage:
This export contains the key results and performance metrics from
each bootstrap validation iteration for both main and kmer_modeling workflows.

STANDARD MODE (default):
- Excludes only predictive_proteins directories to save space
- Copies best cutoff feature tables only (not all cutoffs)
- Includes all phage/strain cluster and feature files
- Includes complete kmer_modeling results
- Includes importance files from best cutoffs

ULTRA-MINIMAL MODE (--minimal flag):
- Excludes predictive_proteins directories  
- Excludes phage/ and strain/ directories entirely
- Excludes kmer_modeling modeling_results (but keeps model_validation)
- From kmer feature_tables, only copies final_feature_table.csv
- Still copies best cutoff feature tables only
- Focuses on core modeling performance and top-level summaries only
"""
    
    if minimal:
        manifest_content += """
ULTRA-MINIMAL MODE NOTES:
- phage/ and strain/ directories completely excluded to save space
- kmer_modeling/modeling/modeling_results/ completely excluded
- kmer_modeling/feature_tables/ reduced to only final_feature_table.csv
- kmer_modeling/model_validation/ is KEPT (contains important predictions)
- This provides the essential files needed for result review with minimal storage
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
        description="Export bootstrap validation results from HPC cluster (main + kmer_modeling)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python export_bootstrap_results_with_kmer.py -i ./bootstrap_validation -o ./export
  python export_bootstrap_results_with_kmer.py --input /scratch/user/bootstrap_validation --output ~/results --minimal

  Standard mode (default): Copies core results, skips predictive_proteins, includes phage/strain features
  Ultra-minimal mode (--minimal): Only copies essential modeling results, skips phage/strain/kmer_modeling_results
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Input directory containing bootstrap_validation results'
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
        help='Ultra-minimal export: skip phage/, strain/, kmer modeling_results (keeps model_validation), only copy kmer final_feature_table.csv'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be copied without actually copying files'
    )
    
    args = parser.parse_args()
    
    if args.dry_run:
        mode_text = "ULTRA-MINIMAL" if args.minimal else "STANDARD"
        print(f"DRY RUN MODE ({mode_text}) - No files will be copied")
        print("")
    
    # Convert to Path objects and export
    source_path = Path(args.input).resolve()
    dest_path = Path(args.output).resolve()
    
    if not args.dry_run:
        export_bootstrap_results(source_path, dest_path, minimal=args.minimal)
    else:
        print(f"Would export from: {source_path}")
        print(f"Would export to: {dest_path}")
        print(f"Minimal mode: {args.minimal}")
        print("Run without --dry-run to perform actual export")


if __name__ == "__main__":
    main()