#!/usr/bin/env python3
"""
Copy strain_median_predictions.csv files from CV directory structure.
Maintains the same directory structure in the output directory.
"""

import os
import shutil
import argparse
from pathlib import Path
import sys

def copy_predictions_files(input_dir, output_dir, target_filename="strain_median_predictions.csv"):
    """
    Copy strain_median_predictions.csv files from CV directory structure.
    
    Args:
        input_dir (str): Input CV directory containing iteration_* subdirectories
        output_dir (str): Output directory to copy files to
        target_filename (str): Target filename to copy (default: strain_median_predictions.csv)
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    if not input_path.exists():
        print(f"Error: Input directory does not exist: {input_dir}")
        return 1
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all iteration directories
    iteration_dirs = [d for d in input_path.iterdir() 
                     if d.is_dir() and d.name.startswith('iteration_')]
    
    if not iteration_dirs:
        print(f"Error: No iteration_* directories found in {input_dir}")
        return 1
    
    print(f"Found {len(iteration_dirs)} iteration directories")
    
    copied_count = 0
    missing_count = 0
    
    for iteration_dir in sorted(iteration_dirs):
        print(f"\nProcessing {iteration_dir.name}...")
        
        # Look for augmentation directory
        augmentation_dir = iteration_dir / "augmentation"
        if not augmentation_dir.exists():
            print(f"  Warning: No augmentation directory found in {iteration_dir.name}")
            continue
        
        # Find all frac_*pct directories
        frac_dirs = [d for d in augmentation_dir.iterdir() 
                    if d.is_dir() and d.name.startswith('frac_') and d.name.endswith('pct')]
        
        if not frac_dirs:
            print(f"  Warning: No frac_*pct directories found in {iteration_dir.name}/augmentation")
            continue
        
        print(f"  Found {len(frac_dirs)} fraction directories")
        
        for frac_dir in sorted(frac_dirs):
            # Construct path to target file
            target_file = frac_dir / "model_validation" / "predict_results" / target_filename
            
            if target_file.exists():
                # Construct output path maintaining directory structure
                relative_path = target_file.relative_to(input_path)
                output_file = output_path / relative_path
                
                # Create output directory structure
                output_file.parent.mkdir(parents=True, exist_ok=True)
                
                # Copy file
                shutil.copy2(target_file, output_file)
                print(f"    ✓ Copied {frac_dir.name}/{target_filename}")
                copied_count += 1
            else:
                print(f"    ✗ Missing {frac_dir.name}/{target_filename}")
                missing_count += 1
    
    print(f"\n" + "="*50)
    print(f"Summary:")
    print(f"  Files copied: {copied_count}")
    print(f"  Files missing: {missing_count}")
    print(f"  Output directory: {output_dir}")
    
    if copied_count > 0:
        print(f"\n✓ Successfully copied {copied_count} prediction files!")
        return 0
    else:
        print(f"\n✗ No files were copied. Check your directory structure.")
        return 1

def main():
    parser = argparse.ArgumentParser(
        description="Copy strain_median_predictions.csv files maintaining directory structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python copy_predictions.py /path/to/cv_dir /path/to/output_dir
  python copy_predictions.py --input /scratch/cv_results --output /scratch/predictions_only
  python copy_predictions.py -i ./ecoli_cv -o ./predictions --filename strain_all_predictions.csv

Directory structure expected:
  input_dir/
  ├── iteration_1/
  │   └── augmentation/
  │       ├── frac_00pct/
  │       │   └── model_validation/predict_results/strain_median_predictions.csv
  │       ├── frac_01pct/
  │       │   └── model_validation/predict_results/strain_median_predictions.csv
  │       └── ...
  ├── iteration_2/
  └── ...
        """)
    
    parser.add_argument('input_dir', nargs='?', 
                       help='Input CV directory containing iteration_* subdirectories')
    parser.add_argument('output_dir', nargs='?',
                       help='Output directory to copy files to')
    parser.add_argument('-i', '--input', dest='input_dir_alt', 
                       help='Input CV directory (alternative to positional arg)')
    parser.add_argument('-o', '--output', dest='output_dir_alt',
                       help='Output directory (alternative to positional arg)')
    parser.add_argument('--filename', default='strain_median_predictions.csv',
                       help='Target filename to copy (default: strain_median_predictions.csv)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be copied without actually copying')
    
    args = parser.parse_args()
    
    # Handle both positional and named arguments
    input_dir = args.input_dir or args.input_dir_alt
    output_dir = args.output_dir or args.output_dir_alt
    
    if not input_dir or not output_dir:
        parser.print_help()
        print("\nError: Both input and output directories are required")
        return 1
    
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Target filename: {args.filename}")
    if args.dry_run:
        print("DRY RUN MODE - No files will actually be copied")
    print()
    
    if args.dry_run:
        # TODO: Implement dry run mode
        print("Dry run mode not yet implemented")
        return 1
    
    return copy_predictions_files(input_dir, output_dir, args.filename)

if __name__ == "__main__":
    sys.exit(main())