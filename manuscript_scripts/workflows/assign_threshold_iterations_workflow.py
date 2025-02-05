import os
import glob
import subprocess
import pandas as pd
import argparse

def get_top_cutoff(method_dir):
    performance_file = os.path.join(method_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
    if os.path.exists(performance_file):
        performance_df = pd.read_csv(performance_file)
        performance_df = performance_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
        return performance_df['cut_off'].iloc[0]
    else:
        print(f"Warning: Performance file not found at {performance_file}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Run assign_threshold_iterations.py for each iteration directory.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing iteration_* subdirectories.")
    parser.add_argument('--fasta_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--select_protein_families', type=str, default=None, help="Path to a file with selected protein families.")
    
    args = parser.parse_args()

    # Find all iteration directories
    iteration_dirs = glob.glob(os.path.join(args.input_dir, 'iteration_*'))

    for iter_dir in sorted(iteration_dirs):
        if os.path.isdir(iter_dir):
            # Determine the cutoff
            cutoff = get_top_cutoff(iter_dir)
            if cutoff is None:
                continue  # Skip this iteration if we can't find the cutoff

            # Construct paths for other arguments
            feature_map = os.path.join(iter_dir, 'strain', 'features', 'selected_features.csv')
            best_hits_tsv = os.path.join(iter_dir, 'model_validation', 'tmp', 'best_hits.tsv')
            phage_feature_table_path = os.path.join(iter_dir, 'phage', 'features', 'feature_table.csv')
            output_dir = os.path.join(iter_dir, 'model_validation', 'threshold_iteration')
            model_dir = os.path.join(iter_dir, 'modeling_results', f'{cutoff}')

            # Check if all required files exist
            if not all(os.path.exists(path) for path in [feature_map, best_hits_tsv, phage_feature_table_path]):
                print(f"Skipping iteration {iter_dir} due to missing files.")
                continue

            final_strain_file = os.path.join(output_dir, 'threshold_1.00/strain_feature_table.csv')
            if os.path.exists(final_strain_file):
                print(f"Skipping iteration {iter_dir} because it has already been processed.")
                continue

            # Construct the command
            command = f"python /usr2/people/anoonan/BRaVE/machine_learning/phage_foundry/ml_training/workflows/assign_threshold_iterations.py --fasta_dir {args.fasta_dir} --feature_map {feature_map} --best_hits_tsv {best_hits_tsv} --phage_feature_table_path {phage_feature_table_path} --output_dir {output_dir} --model_dir {model_dir} --thresholds 0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 --select_protein_families {args.select_protein_families}"

            # Run the command
            try:
                subprocess.run(command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running command for {iter_dir}: {e}")

if __name__ == "__main__":
    main()