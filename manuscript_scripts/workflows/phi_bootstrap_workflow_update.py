import os
import pandas as pd
import shutil
import random
import argparse
from phage_modeling.workflows.protein_family_workflow import run_protein_family_workflow
from phage_modeling.workflows.assign_predict_workflow import assign_predict_workflow

def get_full_strain_list(interaction_matrix, input_strain_dir, strain_column):
    # Read interaction matrix
    interaction_df = pd.read_csv(interaction_matrix)
    strains_in_matrix = interaction_df[strain_column].unique()

    # Get strains from input_strain_dir
    strains_in_dir = [f.split('.')[0] for f in os.listdir(input_strain_dir) if f.endswith('.faa')]

    # Get intersection of strains in interaction matrix and input directory
    full_strain_list = list(set(strains_in_matrix).intersection(set(strains_in_dir)))
    print(f"Found {len(full_strain_list)} strains in the interaction matrix and input directory.")

    return full_strain_list

def split_strains(full_strain_list, iteration, validation_percentage=0.1):
    """
    Splits the full strain list into modeling and validation sets.

    Args:
        full_strain_list (list): List of all available strains.
        iteration (int): The iteration number, used as the seed for reproducibility.
        validation_percentage (float): Percentage of strains to use for validation (default: 10%).

    Returns:
        (list, list): Tuple of (modeling_strains, validation_strains)
    """
    # Set the random seed based on the iteration number
    random.seed(iteration)

    # Randomly shuffle strains and split into 90% modeling and 10% validation
    random.shuffle(full_strain_list)
    split_index = int(len(full_strain_list) * (1 - validation_percentage))

    modeling_strains = full_strain_list[:split_index]
    validation_strains = full_strain_list[split_index:]

    return modeling_strains, validation_strains

def select_best_cutoff(output_dir):
    # Load model performance metrics
    metrics_file = os.path.join(output_dir, "modeling_results/model_performance/model_performance_metrics.csv")
    metrics_df = pd.read_csv(metrics_file)

    metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
    best_cutoff = metrics_df['cut_off'].values[0]

    return best_cutoff

def run_model_validation(
    full_strain_list, 
    input_strain_dir, 
    input_phage_dir, 
    clustering_dir,
    interaction_matrix, 
    output_dir, 
    n_iterations, 
    threads, 
    num_runs_fs, 
    num_runs_modeling, 
    use_dynamic_weights,
    weights_method='log10',
    use_clustering=True,
    cluster_method='hdbscan',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    duplicate_all=False
):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    final_predictions = pd.DataFrame()

    # Check existing iterations
    completed_iterations = []
    for i in range(1, n_iterations + 1):
        iteration_output_dir = os.path.join(output_dir, f'iteration_{i}')
        median_predictions_file = os.path.join(iteration_output_dir, 'model_validation', 'predict_results', 'strain_median_predictions.csv')
        
        # Skip the iteration if median predictions file already exists
        if os.path.exists(median_predictions_file):
            completed_iterations.append(i)
            print(f"Iteration {i} is already complete, skipping.")
            # Load median predictions and append to final_predictions
            median_predictions = pd.read_csv(median_predictions_file)
            median_predictions['iteration'] = i
            final_predictions = pd.concat([final_predictions, median_predictions], ignore_index=True)
            continue

        # Delete the partial directory if it exists to start fresh
        if os.path.exists(iteration_output_dir):
            # Check if modeling was completed
            metrics_file = os.path.join(iteration_output_dir, "modeling_results/model_performance/model_performance_metrics.csv")
            if os.path.exists(metrics_file):
                print(f"Partial directory found for iteration {i}, modeling completed. Will resume with prediction.")
            else:
                print(f"Partial directory found for iteration {i}, modeling incomplete. Resuming modeling...")
                # print(f"Partial directory found for iteration {i}, modeling incomplete. Deleting...")
                # shutil.rmtree(iteration_output_dir)
        
        # Proceed with the iteration
        print(f"Starting iteration {i}...")
        os.makedirs(iteration_output_dir, exist_ok=True)
        modeling_tmp_dir = os.path.join(iteration_output_dir, 'tmp')

        # Split strains into 90% for modeling and 10% for validation
        if clustering_dir:
            modeling_strains_old = os.path.join(clustering_dir, f'iteration_{i}', 'modeling_strains.csv')
            validation_strains_old = os.path.join(clustering_dir, f'iteration_{i}', 'validation_strains.csv')

            modeling_strains_new = os.path.join(iteration_output_dir, 'modeling_strains.csv')
            validation_strains_new = os.path.join(iteration_output_dir, 'validation_strains.csv')

            # If link doesn't exist, create new links
            if not os.path.exists(modeling_strains_new):
                os.symlink(modeling_strains_old, modeling_strains_new)
            if not os.path.exists(validation_strains_new):
                os.symlink(validation_strains_old, validation_strains_new)
        else:
            modeling_strains_path = os.path.join(iteration_output_dir, 'modeling_strains.csv')
            validation_strains_path = os.path.join(iteration_output_dir, 'validation_strains.csv')

            if not os.path.exists(modeling_strains_path) or not os.path.exists(validation_strains_path):
                modeling_strains, validation_strains = split_strains(full_strain_list, iteration=i)

                # Save strain lists for this iteration
                pd.DataFrame(modeling_strains, columns=['strain']).to_csv(os.path.join(iteration_output_dir, 'modeling_strains.csv'), index=False)
                pd.DataFrame(validation_strains, columns=['strain']).to_csv(os.path.join(iteration_output_dir, 'validation_strains.csv'), index=False)
            
            else:
                print("Strain lists already exist. Loading...")
                modeling_strains = pd.read_csv(modeling_strains_path)['strain'].tolist()
                validation_strains = pd.read_csv(validation_strains_path)['strain'].tolist()

        # Step 3: Run protein_family_workflow with modeling strains
        metrics_file = os.path.join(iteration_output_dir, "modeling_results/model_performance/model_performance_metrics.csv")
        if clustering_dir:
            iteration_clustering_dir = os.path.join(clustering_dir, f'iteration_{i}')
        else:
            iteration_clustering_dir = None

        if not os.path.exists(metrics_file):
            run_protein_family_workflow(
                input_path_strain=input_strain_dir,
                input_path_phage=input_phage_dir,
                clustering_dir=iteration_clustering_dir,
                phenotype_matrix=interaction_matrix,
                output_dir=iteration_output_dir,
                tmp_dir=modeling_tmp_dir,
                min_seq_id=0.4,
                coverage=0.8,
                sensitivity=7.5,
                threads=threads,
                phenotype_column='interaction',
                phage_column='phage',
                strain_list=os.path.join(iteration_output_dir, 'modeling_strains.csv'),
                phage_list=interaction_matrix,
                filter_type='strain',
                num_runs_fs=num_runs_fs,
                num_runs_modeling=num_runs_modeling,
                use_dynamic_weights=use_dynamic_weights,
                weights_method=weights_method,
                use_clustering=use_clustering,
                cluster_method=cluster_method,
                n_clusters=n_clusters,
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                cluster_selection_epsilon=cluster_selection_epsilon,
                check_feature_presence=check_feature_presence,
                max_ram=40
            )
        else:
            print("Modeling results already exist. Skipping protein_family_workflow.")

        # Step 4: Get the cutoff with the highest MCC
        best_cutoff = select_best_cutoff(iteration_output_dir)
        model_dir = os.path.join(iteration_output_dir, f'modeling_results', f'{best_cutoff}')

        # Step 5: Predict interactions for validation strains
        validation_output_dir = os.path.join(iteration_output_dir, 'model_validation')
        os.makedirs(validation_output_dir, exist_ok=True)
        validation_tmp_dir = os.path.join(validation_output_dir, 'tmp')

        # check for modified AA directory
        modified_aa_dir = os.path.join(iteration_output_dir, 'strain', 'modified_AAs', 'strain')
        if os.path.exists(modified_aa_dir):
            input_strain_dir = modified_aa_dir
            print("Using modified AA directory for prediction.")

        select_feature_table = os.path.join(iteration_output_dir, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_{best_cutoff}.csv')

        # if clustering_dir:
        #     mmseqs_db = os.path.join(clustering_dir, 'tmp', 'strain', 'mmseqs_db')
        #     cluster_tsv = os.path.join(clustering_dir, 'strain', 'best_hits.tsv')
        #     feature_map = os.path.join(clustering_dir, 'strain', 'features', 'selected_features.csv')
        #     phage_feature_table_path = os.path.join(clustering_dir, 'phage', 'features', 'feature_table.csv')
        # else:
        #     mmseqs_db = os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db')
        #     cluster_tsv = os.path.join(iteration_output_dir, 'strain', 'best_hits.tsv')
        #     feature_map = os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv')

        assign_predict_workflow(
            input_dir=input_strain_dir,
            genome_list=os.path.join(iteration_output_dir, 'validation_strains.csv'),
            mmseqs_db=os.path.join(iteration_output_dir, 'tmp', 'strain', 'mmseqs_db'),
            clusters_tsv=os.path.join(iteration_output_dir, 'strain', 'clusters.tsv'),
            feature_map=os.path.join(iteration_output_dir, 'strain', 'features', 'selected_features.csv'),
            tmp_dir=validation_tmp_dir,
            suffix='faa',
            model_dir=model_dir,
            feature_table=select_feature_table,
            phage_feature_table_path=os.path.join(iteration_output_dir, 'phage', 'features', 'feature_table.csv'),
            output_dir=validation_output_dir,
            threads=threads,
            genome_type='strain',
            sensitivity=7.5,
            coverage=0.8,
            min_seq_id=0.4,
            duplicate_all=duplicate_all
        )

        # Step 6: Concatenate median predictions for this iteration
        median_predictions = pd.read_csv(os.path.join(validation_output_dir, 'predict_results', 'strain_median_predictions.csv'))
        median_predictions['iteration'] = i
        final_predictions = pd.concat([final_predictions, median_predictions], ignore_index=True)

    # Save final concatenated predictions
    final_predictions.to_csv(os.path.join(output_dir, 'final_predictions.csv'), index=False)
    print("Completed all iterations. Final predictions saved.")

# Main function to handle arguments
def main():
    parser = argparse.ArgumentParser(description="Bootstrap model construction and validation.")
    parser.add_argument('--input_strain_dir', type=str, required=True, help="Directory containing strain FASTA files.")
    parser.add_argument('--input_phage_dir', type=str, required=True, help="Directory containing phage FASTA files.")
    parser.add_argument('--interaction_matrix', type=str, required=True, help="Path to the interaction matrix.")
    parser.add_argument('--clustering_dir', type=str, help="Directory containing strain clustering results.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--n_iterations', type=int, default=10, help="Number of iterations.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads to use.")
    parser.add_argument('--strain_column', type=str, default='strain', help="Column in the interaction matrix containing strain names.")
    parser.add_argument('--num_runs_fs', type=int, default=25, help="Number of runs for feature selection.")
    parser.add_argument('--num_runs_modeling', type=int, default=50, help="Number of runs for modeling.")
    parser.add_argument('--use_dynamic_weights', action='store_true', help="Use dynamic weights for feature selection.")
    parser.add_argument('--weights_method', type=str, default='log10', choices=['log10', 'inverse_frequency', 'balanced'], help="Method for calculating dynamic weights.")
    parser.add_argument('--use_clustering', action='store_true', help="Use clustering results for feature selection.")
    parser.add_argument('--cluster_method', type=str, default='hdbscan', choices=['hdbscan', 'hierarchical'], help="Clustering method to use.")
    parser.add_argument('--n_clusters', type=int, default=20, help="Number of clusters for clustering feature selection.")
    parser.add_argument('--min_cluster_size', type=int, default=2, help="Minimum cluster size for clustering feature selection.")
    parser.add_argument('--min_samples', type=int, help="Minimum number of samples for clustering feature selection.")
    parser.add_argument('--cluster_selection_epsilon', type=float, default=0.0, help="Epsilon value for clustering feature selection.")
    parser.add_argument('--check_feature_presence', action='store_true', help="Check for feature presence for train-test splits.")
    parser.add_argument('--duplicate_all', action='store_true', help="Duplicate all genomes in the feature table for predictions.")

    args = parser.parse_args()

    # Step 1: Get the full list of strains to use
    full_strain_list = get_full_strain_list(args.interaction_matrix, args.input_strain_dir, args.strain_column)

    # Step 2: Run model construction and validation for each iteration
    run_model_validation(
        full_strain_list=full_strain_list,
        input_strain_dir=args.input_strain_dir,
        input_phage_dir=args.input_phage_dir,
        clustering_dir=args.clustering_dir,
        interaction_matrix=args.interaction_matrix,
        output_dir=args.output_dir,
        n_iterations=args.n_iterations,
        threads=args.threads,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        duplicate_all=args.duplicate_all
    )

if __name__ == "__main__":
    main()

# python /usr2/people/anoonan/BRaVE/machine_learning/phage_foundry/ml_training/workflows/phi_bootstrap_workflow_update.py --input_strain_dir ~/BRaVE/resources/genome_data/klebsiella/beatriz/kleb_genomes/genome_AAs/ --input_phage_dir ~/BRaVE/resources/genome_data/klebsiella/beatriz/kleb_phage/genome_AAs/ --clustering_dir ~/BRaVE/machine_learning/phage_modeling_manuscript/kleb_Bea/full_workflow_bootstrap --interaction_matrix ~/BRaVE/resources/genome_data/klebsiella/beatriz/interactions/kleb_beatriz_interaction_matrix_full.csv --output_dir kleb_Bea/full_workflow_bootstrap_weighted_log10 --n_iterations 20 --threads 4 --strain_column strain --num_runs_fs 25 --num_runs_modeling 50 --use_dynamic_weights --use_clustering --cluster_method hierarchical --n_clusters 20 --duplicate_all
# python /usr2/people/anoonan/BRaVE/machine_learning/phage_foundry/ml_training/workflows/phi_bootstrap_workflow_update.py --input_strain_dir ~/BRaVE/resources/genome_data/klebsiella/demi/kleb_genomes/genome_AAs/ --input_phage_dir ~/BRaVE/resources/genome_data/klebsiella/demi/kleb_phage/genome_AAs/ --interaction_matrix ~/BRaVE/resources/genome_data/klebsiella/demi/interactions/phage_host_interactions_set1_long.csv --clustering_dir ~/BRaVE/machine_learning/phage_modeling_manuscript/kleb_Dimi1/full_workflow_bootstrap --output_dir kleb_Dimi1/full_workflow_bootstrap_weighted_log10 --n_iterations 20 --threads 4 --strain_column strain --num_runs_fs 25 --num_runs_modeling 50 --duplicate_all --use_dynamic_weights
