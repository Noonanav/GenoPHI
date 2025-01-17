import os
import pandas as pd
import logging
import subprocess
from argparse import ArgumentParser
from phage_modeling.mmseqs2_clustering import create_mmseqs_database, load_strains, create_contig_to_genome_dict, select_best_hits
from phage_modeling.workflows.prediction_workflow import generate_full_feature_table, predict_interactions, calculate_mean_predictions, load_model

def map_features(best_hits_tsv, feature_map, output_dir, genome_contig_mapping, genome, genome_type):
    """
    Maps the features for each new genome based on cluster assignments.

    Args:
        best_hits_tsv (str): Path to the best hits TSV file (output of select_best_hits).
        feature_map (str): Path to the feature mapping CSV file (selected_features.csv).
        output_dir (str): Directory to save the final feature table.
        genome_contig_mapping (dict): Dictionary mapping contigs to genomes.
        genome (str): Genome name of the current genome.
    """
    if not os.path.exists(best_hits_tsv):
        logging.error("Best hits TSV file does not exist: %s", best_hits_tsv)
        return

    try:
        # Load the best hits and feature mapping
        best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster'])
        feature_mapping = pd.read_csv(feature_map)

        # Ensure the 'Cluster' columns in both DataFrames are of the same type
        best_hits_df['Cluster'] = best_hits_df['Cluster'].astype(str)
        feature_mapping['Cluster_Label'] = feature_mapping['Cluster_Label'].astype(str)
        
    except Exception as e:
        logging.error("Error reading input files: %s", e)
        return

    logging.info(f"Mapping features for {genome_type} '{genome}'")

    # Convert the genome_contig_mapping dictionary to a DataFrame for merging
    genome_contig_mapping_df = pd.DataFrame(list(genome_contig_mapping.items()), columns=['contig_id', genome_type])

    # Merge best_hits_df with feature_mapping on the 'Cluster' column
    merged_df = best_hits_df.merge(feature_mapping, left_on='Cluster', right_on='Cluster_Label')

    # Merge with genome_contig_mapping_df to get genome information
    merged_df = merged_df.merge(genome_contig_mapping_df, left_on='Query', right_on='contig_id')

    # Create the binary feature presence table using the genome_type as the index
    feature_presence = merged_df.pivot_table(index=genome_type, columns='Feature', aggfunc='size', fill_value=0)
    feature_presence = (feature_presence > 0).astype(int)

    # Ensure all features are represented
    all_features = feature_mapping['Feature'].unique()
    for feature in all_features:
        if feature not in feature_presence.columns:
            feature_presence[feature] = 0

    # Reindex to ensure the presence of all features
    feature_presence = feature_presence.reindex(columns=all_features, fill_value=0).reset_index()

    return feature_presence

def assign_sequences_to_existing_clusters(query_db, target_db, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv):
    """
    Assigns sequences from the query database to clusters in the existing target database.

    Args:
        query_db (str): Path to the query MMseqs2 database (e.g., validation strains).
        target_db (str): Path to the target MMseqs2 database (e.g., strain feature database).
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        sensitivity (float): Sensitivity for assignment.
        threads (int): Number of threads to use.
        clusters_tsv (str): Path to the clusters TSV file.
    """
    logging.info("Assigning query sequences to existing target clusters...")

    result_db = os.path.join(tmp_dir, "result_db")
    # Run mmseqs2 search from query_db to target_db
    search_command = (
        f"mmseqs search {query_db} {target_db} {result_db} {tmp_dir} "
        f"-c {coverage} --min-seq-id {min_seq_id} -s {sensitivity} "
        f"--threads {threads} -v 3"
    )
    try:
        subprocess.run(search_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during MMseqs search: {e}")
        return None  # Early exit on error

    logging.info("Sequence assignment to existing clusters completed successfully.")

    assigned_tsv = os.path.join(output_dir, "assigned_clusters.tsv")
    # Create a TSV from query_db to target_db results
    createtsv_command = f"mmseqs createtsv {query_db} {target_db} {result_db} {assigned_tsv} --threads {threads} -v 3"
    
    try:
        subprocess.run(createtsv_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during MMseqs createtsv: {e}")
        return None  # Early exit on error
        
    logging.info(f"Assigned clusters TSV saved to {assigned_tsv}")

    best_hits_tsv = os.path.join(output_dir, "best_hits.tsv")
    select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)

    return best_hits_tsv

def process_new_genomes(input_dir, mmseqs_db, suffix, tmp_dir, output_dir, feature_map, clusters_tsv, genome_type, genomes=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4):
    """
    Processes genomes to assign them to existing clusters and generate feature tables.
    """
    os.makedirs(output_dir, exist_ok=True)
    feature_table_dir = os.path.join(output_dir, "feature_tables")
    os.makedirs(feature_table_dir, exist_ok=True)

    all_genomes_feature_tables = []

    if genomes is None:
        genomes = ['.'.join(f.split('.')[:-1]) for f in os.listdir(input_dir) if f.endswith(suffix)]

    for genome in genomes:
        genome_tmp_dir = os.path.join(tmp_dir, genome)
        os.makedirs(genome_tmp_dir, exist_ok=True)

        logging.info(f"Processing {genome_type} {genome}...")

        try:
            query_db = os.path.join(genome_tmp_dir, 'query_db')
            # Attempt to create MMseqs2 database for the current genome
            fasta_files = create_mmseqs_database(input_dir, query_db, suffix, 'directory', [genome], threads)
            
            if not fasta_files:
                logging.warning(f"No FASTA files found for {genome_type} '{genome}' with suffix '{suffix}'. Skipping...")
                continue  # Skip to the next genome if no FASTA files are found

            # Run the assignment and best-hit selection
            best_hits_tsv = assign_sequences_to_existing_clusters(query_db, mmseqs_db, genome_tmp_dir, genome_tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv)

            genome_contig_mapping, _ = create_contig_to_genome_dict(fasta_files, 'directory')

            # Get the feature presence DataFrame for the current genome
            feature_presence = map_features(best_hits_tsv, feature_map, output_dir, genome_contig_mapping, genome, genome_type)

            if feature_presence is not None:
                # Save the feature table for the current genome
                output_path = os.path.join(feature_table_dir, f'{genome}_feature_table.csv')
                feature_presence.to_csv(output_path, index=False)
                logging.info(f"Feature table for {genome_type} '{genome}' saved to {output_path}")
                
                # Append the feature table to the list
                all_genomes_feature_tables.append(feature_presence)

        except FileNotFoundError as e:
            logging.error(f"Error processing {genome_type} {genome}: {e}. Skipping...")
            continue  # Skip to the next genome if an error occurs

        except Exception as e:
            logging.error(f"Unexpected error while processing {genome_type} {genome}: {e}. Skipping...")
            continue  # Catch any other unexpected errors and continue with the next genome

    logging.info("Genome processing complete.")

    # Combine all genome feature tables into one DataFrame
    if all_genomes_feature_tables:
        combined_feature_table = pd.concat(all_genomes_feature_tables)
        combined_output_path = os.path.join(output_dir, 'combined_feature_table.csv')
        combined_feature_table.to_csv(combined_output_path, index=False)
        logging.info(f"Combined feature table saved to {combined_output_path}")

def run_assign_and_predict_workflow(input_dir, mmseqs_db, clusters_tsv, feature_map, tmp_dir, suffix, genome_list, genome_type, genome_column, model_dir, phage_feature_table, output_dir, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4):
    """
    Combines feature assignment with interaction prediction.
    """
    # Step 1: Feature assignment
    if genome_column is None:
        genome_column = genome_type

    if genome_list and os.path.exists(genome_list):
        genomes = load_strains(genome_list, genome_column)
    else:
        genomes = None  # Process all genomes in input_dir if no list is provided

    # Process genomes (phages or strains)
    process_new_genomes(
        input_dir=input_dir,
        genomes=genomes,
        mmseqs_db=mmseqs_db,
        suffix=suffix,
        tmp_dir=tmp_dir,
        output_dir=output_dir,
        feature_map=feature_map,
        clusters_tsv=clusters_tsv,
        genome_type=genome_type,  # Pass genome_type to the function
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        threads=threads
    )

    # Step 2: Prediction
    full_predictions_df = pd.DataFrame()
    phage_feature_table = pd.read_csv(phage_feature_table)

    strain_files = [f for f in os.listdir(output_dir) if f.endswith('_feature_table.csv')]
    for strain_file in strain_files:
        strain_table = pd.read_csv(os.path.join(output_dir, strain_file))
        prediction_feature_table = generate_full_feature_table(strain_table, phage_feature_table)

        all_predictions_df = predict_interactions(model_dir, prediction_feature_table)
        full_predictions_df = pd.concat([full_predictions_df, all_predictions_df])

    # Save all predictions
    full_predictions_df.to_csv(os.path.join(output_dir, 'all_predictions.csv'), index=False)
    mean_predictions_df = calculate_mean_predictions(full_predictions_df)
    mean_predictions_df.to_csv(os.path.join(output_dir, 'mean_predictions.csv'), index=False)

# Main function for CLI
def main():
    parser = ArgumentParser(description="Assign new genes to existing clusters and predict interactions.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing new genome FASTA files.")
    parser.add_argument('--genome_list', type=str, help="CSV file with genome names.")
    parser.add_argument('--genome_type', type=str, choices=['strain', 'phage'], default='strain', help="Type of genome to process.")
    parser.add_argument('--genome_column', type=str, help="Column name for genome identifiers in genome_list.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the existing MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map (selected_features.csv).")
    parser.add_argument('--suffix', type=str, default="faa", help="Suffix for FASTA files.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory with trained models.")
    parser.add_argument('--phage_feature_table', type=str, required=True, help="Path to the phage feature table.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.6, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")

    args = parser.parse_args()

    # Call the main workflow function
    run_assign_and_predict_workflow(
        input_dir=args.input_dir,
        genome_list=args.genome_list,
        genome_type=args.genome_type,
        genome_column=args.genome_column,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        feature_map=args.feature_map,
        tmp_dir=args.tmp_dir,
        suffix=args.suffix,
        model_dir=args.model_dir,
        phage_feature_table=args.phage_feature_table,
        output_dir=args.output_dir,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads
    )


if __name__ == "__main__":
    main()
