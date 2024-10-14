import os
import subprocess
import pandas as pd
import logging
from argparse import ArgumentParser
from Bio import SeqIO
from phage_modeling.mmseqs2_clustering import create_mmseqs_database, assign_sequences_to_clusters, load_strains, create_contig_to_genome_dict, select_best_hits

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

    logging.info("Mapping features for genome '%s'", genome)

    # Convert the genome_contig_mapping dictionary to a DataFrame for merging
    genome_contig_mapping_df = pd.DataFrame(list(genome_contig_mapping.items()), columns=['contig_id', 'genome'])

    # Merge best_hits_df with feature_mapping on the 'Cluster' column
    merged_df = best_hits_df.merge(feature_mapping, left_on='Cluster', right_on='Cluster_Label')

    # Merge with genome_contig_mapping_df to get genome information
    merged_df = merged_df.merge(genome_contig_mapping_df, left_on='Query', right_on='contig_id')

    # Create the binary feature presence table using the genome as the index
    feature_presence = merged_df.pivot_table(index='genome', columns='Feature', aggfunc='size', fill_value=0)
    feature_presence = (feature_presence > 0).astype(int)

    # Ensure all features are represented
    all_features = feature_mapping['Feature'].unique()
    for feature in all_features:
        if feature not in feature_presence.columns:
            feature_presence[feature] = 0

    # Reindex to ensure the presence of all features
    feature_presence = feature_presence.reindex(columns=all_features, fill_value=0).reset_index()

    return feature_presence

def process_new_genomes(input_dir, mmseqs_db, suffix, tmp_dir, output_dir, feature_map, clusters_tsv, genome_type, genomes=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4):
    """
    Processes genomes to assign them to existing clusters and generate feature tables with genome column named after the genome_type.
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

            result_db = os.path.join(genome_tmp_dir, 'result_db')
            assign_sequences_to_clusters(query_db, genome_tmp_dir, genome_tmp_dir, coverage, min_seq_id, sensitivity, threads, mmseqs_db)

            assigned_tsv = os.path.join(genome_tmp_dir, 'assigned_clusters.tsv')
            best_hits_tsv = os.path.join(genome_tmp_dir, 'best_hits.tsv')

            select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)

            genome_contig_mapping, _ = create_contig_to_genome_dict(fasta_files, 'directory')

            # Get the feature presence DataFrame for the current genome
            feature_presence = map_features(best_hits_tsv, feature_map, output_dir, genome_contig_mapping, genome, genome_type)

            if feature_presence is not None:
                # Rename the genome column to the genome_type
                feature_presence.rename(columns={'genome': genome_type}, inplace=True)
                
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

    # Combine all genome feature tables into one DataFrame
    if all_genomes_feature_tables:
        combined_feature_table = pd.concat(all_genomes_feature_tables)
        # Rename the genome column in the combined table to genome_type
        combined_feature_table.rename(columns={'genome': genome_type}, inplace=True)
        combined_output_path = os.path.join(output_dir, 'combined_feature_table.csv')
        combined_feature_table.to_csv(combined_output_path, index=False)
        logging.info(f"Combined feature table saved to {combined_output_path}")

    logging.info("Genome processing complete.")

def run_assign_features_workflow(input_dir, mmseqs_db, clusters_tsv, output_dir, feature_map, tmp_dir, suffix="faa", genome_list=None, genome_type='strain', genome_column=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4):
    """
    Wrapper function to run the full feature assignment workflow.

    Args:
        input_dir (str): Directory containing new genome FASTA files.
        mmseqs_db (str): Path to the existing MMseqs2 database.
        clusters_tsv (str): Path to the clusters TSV file.
        output_dir (str): Directory to save final feature tables.
        feature_map (str): Path to the feature map (selected_features.csv).
        tmp_dir (str): Temporary directory for intermediate files.
        suffix (str): Suffix for FASTA files.
        genome_list (str): Path to a CSV file with genome names.
        genome_type (str): Type of genomes being processed ('phage' or 'strain'). Default is 'strain'.
        genome_column (str): Column name for genome identifiers in the genome_list. Defaults to 'strain' or 'phage' depending on genome_type.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
    """
    # Determine default column based on genome type
    if genome_column is None:
        genome_column = genome_type

    # Load genomes from CSV if provided
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
        genome_type=genome_type,  # <-- Pass the genome_type argument here
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        threads=threads
    )

def main():
    # Set up argument parser
    parser = ArgumentParser(description="Assign new genes to existing clusters and generate feature tables.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing new genome FASTA files.")
    parser.add_argument('--genome_list', type=str, help="CSV file with genome names.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the existing MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save final feature tables.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map (selected_features.csv).")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--suffix', type=str, default="faa", help="Suffix for FASTA files.")
    parser.add_argument('--genome_type', type=str, choices=['strain', 'phage'], default='strain', help="Type of genome to process ('strain' or 'phage'). Default is 'strain'.")
    parser.add_argument('--genome_column', type=str, help="Column name for genome identifiers in genome_list.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.6, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")

    args = parser.parse_args()

    # Run the feature assignment workflow
    run_assign_features_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        output_dir=args.output_dir,
        feature_map=args.feature_map,
        tmp_dir=args.tmp_dir,
        suffix=args.suffix,
        genome_list=args.genome_list,
        genome_type=args.genome_type,
        genome_column=args.genome_column,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
