import os
import subprocess
import pandas as pd
import logging
from argparse import ArgumentParser
from Bio import SeqIO
from phage_modeling.mmseqs2_clustering import create_mmseqs_database, assign_sequences_to_clusters, load_strains, create_contig_to_genome_dict, select_best_hits

def map_features(best_hits_tsv, feature_map, output_dir, strain_contig_mapping, strain):
    """
    Maps the features for each new genome based on cluster assignments.

    Args:
        best_hits_tsv (str): Path to the best hits TSV file (output of select_best_hits).
        feature_map (str): Path to the feature mapping CSV file (selected_features.csv).
        output_dir (str): Directory to save the final feature table.
        strain_contig_mapping (DataFrame): DataFrame mapping contigs to genomes.
        strain (str): Strain name of the current genome.
    """
    if not os.path.exists(best_hits_tsv):
        logging.error("Best hits TSV file does not exist: %s", best_hits_tsv)
        return

    try:
        best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster'])
        feature_mapping = pd.read_csv(feature_map)
    except Exception as e:
        logging.error("Error reading input files: %s", e)
        return

    logging.info("Mapping features for strain '%s'", strain)

    # Merge best_hits_df with feature_mapping
    merged_df = best_hits_df.merge(feature_mapping, left_on='Cluster', right_on='Cluster_Label')
    merged_df = merged_df.merge(strain_contig_mapping, left_on='Query', right_on='contig_id')

    # Create the binary feature presence table
    feature_presence = merged_df.pivot_table(index='strain', columns='Feature', aggfunc='size', fill_value=0)
    feature_presence = (feature_presence > 0).astype(int)

    # Ensure all features are represented
    all_features = feature_mapping['Feature'].unique()
    for feature in all_features:
        if feature not in feature_presence.columns:
            feature_presence[feature] = 0

    # Reindex to ensure the presence of all features
    feature_presence = feature_presence.reindex(columns=all_features, fill_value=0).reset_index()

    output_path = os.path.join(output_dir, f'{strain}_feature_table.csv')
    feature_presence.to_csv(output_path, index=False)
    logging.info("Feature table for strain '%s' saved to %s", strain, output_path)

def process_new_genomes(input_dir, mmseqs_db, suffix, tmp_dir, output_dir, feature_map, clusters_tsv, strains=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4):
    """
    Main function to process new genomes and assign them to existing clusters and features.

    Args:
        input_dir (str): Directory containing new genome FASTA files.
        mmseqs_db (str): Path to the existing MMseqs2 database for searching.
        suffix (str): Suffix for FASTA files.
        tmp_dir (str): Temporary directory for intermediate files.
        output_dir (str): Directory to save final feature tables.
        feature_map (str): Path to the feature map (selected_features.csv).
        clusters_tsv (str): Path to the clusters TSV file.
        strains (list): List of strain names to process (if provided).
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
    """
    os.makedirs(output_dir, exist_ok=True)

    # If strains not provided, infer from files in the input directory
    if strains is None:
        strains = ['.'.join(f.split('.')[:-1]) for f in os.listdir(input_dir) if f.endswith(suffix)]

    for strain in strains:
        strain_tmp_dir = os.path.join(tmp_dir, strain)
        os.makedirs(strain_tmp_dir, exist_ok=True)

        query_db = os.path.join(strain_tmp_dir, 'query_db')
        fasta_files = create_mmseqs_database(input_dir, query_db, suffix, 'directory', [strain], threads)

        result_db = os.path.join(strain_tmp_dir, 'result_db')
        assign_sequences_to_clusters(query_db, strain_tmp_dir, strain_tmp_dir, coverage, min_seq_id, sensitivity, threads, mmseqs_db)

        assigned_tsv = os.path.join(strain_tmp_dir, 'assigned_clusters.tsv')
        best_hits_tsv = os.path.join(strain_tmp_dir, 'best_hits.tsv')
        
        # Use select_best_hits to select the best hit for each query
        select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)

        strain_contig_mapping, _ = create_contig_to_genome_dict(fasta_files, 'directory')

        # Map the features based on the best hits
        map_features(best_hits_tsv, feature_map, output_dir, strain_contig_mapping, strain)

def run_assign_features_workflow(input_dir, mmseqs_db, clusters_tsv, output_dir, feature_map, tmp_dir, suffix="faa", strains_csv=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4):
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
        strains_csv (str): Path to a CSV file with strain names.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
    """
    # Load strains from CSV if provided
    if strains_csv:
        strains = load_strains(strains_csv, 'strain')
    else:
        strains = None  # Process all strains in input_dir if not provided

    process_new_genomes(
        input_dir=input_dir,
        strains=strains,
        mmseqs_db=mmseqs_db,
        suffix=suffix,
        tmp_dir=tmp_dir,
        output_dir=output_dir,
        feature_map=feature_map,
        clusters_tsv=clusters_tsv,
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        threads=threads
    )


def main():
    # Set up argument parser
    parser = ArgumentParser(description="Assign new genes to existing clusters and generate feature tables.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing new genome FASTA files.")
    parser.add_argument('--strains_csv', type=str, help="CSV file with strain names.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the existing MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save final feature tables.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map (selected_features.csv).")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--suffix', type=str, default="faa", help="Suffix for FASTA files.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.6, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")

    args = parser.parse_args()

    # Load strains from CSV if provided, otherwise process all strains in the input directory
    if args.strains_csv:
        strains = load_strains(args.strains_csv, 'strain')
    else:
        strains = None  # This will make the script process all files in the directory

    # Process each strain
    process_new_genomes(
        input_dir=args.input_dir,
        strains=strains,
        mmseqs_db=args.mmseqs_db,
        suffix=args.suffix,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        feature_map=args.feature_map,
        clusters_tsv=args.clusters_tsv,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
