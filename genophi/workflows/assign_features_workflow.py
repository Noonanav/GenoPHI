import os
import subprocess
import pandas as pd
import logging
from argparse import ArgumentParser
from Bio import SeqIO
from collections import defaultdict
import warnings

# Suppress DataFrame fragmentation warning
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from genophi.mmseqs2_clustering import create_mmseqs_database, load_strains, create_contig_to_genome_dict, select_best_hits

def detect_and_modify_duplicates(input_dir, output_dir, suffix='faa', strains_to_process=None, duplicate_all=False):
    """
    Detect and resolve duplicate protein IDs by prefixing them with strain names.

    Args:
        input_dir (str): Directory containing FASTA files.
        output_dir (str): Directory to save modified files.
        suffix (str): File suffix for FASTA files (default is 'faa').
        strains_to_process (list or None): List of strains to process; if None, processes all files.
        duplicate_all (bool): Process all genomes even if duplicates are found.

    Returns:
        str: Path to the directory containing modified FASTA files if duplicates found, else the input directory.
    """
    logging.info("Detecting and resolving duplicate protein IDs...")
    duplicate_found = False
    protein_id_tracker = defaultdict(set)
    modified_dir = os.path.join(output_dir, 'modified_AAs')
    
    # Only create the directory if a duplicate is found
    for file_name in os.listdir(input_dir):
        if file_name.endswith(suffix):
            strain_name = file_name.replace(f".{suffix}", "")
            if not duplicate_all and strains_to_process and strain_name not in strains_to_process:
                continue
            
            file_path = os.path.join(input_dir, file_name)
            modified_file_path = os.path.join(modified_dir, file_name)

            for record in SeqIO.parse(file_path, "fasta"):
                if record.id in protein_id_tracker:
                    duplicate_found = True
                    logging.warning(f"Duplicate protein ID found: {record.id} in strain {strain_name}")
                    break  # Break to avoid unnecessary processing if duplicate found
            
            if duplicate_found:
                break  # No need to continue checking if a duplicate has been found

    if duplicate_found:
        os.makedirs(modified_dir, exist_ok=True)
        logging.info("Duplicate IDs detected. Creating modified files...")
        for file_name in os.listdir(input_dir):
            if file_name.endswith(suffix):
                strain_name = file_name.replace(f".{suffix}", "")
                if not duplicate_all and strains_to_process and strain_name not in strains_to_process:
                    continue
                
                file_path = os.path.join(input_dir, file_name)
                modified_file_path = os.path.join(modified_dir, file_name)

                with open(modified_file_path, "w") as modified_file:
                    for record in SeqIO.parse(file_path, "fasta"):
                        record.id = f"{strain_name}::{record.id}"
                        record.description = ""
                        SeqIO.write(record, modified_file, "fasta")
        logging.info("Resolved duplicate IDs and saved modified files.")
    else:
        logging.info("No duplicate protein IDs found. Using original files.")
        modified_dir = input_dir  # Use the original input directory if no duplicates

    return modified_dir

def assign_sequences_to_clusters(db_name, target_db, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv, clear_tmp):
    """
    Assigns sequences to existing clusters using MMseqs2 search and creates a best hits TSV file.

    Args:
        db_name (str): Path to the MMseqs2 database.
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        sensitivity (float): Sensitivity for assignment.
        threads (int): Number of threads to use.
        clusters_tsv (str): Path to the clusters TSV file.

    Returns:
        str: Path to the best hits TSV file.
    """
    logging.info("Assigning sequences to clusters...")
    
    result_db = os.path.join(tmp_dir, "result_db")
    search_command = (
        f"mmseqs search {db_name} {target_db} {result_db} {tmp_dir} "
        f"-c {coverage} --min-seq-id {min_seq_id} -s {sensitivity} "
        f"--threads {threads} -v 3"
    )
    subprocess.run(search_command, shell=True, check=True)
    logging.info("Sequence assignment completed successfully.")
    
    assigned_tsv = os.path.join(output_dir, "assigned_clusters.tsv")
    createtsv_command = f"mmseqs createtsv {db_name} {target_db} {result_db} {assigned_tsv} --threads {threads} -v 3"
    subprocess.run(createtsv_command, shell=True, check=True)
    logging.info(f"Assigned clusters TSV saved to {assigned_tsv}")
    
    best_hits_tsv = os.path.join(output_dir, "best_hits.tsv")
    select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)

    # Delete the intermediate clustering files
    if clear_tmp:
        logging.info("Clearing temporary files...")
        for file in os.listdir(tmp_dir):
            file_path = os.path.join(tmp_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
    
    return best_hits_tsv

def map_features(best_hits_tsv, feature_map, genome_contig_mapping, genome_type):
    """
    Maps features for all genomes at once based on cluster assignments.

    Args:
        best_hits_tsv (str): Path to the best hits TSV file.
        feature_map (str): Path to the feature mapping CSV file.
        genome_contig_mapping (dict): Dictionary mapping contigs to genomes.
        genome_type (str): Type of genome ('strain' or 'phage').

    Returns:
        pd.DataFrame: Combined feature presence table for all genomes.
    """
    if not os.path.exists(best_hits_tsv):
        logging.error(f"Best hits TSV file does not exist: {best_hits_tsv}")
        return None

    try:
        best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster'])
        feature_mapping = pd.read_csv(feature_map)
        
        best_hits_df['Cluster'] = best_hits_df['Cluster'].astype(str)
        feature_mapping['Cluster_Label'] = feature_mapping['Cluster_Label'].astype(str)
        
    except Exception as e:
        logging.error(f"Error reading input files: {e}")
        return None

    logging.info(f"Mapping features for all {genome_type}s...")

    genome_contig_mapping_df = pd.DataFrame(list(genome_contig_mapping.items()), columns=['contig_id', genome_type])

    merged_df = best_hits_df.merge(feature_mapping, left_on='Cluster', right_on='Cluster_Label')
    merged_df['Query'] = merged_df['Query'].astype(str)
    merged_df = merged_df.merge(genome_contig_mapping_df, left_on='Query', right_on='contig_id')

    feature_presence = merged_df.pivot_table(index=genome_type, columns='Feature', aggfunc='size', fill_value=0)
    feature_presence = (feature_presence > 0).astype(int)

    all_features = feature_mapping['Feature'].unique()
    for feature in all_features:
        if feature not in feature_presence.columns:
            feature_presence[feature] = 0

    feature_presence = feature_presence.reindex(columns=all_features, fill_value=0).reset_index()

    return feature_presence

def run_assign_features_workflow(input_dir, mmseqs_db, tmp_dir, output_dir, feature_map, clusters_tsv, genome_type, genome_list=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4, suffix='faa', duplicate_all=False):
    """
    Process all genomes in the input directory at once, with optional list of genomes to process.

    Args:
        input_dir (str): Directory containing genome FASTA files.
        mmseqs_db (str): Path to the existing MMseqs2 database.
        tmp_dir (str): Temporary directory for intermediate files.
        output_dir (str): Directory to save results.
        feature_map (str): Path to the feature mapping CSV file.
        clusters_tsv (str): Path to the clusters TSV file.
        genome_type (str): Type of genomes ('strain' or 'phage').
        genome_list (str or None): Path to a file with strain names or None for all.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
        suffix (str): Suffix for FASTA files.
        duplicate_all (bool): Process all genomes even if duplicates are found.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Load strains if a list is provided
    strains_to_process = None
    if genome_list and os.path.exists(genome_list):
        strains_to_process = load_strains(genome_list, genome_type)

    # Use detect_and_modify_duplicates with a check for duplicates
    dir_to_use = detect_and_modify_duplicates(input_dir, output_dir, suffix, strains_to_process, duplicate_all=duplicate_all)

    logging.info(f"Creating a combined MMseqs2 database for all {genome_type}s...")
    combined_db = os.path.join(tmp_dir, "combined_db")
    fasta_files = create_mmseqs_database(dir_to_use, combined_db, suffix, 'directory', strains_to_process, threads)

    if not fasta_files:
        logging.error("No FASTA files found for processing.")
        return

    logging.info(f"Assigning {genome_type} sequences to existing clusters...")
    result_db = os.path.join(tmp_dir, "result_db")
    assign_sequences_to_clusters(combined_db, mmseqs_db, tmp_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv, clear_tmp=False)

    assigned_tsv = os.path.join(tmp_dir, 'assigned_clusters.tsv')
    best_hits_tsv = os.path.join(tmp_dir, 'best_hits.tsv')
    select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)

    genome_contig_mapping, _ = create_contig_to_genome_dict(fasta_files, 'directory')

    logging.info(f"Generating feature tables for all {genome_type}s...")
    feature_presence = map_features(best_hits_tsv, feature_map, genome_contig_mapping, genome_type)

    if feature_presence is not None:
        combined_feature_table_path = os.path.join(output_dir, f'{genome_type}_combined_feature_table.csv')
        feature_presence.to_csv(combined_feature_table_path, index=False)
        logging.info(f"Combined feature table saved to {combined_feature_table_path}")

def main():
    parser = ArgumentParser(description="Process all genomes and generate combined feature tables.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the existing MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map (selected_features.csv).")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--genome_type', type=str, choices=['strain', 'phage'], default='strain', help="Type of genome to process.")
    parser.add_argument('--genome_list', type=str, help="Path to file with list of genomes to process.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.6, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")
    parser.add_argument('--suffix', type=str, default='faa', help="Suffix for FASTA files.")
    parser.add_argument('--duplicate_all', action='store_true', help="Process all genomes even if duplicates are found.")

    args = parser.parse_args()

    run_assign_features_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        feature_map=args.feature_map,
        clusters_tsv=args.clusters_tsv,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        duplicate_all=args.duplicate_all
    )

if __name__ == "__main__":
    main()