import os
import pandas as pd
import subprocess
import logging
import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from collections import defaultdict
from sklearn.cluster import AgglomerativeClustering

# Set thread limits for libraries to prevent oversubscription
os.environ['OMP_NUM_THREADS'] = '12'
os.environ['OPENBLAS_NUM_THREADS'] = '12'
os.environ['MKL_NUM_THREADS'] = '12'
os.environ['VECLIB_MAXIMUM_THREADS'] = '12'
os.environ['NUMEXPR_NUM_THREADS'] = '12'
os.environ['NUMBA_NUM_THREADS'] = '12'

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Utility Functions
def load_strains(strain_list, strain_column):
    """
    Loads a list of strains from a CSV file and returns the unique strains.

    Args:
        strain_list (str): Path to the strain list CSV file.
        strain_column (str): Name of the column that contains strain names.

    Returns:
        list: List of unique strain names.

    Raises:
        ValueError: If the specified column is not found in the file.
        Exception: If an error occurs while reading the file.
    """
    try:
        strains_df = pd.read_csv(strain_list)
        if strain_column not in strains_df.columns:
            raise ValueError(f"Column '{strain_column}' not found in strain list file.")
        return list(strains_df[strain_column].unique())
    except Exception as e:
        logging.error(f"Error loading strain list: {e}")
        raise
    
def detect_duplicate_ids(input_path, suffix='faa', strains_to_process=None, input_type='directory'):
    """
    Detects duplicate protein IDs across relevant .faa files in the input directory.
    
    Args:
        input_path (str): Directory containing .faa files.
        suffix (str): Suffix for the input FASTA files (default is 'faa').
        strains_to_process (list or None): List of strains to process; if None, processes all .faa files.
    
    Returns:
        bool: True if duplicates are found, False otherwise.
    """
    logging.info("Detecting duplicate protein IDs...")
    duplicate_found = False
    protein_id_tracker = defaultdict(set)

    if input_type == 'directory':
        file_list = os.listdir(input_path)
    elif input_type == 'file':
        file_list = [input_path]
    else:
        logging.error(f"Invalid input type: {input_type}")
        return duplicate_found

    # Identify and track protein IDs in specified strains
    for file_name in file_list:
        if file_name.endswith(suffix):
            strain_name = file_name.replace(f".{suffix}", "")
            if strains_to_process and strain_name not in strains_to_process:
                continue  # Skip strains not in the strain list
            
            file_path = os.path.join(input_path, file_name)
            for record in SeqIO.parse(file_path, "fasta"):
                if record.id in protein_id_tracker:
                    duplicate_found = True
                    logging.warning(f"Duplicate protein ID found: {record.id} in strain {strain_name}")
                    break
                protein_id_tracker[record.id].add(strain_name)
    
    return duplicate_found

def modify_duplicate_ids(input_path, output_dir, suffix='faa', strains_to_process=None, strain_column='strain'):
    """
    Detects duplicate protein IDs and modifies all protein IDs in relevant .faa files 
    by prefixing them with genome names to ensure uniqueness.

    Args:
        input_path (str): Directory containing .faa files.
        suffix (str): Suffix for the input FASTA files (default is 'faa').
        strains_to_process (list or None): List of strains to process; if None, processes all .faa files.
    
    Returns:
        list: Path to modified .faa files.
    """
    logging.info(f"Duplicate protein IDs found; modifying protein IDs and saving to {output_dir}/modified_AAs/{strain_column}")
    modified_file_dir = os.path.join(output_dir, 'modified_AAs', strain_column)
    os.makedirs(modified_file_dir, exist_ok=True)

    for file_name in os.listdir(input_path):
        if file_name.endswith(suffix):
            strain_name = file_name.replace(f".{suffix}", "")
            if strains_to_process and strain_name not in strains_to_process:
                continue
            
            file_path = os.path.join(input_path, file_name)
            modified_file_path = os.path.join(modified_file_dir, file_name)
            
            with open(modified_file_path, "w") as modified_file:
                for record in SeqIO.parse(file_path, "fasta"):
                    # Update ID with <genome_id>::<protein_ID> format
                    record.id = f"{strain_name}::{record.id}"
                    record.description = ""  # Clear description to avoid duplication
                    SeqIO.write(record, modified_file, "fasta")

            logging.info(f"Modified protein IDs in file: {modified_file_path}")
    
    return modified_file_dir

def create_mmseqs_database(input_path, db_name, suffix, input_type, strains, threads):
    """
    Creates an MMseqs2 database from FASTA files.

    Args:
        input_path (str): Path to the input directory or file.
        db_name (str): Name of the output MMseqs2 database.
        suffix (str): Suffix for FASTA files to be included.
        input_type (str): Type of input, either 'directory' or 'file'.
        strains (list or None): List of strains to filter, or None for no filtering.
        threads (int): Number of threads to use for MMseqs2.

    Returns:
        list: List of processed FASTA file paths.

    Raises:
        FileNotFoundError: If no matching FASTA files are found.
    """
    logging.info("Creating MMseqs2 database...")
    fasta_files = []
    strains = [str(s) for s in strains] if strains else None
    if input_type == 'directory':
        logging.info(f"Searching for FASTA files with suffix '{suffix}' in {input_path}.")
        for fasta in os.listdir(input_path):
            if fasta.endswith(suffix):
                strain_name = fasta.replace(f".{suffix}", "")
                if strains is None or strain_name in strains:
                    fasta_files.append(os.path.join(input_path, fasta))
    else:
        logging.info(f"Using input file: {input_path}")
        fasta_files.append(input_path)
    
    if not fasta_files:
        logging.error("No matching FASTA files found. Exiting.")
        raise FileNotFoundError("No FASTA files to process.")
    
    logging.info(f"Creating MMseqs2 database for {len(fasta_files)} files using xargs")
    temp_file = os.path.join(os.path.dirname(db_name), "file_list.tmp")
    with open(temp_file, 'w') as f:
        for fasta_file in fasta_files:
            f.write(f"{fasta_file}\0")

    createdb_command = f"cat {temp_file} | xargs -0 cat | mmseqs createdb stdin {db_name} -v 3"
    subprocess.run(createdb_command, shell=True, check=True)
    os.remove(temp_file)  # Clean up
    
    return fasta_files

def validate_checkpoint_file(filepath, min_size=1, file_type='general'):
    """
    Validates if a checkpoint file exists and meets basic criteria.
    
    Args:
        filepath (str): Path to the file to validate
        min_size (int): Minimum file size in bytes
        file_type (str): Type of file for specific validation ('tsv', 'csv', 'mmseqs_db')
    
    Returns:
        bool: True if file is valid, False otherwise
    """
    if not os.path.exists(filepath):
        return False
    
    # Check file size
    if os.path.getsize(filepath) < min_size:
        logging.warning(f"File {filepath} exists but is too small ({os.path.getsize(filepath)} bytes)")
        return False
    
    # Type-specific validation
    if file_type == 'mmseqs_db':
        # For MMseqs database, check for .dbtype file
        dbtype_file = f"{filepath}.dbtype"
        if not os.path.exists(dbtype_file):
            logging.warning(f"MMseqs database {filepath} missing .dbtype file")
            return False
        return True
    
    elif file_type in ['tsv', 'csv']:
        # Basic structure validation for tabular files
        try:
            delimiter = '\t' if file_type == 'tsv' else ','
            # Just read first few lines to check structure
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                if not first_line:
                    logging.warning(f"File {filepath} appears to be empty")
                    return False
                
                # Check if it's properly delimited
                parts = first_line.split(delimiter)
                if len(parts) < 2:
                    logging.warning(f"File {filepath} doesn't appear to be properly {file_type} formatted")
                    return False
            return True
        except Exception as e:
            logging.warning(f"Error validating {filepath}: {e}")
            return False
    
    # For general files, just check existence and size
    return True

def cleanup_partial_results(output_dir, tmp_dir, stage):
    """
    Removes incomplete intermediate files for a given stage.
    
    Args:
        output_dir (str): Output directory path
        tmp_dir (str): Temporary directory path  
        stage (str): Stage to clean ('clustering', 'assignment', 'matrix')
    """
    logging.info(f"Cleaning up partial results for stage: {stage}")
    
    try:
        if stage == 'clustering':
            # Remove clustering outputs but keep database
            files_to_remove = [
                os.path.join(output_dir, "clusters"),
                os.path.join(output_dir, "clusters.tsv"),
            ]
            # Also clean any clustering intermediates in tmp
            for file in os.listdir(tmp_dir):
                if 'cluster' in file.lower():
                    files_to_remove.append(os.path.join(tmp_dir, file))
                    
        elif stage == 'assignment':
            # Remove assignment outputs but keep clustering results
            files_to_remove = [
                os.path.join(output_dir, "assigned_clusters.tsv"),
                os.path.join(output_dir, "best_hits.tsv"),
                os.path.join(tmp_dir, "result_db"),
            ]
            # Clean result_db related files
            for file in os.listdir(tmp_dir):
                if file.startswith('result_db'):
                    files_to_remove.append(os.path.join(tmp_dir, file))
                    
        elif stage == 'matrix':
            # Remove matrix output but keep assignment results
            files_to_remove = [
                os.path.join(output_dir, "presence_absence_matrix.csv"),
            ]
            
        # Remove files/directories
        for item in files_to_remove:
            try:
                if os.path.isfile(item):
                    os.remove(item)
                    logging.debug(f"Removed file: {item}")
                elif os.path.isdir(item):
                    import shutil
                    shutil.rmtree(item)
                    logging.debug(f"Removed directory: {item}")
            except FileNotFoundError:
                pass  # File already doesn't exist
            except Exception as e:
                logging.warning(f"Could not remove {item}: {e}")
                
        logging.info(f"Cleanup completed for stage: {stage}")
        
    except Exception as e:
        logging.error(f"Error during cleanup of stage {stage}: {e}")

def create_contig_to_genome_dict(fasta_files, input_type, suffix='faa'):
    """
    Creates a mapping from contigs to genomes based on input FASTA files.

    Args:
        fasta_files (list): List of FASTA file paths.
        input_type (str): Type of input, either 'directory' or 'file'.

    Returns:
        tuple: Dictionary mapping contig IDs to genome names, and a list of genome names.
    """
    contig_to_genome = {}
    genome_list = []
    logging.info("Creating contig to genome dictionary...")
    
    if input_type == 'directory':
        for fasta in fasta_files:
            genome_name = os.path.basename(fasta).replace(f".{suffix}", "")
            genome_list.append(genome_name)
            for record in SeqIO.parse(fasta, "fasta"):
                contig_to_genome[record.id] = genome_name
    else:
        for record in SeqIO.parse(fasta_files[0], "fasta"):
            protein_id = record.id
            genome_name = '_'.join(protein_id.split(' # ')[0].split('_')[:-1])
            genome_list.append(genome_name)
            contig_to_genome[protein_id] = genome_name
        genome_list = list(set(genome_list))
    
    logging.info(f"Created contig to genome dictionary with {len(contig_to_genome)} entries.")
    return contig_to_genome, genome_list

def run_mmseqs_cluster(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads):
    """
    Runs MMseqs2 clustering and creates a clusters TSV file.

    Args:
        db_name (str): Path to the MMseqs2 database.
        output_dir (str): Directory to save clustering results.
        tmp_dir (str): Temporary directory for intermediate files.
        coverage (float): Minimum coverage for clustering.
        min_seq_id (float): Minimum sequence identity for clustering.
        sensitivity (float): Sensitivity for clustering.
        threads (int): Number of threads to use for clustering.

    Returns:
        str: Path to the clusters TSV file.
    """
    logging.info("Running MMseqs2 clustering...")
    cluster_output = os.path.join(output_dir, "clusters")
    cluster_command = (
        f"mmseqs cluster {db_name} {cluster_output} {tmp_dir} "
        f"-c {coverage} --min-seq-id {min_seq_id} -s {sensitivity} "
        f"--threads {threads} -v 3"
    )
    subprocess.run(cluster_command, shell=True, check=True)
    logging.info("Clustering completed successfully.")
    
    clusters_tsv = os.path.join(output_dir, "clusters.tsv")
    createtsv_command = f"mmseqs createtsv {db_name} {db_name} {cluster_output} {clusters_tsv} --threads {threads} -v 3"
    subprocess.run(createtsv_command, shell=True, check=True)
    logging.info(f"Clusters TSV saved to {clusters_tsv}")
    
    return clusters_tsv

def assign_sequences_to_clusters(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv, clear_tmp):
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
    assigned_tsv = os.path.join(output_dir, "assigned_clusters.tsv")
    best_hits_tsv = os.path.join(output_dir, "best_hits.tsv")
    
    if validate_checkpoint_file(best_hits_tsv, file_type='tsv'):
        logging.info("Found existing best_hits.tsv, skipping assignment")
        return best_hits_tsv
        
    if validate_checkpoint_file(assigned_tsv, file_type='tsv'):
        logging.info("Found existing assigned_clusters.tsv, processing to best hits") 
        select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)
        return best_hits_tsv
    
    result_db = os.path.join(tmp_dir, "result_db")
    search_command = (
        f"mmseqs search {db_name} {db_name} {result_db} {tmp_dir} "
        f"-c {coverage} --min-seq-id {min_seq_id} -s {sensitivity} "
        f"--threads {threads} -v 3"
    )
    subprocess.run(search_command, shell=True, check=True)
    logging.info("Sequence assignment completed successfully.")
    
    assigned_tsv = os.path.join(output_dir, "assigned_clusters.tsv")
    createtsv_command = f"mmseqs createtsv {db_name} {db_name} {result_db} {assigned_tsv} --threads {threads} -v 3"
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

def select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv):
    """
    Selects the best hit for each query and saves the results.

    Args:
        assigned_tsv (str): Path to the assigned clusters TSV file.
        best_hits_tsv (str): Path to save the best hits TSV file.
        clusters_tsv (str): Path to the clusters TSV file.
    """
    logging.info("Selecting best hits...")

    # Read the assigned clusters file and ensure the correct column names
    dtype_dict = {
        'Query': 'string', 'Target': 'string', 'Score': 'float32',
        'SeqIdentity': 'float32', 'E-value': 'float32',
        'qStartPos': 'int32', 'qEndPos': 'int32', 'qLen': 'int32',
        'tStartPos': 'int32', 'tEndPos': 'int32', 'tLen': 'int32'
    }

    assigned_df = pd.read_csv(assigned_tsv, sep='\t', header=None, 
                            names=list(dtype_dict.keys()), dtype=dtype_dict)

    # Confirm assigned_df is not empty
    if assigned_df.empty:
        logging.error(f"No data found in assigned clusters file: {assigned_tsv}")
        return

    # Sort by Query, then by SeqIdentity and E-value to find the best hits
    assigned_df = assigned_df.sort_values(by=['Query', 'SeqIdentity', 'E-value'], ascending=[True, False, True])

    # Drop duplicates to keep the best hit for each Query
    best_hits_df = assigned_df.drop_duplicates(subset=['Query'], keep='first')

    # Read the clusters TSV and ensure the columns are correctly named
    clusters_df = pd.read_csv(clusters_tsv, sep='\t', header=None, names=['Cluster', 'Contig'])

    # Confirm clusters_df is not empty
    if clusters_df.empty:
        logging.error(f"No data found in clusters file: {clusters_tsv}")
        return

    # Ensure both columns being merged on are strings to avoid type mismatch
    best_hits_df['Target'] = best_hits_df['Target'].astype(str)
    clusters_df['Contig'] = clusters_df['Contig'].astype(str)

    # Merge best hits with clusters to map the target sequences to clusters
    best_hits_df = pd.merge(best_hits_df, clusters_df, left_on='Target', right_on='Contig', how='left')

    # Check for any missing clusters
    missing_clusters = best_hits_df['Cluster'].isna().sum()
    if missing_clusters > 0:
        logging.warning(f"{missing_clusters} entries in 'Target' did not match any 'Contig' in clusters file.")
    
    # Save the final result with only the Query and Cluster columns
    best_hits_df[['Query', 'Cluster']].to_csv(best_hits_tsv, sep='\t', index=False, header=False)
    logging.info(f"Best hits saved to {best_hits_tsv}")

def generate_presence_absence_matrix(best_hits_tsv, output_csv_path, contig_to_genome, genome_list):
    """
    Generates a presence-absence matrix based on best hits and cluster assignments.

    Args:
        best_hits_tsv (str): Path to the best hits TSV file.
        output_csv_path (str): Path to save the presence-absence matrix.
        contig_to_genome (dict): Dictionary mapping contigs to genomes.
        genome_list (list): List of genome names.
    """
    logging.info("Generating presence-absence matrix...")
    
    cluster_dict = {}
    with open(best_hits_tsv, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            contig_id, cluster_id = parts[0], parts[1]
            genome = contig_to_genome.get(contig_id)
            if not genome:
                logging.warning(f"Contig ID '{contig_id}' not found in mapping.")
                continue
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = set()
            cluster_dict[cluster_id].add(genome)

    presence_absence_data = {'Genome': genome_list}
    for cluster_id, genomes in cluster_dict.items():
        presence_absence_data[cluster_id] = [1 if genome in genomes else 0 for genome in genome_list]

    presence_absence_df = pd.DataFrame(presence_absence_data)
    presence_absence_df.to_csv(output_csv_path, index=False)
    logging.info(f"Presence-absence matrix saved to {output_csv_path}")

def compare_cluster_and_search_results(clusters_tsv, best_hits_tsv, output_dir):
    """
    Compares the original clustering results to the best hit assignments and saves the comparison.

    Args:
        clusters_tsv (str): Path to the original clusters TSV file.
        best_hits_tsv (str): Path to the best hits TSV file.
        output_dir (str): Directory to save the comparison results.
    """
    logging.info("Comparing clusters...")
    
    clusters_df = pd.read_csv(clusters_tsv, sep='\t', header=None, names=['Cluster_original', 'Contig'])
    best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster_assigned'])
    
    comparison_df = pd.merge(clusters_df, best_hits_df, left_on='Contig', right_on='Query', how='outer', suffixes=('_original', '_assigned'))
    comparison_output_path = os.path.join(output_dir, "cluster_assignment_comparison.csv")
    comparison_df.to_csv(comparison_output_path, index=False)
    logging.info(f"Cluster assignment comparison saved to {comparison_output_path}")

# Additional script functions
def filter_presence_absence(presence_absence, select, select_column):
    """
    Filters the presence-absence matrix based on a selected strain list and removes clusters
    present in only one genome.

    Args:
        presence_absence (DataFrame): The presence-absence matrix.
        select (str): Path to the file containing selected strains.
        select_column (str): Column in the selection file that contains strain names.

    Returns:
        DataFrame: Filtered presence-absence matrix.
    """
    logging.info("Filtering presence-absence table...")

    # Filter by strain list if provided
    if select:
        select_df = pd.read_csv(select)
        select_list = list(set(select_df[select_column].tolist()))
        presence_absence = presence_absence[presence_absence['Genome'].isin(select_list)]
        logging.info(f"Filtered to {len(select_list)} selected genomes.")

    # Remove clusters present in only one genome
    cluster_counts = presence_absence.iloc[:, 1:].sum(axis=0)  # Sum across rows for each cluster
    valid_clusters = cluster_counts[cluster_counts > 1].index  # Keep clusters with counts > 1
    all_cols = len(presence_absence.columns)
    presence_absence = presence_absence[['Genome'] + list(valid_clusters)]
    logging.info(f"Removed {all_cols - len(valid_clusters) - 1} clusters present in only one genome.")

    # Remove columns with all zeros
    all_cols = len(presence_absence.columns)
    presence_absence = presence_absence.loc[:, (presence_absence != 0).any(axis=0)]
    logging.info(f"Removed {all_cols - len(presence_absence.columns)} columns with all zeros.")

    return presence_absence

def get_genome_assignments_tables(presence_absence, genome_column_name):
    """
    Generates genome assignments from the presence-absence matrix.

    Args:
        presence_absence (DataFrame): The presence-absence matrix.
        genome_column_name (str): The column name that contains genome information (e.g., 'strain' or 'phage').

    Returns:
        DataFrame: Genome assignments in long format.
    """
    logging.info("Getting genome assignments...")
    presence_absence.rename(columns={'Genome': genome_column_name}, inplace=True)
    genome_assignments = presence_absence.melt(id_vars=genome_column_name, var_name="Cluster_Label", value_name="Presence")
    genome_assignments = genome_assignments[genome_assignments['Presence'] == 1]
    return genome_assignments.drop(columns=["Presence"])

def feature_selection_optimized(presence_absence, source, genome_column_name):
    """
    Optimizes feature selection by identifying perfect co-occurrence of features using hashing.

    Args:
        presence_absence (DataFrame): The presence-absence matrix.
        source (str): A prefix for naming the selected features.
        genome_column_name (str): The column name that contains genome information (e.g., 'strain' or 'phage').

    Returns:
        DataFrame: Optimized feature selection results.
    """
    logging.info("Optimizing feature selection using hashing...")

    # Set index using the genome_column_name
    presence_absence.set_index(genome_column_name, inplace=True)

    # Ensure binary presence-absence format
    presence_absence = presence_absence.applymap(lambda x: 1 if x > 0 else 0)

    # Compute hashes for each column
    logging.info("Hashing columns to identify identical patterns...")
    column_hashes = presence_absence.apply(lambda col: hash(tuple(col)), axis=0)

    # Group columns by hash
    logging.info("Grouping columns by hash...")
    hash_to_columns = {}
    for col, col_hash in column_hashes.items():
        hash_to_columns.setdefault(col_hash, []).append(col)

    # Identify unique clusters based on hash groups
    unique_clusters = list(hash_to_columns.values())
    logging.info(f"Identified {len(unique_clusters)} unique clusters.")

    # Prepare output DataFrame
    logging.info("Preparing output DataFrame...")
    data = [
        (f"{source[0]}c_{idx}", cluster)
        for idx, cluster_group in enumerate(unique_clusters)
        for cluster in cluster_group
    ]
    logging.info(f"Feature selection completed with {len(unique_clusters)} unique clusters.")

    selected_features = pd.DataFrame(data, columns=["Feature", "Cluster_Label"])

    return selected_features

def feature_assignment(genome_assignments, selected_features, genome_column_name):
    """
    Assigns features to genomes based on the selected features.

    Args:
        genome_assignments (DataFrame): Genome assignments.
        selected_features (DataFrame): Selected features.
        genome_column_name (str): The column name that contains genome information (e.g., 'strain' or 'phage').

    Returns:
        tuple: DataFrame of feature assignments and feature table in wide format.
    """
    logging.info("Assigning features to genomes...")
    
    # Merge genome assignments with selected features
    assignment_df = genome_assignments.merge(selected_features, on="Cluster_Label", how="inner")
    assignment_df = assignment_df.drop(columns=["Cluster_Label"]).drop_duplicates()

    # Create the feature table in wide format using pivot_table
    feature_table = assignment_df.pivot_table(index=genome_column_name, columns="Feature", aggfunc="size", fill_value=0)
    
    # Reset the index to turn the genome_column_name (strain/phage) back into a regular column
    feature_table = feature_table.reset_index()

    return assignment_df, feature_table


def run_clustering_workflow(input_path, output_dir, tmp_dir="tmp", min_seq_id=0.6, coverage=0.8, sensitivity=7.5, suffix='faa', threads=4, strain_list='none', strain_column='strain', compare=False, bootstrapping=False, clear_tmp=False, force_restart=False):
    """
    Runs a full MMseqs2 clustering workflow including presence-absence matrix generation.
    
    This function processes input FASTA files to create a database using MMseqs2, runs clustering, 
    and assigns sequences to clusters. It then generates a presence-absence matrix for the identified clusters.
    
    Args:
        input_path (str): Path to the input directory or file containing FASTA sequences.
        output_dir (str): Directory to save clustering and presence-absence results.
        tmp_dir (str): Temporary directory for intermediate files.
        min_seq_id (float): Minimum sequence identity for clustering.
        coverage (float): Minimum coverage required for clustering.
        sensitivity (float): Sensitivity level for clustering.
        suffix (str): Suffix for input FASTA files to include in the database.
        threads (int): Number of threads to use for MMseqs2.
        strain_list (str): Path to a strain list file, or 'none' to process all strains.
        strain_column (str): Column name in the strain list file that contains strain names.
        compare (bool): Whether to compare the original clusters with assigned clusters.
        bootstrapping (bool): Whether this is part of bootstrapping workflow.
        clear_tmp (bool): Whether to clear the temporary directory after processing.
        force_restart (bool): Whether to bypass all checkpoints and restart from beginning.
    
    Raises:
        FileNotFoundError: If no FASTA files are found in the input path.
        ValueError: If there is an issue with loading the strain list.
    """

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Define paths
    presence_absence_csv = os.path.join(output_dir, "presence_absence_matrix.csv")
    clusters_tsv = os.path.join(output_dir, "clusters.tsv")
    best_hits_tsv = os.path.join(output_dir, "best_hits.tsv")
    db_name = os.path.join(tmp_dir, "mmseqs_db")
    
    # Quick exit if final output exists
    if not force_restart and validate_checkpoint_file(presence_absence_csv, file_type='csv'):
        logging.info("Found complete presence-absence matrix. Skipping entire workflow.")
        if compare:
            compare_cluster_and_search_results(clusters_tsv, best_hits_tsv, output_dir)
        return

    # Validate checkpoints and clean inconsistent states
    database_valid = not force_restart and validate_checkpoint_file(db_name, file_type='mmseqs_db')
    clustering_valid = not force_restart and validate_checkpoint_file(clusters_tsv, file_type='tsv')
    assignment_valid = not force_restart and validate_checkpoint_file(best_hits_tsv, file_type='tsv')
    
    # Clean inconsistent states
    if assignment_valid and not clustering_valid:
        logging.info("Assignment exists without clustering - cleaning assignment")
        cleanup_partial_results(output_dir, tmp_dir, 'assignment')
        assignment_valid = False
        
    if clustering_valid and not database_valid:
        logging.info("Clustering exists without database - cleaning clustering")  
        cleanup_partial_results(output_dir, tmp_dir, 'clustering')
        clustering_valid = False
        
    # Log what will be reused
    if database_valid:
        logging.info("Reusing existing database")
    if clustering_valid:
        logging.info("Reusing existing clustering")
    if assignment_valid:
        logging.info("Reusing existing assignment")

    # Prepare workflow variables (always needed)
    input_type = 'directory' if os.path.isdir(input_path) else 'file'
    strain_list_value = None if strain_list == 'none' else strain_list
    strains_to_process = load_strains(strain_list_value, strain_column) if strain_list_value else None

    # Handle duplicate IDs (this may modify input_path)
    # Check if modified directory already exists to avoid unnecessary work
    expected_modified_dir = os.path.join(output_dir, 'modified_AAs', strain_column)
    duplicate_found = detect_duplicate_ids(input_path, suffix, strains_to_process, input_type)
    
    if duplicate_found:
        if input_type == 'directory':
            if os.path.exists(expected_modified_dir) and os.listdir(expected_modified_dir):
                logging.info("Using existing modified AA files")
                input_path = expected_modified_dir
            else:
                if bootstrapping:
                    input_path = modify_duplicate_ids(input_path, output_dir, suffix, None, strain_column)
                else:
                    input_path = modify_duplicate_ids(input_path, output_dir, suffix, strains_to_process, strain_column)
        else:
            logging.error("Duplicate protein IDs found in input file; please modify the IDs manually.")
            return

    # Stage 1: Database creation
    if not database_valid:
        fasta_files = create_mmseqs_database(input_path, db_name, suffix, input_type, strains_to_process, threads)
    else:
        # Build fasta_files list for later stages (using same logic as create_mmseqs_database)
        fasta_files = []
        strains = [str(s) for s in strains_to_process] if strains_to_process else None
        if input_type == 'directory':
            for fasta in os.listdir(input_path):
                if fasta.endswith(suffix):
                    strain_name = fasta.replace(f".{suffix}", "")
                    if strains is None or strain_name in strains:
                        fasta_files.append(os.path.join(input_path, fasta))
        else:
            fasta_files.append(input_path)

    # Create contig mapping (always needed for matrix generation)
    contig_to_genome, genome_list = create_contig_to_genome_dict(fasta_files, input_type, suffix)

    # Stage 2: Clustering  
    if not clustering_valid:
        clusters_tsv = run_mmseqs_cluster(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads)

    # Stage 3: Assignment
    if not assignment_valid:
        best_hits_tsv = assign_sequences_to_clusters(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv, clear_tmp)

    # Stage 4: Matrix generation (only if not already valid)
    if not validate_checkpoint_file(presence_absence_csv, file_type='csv'):
        generate_presence_absence_matrix(best_hits_tsv, presence_absence_csv, contig_to_genome, genome_list)
    else:
        logging.info("Reusing existing presence-absence matrix")

    if compare:
        compare_cluster_and_search_results(clusters_tsv, best_hits_tsv, output_dir)

def run_feature_assignment(input_file, output_dir, source='strain', select='none', select_column='strain', input_type='directory', max_ram=8, threads=4):
    """
    Runs the feature assignment workflow from the presence-absence matrix.
    
    This function processes a presence-absence matrix, assigns features to genomes based on selected features, 
    and outputs both feature assignments and a feature table. It allows for optional filtering based on a strain list.
    
    Args:
        input_file (str): Path to the presence-absence matrix CSV file.
        output_dir (str): Directory to save the results (selected features, assignments, and feature table).
        source (str): Prefix for naming the selected features (e.g., 'strain' or 'phage').
        select (str): Path to a strain list file for filtering, or 'none' to skip filtering.
        select_column (str): Column name in the strain list file that contains strain names.
    
    Raises:
        FileNotFoundError: If the input presence-absence matrix file is not found.
        ValueError: If there is an issue with the strain list during filtering.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load the presence-absence matrix
    presence_absence = pd.read_csv(input_file)

    # Convert 'none' to None for filtering purposes
    select_value = None if select == 'none' else select

    # If a strain list is provided, filter the presence-absence matrix
    if select_value:
        if input_type == 'file':
            if select_column in presence_absence.columns:
                logging.info(f"Filtering presence-absence matrix based on selected {select_column}...")
                presence_absence = filter_presence_absence(presence_absence, select_value, select_column)
            elif 'Genome' in presence_absence.columns:
                logging.info(f"Renaming 'Genome' column to {select_column} and filtering presence-absence matrix...")
                # Rename 'Genome' column to the selected column name
                presence_absence = presence_absence.rename(columns={'Genome': select_column})
                presence_absence = filter_presence_absence(presence_absence, select_value, select_column)
            else:
                logging.error(f"Column '{select_column}' not found in presence-absence matrix.")
                raise ValueError(f"Column '{select_column}' not found in presence-absence matrix.")
        else:
            presence_absence = filter_presence_absence(presence_absence, select_value, select_column)

    genome_column_name = source

    # Assign features to genomes
    genome_assignments = get_genome_assignments_tables(presence_absence, genome_column_name)

    # Pass the correct genome column name to feature_selection_optimized
    selected_features = feature_selection_optimized(presence_absence, source, genome_column_name)

    # Save the selected features and feature assignments
    selected_features_path = os.path.join(output_dir, 'selected_features.csv')
    selected_features.to_csv(selected_features_path, index=False)

    feature_assignments, feature_table = feature_assignment(genome_assignments, selected_features, genome_column_name)
    feature_assignments.to_csv(os.path.join(output_dir, 'feature_assignments.csv'), index=False)
    feature_table.to_csv(os.path.join(output_dir, 'feature_table.csv'), index=False)

def cluster_and_filter_features(
    feature_table,
    feature_type='strain', 
    sample_column='strain',
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    output_dir=None
):
    """
    Clusters samples by feature content and filters features based on cluster presence.
    This is for pre-processing, separate from train-test split clustering.
    
    Returns:
        filtered_feature_table (DataFrame): Feature table with cluster-filtered features
        cluster_assignments (DataFrame): Sample-to-cluster assignments (optional)
    """
    if not use_feature_clustering:
        return feature_table, None
        
    # Extract feature columns for clustering
    feature_columns = [col for col in feature_table.columns if col.startswith(f"{feature_type[0]}c_")]
    
    if not feature_columns:
        logging.warning(f"No {feature_type} feature columns found for clustering.")
        return feature_table, None
    
    # Prepare clustering data
    clustering_data = feature_table[[sample_column] + feature_columns].drop_duplicates()
    
    # Handle edge case where we have fewer samples than requested clusters
    n_samples = len(clustering_data)
    if n_samples < 2:
        logging.warning(f"Insufficient samples ({n_samples}) for clustering. Skipping feature clustering.")
        return feature_table, None
    
    actual_n_clusters = min(feature_n_clusters, n_samples - 1)
    if actual_n_clusters != feature_n_clusters:
        logging.warning(f"Reduced clusters from {feature_n_clusters} to {actual_n_clusters} due to sample size")
    
    # Perform hierarchical clustering
    clusterer = AgglomerativeClustering(n_clusters=actual_n_clusters)
    cluster_labels = clusterer.fit_predict(clustering_data[feature_columns])
    
    # Create cluster assignments
    clustering_data['cluster'] = cluster_labels
    cluster_assignments = clustering_data[[sample_column, 'cluster']]
    
    # Merge back to full feature table
    feature_table_with_clusters = feature_table.merge(cluster_assignments, on=sample_column, how='left')
    
    # Filter features by cluster presence
    if feature_min_cluster_presence > 1:
        # Group by cluster and see which features are present
        cluster_feature_presence = feature_table_with_clusters.groupby('cluster')[feature_columns].apply(
            lambda group: (group > 0).any()
        )
        
        # Count how many clusters each feature appears in
        feature_cluster_counts = cluster_feature_presence.sum(axis=0)
        
        # Keep features that appear in at least feature_min_cluster_presence clusters
        valid_features = feature_cluster_counts[feature_cluster_counts >= feature_min_cluster_presence].index.tolist()
        
        # Filter the feature table
        other_columns = [col for col in feature_table.columns if col not in feature_columns]
        filtered_feature_table = feature_table[other_columns + valid_features]
        
        logging.info(f"Pre-processing cluster filtering: kept {len(valid_features)}/{len(feature_columns)} {feature_type} features")
        logging.info(f"Removed {len(feature_columns) - len(valid_features)} features present in < {feature_min_cluster_presence} clusters")
    else:
        filtered_feature_table = feature_table
    
    # Save cluster assignments if output_dir provided
    if output_dir:
        cluster_file = os.path.join(output_dir, f"{feature_type}_preprocess_clusters.csv")
        cluster_assignments.to_csv(cluster_file, index=False)
        logging.info(f"Pre-processing cluster assignments saved to: {cluster_file}")
    
    return filtered_feature_table, cluster_assignments

def merge_feature_tables(
    strain_features, 
    phenotype_matrix, 
    output_dir, 
    sample_column='strain', 
    phage_features=None, 
    remove_suffix=False, 
    output_file=None,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical', 
    feature_n_clusters=20,
    feature_min_cluster_presence=2
):
    """
    Merges strain (and optionally phage) feature tables with a phenotype matrix.

    Args:
        strain_features (str): Path to the strain feature table.
        phenotype_matrix (str): Path to the phenotype matrix.
        output_dir (str): Directory to save the merged feature table.
        sample_column (str): The column name used as a sample identifier (default is 'strain').
        phage_features (str or None): Path to the phage feature table. If None, only strain features will be merged.
        remove_suffix (bool): Whether to remove suffix from the genome names.
        output_file (str or None): Optional output file name prefix.

    Returns:
        str: Path to the merged feature table.
    """
    logging.info('Starting to merge feature tables')

    # Helper function to read CSV with optional renaming
    def read_csv_with_check(filepath, rename_col=None, new_col=None):
        try:
            logging.info(f'Reading file: {filepath}')
            df = pd.read_csv(filepath)
            if rename_col and rename_col in df.columns:
                df.rename(columns={rename_col: new_col}, inplace=True)
                logging.info(f'Renamed column {rename_col} to {new_col} in {filepath}')
            return df
        except Exception as e:
            logging.error(f'Error reading {filepath}: {e}')
            raise

    # Load strain features
    strain_features_df = read_csv_with_check(strain_features, rename_col='Genome', new_col=sample_column)
    strain_features_df = strain_features_df.astype({col: 'uint8' for col in strain_features_df.columns if col.startswith(('sc_', 'pc_'))})

    if remove_suffix:
        strain_features_df[sample_column] = strain_features_df[sample_column].str.split('.').str[0]
    strain_features_df[sample_column] = strain_features_df[sample_column].astype(str)
    
    # Check if strain features table has actual features
    strain_has_features = len(strain_features_df.columns) > 1 and len(strain_features_df) > 0
    
    # Load phenotype matrix
    phenotype_matrix_df = read_csv_with_check(phenotype_matrix)
    phenotype_matrix_df[sample_column] = phenotype_matrix_df[sample_column].astype(str)

    if phage_features:
        # If phage features are provided, merge strain, phage, and phenotype matrices
        phage_features_df = read_csv_with_check(phage_features, rename_col='Genome', new_col='phage')
        phage_features_df = phage_features_df.astype({col: 'uint8' for col in phage_features_df.columns if col.startswith(('sc_', 'pc_'))})
        
        # Check if phage features table has actual features
        phage_has_features = len(phage_features_df.columns) > 1 and len(phage_features_df) > 0

        if sample_column not in phenotype_matrix_df.columns:
            logging.error(f'The phenotype matrix does not contain the "{sample_column}" column.')
            raise KeyError(f'Missing "{sample_column}" column in phenotype matrix.')
        
        if 'phage' not in phenotype_matrix_df.columns:
            logging.error('The phenotype matrix does not contain the "phage" column.')
            raise KeyError('Missing "phage" column in phenotype matrix.')

        try:
            # Start with phenotype matrix
            feature_table = phenotype_matrix_df
            
            # Only merge strain features if they exist
            if strain_has_features:
                feature_table = feature_table.merge(strain_features_df, on=sample_column, how='inner')
            
            # Only merge phage features if they exist
            if phage_has_features:
                feature_table = feature_table.merge(phage_features_df, on='phage', how='inner')
                
        except Exception as e:
            logging.error(f'Error merging strain and phage tables with phenotype matrix: {e}')
            raise
    else:
        # If no phage features, merge only strain features with phenotype matrix
        if sample_column not in phenotype_matrix_df.columns:
            logging.error(f'The phenotype matrix does not contain the "{sample_column}" column.')
            raise KeyError(f'Missing "{sample_column}" column in phenotype matrix.')

        try:
            # Start with phenotype matrix
            feature_table = phenotype_matrix_df
            
            # Only merge strain features if they exist
            if strain_has_features:
                feature_table = feature_table.merge(strain_features_df, on=sample_column, how='inner')
                
        except Exception as e:
            logging.error(f'Error merging strain features with phenotype matrix: {e}')
            raise

    if use_feature_clustering:
        logging.info("Applying pre-processing feature clustering...")
        
        # Cluster strain features if they exist
        strain_feature_columns = [col for col in feature_table.columns if col.startswith('sc_')]
        if strain_feature_columns:
            feature_table, _ = cluster_and_filter_features(
                feature_table,
                feature_type='strain',
                sample_column=sample_column,
                use_feature_clustering=True,
                feature_cluster_method=feature_cluster_method,
                feature_n_clusters=feature_n_clusters,
                feature_min_cluster_presence=feature_min_cluster_presence,
                output_dir=output_dir
            )
        
        # Cluster phage features if they exist
        phage_feature_columns = [col for col in feature_table.columns if col.startswith('pc_')]
        if phage_feature_columns and 'phage' in feature_table.columns:
            feature_table, _ = cluster_and_filter_features(
                feature_table,
                feature_type='phage',
                sample_column='phage',
                use_feature_clustering=True,
                feature_cluster_method=feature_cluster_method,
                feature_n_clusters=feature_n_clusters,
                feature_min_cluster_presence=feature_min_cluster_presence,
                output_dir=output_dir
            )

    # Determine output filename
    output_filename = f"{output_file}_full_feature_table.csv" if output_file else "full_feature_table.csv"
    feature_table_path = os.path.join(output_dir, output_filename)

    try:
        feature_table.to_csv(feature_table_path, index=False)
        logging.info(f'Successfully saved merged feature table to {feature_table_path}')
    except Exception as e:
        logging.error(f'Error saving merged feature table: {e}')
        raise

    return feature_table_path
