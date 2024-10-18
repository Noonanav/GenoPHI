import os
import pandas as pd
import subprocess
import logging
from Bio import SeqIO

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
    if input_type == 'directory':
        logging.info(f"Searching for FASTA files with suffix '{suffix}' in directory.")
        for fasta in os.listdir(input_path):
            if fasta.endswith(suffix):
                strain_name = os.path.splitext(fasta)[0]
                if strains is None or strain_name in strains:
                    fasta_files.append(os.path.join(input_path, fasta))
    else:
        logging.info(f"Using input file: {input_path}")
        fasta_files.append(input_path)
    
    if not fasta_files:
        logging.error("No matching FASTA files found. Exiting.")
        raise FileNotFoundError("No FASTA files to process.")
    
    fasta_files_str = " ".join(fasta_files)
    createdb_command = f"mmseqs createdb {fasta_files_str} {db_name} -v 3"
    logging.info(f"Running command: {createdb_command}")
    subprocess.run(createdb_command, shell=True, check=True)
    logging.info("MMseqs2 database created successfully.")
    
    return fasta_files

def create_contig_to_genome_dict(fasta_files, input_type):
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
            genome_name = os.path.basename(fasta).split(".")[0]
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

def assign_sequences_to_clusters(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv):
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
    assigned_df = pd.read_csv(assigned_tsv, sep='\t', header=None, names=[
        'Query', 'Target', 'Score', 'SeqIdentity', 'E-value', 'qStartPos', 'qEndPos', 'qLen', 'tStartPos', 'tEndPos', 'tLen'
    ])

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
    Filters the presence-absence matrix based on a selected strain list.

    Args:
        presence_absence (DataFrame): The presence-absence matrix.
        select (str): Path to the file containing selected strains.
        select_column (str): Column in the selection file that contains strain names.

    Returns:
        DataFrame: Filtered presence-absence matrix.
    """
    logging.info("Filtering presence-absence table...")
    select_df = pd.read_csv(select)
    select_list = list(set(select_df[select_column].tolist()))
    presence_absence = presence_absence[presence_absence['Genome'].isin(select_list)]
    
    all_cols = len(presence_absence.columns)
    presence_absence = presence_absence.loc[:, (presence_absence != 0).any(axis=0)]
    logging.info(f"Removed {all_cols - len(presence_absence.columns) - 1} columns with all zeros")
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
    Optimizes feature selection by identifying perfect co-occurrence of features.

    Args:
        presence_absence (DataFrame): The presence-absence matrix.
        source (str): A prefix for naming the selected features.
        genome_column_name (str): The column name that contains genome information (e.g., 'strain' or 'phage').

    Returns:
        DataFrame: Optimized feature selection results.
    """
    logging.info("Optimizing feature selection...")

    # Set index using the genome_column_name instead of 'Genome'
    presence_absence.set_index(genome_column_name, inplace=True)
    
    boolean_matrix = presence_absence.astype(bool)
    perfect_cooccurrence = {col: set(boolean_matrix.columns[boolean_matrix.eq(boolean_matrix[col], axis=0).all()]) for col in boolean_matrix.columns}
    
    unique_clusters = []
    seen = set()
    for cluster in perfect_cooccurrence.values():
        if not cluster.intersection(seen):
            unique_clusters.append(list(cluster))
        seen.update(cluster)
    
    data = [(f"{source[0]}c_{idx}", cluster) for idx, cluster_group in enumerate(unique_clusters) for cluster in cluster_group]
    return pd.DataFrame(data, columns=["Feature", "Cluster_Label"])

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


def run_clustering_workflow(input_path, output_dir, tmp_dir="tmp", min_seq_id=0.6, coverage=0.8, sensitivity=7.5, suffix='faa', threads=4, strain_list='none', strain_column='strain', compare=False):
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
    
    Raises:
        FileNotFoundError: If no FASTA files are found in the input path.
        ValueError: If there is an issue with loading the strain list.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # Convert 'none' to None internally for processing
    strain_list_value = None if strain_list == 'none' else strain_list

    input_type = 'directory' if os.path.isdir(input_path) else 'file'
    strains = load_strains(strain_list_value, strain_column) if strain_list_value else None

    db_name = os.path.join(tmp_dir, "mmseqs_db")
    fasta_files = create_mmseqs_database(input_path, db_name, suffix, input_type, strains, threads)

    contig_to_genome, genome_list = create_contig_to_genome_dict(fasta_files, input_type)

    clusters_tsv = run_mmseqs_cluster(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads)

    best_hits_tsv = assign_sequences_to_clusters(db_name, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv)

    presence_absence_csv = os.path.join(output_dir, "presence_absence_matrix.csv")
    generate_presence_absence_matrix(best_hits_tsv, presence_absence_csv, contig_to_genome, genome_list)

    if compare:
        compare_cluster_and_search_results(clusters_tsv, best_hits_tsv, output_dir)

def run_feature_assignment(input_file, output_dir, source='bacteria', select='none', select_column='strain'):
    """
    Runs the feature assignment workflow from the presence-absence matrix.
    
    This function processes a presence-absence matrix, assigns features to genomes based on selected features, 
    and outputs both feature assignments and a feature table. It allows for optional filtering based on a strain list.
    
    Args:
        input_file (str): Path to the presence-absence matrix CSV file.
        output_dir (str): Directory to save the results (selected features, assignments, and feature table).
        source (str): Prefix for naming the selected features (e.g., 'bacteria' or 'phage').
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


def merge_feature_tables(strain_features, phenotype_matrix, output_dir, sample_column='strain', phage_features=None, remove_suffix=False, output_file=None):
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
    if remove_suffix:
        strain_features_df[sample_column] = strain_features_df[sample_column].str.split('.').str[0]
    
    # Load phenotype matrix
    phenotype_matrix_df = read_csv_with_check(phenotype_matrix)

    if phage_features:
        # If phage features are provided, merge strain, phage, and phenotype matrices
        phage_features_df = read_csv_with_check(phage_features, rename_col='Genome', new_col='phage')

        if sample_column not in phenotype_matrix_df.columns:
            logging.error(f'The phenotype matrix does not contain the "{sample_column}" column.')
            raise KeyError(f'Missing "{sample_column}" column in phenotype matrix.')
        
        if 'phage' not in phenotype_matrix_df.columns:
            logging.error('The phenotype matrix does not contain the "phage" column.')
            raise KeyError('Missing "phage" column in phenotype matrix.')

        try:
            feature_table = phenotype_matrix_df.merge(strain_features_df, on=sample_column, how='inner')
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
            feature_table = phenotype_matrix_df.merge(strain_features_df, on=sample_column, how='inner')
        except Exception as e:
            logging.error(f'Error merging strain features with phenotype matrix: {e}')
            raise

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
