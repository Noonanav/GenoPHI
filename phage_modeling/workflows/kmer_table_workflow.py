import argparse
import logging
import os
import time
import psutil
import pandas as pd
import numpy as np
from tqdm import tqdm
import gc
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from phage_modeling.mmseqs2_clustering import merge_feature_tables
from phage_modeling.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table

# Set up logging
def setup_logging(output_dir, log_filename="kmer_modeling_workflow.log"):
    """
    Set up logging to both console and file if logging is not already configured.

    Args:
        output_dir (str): Directory where the log file will be saved.
        log_filename (str): Name of the log file. Default is "workflow.log".
    """
    if not logging.getLogger().hasHandlers():
        os.makedirs(output_dir, exist_ok=True)
        log_file = os.path.join(output_dir, log_filename)

        # Configure root logger
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file, mode='w'),  # Overwrite log file
                logging.StreamHandler()
            ]
        )
        logging.info("Logging initialized. Logs will be written to: %s", log_file)
    else:
        logging.info("Logging is already configured by the calling workflow.")

# Generate a report
def write_report(output_dir, start_time, end_time, max_ram_usage, avg_cpu_usage, peak_cpu_usage, input_genomes, select_kmers, features):
    report_file = os.path.join(output_dir, "workflow_report.txt")
    with open(report_file, "w") as report:
        report.write("Workflow Report\n")
        report.write("=" * 40 + "\n")
        report.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        report.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        report.write(f"Total Runtime: {end_time - start_time:.2f} seconds\n")
        report.write(f"Max RAM Usage: {max_ram_usage / (1024 ** 3):.2f} GB\n")
        report.write(f"Average CPU Usage: {avg_cpu_usage:.2f}%\n")
        report.write(f"Peak CPU Usage: {peak_cpu_usage:.2f}%\n")
        report.write(f"Input Genomes: {input_genomes}\n")
        report.write(f"K-mers Identified: {select_kmers}\n")
        report.write(f"Features in Final Table: {features}\n")
    logging.info(f"Report saved to: {report_file}")

def split_into_kmers(sequence, k, id):
    """Splits an amino acid sequence into k-mers of length k."""
    logging.info(f"Splitting {id} sequence into k-mers of length {k}.")
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

def generate_kmer_dict(seqrecords, k):
    """Generates a dictionary of k-mers from SeqRecord objects."""
    logging.info(f"Generating k-mer dictionary for k = {k}.")
    kmers_dict = {}
    for seqrecord in seqrecords:
        seq = str(seqrecord.seq).replace('*', '')  # Remove stop codons (*)
        id = seqrecord.id
        kmers = split_into_kmers(seq, k, id)
        for kmer in kmers:
            if kmer in kmers_dict:
                kmers_dict[kmer].append(seqrecord.id)
            else:
                kmers_dict[kmer] = [seqrecord.id]
    logging.info(f"Generated k-mer dictionary with {len(kmers_dict)} unique k-mers.")
    return kmers_dict

def generate_kmer_df(seqrecords, k):
    """Generates a DataFrame of k-mers and their corresponding protein IDs."""
    logging.info(f"Generating k-mer DataFrame for k = {k}.")
    kmers_dict = generate_kmer_dict(seqrecords, k)
    kmer_df = pd.DataFrame({'protein_ID': list(kmers_dict.values())}, index=kmers_dict.keys())
    kmer_df = kmer_df.reset_index().rename(columns={'index': 'kmer'})
    kmer_df = kmer_df.explode('protein_ID')  # Expand protein IDs to separate rows
    logging.info(f"Generated k-mer DataFrame with {len(kmer_df)} rows.")
    return kmer_df

def construct_feature_table(fasta_file, protein_csv, k, id_col, one_gene, output_dir, output_name, k_range=False, ignore_families=False):
    """Constructs a feature table from k-mers and merges it with protein feature data."""
    logging.info(f"Starting feature table construction for {fasta_file} with k = {k} (range: {k_range}).")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load sequences
    seqrecords = list(SeqIO.parse(fasta_file, 'fasta'))
    if not seqrecords:
        logging.error(f"No sequences found in {fasta_file}.")
        return

    # Generate k-mer data
    full_kmer_df = pd.DataFrame()
    if k_range:
        for k_len in range(3, k + 1):
            logging.info(f"Generating k-mers for k = {k_len}.")
            kmers_df_temp = generate_kmer_df(seqrecords, k_len)
            kmers_df_temp['k'] = k_len
            full_kmer_df = pd.concat([full_kmer_df, kmers_df_temp], ignore_index=True)
            del kmers_df_temp  # Clear temporary kmer DataFrame
            gc.collect()
    else:
        full_kmer_df = generate_kmer_df(seqrecords, k)
        full_kmer_df['k'] = k

    # Load the protein feature data
    logging.info(f"Loading protein features from {protein_csv}.")
    gene_data = pd.read_csv(protein_csv)
    
    if 'protein_ID' not in gene_data.columns:
        logging.error(f"'protein_ID' column not found in {protein_csv}.")
        raise KeyError(f"'protein_ID' column not found in {protein_csv}")

    # Merge k-mer data with gene features
    logging.info("Merging k-mer data with protein feature data.")
    full_kmer_df_merged = full_kmer_df.merge(gene_data, on="protein_ID", how='left')
    del full_kmer_df  # Clear full k-mer DataFrame
    gc.collect()

    # Aggregate k-mer data by feature and cluster
    logging.info("Aggregating k-mer data by feature and cluster.")
    full_kmer_df_counts = full_kmer_df_merged.copy()
    
    if ignore_families:
        logging.info("Ignoring protein families when aggregating k-mers.")

    group_cols = ['kmer'] if ignore_families else ['kmer', 'Feature', 'cluster']
    full_kmer_df_counts['counts'] = full_kmer_df_counts.groupby(group_cols)['protein_ID'].transform('nunique')
    full_kmer_df_counts['k_length'] = full_kmer_df_counts['kmer'].apply(len)

    # Filter unique k-mers based on 'one_gene' parameter
    if one_gene:
        logging.info("Including features with 1 gene.")
        full_kmer_df_counts = full_kmer_df_counts[full_kmer_df_counts['counts'] >= 1]
    else:
        logging.info("Filtering out features with only 1 gene.")
        full_kmer_df_counts = full_kmer_df_counts[full_kmer_df_counts['counts'] > 1]

    # Create kmer_id by concatenating cluster and kmer, and filter by k length
    logging.info("Creating kmer IDs and filtering by k length.")
    if ignore_families:
        full_kmer_df_counts['kmer_id'] = full_kmer_df_counts['kmer']
    else:
        full_kmer_df_counts['kmer_id'] = full_kmer_df_counts['cluster'] + '_' + full_kmer_df_counts['kmer']
    full_kmer_df_counts = full_kmer_df_counts.drop_duplicates(subset=[id_col, 'kmer_id'])

    # Construct presence-absence matrix
    logging.info("Constructing final presence-absence matrix.")
    feature_table = full_kmer_df_counts[[id_col, 'kmer_id']].copy()
    del full_kmer_df_counts  # Clear counts DataFrame
    gc.collect()

    feature_table['presence'] = 1
    feature_table = feature_table.pivot_table(index=id_col, columns='kmer_id', values='presence', fill_value=0)

    # Convert matrix to binary (1 for presence, 0 for absence)
    feature_table = feature_table.applymap(lambda x: 1 if x > 0 else 0).reset_index()

    # Save the presence-absence matrix as a CSV
    feature_table_output = os.path.join(output_dir, f"{output_name}_feature_table.csv")
    feature_table.to_csv(feature_table_output, index=False)
    logging.info(f"Presence-absence matrix saved to {feature_table_output}")

    del feature_table  # Clear the feature table
    gc.collect()

    return feature_table_output

def get_genome_assignments_tables(presence_absence, genome_column_name, output_dir, prefix=None):
    """Generates genome assignments from the presence-absence matrix."""
    logging.info("Getting genome assignments...")
    presence_absence.rename(columns={'Genome': genome_column_name}, inplace=True)
    genome_assignments = presence_absence.melt(id_vars=genome_column_name, var_name="Cluster_Label", value_name="Presence")
    genome_assignments = genome_assignments[genome_assignments['Presence'] == 1]

    if prefix:
        genome_assignments_output = os.path.join(output_dir, f"{prefix}_genome_assignments.csv")
    else:
        genome_assignments_output = os.path.join(output_dir, "genome_assignments.csv")
    genome_assignments.to_csv(genome_assignments_output, index=False)
    logging.info(f"Genome assignments saved to {genome_assignments_output}")

    return genome_assignments.drop(columns=["Presence"])

def feature_selection_optimized(presence_absence, source, genome_column_name, output_dir=None, prefix=None):
    """
    Optimizes feature selection by identifying perfect co-occurrence of features using hashing for performance.

    Args:
        presence_absence (DataFrame): The presence-absence matrix.
        source (str): A prefix for naming the selected features.
        genome_column_name (str): The column name that contains genome information.
        output_dir (str, optional): Directory to save the selected features CSV.
        prefix (str, optional): Prefix for the output filename.

    Returns:
        DataFrame: Optimized feature selection results.
    """
    logging.info("Optimizing feature selection using hashing...")

    # Set index using the genome_column_name
    presence_absence.set_index(genome_column_name, inplace=True)

    # Ensure binary presence-absence format
    presence_absence = presence_absence.applymap(lambda x: 1 if x > 0 else 0)

    # Hash each column to identify identical patterns
    column_hashes = presence_absence.apply(lambda col: hash(tuple(col)), axis=0)

    # Group columns by hash
    hash_to_columns = {}
    for col, col_hash in column_hashes.items():
        hash_to_columns.setdefault(col_hash, []).append(col)

    # Create unique clusters
    unique_clusters = list(hash_to_columns.values())

    # Create the selected features DataFrame
    data = [
        (f"{source[0]}c_{idx}", cluster)
        for idx, cluster_group in enumerate(unique_clusters)
        for cluster in cluster_group
    ]
    selected_features = pd.DataFrame(data, columns=["Feature", "Cluster_Label"])

    # Optionally save to CSV if output_dir is provided
    if output_dir:
        selected_features_output = os.path.join(output_dir, f"{prefix}_selected_features.csv" if prefix else "selected_features.csv")
        selected_features.to_csv(selected_features_output, index=False)
        logging.info(f"Selected features saved to {selected_features_output}")

    return selected_features

def feature_assignment(genome_assignments, selected_features, genome_column_name, output_dir, prefix=None, all_genomes=None):
    """
    Assigns features to genomes based on the selected features and ensures all genomes are included in the output,
    even those without predictive features.

    Args:
        genome_assignments (DataFrame): DataFrame of genome assignments.
        selected_features (DataFrame): DataFrame of selected features.
        genome_column_name (str): Column name representing genome IDs.
        output_dir (str): Directory to save outputs.
        prefix (str, optional): Prefix for output files.
        all_genomes (list, optional): List of all genomes to ensure they appear in the output, even if missing features.

    Returns:
        assignment_df (DataFrame): DataFrame of assigned features to genomes.
        feature_table (DataFrame): Final feature table with presence-absence for each genome.
        final_feature_table_output (str): Path to the saved final feature table.
    """
    logging.info("Assigning features to genomes...")
    assignment_df = genome_assignments.merge(selected_features, on="Cluster_Label", how="inner")
    assignment_df = assignment_df.drop(columns=["Cluster_Label"]).drop_duplicates()

    # Create feature table in wide format
    feature_table = assignment_df.pivot_table(index=genome_column_name, columns="Feature", aggfunc=lambda x: 1, fill_value=0).reset_index()

    # Ensure all genomes are included with zeros for missing features
    if all_genomes:
        missing_genomes = set(all_genomes) - set(feature_table[genome_column_name])
        if missing_genomes:
            logging.info(f"Adding {len(missing_genomes)} missing genomes with zero values.")
            missing_df = pd.DataFrame({genome_column_name: list(missing_genomes)})
            for feature in feature_table.columns:
                if feature != genome_column_name:
                    missing_df[feature] = 0
            feature_table = pd.concat([feature_table, missing_df], ignore_index=True)

    # Save the feature assignment and final feature table
    feature_assignment_output = os.path.join(output_dir, "feature_assignment.csv")
    assignment_df.to_csv(feature_assignment_output, index=False)
    logging.info(f"Feature assignment saved to {feature_assignment_output}")

    if prefix:
        final_feature_table_output = os.path.join(output_dir, f"{prefix}_final_feature_table.csv")
    else:
        final_feature_table_output = os.path.join(output_dir, "final_feature_table.csv")
    feature_table.to_csv(final_feature_table_output, index=False)
    logging.info(f"Final feature table saved to {final_feature_table_output}")

    return assignment_df, feature_table, final_feature_table_output

def load_genome_list(file_path, genome_column):
    """Load a list of genomes from a file, ensuring deduplication."""
    if file_path and os.path.exists(file_path):
        genome_df = pd.read_csv(file_path)
        return list(genome_df[genome_column].unique())
    else:
        logging.warning(f"Genome list file {file_path} not found or not provided.")
        return None

def run_kmer_table_workflow(strain_fasta, protein_csv, k, id_col, one_gene, output_dir, k_range=False,
                            phenotype_matrix=None, phage_fasta=None, protein_csv_phage=None, remove_suffix=False,
                            sample_column='strain', phenotype_column='interaction', modeling=False, filter_type='strain',
                            num_features=100, num_runs_fs=10, num_runs_modeling=20, method='rfe', strain_list=None,
                            phage_list=None, threads=4, task_type='classification', max_features='none', ignore_families=False, 
                            max_ram=8, use_shap=False):
    """
    Executes a full workflow for k-mer-based feature table construction, including strain and phage clustering,
    feature selection, phenotype merging, and optional modeling.

    Args:
        strain_fasta (str): Path to the FASTA file containing strain amino acid sequences.
        protein_csv (str): Path to the CSV file containing protein feature data for strains.
        k (int): Length of k-mers to generate from amino acid sequences.
        id_col (str): Column name in the data to use as the genome identifier.
        one_gene (bool): If True, includes features with only one gene; if False, filters out single-gene features.
        output_dir (str): Directory to save all output files, including intermediate and final feature tables.
        k_range (bool, optional): If True, generates k-mers over a range of lengths from 3 to k. Default is False.
        phenotype_matrix (str, optional): Path to the phenotype matrix CSV for merging with feature tables. Default is None.
        phage_fasta (str, optional): Path to the FASTA file containing phage amino acid sequences. Default is None.
        protein_csv_phage (str, optional): Path to the CSV file containing phage protein feature data. If None, uses `protein_csv`.
        remove_suffix (bool, optional): If True, removes suffixes from genome names when merging with phenotype matrix. Default is False.
        sample_column (str, optional): Column name for sample identifiers in the merged phenotype matrix. Default is 'strain'.
        phenotype_column (str, optional): Column name for the phenotype in the phenotype matrix, used in modeling. Default is 'interaction'.
        modeling (bool, optional): If True, runs the modeling workflow after feature table generation. Default is False.
        filter_type (str, optional): Type of feature filtering to apply; options include 'strain', 'phage', or 'none'. Default is 'strain'.
        num_features (int, optional): Number of features to select during feature selection. Default is 100.
        num_runs_fs (int, optional): Number of iterations for feature selection. Default is 10.
        num_runs_modeling (int, optional): Number of modeling iterations per feature table. Default is 20.
        method (str, optional): Feature selection method, such as 'rfe' or 'lasso'. Default is 'rfe'.
        strain_list (str, optional): Path to a list of strains to include in the analysis. Default is None.
        phage_list (str, optional): Path to a list of phages to include in the analysis. Default is None.
        ignore_families (bool, optional): If True, ignores protein families when defining k-mer features. Default is False.
        threads (int, optional): Number of threads to use for parallel processing. Default is 4.
        max_ram (float, optional): Maximum allowable RAM usage in GB for feature selection.

    Returns:
        None. Saves the final feature tables and optional modeling results to `output_dir`.
    """
    os.makedirs(output_dir, exist_ok=True)
    setup_logging(output_dir)

    # Track time, memory, and CPU usage
    start_time = time.time()
    ram_monitor = psutil.Process()
    max_ram_usage = 0
    cpu_percentages = []  # Track CPU usage over time

    # Initialize metrics for the report
    input_genomes = select_kmers = features = 0

    try:
        logging.info("Running the full k-mer feature table workflow...")
        feature_output_dir = os.path.join(output_dir, "feature_tables")
        os.makedirs(feature_output_dir, exist_ok=True)

        # Monitor CPU usage in a separate thread
        def monitor_cpu():
            while True:
                cpu_percentages.append(psutil.cpu_percent(interval=1))
                if stop_monitoring:
                    break

        stop_monitoring = False
        from threading import Thread
        cpu_monitor_thread = Thread(target=monitor_cpu)
        cpu_monitor_thread.start()

        # Step 1: Construct strain feature table
        strain_feature_table_path = construct_feature_table(strain_fasta, protein_csv, k, id_col, one_gene,
                                                            feature_output_dir, "strain", k_range, ignore_families=ignore_families)
        strain_feature_table = pd.read_csv(strain_feature_table_path)
        input_genomes += strain_feature_table[id_col].nunique()

        # Step 2: Get genome assignments for strain
        genome_assignments = get_genome_assignments_tables(strain_feature_table, id_col, feature_output_dir)

        # Step 3: Optimize feature selection for strain
        selected_features = feature_selection_optimized(
            strain_feature_table, "selected", id_col, feature_output_dir
        )

        del strain_feature_table  # Clear strain feature table
        gc.collect()

        # Step 4: Assign features to genomes
        assignment_df, final_feature_table, final_feature_table_output = feature_assignment(
            genome_assignments, selected_features, id_col, feature_output_dir
        )

        select_kmers += len(final_feature_table.columns) - 1  # Exclude genome ID column

        del selected_features, genome_assignments, assignment_df, final_feature_table  # Clear selected features
        gc.collect()

        # Step 5: Optionally, construct phage feature table
        phage_final_feature_table_output = None
        if phage_fasta:
            phage_protein_csv = protein_csv_phage if protein_csv_phage else protein_csv
            phage_feature_table_path = construct_feature_table(phage_fasta, phage_protein_csv, k, 'phage', one_gene,
                                                               feature_output_dir, "phage", k_range, ignore_families=ignore_families)
            phage_feature_table = pd.read_csv(phage_feature_table_path)

            phage_genome_assignments = get_genome_assignments_tables(phage_feature_table, 'phage', feature_output_dir, prefix='phage')

            phage_selected_features = feature_selection_optimized(
                phage_feature_table, "phage_selected", 'phage', feature_output_dir, prefix='phage'
            )

            del phage_feature_table
            gc.collect()

            phage_assignment_df, phage_final_feature_table, phage_final_feature_table_output = feature_assignment(
                phage_genome_assignments, phage_selected_features, 'phage', feature_output_dir, prefix='phage'
            )

            del phage_selected_features, phage_genome_assignments, phage_assignment_df, phage_final_feature_table  # Clear phage selected features
            gc.collect()

        # Step 6: Merge feature tables
        if phenotype_matrix:
            merged_table_path = merge_feature_tables(
                strain_features=final_feature_table_output,
                phenotype_matrix=phenotype_matrix,
                output_dir=output_dir,
                sample_column=sample_column,
                phage_features=phage_final_feature_table_output,
                remove_suffix=remove_suffix
            )
            logging.info(f"Merged feature table saved to: {merged_table_path}")
            final_df = pd.read_csv(merged_table_path)
            features = len(final_df.columns) - 1
            num_rows = len(final_df)
            if num_rows < 500:
                num_features = 50
            elif num_rows < 2000:
                num_features = 100
            else:
                num_features = int(num_rows * 0.05)

            max_features = num_features
        else:
            features = len(final_feature_table.columns) - 1

        # Optional: Run modeling
        if modeling:
            modeling_output_dir = os.path.join(output_dir, "modeling")
            os.makedirs(modeling_output_dir, exist_ok=True)
            run_modeling_workflow_from_feature_table(
                full_feature_table=merged_table_path,
                output_dir=modeling_output_dir,
                threads=threads,
                num_features=num_features,
                filter_type=filter_type,
                num_runs_fs=num_runs_fs,
                num_runs_modeling=num_runs_modeling,
                sample_column=sample_column,
                phenotype_column=phenotype_column,
                method=method,
                task_type=task_type,
                binary_data=True,
                max_features=max_features,
                max_ram=max_ram,
                use_shap=use_shap
            )
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
    finally:
        # Stop CPU monitoring and compute CPU usage stats
        stop_monitoring = True
        cpu_monitor_thread.join()

        avg_cpu_usage = sum(cpu_percentages) / len(cpu_percentages) if cpu_percentages else 0
        peak_cpu_usage = max(cpu_percentages) if cpu_percentages else 0

        # Track end time and max RAM usage
        end_time = time.time()
        max_ram_usage = max(max_ram_usage, ram_monitor.memory_info().rss)

        # Write the workflow report
        write_report(output_dir, start_time, end_time, max_ram_usage, avg_cpu_usage, peak_cpu_usage, input_genomes, select_kmers, features)

# Command-line interface

def main():
    """
    Command-line interface for generating k-mer feature tables and running the full workflow.
    """
    parser = argparse.ArgumentParser(description='Generates k-mer feature tables and final presence-absence matrix from AA sequences.')

    # Input data
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-i', '--strain_fasta', required=True, help="Path to the FASTA file containing strain AA sequences.")
    input_group.add_argument('-ip', '--phage_fasta', default=None, help="Path to the FASTA file containing phage AA sequences.")
    input_group.add_argument('-p', '--protein_csv', required=True, help="Path to the CSV file containing protein feature data for strain.")
    input_group.add_argument('--protein_csv_phage', default=None, help="Path to the CSV file containing protein feature data for phage.")

    # Optional input parameters
    optional_input_group = parser.add_argument_group('Optional input parameters')
    optional_input_group.add_argument('--k', type=int, required=True, help="k-mer length.")
    optional_input_group.add_argument('--id_col', default="strain", help="Column name for genome ID.")
    optional_input_group.add_argument('--one_gene', action='store_true', help="Include features with 1 gene.")
    optional_input_group.add_argument('--k_range', action='store_true', help="Generate range of k-mer lengths from 3 to k.")
    optional_input_group.add_argument('--phenotype_matrix', default=None, help="Path to the phenotype matrix CSV for merging. (optional)")
    optional_input_group.add_argument('--remove_suffix', action='store_true', help="Remove suffix from genome names when merging.")
    optional_input_group.add_argument('--sample_column', default='strain', help="Sample identifier column name (default: 'strain').")
    optional_input_group.add_argument('--phenotype_column', default='interaction', help="Phenotype column name in the phenotype matrix.")
    optional_input_group.add_argument('--strain_list', default=None, help="Full list of strains to include in the analysis.")
    optional_input_group.add_argument('--phage_list', default=None, help="Full list of phages to include in the analysis.")
    optional_input_group.add_argument('--ignore_families', action='store_true', help="Ignore protein families when defining k-mer features")

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('-o', '--output_dir', required=True, help="Directory to save all output files.")

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--modeling', action='store_true', help="Run modeling workflow after feature table generation.")
    fs_modeling_group.add_argument('--filter_type', default='strain', help="Type of feature filtering to apply (default: 'strain').")
    fs_modeling_group.add_argument('--num_features', type=int, default=100, help="Number of features to select for modeling (default: 100).")
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help="Number of runs for feature selection (default: 10).")
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=20, help="Number of runs for modeling (default: 20).")
    fs_modeling_group.add_argument('--method', default='rfe', help="Feature selection method to use (default: 'rfe').")
    fs_modeling_group.add_argument('--task_type', default='classification', choices=['classification', 'regression'], help="Specify 'classification' or 'regression' task.")
    fs_modeling_group.add_argument('--max_features', default='none', help='Maximum number of features to include in the feature tables.')
    fs_modeling_group.add_argument('--use_shap', action='store_true', help='If True, calculates SHAP feature importance values.')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help="Number of threads to use (default: 4).")
    general_group.add_argument('--max_ram', type=float, default=8, help="Maximum RAM usage in GB for feature selection (default: 8).")

    args = parser.parse_args()

    # Call the full workflow function
    run_kmer_table_workflow(
        strain_fasta=args.strain_fasta,
        protein_csv=args.protein_csv,
        k=args.k,
        id_col=args.id_col,
        one_gene=args.one_gene,
        output_dir=args.output_dir,
        k_range=args.k_range,
        phenotype_matrix=args.phenotype_matrix,
        phage_fasta=args.phage_fasta,
        protein_csv_phage=args.protein_csv_phage,
        remove_suffix=args.remove_suffix,
        sample_column=args.sample_column,
        modeling=args.modeling,
        phenotype_column=args.phenotype_column,
        filter_type=args.filter_type,
        num_features=args.num_features,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        method=args.method,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        threads=args.threads,
        task_type=args.task_type,
        max_features=args.max_features,
        ignore_families=args.ignore_families,
        max_ram=args.max_ram,
        use_shap=args.use_shap
    )

if __name__ == "__main__":
    main()