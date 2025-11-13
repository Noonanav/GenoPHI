import os
import subprocess
import pandas as pd
import logging
from argparse import ArgumentParser
from Bio import SeqIO
from collections import defaultdict
import warnings
from tqdm import tqdm
from tqdm.auto import tqdm as tqdm_auto
import multiprocessing
from functools import partial

# Suppress DataFrame fragmentation warning
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from genophi.mmseqs2_clustering import create_mmseqs_database, load_strains, create_contig_to_genome_dict, select_best_hits

def load_aa_sequences_from_files(fasta_files):
    """
    Loads amino acid sequences from multiple FASTA files into a DataFrame.

    Parameters:
    fasta_files (list): List of paths to FASTA files to load.

    Returns:
    DataFrame: A DataFrame with 'protein_ID' and 'sequence' columns.
    """
    logging.info(f"Loading amino acid sequences from {len(fasta_files)} FASTA files")
    
    all_records = []
    loaded_ids = set()  # Track IDs to avoid duplicates
    
    # Process each FASTA file
    for fasta_file in tqdm(fasta_files, desc="Loading FASTA files"):
        try:
            # Load records from this file
            file_records = list(SeqIO.parse(fasta_file, 'fasta'))
            logging.info(f"Loaded {len(file_records)} records from {fasta_file}")
            
            # Add only new records to avoid duplicates
            for record in file_records:
                if record.id not in loaded_ids:
                    all_records.append(record)
                    loaded_ids.add(record.id)
        except Exception as e:
            logging.error(f"Error loading sequences from {fasta_file}: {str(e)}")
    
    # Create the DataFrame with consistent-length lists
    aa_sequences_df = pd.DataFrame({
        'protein_ID': [record.id for record in all_records],
        'sequence': [str(record.seq) for record in all_records]
    })
    
    logging.info(f"Loaded {len(aa_sequences_df)} unique sequences from all files.")
    return aa_sequences_df

def detect_and_modify_duplicates(input_dir, output_dir, suffix='faa', strains_to_process=None):
    """
    Detect and resolve duplicate protein IDs by prefixing them with strain names.

    Args:
        input_dir (str): Directory containing FASTA files.
        output_dir (str): Directory to save modified files.
        suffix (str): File suffix for FASTA files (default is 'faa').
        strains_to_process (list or None): List of strains to process; if None, processes all files.

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
            if strains_to_process and strain_name not in strains_to_process:
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
                if strains_to_process and strain_name not in strains_to_process:
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
    Reuses existing output if available to avoid unnecessary recomputation.

    Args:
        db_name (str): Path to the MMseqs2 database.
        target_db (str): Path to the target MMseqs2 database.
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        sensitivity (float): Sensitivity for assignment.
        threads (int): Number of threads to use.
        clusters_tsv (str): Path to the clusters TSV file.
        clear_tmp (bool): Whether to clear temporary files.

    Returns:
        str: Path to the best hits TSV file.
    """
    result_db = os.path.join(tmp_dir, "result_db")
    assigned_tsv = os.path.join(output_dir, "assigned_clusters.tsv")
    best_hits_tsv = os.path.join(output_dir, "best_hits.tsv")
    
    # Check if best hits TSV already exists
    if os.path.exists(best_hits_tsv):
        logging.info(f"Found existing best hits file: {best_hits_tsv}. Reusing it.")
        return best_hits_tsv
        
    # Check if assigned clusters TSV already exists
    if os.path.exists(assigned_tsv):
        logging.info(f"Found existing assigned clusters file: {assigned_tsv}. Reusing it.")
    else:
        # Check if result database already exists
        result_db_files_exist = all(os.path.exists(f"{result_db}.{ext}") for ext in ["dbtype", "index"])
        
        if result_db_files_exist:
            logging.info(f"Found existing result database: {result_db}. Reusing it.")
        else:
            logging.info("Assigning sequences to clusters...")
            search_command = (
                f"mmseqs search {db_name} {target_db} {result_db} {tmp_dir} "
                f"-c {coverage} --min-seq-id {min_seq_id} -s {sensitivity} "
                f"--threads {threads} -v 3"
            )
            subprocess.run(search_command, shell=True, check=True)
            logging.info("Sequence assignment completed successfully.")
        
        # Create TSV from the result database
        createtsv_command = f"mmseqs createtsv {db_name} {target_db} {result_db} {assigned_tsv} --threads {threads} -v 3"
        subprocess.run(createtsv_command, shell=True, check=True)
        logging.info(f"Assigned clusters TSV saved to {assigned_tsv}")
    
    # Generate best hits file if it doesn't exist
    select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)
    logging.info(f"Best hits saved to {best_hits_tsv}")

    # Delete the intermediate clustering files
    if clear_tmp:
        logging.info("Clearing temporary files...")
        for file in os.listdir(tmp_dir):
            file_path = os.path.join(tmp_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
    
    return best_hits_tsv

def process_row_for_kmer_matching(row_tuple, protein_sequences, kmer_mapping):
    """Process a single row for kmer matching, suitable for multiprocessing."""
    try:
        # Unpack the tuple
        idx, row = row_tuple
        
        protein_id = row['Query']
        sequence = protein_sequences.get(protein_id, "")

        if protein_id not in protein_sequences:
            print(f"Sequence not found for {protein_id}")
        
        # Just check for the specific kmer in this row
        kmer = row['kmer']
        
        # Check if this kmer exists in the sequence
        matching_count = 1 if kmer in sequence else 0
        return (idx, matching_count)
    except Exception as e:
        print(f"Error in process_row: {e}")
        return (idx, 0)

def map_features_with_kmers_and_sequences(best_hits_tsv, feature_map, filtered_kmers, 
                                          genome_contig_mapping, genome_type, 
                                          fasta_files, aa_sequence_file, threshold, threads=None):
    """
    Maps features for all genomes based on cluster assignments and kmer presence in sequences.

    Args:
        best_hits_tsv (str): Path to the best hits TSV file.
        feature_map (str): Path to the feature mapping CSV file.
        filtered_kmers (str): Path to the filtered kmers CSV file.
        genome_contig_mapping (dict): Dictionary mapping contigs to genomes.
        genome_type (str): Type of genome ('strain' or 'phage').
        fasta_files (list): List of FASTA files containing protein sequences.
        aa_sequence_file (str): Path to the reference FASTA file (used if fasta_files is empty).
        threshold (float): Minimum percentage of kmers per feature that need to match.
        threads (int, optional): Number of threads to use for parallel processing.

    Returns:
        pd.DataFrame: Combined feature presence table for all genomes.
    """
    if not os.path.exists(best_hits_tsv) or not os.path.exists(filtered_kmers):
        logging.error(f"Required input files do not exist: {best_hits_tsv}, {filtered_kmers}")
        return None

    try:
        # Load necessary files
        best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster'])
        
        # Filter out rows with NaN clusters before merging
        best_hits_df = best_hits_df.dropna(subset=['Cluster'])
        logging.info(f"After filtering NaN clusters: {len(best_hits_df)} rows remain")
        
        feature_mapping = pd.read_csv(feature_map)
        kmer_mapping = pd.read_csv(filtered_kmers)

        # Print column information for debugging
        logging.info(f"Feature mapping columns: {feature_mapping.columns}")
        logging.info(f"Kmer mapping columns: {kmer_mapping.columns}")
        logging.info(f"Kmer mapping unique features: {kmer_mapping['Feature'].nunique()}")
        logging.info(f"Kmer mapping unique kmers: {kmer_mapping['kmer'].nunique()}")
        logging.info(f"Kmer mapping:")
        logging.info(kmer_mapping.head(10))
        logging.info(f"Best hits columns: {best_hits_df.columns}")
        logging.info(kmer_mapping.head(10))

        # Convert data types explicitly
        best_hits_df['Cluster'] = best_hits_df['Cluster'].astype(str)
        feature_mapping['Cluster_Label'] = feature_mapping['Cluster_Label'].astype(str)
        kmer_mapping['protein_family'] = kmer_mapping['protein_family'].astype(str)

        logging.info(f"Unique protein families in kmer mapping: {kmer_mapping['protein_family'].nunique()}")
        logging.info(f"Unique kmers in kmer mapping: {kmer_mapping['kmer'].nunique()}")

        # First try to load sequences from fasta_files, then fall back to aa_sequence_file
        if fasta_files:
            # Load sequences from all FASTA files in fasta_files
            aa_sequences_df = load_aa_sequences_from_files(fasta_files)
            logging.info(f"Loaded sequences from {len(fasta_files)} FASTA files")
        else:
            # Fall back to loading from the reference file if fasta_files is empty
            aa_sequences_df = load_aa_sequences(aa_sequence_file)
            logging.info(f"Loaded sequences from reference file: {aa_sequence_file}")

        print('AA sequences:')
        print(aa_sequences_df.head(10))
        protein_sequences = aa_sequences_df.set_index('protein_ID')['sequence'].to_dict()
        print(f"Loaded {len(protein_sequences)} protein sequences")

        logging.info(f"Mapping features with kmers for {genome_type}s with threshold: {threshold}...")
        
        try:
            genome_contig_mapping_df = pd.DataFrame(list(genome_contig_mapping.items()), columns=['contig_id', genome_type])
            print('Genome contig mapping:')
            print(genome_contig_mapping_df.head(10))
            merged_df = best_hits_df.merge(genome_contig_mapping_df, left_on='Query', right_on='contig_id')
            # print('Merged data:')
            # print(merged_df.head(10))
            logging.info(f"First merge successful: {len(merged_df)} rows")
        except Exception as e:
            logging.error(f"Error in first merge: {e}")
            return None
            
        try:
            merged_df = merged_df.merge(kmer_mapping, left_on='Cluster', right_on='protein_family', how='inner')
            logging.info(f"Unique protein families in merged_df: {merged_df['protein_family'].nunique()}")
            logging.info(f"Unique kmers in merged_df: {merged_df['kmer'].nunique()}")
            # print('Merged data:')
            # print(merged_df.head(10))
            logging.info(f"Second merge successful: {len(merged_df)} rows")
        except Exception as e:
            logging.error(f"Error in second merge: {e}")
            return None

        print('Merged data:')
        print(merged_df.head(10))

        # Starting the post-merge operations with detailed debugging
        logging.info("Starting kmer matching calculation...")
        
        # Debug: Check if 'kmer' column exists in merged_df
        if 'kmer' not in merged_df.columns:
            logging.error("'kmer' column not found in merged_df after the third merge!")
            logging.info(f"Available columns: {merged_df.columns.tolist()}")
            return None
        
        try:
            # Determine the number of threads to use
            if threads is None:
                num_cores = 1
            else:
                num_cores = threads
                
            total_rows = len(merged_df)
            
            # If there are very few rows or only one core, use the non-parallel approach
            if total_rows < 1000 or num_cores == 1:
                logging.info(f"Using non-parallel approach for {total_rows} rows")
                # Add the tqdm progress bar to pandas
                tqdm_auto.pandas(desc="Matching kmers")
                # Use progress_apply instead of apply
                merged_df['Matching_Kmers'] = merged_df.progress_apply(
                    lambda row: process_row_for_kmer_matching((row.name, row), protein_sequences, kmer_mapping)[1], 
                    axis=1
                )
            else:
                logging.info(f"Using parallel approach with {num_cores} cores for {total_rows} rows")
                
                # Create a partial function with the constant arguments
                process_func = partial(
                    process_row_for_kmer_matching,  # Use the module-level function
                    protein_sequences=protein_sequences,
                    kmer_mapping=kmer_mapping
                )
                
                # Create a pool of workers
                with multiprocessing.Pool(processes=num_cores) as pool:
                    # Calculate chunk size for good parallelism
                    chunk_size = max(1, min(100, total_rows // (num_cores * 4)))
                    
                    # Prepare the data as a list of (index, row) tuples
                    row_tuples = list(merged_df.iterrows())
                    
                    # Execute the parallel processing with a progress bar
                    results = list(tqdm(
                        pool.imap(process_func, row_tuples, chunksize=chunk_size),
                        total=total_rows,
                        desc="Matching kmers"
                    ))
                
                # Process the results
                result_dict = dict(results)
                
                # Create a Series with the correct indices and assign back to DataFrame
                merged_df['Matching_Kmers'] = pd.Series(result_dict)
            
            logging.info("Successfully calculated Matching_Kmers")
            logging.info(f"Unique protein families in merged_df: {merged_df['protein_family'].nunique()}")
            logging.info(f"Unique features in merged_df: {merged_df['Feature'].nunique()}")
            logging.info(f"Unique kmers in merged_df: {merged_df['kmer'].nunique()}")

            print(merged_df.head(10))
            
            # Check for missing values before groupby
            missing_genome = merged_df[genome_type].isna().sum()
            missing_feature = merged_df['Feature'].isna().sum()
            if missing_genome > 0 or missing_feature > 0:
                logging.error(f"Found {missing_genome} missing {genome_type} values and {missing_feature} missing Feature values before groupby!")
                merged_df = merged_df.dropna(subset=[genome_type, 'Feature'])
                logging.info(f"After filtering NaN values: {len(merged_df)} rows remain")
            logging.info("Checked for missing values before groupby")
            
            # Group by genome and feature
            logging.info("Starting groupby operation...")
            try:
                kmer_counts = merged_df.groupby([genome_type, 'Feature', 'cluster', 'kmer']).agg(
                    Matching_Kmers=('Matching_Kmers', 'sum')
                ).reset_index()

                logging.info(f"Counted identified kmers")
                logging.info(f"kmer_counts:")
                logging.info(kmer_counts.head(10))

                kmer_counts['Matching_Kmers'] = [1 if x > 0 else 0 for x in kmer_counts['Matching_Kmers']]

                logging.info(f"kmer_counts:")
                logging.info(kmer_counts.head(10))


                kmer_counts = kmer_counts.groupby([genome_type, 'Feature']).agg(
                    Total_Kmers=('cluster', 'nunique'),
                    Matching_Kmers=('Matching_Kmers', 'sum')
                ).reset_index()
                
                logging.info(f"Groupby successful, resulting in {len(kmer_counts)} rows")
                
                # Apply threshold
                logging.info("Applying threshold to kmer counts...")
                kmer_counts['Kmer_Percentage'] = kmer_counts['Matching_Kmers'] / kmer_counts['Total_Kmers']
                kmer_counts['Meets_Threshold'] = kmer_counts['Kmer_Percentage'] >= threshold if threshold <= 1 else kmer_counts['Matching_Kmers'] >= threshold

                print('Kmer counts:')
                print(kmer_counts.head(10))

                # In map_features_with_kmers_and_sequences function, after calculating Kmer_Percentage
                kmer_distribution = kmer_counts['Kmer_Percentage'].describe()
                logging.info(f"Kmer percentage distribution: {kmer_distribution}")

                # Also log counts of features meeting each threshold
                for test_threshold in [0.001, 0.01, 0.1, 0.2, 0.5, 0.8]:
                    meets_count = (kmer_counts['Kmer_Percentage'] >= test_threshold).sum()
                    total_count = len(kmer_counts)
                    logging.info(f"Threshold {test_threshold}: {meets_count}/{total_count} features meet threshold ({meets_count/total_count:.2%})")
                
                # Create pivot table
                logging.info("Creating pivot table...")
                feature_presence = kmer_counts.pivot(index=genome_type, columns='Feature', values='Meets_Threshold').fillna(0).astype(int)
                logging.info(f"Pivot successful, resulting in a table with {feature_presence.shape[0]} rows and {feature_presence.shape[1]} columns")
                
                # Ensure all features are included
                all_features = feature_mapping['Feature'].unique()
                for feature in all_features:
                    if feature not in feature_presence.columns:
                        feature_presence[feature] = 0
                
                feature_presence = feature_presence.reindex(columns=all_features, fill_value=0).reset_index()
                logging.info("Successfully created feature presence table")
                
                return feature_presence
            except Exception as e:
                logging.error(f"Error in groupby/pivot operations: {e}")
                import traceback
                logging.error(traceback.format_exc())
                return None
        except Exception as e:
            logging.error(f"Error calculating matching kmers: {e}")
            import traceback
            logging.error(traceback.format_exc())
            return None
    except Exception as e:
        logging.error(f"Error in map_features_with_kmers_and_sequences: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return None

def run_kmer_assign_features_workflow(input_dir, mmseqs_db, tmp_dir, output_dir, feature_map, filtered_kmers, aa_sequence_file, clusters_tsv, genome_type, genome_list=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4, suffix='faa', threshold=0.5, reuse_existing=True):
    """
    Process all genomes in the input directory at once, with optional list of genomes to process.
    Can reuse existing output files to avoid redoing computation.

    Args:
        input_dir (str): Directory containing genome FASTA files.
        mmseqs_db (str): Path to the existing MMseqs2 database.
        tmp_dir (str): Temporary directory for intermediate files.
        output_dir (str): Directory to save results.
        feature_map (str): Path to the feature mapping CSV file.
        filtered_kmers (str): Path to the filtered kmers CSV file.
        aa_sequence_file (str): Path to the FASTA file containing amino acid sequences.
        clusters_tsv (str): Path to the clusters TSV file.
        genome_type (str): Type of genomes ('strain' or 'phage').
        genome_list (str or None): Path to a file with strain names or None for all.
        sensitivity (float): Sensitivity for MMseqs2 search.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        threads (int): Number of threads for MMseqs2.
        suffix (str): Suffix for FASTA files.
        threshold (float): Threshold for kmer matching percentage.
        reuse_existing (bool): Whether to reuse existing output files.
    """
    os.makedirs(output_dir, exist_ok=True)
    if not os.path.islink(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)
    
    # Check if final output already exists
    combined_feature_table_path = os.path.join(output_dir, f'{genome_type}_combined_feature_table.csv')
    
    if reuse_existing and os.path.exists(combined_feature_table_path):
        logging.info(f"Found existing feature table: {combined_feature_table_path}. Reusing it.")
        return
    
    # Load strains if a list is provided
    strains_to_process = None
    if genome_list and os.path.exists(genome_list):
        strains_to_process = load_strains(genome_list, genome_type)
        if strains_to_process:
            strains_to_process = [str(s) for s in strains_to_process]
            logging.info(f"Loaded {len(strains_to_process)} strains from {genome_list}")

    # Use detect_and_modify_duplicates with a check for duplicates
    dir_to_use = detect_and_modify_duplicates(input_dir, output_dir, suffix, strains_to_process)

    # Check for existing combined database
    combined_db = os.path.join(tmp_dir, "combined_db")
    combined_db_exists = all(os.path.exists(f"{combined_db}.{ext}") for ext in ["dbtype", "index"])
    
    if reuse_existing and combined_db_exists:
        logging.info(f"Found existing combined database: {combined_db}. Reusing it.")
        if os.path.exists(os.path.join(dir_to_use, 'strain')):
            dir_to_use = os.path.join(dir_to_use, 'strain')
            logging.info(f"Using modified FASTA files from {dir_to_use} due to duplicate protein IDs.")
        # We need to reconstruct fasta_files list to match what create_mmseqs_database would return
        fasta_files = []
        
        logging.info(f"Searching for FASTA files with suffix '{suffix}' in directory {dir_to_use}")
        for file_name in os.listdir(dir_to_use):
            if file_name.endswith(suffix):
                strain_name = file_name.replace(f".{suffix}", "")
                if strains_to_process is None or strain_name in strains_to_process:
                    fasta_files.append(os.path.join(dir_to_use, file_name))
        
        logging.info(f"Found {len(fasta_files)} FASTA files for processing")
        
        # Debug: print some of the FASTA files to verify paths
        if fasta_files:
            sample_files = fasta_files[:3] if len(fasta_files) > 3 else fasta_files
            logging.info(f"Sample FASTA files: {sample_files}")
    else:
        logging.info(f"Creating a combined MMseqs2 database for all {genome_type}s...")
        fasta_files = create_mmseqs_database(dir_to_use, combined_db, suffix, 'directory', strains_to_process, threads)

    if not fasta_files:
        logging.error("No FASTA files found for processing.")
        return

    # Check for existing best hits file
    best_hits_tsv = os.path.join(output_dir, 'best_hits.tsv')
    if not os.path.exists(best_hits_tsv):
        best_hits_tsv = os.path.join(tmp_dir, 'best_hits.tsv')
        
    if reuse_existing and os.path.exists(best_hits_tsv):
        logging.info(f"Found existing best hits file: {best_hits_tsv}. Reusing it.")
    else:
        logging.info(f"Assigning {genome_type} sequences to existing clusters...")
        result_db = os.path.join(tmp_dir, "result_db")
        best_hits_tsv = assign_sequences_to_clusters(
            combined_db, mmseqs_db, output_dir, tmp_dir, 
            coverage, min_seq_id, sensitivity, threads, 
            clusters_tsv, clear_tmp=False
        )

    # Create genome-to-contig mapping for feature assignment
    logging.info("Creating contig to genome dictionary...")
    try:
        genome_contig_mapping, _ = create_contig_to_genome_dict(fasta_files, 'directory')
        logging.info(f"Created genome contig mapping with {len(genome_contig_mapping)} entries")
    except Exception as e:
        logging.error(f"Error creating contig to genome dictionary: {str(e)}")
        # Add more detailed error information
        import traceback
        logging.error(traceback.format_exc())
        return

    logging.info(f"Generating feature tables for all {genome_type}s...")
    feature_presence = map_features_with_kmers_and_sequences(
        best_hits_tsv, feature_map, filtered_kmers, 
        genome_contig_mapping, genome_type, 
        fasta_files, aa_sequence_file, threshold, threads
    )

    if feature_presence is not None:
        feature_presence.to_csv(combined_feature_table_path, index=False)
        logging.info(f"Combined feature table saved to {combined_feature_table_path}")
    else:
        logging.error(f"Failed to generate feature presence table for {genome_type}s.")

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
    parser.add_argument('--filtered_kmers', type=str, required=True, help="Path to the filtered kmers CSV file.")
    parser.add_argument('--aa_sequence_file', type=str, required=True, help="Path to the FASTA file containing amino acid sequences.")
    parser.add_argument('--threshold', type=float, default=0.001, help="Threshold for kmer matching percentage.")
    parser.add_argument('--reuse_existing', action='store_true', help="Reuse existing output files if available.")

    args = parser.parse_args()

    run_kmer_assign_features_workflow(
        input_dir=args.input_dir,
        mmseqs_db=args.mmseqs_db,
        tmp_dir=args.tmp_dir,
        output_dir=args.output_dir,
        feature_map=args.feature_map,
        filtered_kmers=args.filtered_kmers,
        aa_sequence_file=args.aa_sequence_file,
        clusters_tsv=args.clusters_tsv,
        genome_type=args.genome_type,
        genome_list=args.genome_list,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        suffix=args.suffix,
        threshold=args.threshold,
        reuse_existing=args.reuse_existing
    )
    
if __name__ == "__main__":
    main()