import os
import pandas as pd
import argparse
import logging
import csv
import gc
import time
import psutil
from Bio import SeqIO
from collections import defaultdict
from genophi.workflows.kmer_table_workflow import run_kmer_table_workflow


def setup_logging(output_dir, log_filename="kmer_workflow.log"):
    """
    Set up logging to both console and file if logging is not already configured.

    Args:
        output_dir (str): Directory where the log file will be saved.
        log_filename (str): Name of the log file. Default is "kmer_workflow.log".
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


def write_csv_log(output_dir, data):
    """
    Writes a CSV log file with variable names and their values.

    Args:
        output_dir (str): Directory where the CSV log file will be saved.
        data (dict): Dictionary of variable names and their values.
    """
    csv_file = os.path.join(output_dir, 'kmer_workflow_report.csv')
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Variable', 'Value'])
        for key, value in data.items():
            writer.writerow([key, value])
    logging.info(f"CSV log saved to {csv_file}.")


def write_report(output_dir, start_time, end_time, ram_usage, avg_cpu_usage, max_cpu_usage, 
                 input_genomes, kmer_features, strain_feature_count, phage_feature_count=None):
    """
    Writes a detailed workflow report to a text file.

    Args:
        output_dir (str): Directory where the report file will be saved.
        start_time (float): Workflow start time.
        end_time (float): Workflow end time.
        ram_usage (int): Maximum RAM usage in bytes.
        avg_cpu_usage (float): Average CPU usage during workflow.
        max_cpu_usage (float): Maximum CPU usage during workflow.
        input_genomes (int): Number of input genomes.
        kmer_features (int): Number of k-mer features generated.
        strain_feature_count (int): Number of strain features.
        phage_feature_count (int, optional): Number of phage features.
    """
    report_file = os.path.join(output_dir, "kmer_workflow_report.txt")
    with open(report_file, "w") as report:
        report.write("K-mer Workflow Report\n")
        report.write("=" * 40 + "\n")
        report.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        report.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        report.write(f"Total Runtime: {end_time - start_time:.2f} seconds\n")
        report.write(f"Max RAM Usage: {ram_usage / (1024 ** 3):.2f} GB\n")
        report.write(f"Average CPU Usage: {avg_cpu_usage:.2f}%\n")
        report.write(f"Max CPU Usage: {max_cpu_usage:.2f}%\n")
        report.write(f"Input Genomes: {input_genomes}\n")
        report.write(f"Total strain features: {strain_feature_count}\n")
        if phage_feature_count:
            report.write(f"Total phage features: {phage_feature_count}\n")
        report.write(f"K-mer Features Generated: {kmer_features}\n")
    logging.info(f"Report saved to: {report_file}")


def detect_duplicate_ids(input_path, suffix='faa', strains_to_process=None, input_type='directory'):
    """
    Detects duplicate protein IDs across relevant .faa files in the input directory.
    
    Args:
        input_path (str): Directory containing .faa files.
        suffix (str): Suffix for the input FASTA files (default is 'faa').
        strains_to_process (list or None): List of strains to process; if None, processes all .faa files.
        input_type (str): Type of input ('directory' or 'file').
    
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
        output_dir (str): Output directory for modified files.
        suffix (str): Suffix for the input FASTA files (default is 'faa').
        strains_to_process (list or None): List of strains to process; if None, processes all .faa files.
        strain_column (str): Column name for strain identifier.
    
    Returns:
        str: Path to directory containing modified .faa files.
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


def get_full_strain_list(interaction_matrix, input_strain_dir, strain_column):
    """
    Get the intersection of strains present in both the interaction matrix and input directory.
    
    Args:
        interaction_matrix (str): Path to interaction matrix CSV file.
        input_strain_dir (str): Directory containing strain FASTA files.
        strain_column (str): Column name for strain identifiers in interaction matrix.
        
    Returns:
        list: List of strain names present in both matrix and directory.
    """
    logging.info(f'Reading interaction matrix: {interaction_matrix}')
    interaction_df = pd.read_csv(interaction_matrix)
    logging.info(f'Interaction matrix shape: {interaction_df.shape}')
    logging.info(f'Interaction matrix columns: {list(interaction_df.columns)}')
    
    if strain_column not in interaction_df.columns:
        logging.error(f'ERROR: Column \'{strain_column}\' not found in interaction matrix')
        logging.info(f'Available columns: {list(interaction_df.columns)}')
        return []
    
    strains_in_matrix = interaction_df[strain_column].unique()
    strains_in_matrix = [str(s) for s in strains_in_matrix]
    logging.info(f'Found {len(strains_in_matrix)} unique strains in interaction matrix')
    logging.info(f'First 10 strains from matrix: {list(strains_in_matrix[:10])}')
    
    logging.info(f'Reading strain directory: {input_strain_dir}')
    strain_files = [f for f in os.listdir(input_strain_dir) if f.endswith('.faa')]
    strains_in_dir = ['.'.join(f.split('.')[:-1]) for f in strain_files]
    logging.info(f'Found {len(strains_in_dir)} strain files in directory')
    logging.info(f'First 10 strains from directory: {strains_in_dir[:10]}')
    
    full_strain_list = list(set(strains_in_matrix).intersection(set(strains_in_dir)))
    logging.info(f'Intersection: {len(full_strain_list)} strains found in both matrix and directory')
    
    if len(full_strain_list) == 0:
        logging.error('ERROR: No strains found in both interaction matrix and input directory!')
        logging.error('This might be due to naming mismatches between the files and matrix.')
    
    return full_strain_list


def get_full_phage_list(interaction_matrix, input_phage_dir):
    """
    Get the intersection of phages present in both the interaction matrix and input directory.
    
    Args:
        interaction_matrix (str): Path to interaction matrix CSV file.
        input_phage_dir (str): Directory containing phage FASTA files.
        
    Returns:
        list: List of phage names present in both matrix and directory.
    """
    logging.info(f'Reading interaction matrix for phages: {interaction_matrix}')
    interaction_df = pd.read_csv(interaction_matrix)
    
    # Assume phage column is named 'phage' - adjust if different
    phage_column = 'phage'
    if phage_column not in interaction_df.columns:
        logging.warning(f'WARNING: Column \'{phage_column}\' not found in interaction matrix')
        # Try to find a phage-related column
        potential_cols = [col for col in interaction_df.columns if 'phage' in col.lower()]
        if potential_cols:
            phage_column = potential_cols[0]
            logging.info(f'Using column: {phage_column}')
        else:
            logging.info('No phage column found, using all phage files from directory')
            phage_files = [f for f in os.listdir(input_phage_dir) if f.endswith('.faa')]
            return ['.'.join(f.split('.')[:-1]) for f in phage_files]
    
    phages_in_matrix = interaction_df[phage_column].unique()
    phages_in_matrix = [str(s) for s in phages_in_matrix]
    logging.info(f'Found {len(phages_in_matrix)} unique phages in interaction matrix')
    
    logging.info(f'Reading phage directory: {input_phage_dir}')
    phage_files = [f for f in os.listdir(input_phage_dir) if f.endswith('.faa')]
    phages_in_dir = ['.'.join(f.split('.')[:-1]) for f in phage_files]
    logging.info(f'Found {len(phages_in_dir)} phage files in directory')
    
    full_phage_list = list(set(phages_in_matrix).intersection(set(phages_in_dir)))
    logging.info(f'Intersection: {len(full_phage_list)} phages found in both matrix and directory')
    
    return full_phage_list


def create_filtered_fasta(input_dir, genome_list, output_file, suffix='faa'):
    """
    Create FASTA file with only specified genomes.
    
    Args:
        input_dir (str): Directory containing FASTA files.
        genome_list (list): List of genome names to include.
        output_file (str): Path to output FASTA file.
        suffix (str): File extension for FASTA files.
    """
    logging.info(f'Creating filtered FASTA from {len(genome_list)} genomes')
    
    with open(output_file, 'w') as outfile:
        for genome in genome_list:
            genome_file = os.path.join(input_dir, f'{genome}.{suffix}')
            if os.path.exists(genome_file):
                for record in SeqIO.parse(genome_file, 'fasta'):
                    SeqIO.write(record, outfile, 'fasta')
            else:
                logging.warning(f'FASTA file not found for genome {genome}: {genome_file}')
    
    logging.info(f'Filtered FASTA created at {output_file}')


def create_filtered_protein_csv(input_csv, genome_list, output_csv, genome_column):
    """
    Create filtered protein CSV with only specified genomes.
    
    Args:
        input_csv (str): Path to input protein CSV file.
        genome_list (list): List of genome names to include.
        output_csv (str): Path to output CSV file.
        genome_column (str): Column name for genome identifiers.
    """
    logging.info('Creating filtered protein CSV')
    
    df = pd.read_csv(input_csv)
    filtered_df = df[df[genome_column].isin(genome_list)]
    filtered_df.to_csv(output_csv, index=False)
    
    logging.info(f'Filtered protein CSV created at {output_csv} with {len(filtered_df)} proteins')



def generate_filtered_protein_mapping_csv(fasta_dir, output_csv, genome_col_name, genome_list, file_extension='.faa'):
    """
    Generate protein-to-genome mapping CSV from FASTA directory, filtering by genome list.
    
    Args:
        fasta_dir (str): Directory containing FASTA files.
        output_csv (str): Path to output CSV file.
        genome_col_name (str): Name for genome column ('strain' or 'phage').
        genome_list (list): List of genome names to include.
        file_extension (str): File extension for FASTA files.
        
    Returns:
        str: Path to generated CSV file.
    """
    data = []
    
    for genome_name in genome_list:  # Only process genomes in the list
        fasta_file = f'{genome_name}{file_extension}'
        fasta_path = os.path.join(fasta_dir, fasta_file)
        
        if os.path.exists(fasta_path):
            for record in SeqIO.parse(fasta_path, 'fasta'):
                data.append({
                    'protein_ID': record.id,
                    genome_col_name: genome_name
                })
        else:
            logging.warning(f'FASTA file not found for genome {genome_name}: {fasta_path}')
    
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    logging.info(f'Generated filtered {output_csv} with {len(df)} protein mappings from {len(genome_list)} genomes')
    return output_csv


def run_kmer_workflow(
    input_strain_dir, 
    output_dir, 
    phenotype_matrix, 
    input_phage_dir=None,
    k=5,
    k_range=False,
    one_gene=False,
    suffix='faa',
    threads=4, 
    strain_list=None, 
    phage_list=None, 
    strain_column='strain', 
    phage_column='phage', 
    sample_column='strain', 
    phenotype_column='interaction', 
    num_features=100, 
    filter_type='strain', 
    num_runs_fs=10, 
    num_runs_modeling=10,
    method='rfe',
    task_type='classification', 
    max_features='none', 
    max_ram=8, 
    use_dynamic_weights=False, 
    weights_method='log10',
    use_clustering=False,
    cluster_method='hierarchical',
    n_clusters=20,
    min_cluster_size=5,
    min_samples=None,
    cluster_selection_epsilon=0.0,
    check_feature_presence=False,
    filter_by_cluster_presence=False,
    min_cluster_presence=2, 
    use_shap=False,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2,
    remove_suffix=False,
    run_predictive_analysis=False
):
    """
    Complete k-mer workflow: K-mer feature generation, feature selection, modeling, and analysis.
    
    This workflow is the k-mer equivalent of the protein family workflow, taking FASTA directories
    as input instead of requiring MMSeqs2 clustering.
    """
    
    os.makedirs(output_dir, exist_ok=True)
    setup_logging(output_dir)

    # Collect inputs for CSV log
    inputs = locals()

    # Track time and resource usage
    start_time = time.time()
    ram_monitor = psutil.Process()
    cpu_usage_points = []
    max_ram_usage = 0

    # Initialize counters for report
    input_genomes = 0
    kmer_features = 0
    strain_feature_count = 0
    phage_feature_count = None

    try:
        logging.info("=== Starting K-mer Feature Generation Workflow ===")
        
        # Validate inputs
        if not os.path.exists(input_strain_dir):
            raise FileNotFoundError(f"Strain input directory not found: {input_strain_dir}")
        if not os.path.exists(phenotype_matrix):
            raise FileNotFoundError(f"Phenotype matrix not found: {phenotype_matrix}")
        if input_phage_dir and not os.path.exists(input_phage_dir):
            raise FileNotFoundError(f"Phage input directory not found: {input_phage_dir}")

        # Step 1: Get full strain and phage lists
        logging.info("Step 1: Identifying available strains and phages...")
        
        full_strain_list = get_full_strain_list(phenotype_matrix, input_strain_dir, strain_column)
        full_phage_list = []
        if input_phage_dir:
            full_phage_list = get_full_phage_list(phenotype_matrix, input_phage_dir)

        if not full_strain_list:
            raise ValueError("No valid strains found in both interaction matrix and input directory")
        if input_phage_dir and not full_phage_list:
            raise ValueError("No valid phages found in both interaction matrix and input directory")

        # Count input genomes
        input_genomes = len(full_strain_list)
        if full_phage_list:
            input_genomes += len(full_phage_list)
        logging.info(f"Total input genomes: {input_genomes}")

        # Step 2: Handle duplicate protein IDs
        logging.info("Step 2: Checking for duplicate protein IDs...")
        
        strain_dir_to_use = input_strain_dir
        phage_dir_to_use = input_phage_dir
        
        # Check strain duplicates
        strain_duplicate_found = detect_duplicate_ids(input_strain_dir, suffix, full_strain_list, 'directory')
        if strain_duplicate_found:
            logging.info('Duplicate protein IDs found in strain directory; modifying all protein IDs.')
            strain_dir_to_use = modify_duplicate_ids(input_strain_dir, output_dir, suffix, full_strain_list, strain_column)

        # Check phage duplicates if phage directory provided
        if input_phage_dir:
            phage_duplicate_found = detect_duplicate_ids(input_phage_dir, suffix, full_phage_list, 'directory')
            if phage_duplicate_found:
                logging.info('Duplicate protein IDs found in phage directory; modifying all protein IDs.')
                phage_dir_to_use = modify_duplicate_ids(input_phage_dir, output_dir, suffix, full_phage_list, phage_column)

        # Step 3: Prepare input files for k-mer workflow
        logging.info("Step 3: Preparing input files for k-mer feature generation...")
        
        # Create filtered FASTA files and protein mapping CSVs
        strain_fasta = os.path.join(output_dir, 'strain_combined.faa')
        strain_csv = os.path.join(output_dir, 'strain_proteins.csv')
        
        create_filtered_fasta(strain_dir_to_use, full_strain_list, strain_fasta, suffix)
        generate_filtered_protein_mapping_csv(strain_dir_to_use, strain_csv, strain_column, full_strain_list, f'.{suffix}')
        
        # Prepare phage files if provided
        phage_fasta = None
        phage_csv = None
        if input_phage_dir:
            phage_fasta = os.path.join(output_dir, 'phage_combined.faa')
            phage_csv = os.path.join(output_dir, 'phage_proteins.csv')
            
            create_filtered_fasta(phage_dir_to_use, full_phage_list, phage_fasta, suffix)
            generate_filtered_protein_mapping_csv(phage_dir_to_use, phage_csv, phage_column, full_phage_list, f'.{suffix}')

        # Monitor resource usage
        max_ram_usage = max(max_ram_usage, ram_monitor.memory_info().rss)
        cpu_usage_points.append(psutil.cpu_percent(interval=None))

        # Step 4: Run k-mer table workflow with modeling
        logging.info("Step 4: Running k-mer table workflow with feature selection and modeling...")
        
        merged_feature_table = run_kmer_table_workflow(
            strain_fasta=strain_fasta,
            protein_csv=strain_csv,
            k=k,
            id_col=strain_column,
            one_gene=one_gene,
            output_dir=output_dir,
            k_range=k_range,
            phenotype_matrix=phenotype_matrix,
            phage_fasta=phage_fasta,
            protein_csv_phage=phage_csv,
            remove_suffix=remove_suffix,
            sample_column=sample_column,
            phenotype_column=phenotype_column,
            modeling=True,  # Always run modeling in this workflow
            filter_type=filter_type,  # Essential for train-test split functionality
            num_features=num_features,
            num_runs_fs=num_runs_fs,
            num_runs_modeling=num_runs_modeling,
            method=method,
            strain_list=strain_list,  # Pass through strain list parameter
            phage_list=phage_list,   # Pass through phage list parameter
            threads=threads,
            task_type=task_type,
            max_features=max_features,
            ignore_families=True,  # Ignore families for k-mer workflow
            max_ram=max_ram,
            use_shap=use_shap,
            use_clustering=use_clustering,
            cluster_method=cluster_method,
            n_clusters=n_clusters,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cluster_selection_epsilon=cluster_selection_epsilon,
            use_dynamic_weights=use_dynamic_weights,
            weights_method=weights_method,
            check_feature_presence=check_feature_presence,
            filter_by_cluster_presence=filter_by_cluster_presence,
            min_cluster_presence=min_cluster_presence,
            use_feature_clustering=use_feature_clustering,
            feature_cluster_method=feature_cluster_method,
            feature_n_clusters=feature_n_clusters,
            feature_min_cluster_presence=feature_min_cluster_presence
        )

        # Get feature counts for reporting
        if os.path.exists(merged_feature_table):
            feature_df = pd.read_csv(merged_feature_table)
            kmer_features = len(feature_df.columns) - 3  # Subtract non-feature columns
            
            # Count strain vs phage features
            strain_features = [col for col in feature_df.columns if col.startswith('sc_')]
            phage_features = [col for col in feature_df.columns if col.startswith('pc_')]
            strain_feature_count = len(strain_features)
            if phage_features:
                phage_feature_count = len(phage_features)
            
            logging.info(f"Generated {kmer_features} total k-mer features")
            logging.info(f"  - {strain_feature_count} strain features")
            if phage_feature_count:
                logging.info(f"  - {phage_feature_count} phage features")
            
            del feature_df
            gc.collect()
        else:
            logging.warning("Merged feature table not found, cannot count features")

        # Monitor resource usage
        max_ram_usage = max(max_ram_usage, ram_monitor.memory_info().rss)
        cpu_usage_points.append(psutil.cpu_percent(interval=None))

        # Step 5: Optional predictive analysis workflow
        if run_predictive_analysis:
            logging.info("Step 5: Running predictive k-mer analysis...")
            
            # Get the best performing cutoff from modeling results
            metrics_file = os.path.join(output_dir, 'modeling', 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
            if os.path.exists(metrics_file):
                performance_df = pd.read_csv(metrics_file)
                top_cutoff = performance_df.iloc[0]['cut_off'].split('_')[-1]
                
                logging.info(f"Best performing cutoff: {top_cutoff}")
                
                # Set up paths for predictive analysis
                feature_file_path = os.path.join(output_dir, 'modeling', 'feature_selection', 'filtered_feature_tables', f'select_feature_table_cutoff_{top_cutoff}.csv')
                selected_features_path = os.path.join(output_dir, 'feature_tables', 'selected_features.csv')
                model_dir = os.path.join(output_dir, 'modeling', 'modeling_results', f'cutoff_{top_cutoff}')
                
                # This step could include k-mer specific predictive analysis
                # For now, we log that the setup is complete
                logging.info(f"Predictive analysis files prepared:")
                logging.info(f"  - Feature file: {feature_file_path}")
                logging.info(f"  - Selected features: {selected_features_path}")
                logging.info(f"  - Model directory: {model_dir}")
                
            else:
                logging.warning("Model performance metrics not found, skipping predictive analysis")

        logging.info("=== K-mer Workflow Completed Successfully ===")
        
        return merged_feature_table

    except Exception as e:
        logging.error(f"An error occurred in k-mer workflow: {e}")
        raise
    finally:
        end_time = time.time()
        max_ram_usage = max(max_ram_usage, ram_monitor.memory_info().rss)
        cpu_usage_points.append(psutil.cpu_percent(interval=None))

        avg_cpu_usage = sum(cpu_usage_points) / len(cpu_usage_points) if cpu_usage_points else 0
        max_cpu_usage = max(cpu_usage_points) if cpu_usage_points else 0

        # Write report and logs
        write_report(output_dir, start_time, end_time, max_ram_usage, avg_cpu_usage, 
                    max_cpu_usage, input_genomes, kmer_features, strain_feature_count, phage_feature_count)

        # Save CSV log with all parameters and report values
        inputs.update({
            'start_time': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time)),
            'end_time': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time)),
            'runtime_seconds': end_time - start_time,
            'ram_usage_gb': max_ram_usage / (1024 ** 3),
            'avg_cpu_usage': avg_cpu_usage,
            'max_cpu_usage': max_cpu_usage,
            'input_genomes': input_genomes,
            'kmer_features': kmer_features,
            'strain_feature_count': strain_feature_count,
            'phage_feature_count': phage_feature_count
        })
        write_csv_log(output_dir, inputs)


# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Complete k-mer workflow: Feature generation, selection, modeling, and analysis.')
    
    # Input data
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-is', '--input_strain_dir', type=str, required=True, 
                            help='Directory containing strain FASTA files.')
    input_group.add_argument('-ip', '--input_phage_dir', type=str, 
                            help='Directory containing phage FASTA files (optional).')
    input_group.add_argument('-pm', '--phenotype_matrix', type=str, required=True, 
                            help='Path to the phenotype matrix.')

    # K-mer specific parameters
    kmer_group = parser.add_argument_group('K-mer parameters')
    kmer_group.add_argument('--k', type=int, default=5, help='K-mer length (default: 5).')
    kmer_group.add_argument('--k_range', action='store_true', 
                           help='Generate k-mers from length 3 to k.')
    kmer_group.add_argument('--one_gene', action='store_true',
                           help='Include features with only one gene.')

    # Optional input arguments
    optional_input_group = parser.add_argument_group('Optional input arguments')
    optional_input_group.add_argument('--suffix', type=str, default='faa', 
                                     help='File extension for FASTA files (default: faa).')
    optional_input_group.add_argument('--strain_list', type=str, 
                                     help='Path to a strain list file for filtering.')
    optional_input_group.add_argument('--phage_list', type=str, 
                                     help='Path to a phage list file for filtering.')
    optional_input_group.add_argument('--strain_column', type=str, default='strain', 
                                     help='Column name for strain identifiers (default: strain).')
    optional_input_group.add_argument('--phage_column', type=str, default='phage', 
                                     help='Column name for phage identifiers (default: phage).')
    optional_input_group.add_argument('--sample_column', type=str, default='strain', 
                                     help='Column name for sample identifiers (default: strain).')
    optional_input_group.add_argument('--phenotype_column', type=str, default='interaction', 
                                     help='Column name for phenotype data (default: interaction).')
    optional_input_group.add_argument('--remove_suffix', action='store_true',
                                     help='Remove suffix from genome names when merging.')
    optional_input_group.add_argument('--use_shap', action='store_true', 
                                     help='Use SHAP values for analysis (default: False).')
    optional_input_group.add_argument('--run_predictive_analysis', action='store_true',
                                     help='Run predictive analysis workflow.')

    # Feature selection and modeling parameters  
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--filter_type', type=str, default='strain', 
                                  help="Filter type for train-test split ('none', 'strain', 'phage') (default: strain).")
    fs_modeling_group.add_argument('--method', type=str, default='rfe', 
                                  choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                                  help="Feature selection method (default: rfe).")
    fs_modeling_group.add_argument('--num_features', type=int, default=100, 
                                  help='Number of features to select (default: 100).')
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, 
                                  help='Number of feature selection iterations (default: 10).')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=10, 
                                  help='Number of runs per feature table for modeling (default: 10).')
    fs_modeling_group.add_argument('--task_type', type=str, default='classification', 
                                  choices=['classification', 'regression'], 
                                  help="Task type for modeling (default: classification).")
    fs_modeling_group.add_argument('--max_features', default='none', 
                                  help='Maximum number of features to include.')

    # Advanced feature selection parameters
    advanced_group = parser.add_argument_group('Advanced feature selection')
    advanced_group.add_argument('--use_dynamic_weights', action='store_true', 
                               help='Use dynamic weights for feature selection.')
    advanced_group.add_argument('--weights_method', type=str, default='log10', 
                               choices=['log10', 'inverse_frequency', 'balanced'], 
                               help='Method for calculating dynamic weights (default: log10).')
    advanced_group.add_argument('--use_clustering', action='store_true',
                               help='Use clustering for feature selection.')
    advanced_group.add_argument('--cluster_method', type=str, default='hdbscan', 
                               choices=['hdbscan', 'hierarchical'], 
                               help='Clustering method (default: hdbscan).')
    advanced_group.add_argument('--n_clusters', type=int, default=20, 
                               help='Number of clusters for hierarchical clustering (default: 20).')
    advanced_group.add_argument('--min_cluster_size', type=int, default=5, 
                               help='Minimum cluster size for HDBSCAN (default: 5).')
    advanced_group.add_argument('--min_samples', type=int, 
                               help='Minimum number of samples for HDBSCAN.')
    advanced_group.add_argument('--cluster_selection_epsilon', type=float, default=0.0, 
                               help='Epsilon value for HDBSCAN (default: 0.0).')
    advanced_group.add_argument('--check_feature_presence', action='store_true', 
                               help='Check for feature presence during train-test split.')
    advanced_group.add_argument('--filter_by_cluster_presence', action='store_true', 
                               help='Filter features by cluster presence.')
    advanced_group.add_argument('--min_cluster_presence', type=int, default=2, 
                               help='Minimum cluster presence for features (default: 2).')

    # Feature clustering parameters
    feature_clustering_group = parser.add_argument_group('Feature clustering (pre-processing)')
    feature_clustering_group.add_argument('--use_feature_clustering', action='store_true', 
                                        help='Enable pre-processing cluster-based feature filtering')
    feature_clustering_group.add_argument('--feature_cluster_method', default='hierarchical', 
                                        choices=['hierarchical'], help='Pre-processing clustering method')
    feature_clustering_group.add_argument('--feature_n_clusters', type=int, default=20, 
                                        help='Number of clusters for pre-processing (default: 20)')
    feature_clustering_group.add_argument('--feature_min_cluster_presence', type=int, default=2, 
                                        help='Min cluster presence for pre-processing (default: 2)')

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('-o', '--output_dir', type=str, required=True, help='Output directory.')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads (default: 4).')
    general_group.add_argument('--max_ram', type=float, default=8, 
                              help='Maximum RAM usage in GB (default: 8).')

    args = parser.parse_args()

    # Run the full k-mer workflow
    run_kmer_workflow(
        input_strain_dir=args.input_strain_dir,
        input_phage_dir=args.input_phage_dir,
        phenotype_matrix=args.phenotype_matrix,
        output_dir=args.output_dir,
        k=args.k,
        k_range=args.k_range,
        one_gene=args.one_gene,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        method=args.method,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram,
        use_shap=args.use_shap,
        use_dynamic_weights=args.use_dynamic_weights,
        weights_method=args.weights_method,
        use_clustering=args.use_clustering,
        cluster_method=args.cluster_method,
        n_clusters=args.n_clusters,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
        cluster_selection_epsilon=args.cluster_selection_epsilon,
        check_feature_presence=args.check_feature_presence,
        filter_by_cluster_presence=args.filter_by_cluster_presence,
        min_cluster_presence=args.min_cluster_presence,
        use_feature_clustering=args.use_feature_clustering,
        feature_cluster_method=args.feature_cluster_method,
        feature_n_clusters=args.feature_n_clusters,
        feature_min_cluster_presence=args.feature_min_cluster_presence,
        remove_suffix=args.remove_suffix,
        run_predictive_analysis=args.run_predictive_analysis
    )


if __name__ == "__main__":
    main()