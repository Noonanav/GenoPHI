import os
import argparse
import time
import logging
import psutil
import numpy as np
import pandas as pd
from genophi.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables

# Configure logging
def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "feature_clustering_workflow.log")
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)

# Generate report file
def write_report(output_dir, start_time, end_time, ram_usage, avg_cpu_usage, max_cpu_usage, input_genomes, output_genomes, protein_families, features):
    report_file = os.path.join(output_dir, "feature_clustering_workflow_report.txt")
    with open(report_file, "w") as report:
        report.write("Workflow Report\n")
        report.write("=" * 40 + "\n")
        report.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        report.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        report.write(f"Total Runtime: {end_time - start_time:.2f} seconds\n")
        report.write(f"Max RAM Usage: {ram_usage / (1024 ** 3):.2f} GB\n")
        report.write(f"Average CPU Usage: {avg_cpu_usage:.2f}%\n")
        report.write(f"Max CPU Usage: {max_cpu_usage:.2f}%\n")
        report.write(f"Input Genomes: {input_genomes}\n")
        report.write(f"Output Genomes: {output_genomes}\n")
        report.write(f"Protein Families Identified: {protein_families}\n")
        report.write(f"Features in Final Table: {features}\n")
    logging.info(f"Report saved to: {report_file}")

# Main workflow function
def run_full_feature_workflow(
    input_path_strain, 
    output_dir, 
    phenotype_matrix, 
    tmp_dir="tmp", 
    input_path_phage=None, 
    min_seq_id=0.6, 
    coverage=0.8, 
    sensitivity=7.5, 
    suffix='faa', 
    threads=4, 
    strain_list=None, 
    strain_column='strain', 
    phage_list=None, 
    phage_column='phage', 
    compare=False, 
    source_strain='strain', 
    source_phage='phage', 
    max_ram=8,
    use_feature_clustering=False,
    feature_cluster_method='hierarchical',
    feature_n_clusters=20,
    feature_min_cluster_presence=2
):
    """
    Combines MMseqs2 clustering, feature assignment for strain (and optionally phage) genomes, 
    and merges feature tables with the phenotype matrix.
    
    Args:
        max_ram (float): Maximum allowable RAM usage in GB for feature selection.
    """
    # Set up logging
    os.makedirs(output_dir, exist_ok=True)
    setup_logging(output_dir)

    # Track time and resource usage
    start_time = time.time()
    ram_monitor = psutil.Process()
    cpu_usage_points = []  # Track CPU usage at different stages
    max_ram_usage = 0

    # Counters for the report
    input_genomes = output_genomes = protein_families = features = 0

    try:
        # Define separate temporary directories for strain and phage
        strain_tmp_dir = os.path.join(tmp_dir, "strain")

        # Run clustering and feature assignment for strain
        logging.info("Running clustering workflow for strain genomes...")
        cpu_usage_points.append(psutil.cpu_percent(interval=None))
        strain_output_dir = os.path.join(output_dir, "strain")
        run_clustering_workflow(input_path_strain, strain_output_dir, strain_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list or 'none', strain_column, compare)

        presence_absence_strain = os.path.join(strain_output_dir, "presence_absence_matrix.csv")
        feature_output_dir_strain = os.path.join(strain_output_dir, "features")

        logging.info("Running feature assignment workflow for strain genomes...")
        run_feature_assignment(
            input_file=presence_absence_strain, 
            output_dir=feature_output_dir_strain, 
            source=source_strain, 
            select=strain_list or 'none', 
            select_column=strain_column,
            max_ram=max_ram,  # Pass max_ram here
            threads=threads
        )

        # Count genomes and protein families for strain
        strain_matrix = os.path.join(strain_output_dir, "presence_absence_matrix.csv")
        strain_df = pd.read_csv(strain_matrix)
        input_genomes += len(strain_df['Genome'].unique())
        protein_families += len(strain_df.columns) - 1  # Exclude 'Genome'

        if input_path_phage:
            # Run clustering and feature assignment for phage
            logging.info("Running clustering workflow for phage genomes...")
            cpu_usage_points.append(psutil.cpu_percent(interval=None))
            phage_output_dir = os.path.join(output_dir, "phage")
            phage_tmp_dir = os.path.join(tmp_dir, "phage")

            run_clustering_workflow(input_path_phage, phage_output_dir, phage_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, phage_list or 'none', phage_column, compare)

            presence_absence_phage = os.path.join(phage_output_dir, "presence_absence_matrix.csv")
            feature_output_dir_phage = os.path.join(phage_output_dir, "features")

            logging.info("Running feature assignment workflow for phage genomes...")
            run_feature_assignment(
                input_file=presence_absence_phage, 
                output_dir=feature_output_dir_phage, 
                source=source_phage, 
                select=phage_list or 'none', 
                select_column=phage_column,
                max_ram=max_ram,  # Pass max_ram here
                threads=threads
            )

            # Merge strain and phage feature tables
            logging.info("Merging feature tables for strain and phage genomes...")
            strain_features = os.path.join(feature_output_dir_strain, "feature_table.csv")
            phage_features = os.path.join(feature_output_dir_phage, "feature_table.csv")

            merged_output_dir = os.path.join(output_dir, "merged")
            os.makedirs(merged_output_dir, exist_ok=True)

            merge_feature_tables(
                strain_features=strain_features, 
                phenotype_matrix=phenotype_matrix,
                output_dir=merged_output_dir,
                sample_column=strain_column,
                phage_features=phage_features,
                remove_suffix=False,
                use_feature_clustering=use_feature_clustering,
                feature_cluster_method=feature_cluster_method,
                feature_n_clusters=feature_n_clusters,
                feature_min_cluster_presence=feature_min_cluster_presence
            )
            logging.info(f"Merged feature table saved in: {merged_output_dir}")
        else:
            # Only strain data: merge with phenotype_matrix
            logging.info("Merging strain features with the phenotype matrix...")
            strain_features = os.path.join(feature_output_dir_strain, "feature_table.csv")

            merge_feature_tables(
                strain_features=strain_features,
                phenotype_matrix=phenotype_matrix,
                output_dir=output_dir,
                sample_column=strain_column,
                remove_suffix=False,
                use_feature_clustering=use_feature_clustering,
                feature_cluster_method=feature_cluster_method,
                feature_n_clusters=feature_n_clusters,
                feature_min_cluster_presence=feature_min_cluster_presence
            )
            logging.info(f"Strain feature table merged with phenotype matrix and saved at: {output_dir}")

        # Track final output genomes and features
        final_table = os.path.join(output_dir, "full_feature_table.csv")
        final_df = pd.read_csv(final_table)
        output_genomes = len(final_df[strain_column].unique())
        features = len(final_df.columns) - 1  # Exclude strain_column

    except Exception as e:
        logging.error(f"An error occurred during workflow execution: {e}")
        raise
    finally:
        end_time = time.time()
        max_ram_usage = max(max_ram_usage, ram_monitor.memory_info().rss)
        cpu_usage_points.append(psutil.cpu_percent(interval=None))

        # Calculate average and maximum CPU usage
        avg_cpu_usage = sum(cpu_usage_points) / len(cpu_usage_points)
        max_cpu_usage = max(cpu_usage_points)

        # Write report
        write_report(output_dir, start_time, end_time, max_ram_usage, avg_cpu_usage, max_cpu_usage, input_genomes, output_genomes, protein_families, features)

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run full feature table generation and merging workflow.')
    parser.add_argument('-ih', '--input_strain', type=str, required=True, help='Input path for strain clustering (directory or file).')
    parser.add_argument('-ip', '--input_phage', type=str, help='Input path for phage clustering (directory or file). Optional; if not provided, only strain data will be used.')
    parser.add_argument('-pm', '--phenotype_matrix', type=str, required=True, help='Path to the phenotype matrix.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory to save results.')
    parser.add_argument('--tmp', type=str, default="tmp", help='Temporary directory for intermediate files.')
    parser.add_argument('--min_seq_id', type=float, default=0.6, help='Minimum sequence identity for clustering.')
    parser.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering.')
    parser.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering.')
    parser.add_argument('--suffix', type=str, default='faa', help='Suffix for input FASTA files.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use.')
    parser.add_argument('--strain_list', type=str, help='Path to a strain list file for filtering.')
    parser.add_argument('--strain_column', type=str, default='strain', help='Column in the strain list containing strain names.')
    parser.add_argument('--phage_list', type=str, help='Path to a phage list file for filtering.')
    parser.add_argument('--phage_column', type=str, default='phage', help='Column in the phage list containing phage names.')
    parser.add_argument('--compare', action='store_true', help='Compare original clusters with assigned clusters.')
    parser.add_argument('--source_strain', type=str, default='strain', help='Prefix for naming selected features for strain in the assignment step.')
    parser.add_argument('--source_phage', type=str, default='phage', help='Prefix for naming selected features for phage in the assignment step.')
    parser.add_argument('--max_ram', type=float, default=8, help='Maximum RAM usage in GB for feature selection.')
    parser.add_argument('--use_feature_clustering', action='store_true', help='Enable pre-processing cluster-based feature filtering')
    parser.add_argument('--feature_cluster_method', default='hierarchical', choices=['hierarchical'], help='Pre-processing clustering method')
    parser.add_argument('--feature_n_clusters', type=int, default=20, help='Number of clusters for pre-processing feature clustering')
    parser.add_argument('--feature_min_cluster_presence', type=int, default=2, help='Min clusters a feature must appear in during pre-processing')

    args = parser.parse_args()

    # Run the full feature workflow
    run_full_feature_workflow(
        input_path_strain=args.input_strain,
        input_path_phage=args.input_phage,  # Optional; may be None if not provided
        phenotype_matrix=args.phenotype_matrix,
        output_dir=args.output,
        tmp_dir=args.tmp,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        strain_column=args.strain_column,
        phage_list=args.phage_list,
        phage_column=args.phage_column,
        compare=args.compare,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        max_ram=args.max_ram,
        use_feature_clustering=args.use_feature_clustering,
        feature_cluster_method=args.feature_cluster_method,
        feature_n_clusters=args.feature_n_clusters,
        feature_min_cluster_presence=args.feature_min_cluster_presence
    )

if __name__ == "__main__":
    main()
