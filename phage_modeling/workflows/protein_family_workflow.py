import os
import pandas as pd
import argparse
import logging
import csv
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables
from phage_modeling.feature_selection import run_feature_selection_iterations, generate_feature_tables
from phage_modeling.select_feature_modeling import run_experiments
from phage_modeling.workflows.feature_annotations_workflow import run_predictive_proteins_workflow
import time
import psutil

def setup_logging(output_dir):
    """
    Sets up logging to both a file and console.

    Args:
        output_dir (str): Directory where the log file will be saved.
    """
    log_file = os.path.join(output_dir, 'workflow.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logging.info("Logging setup complete.")

def write_csv_log(output_dir, data):
    """
    Writes a CSV log file with variable names and their values.

    Args:
        output_dir (str): Directory where the CSV log file will be saved.
        data (dict): Dictionary of variable names and their values.
    """
    csv_file = os.path.join(output_dir, 'workflow_report.csv')
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Variable', 'Value'])
        for key, value in data.items():
            writer.writerow([key, value])
    logging.info(f"CSV log saved to {csv_file}.")

def write_report(output_dir, start_time, end_time, ram_usage, avg_cpu_usage, max_cpu_usage, input_genomes, output_genomes, protein_families, features):
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
        output_genomes (int): Number of output genomes.
        protein_families (int): Number of protein families identified.
        features (int): Number of features in the final table.
    """
    report_file = os.path.join(output_dir, "workflow_report.txt")
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

def run_protein_family_workflow(input_path_strain, output_dir, phenotype_matrix, tmp_dir="tmp", 
                      input_path_phage=None, min_seq_id=0.6, coverage=0.8, sensitivity=7.5, 
                      suffix='faa', threads=4, strain_list='none', phage_list='none', 
                      strain_column='strain', phage_column='phage', compare=False, 
                      source_strain='strain', source_phage='phage', num_features=100, 
                      filter_type='none', num_runs_fs=10, num_runs_modeling=10, 
                      sample_column='strain', phenotype_column=None, method='rfe',
                      annotation_table_path=None, protein_id_col="protein_ID",
                      task_type='classification', max_features='none', max_ram=8):
    """
    Complete workflow: Feature table generation, feature selection, modeling, and predictive proteins extraction.
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
    input_genomes = output_genomes = protein_families = features = 0

    try:
        logging.info("Step 1: Running feature table generation for strain...")
        strain_output_dir = os.path.join(output_dir, "strain")
        strain_tmp_dir = os.path.join(output_dir, "tmp", "strain")
        run_clustering_workflow(input_path_strain, strain_output_dir, strain_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list, strain_column, compare)
        run_feature_assignment(
            os.path.join(strain_output_dir, "presence_absence_matrix.csv"), 
            os.path.join(strain_output_dir, "features"), 
            source=source_strain, 
            select=strain_list, 
            select_column=strain_column,
            max_ram=max_ram
        )

        strain_matrix = os.path.join(strain_output_dir, "presence_absence_matrix.csv")
        strain_df = pd.read_csv(strain_matrix)
        input_genomes += len(strain_df['Genome'].unique())
        protein_families += len(strain_df.columns) - 1

        if input_path_phage:
            logging.info("Step 1 (continued): Running feature table generation for phage...")
            phage_output_dir = os.path.join(output_dir, "phage")
            phage_tmp_dir = os.path.join(output_dir, "tmp", "phage")
            run_clustering_workflow(input_path_phage, phage_output_dir, phage_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, phage_list, phage_column, compare)
            run_feature_assignment(
                os.path.join(phage_output_dir, "presence_absence_matrix.csv"), 
                os.path.join(phage_output_dir, "features"), 
                source=source_phage, 
                select=phage_list, 
                select_column=phage_column,
                max_ram=max_ram
            )

            merged_output_dir = os.path.join(output_dir, "merged")
            os.makedirs(merged_output_dir, exist_ok=True)
            feature_selection_input = merge_feature_tables(
                strain_features=os.path.join(strain_output_dir, "features", "feature_table.csv"),
                phenotype_matrix=phenotype_matrix,
                output_dir=merged_output_dir,
                sample_column=sample_column,
                phage_features=os.path.join(phage_output_dir, "features", "feature_table.csv"),
                remove_suffix=False
            )
            logging.info(f"Merged feature table saved in: {merged_output_dir}")
        else:
            logging.info("Merging strain feature table with phenotype matrix...")
            strain_features = os.path.join(strain_output_dir, "features", "feature_table.csv")
            feature_selection_input = merge_feature_tables(
                strain_features=strain_features,
                phenotype_matrix=phenotype_matrix,
                output_dir=output_dir,
                sample_column=sample_column,
                remove_suffix=False
            )
            logging.info(f"Strain feature table merged and saved at: {feature_selection_input}")

        logging.info("Step 2: Running feature selection iterations...")
        base_fs_output_dir = os.path.join(output_dir, 'feature_selection')
        run_feature_selection_iterations(
            input_path=feature_selection_input,
            base_output_dir=base_fs_output_dir,
            threads=threads,
            num_features=num_features,
            filter_type=filter_type,
            num_runs=num_runs_fs,
            method=method,
            sample_column=sample_column,
            phenotype_column=phenotype_column,
            task_type=task_type,
            max_ram=max_ram
        )

        logging.info("Generating feature tables from feature selection results...")
        max_features = None if max_features == 'none' else int(max_features)
        filter_table_dir = os.path.join(base_fs_output_dir, 'filtered_feature_tables')
        generate_feature_tables(
            model_testing_dir=base_fs_output_dir,
            full_feature_table_file=feature_selection_input,
            filter_table_dir=filter_table_dir,
            phenotype_column=phenotype_column,
            sample_column=sample_column,
            cut_offs=[3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50],
            binary_data=True,
            max_features=max_features,
            filter_type=filter_type
        )

        logging.info("Step 4: Running modeling experiments...")
        run_experiments(
            input_dir=filter_table_dir,
            base_output_dir=os.path.join(output_dir, 'modeling_results'),
            threads=threads,
            num_runs=num_runs_modeling,
            set_filter=filter_type,
            sample_column=sample_column,
            phenotype_column=phenotype_column,
            task_type=task_type,
            binary_data=True
        )

        logging.info("Step 5: Selecting top-performing cutoff and running predictive proteins workflow...")
        metrics_file = os.path.join(output_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
        performance_df = pd.read_csv(metrics_file)
        top_cutoff = performance_df.iloc[0]['cut_off'].split('_')[-1]

        feature_file_path = os.path.join(output_dir, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_cutoff_{top_cutoff}.csv')
        feature2cluster_path = os.path.join(output_dir, 'strain', 'features', 'selected_features.csv')
        cluster2protein_path = os.path.join(output_dir, 'strain', 'clusters.tsv')
        fasta_dir_or_file = input_path_strain
        modeling_dir = os.path.join(output_dir, 'modeling_results', f'cutoff_{top_cutoff}')
        predictive_proteins_output_dir = os.path.join(output_dir, 'modeling_results', 'model_performance', 'predictive_proteins')
        feature_assignments_path = os.path.join(output_dir, 'strain', 'features', 'feature_assignments.csv')

        run_predictive_proteins_workflow(
            feature_file_path=feature_file_path,
            feature2cluster_path=feature2cluster_path,
            cluster2protein_path=cluster2protein_path,
            fasta_dir_or_file=fasta_dir_or_file,
            modeling_dir=modeling_dir,
            output_dir=predictive_proteins_output_dir,
            output_fasta='predictive_AA_seqs.faa',
            protein_id_col=protein_id_col,
            annotation_table_path=annotation_table_path,  # Optional
            feature_assignments_path=feature_assignments_path,  # Optional
            strain_column='strain',
            feature_type='strain'
        )

        if input_path_phage:
            logging.info("Running predictive proteins workflow for phage...")
            feature2cluster_path = os.path.join(output_dir, 'phage', 'features', 'selected_features.csv')
            cluster2protein_path = os.path.join(output_dir, 'phage', 'clusters.tsv')
            fasta_dir_or_file = input_path_phage
            feature_assignments_path = os.path.join(output_dir, 'phage', 'features', 'feature_assignments.csv')

            run_predictive_proteins_workflow(
                feature_file_path=feature_file_path,
                feature2cluster_path=feature2cluster_path,
                cluster2protein_path=cluster2protein_path,
                fasta_dir_or_file=fasta_dir_or_file,
                modeling_dir=modeling_dir,
                output_dir=predictive_proteins_output_dir,
                output_fasta='predictive_AA_seqs.faa',
                protein_id_col=protein_id_col,
                feature_assignments_path=feature_assignments_path,  # Optional
                strain_column='phage',
                feature_type='phage'
            )
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
    finally:
        end_time = time.time()
        max_ram_usage = max(max_ram_usage, ram_monitor.memory_info().rss)
        cpu_usage_points.append(psutil.cpu_percent(interval=None))

        avg_cpu_usage = sum(cpu_usage_points) / len(cpu_usage_points)
        max_cpu_usage = max(cpu_usage_points)

        # Write report and logs
        write_report(output_dir, start_time, end_time, max_ram_usage, avg_cpu_usage, max_cpu_usage, input_genomes, output_genomes, protein_families, features)

        # Save CSV log with all parameters and report values
        inputs.update({
            'start_time': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time)),
            'end_time': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time)),
            'ram_usage': max_ram_usage / (1024 ** 3),
            'avg_cpu_usage': avg_cpu_usage,
            'max_cpu_usage': max_cpu_usage,
            'input_genomes': input_genomes,
            'output_genomes': output_genomes,
            'protein_families': protein_families,
            'features': features
        })
        write_csv_log(output_dir, inputs)


# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Complete workflow: Feature table generation, feature selection, modeling, and predictive proteins extraction.')
    
    # Input data
    input_group = parser.add_argument_group('Input data')
    input_group.add_argument('-ih', '--input_strain', type=str, required=True, help='Path to the input directory or file for strain clustering.')
    input_group.add_argument('-ip', '--input_phage', type=str, help='Path to the input directory or file for phage clustering. Optional; if not provided, only strain data will be used.')
    input_group.add_argument('-pm', '--phenotype_matrix', type=str, required=True, help='Path to the phenotype matrix.')

    # Optional input arguments
    optional_input_group = parser.add_argument_group('Optional input arguments')
    optional_input_group.add_argument('--suffix', type=str, default='faa', help='Suffix for input FASTA files (default: faa).')
    optional_input_group.add_argument('--strain_list', type=str, default='none', help='Path to a strain list file for filtering (default: none).')
    optional_input_group.add_argument('--phage_list', type=str, default='none', help='Path to a phage list file for filtering (default: none).')
    optional_input_group.add_argument('--strain_column', type=str, default='strain', help='Column in the strain list containing strain names (default: strain).')
    optional_input_group.add_argument('--phage_column', type=str, default='phage', help='Column in the phage list containing phage names (default: phage).')
    optional_input_group.add_argument('--source_strain', type=str, default='strain', help='Prefix for naming selected features for strain in the assignment step (default: strain).')
    optional_input_group.add_argument('--source_phage', type=str, default='phage', help='Prefix for naming selected features for phage in the assignment step (default: phage).')
    optional_input_group.add_argument('--sample_column', type=str, default='strain', help='Column name for the sample identifier (default: strain).')
    optional_input_group.add_argument('--phenotype_column', type=str, default='interaction', help='Column name for the phenotype (optional).')
    optional_input_group.add_argument('--annotation_table_path', type=str, help="Path to an optional annotation table (CSV/TSV).")
    optional_input_group.add_argument('--protein_id_col', type=str, default="protein_ID", help="Column name for protein IDs in the predictive_proteins DataFrame.")

    # Output arguments
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument('-o', '--output', type=str, required=True, help='Output directory to save results.')

    # Clustering parameters
    clustering_group = parser.add_argument_group('Clustering')
    clustering_group.add_argument('--min_seq_id', type=float, default=0.6, help='Minimum sequence identity for clustering (default: 0.6).')
    clustering_group.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering (default: 0.8).')
    clustering_group.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering (default: 7.5).')
    clustering_group.add_argument('--compare', action='store_true', help='Compare original clusters with assigned clusters.')

    # Feature selection and modeling parameters
    fs_modeling_group = parser.add_argument_group('Feature selection and modeling')
    fs_modeling_group.add_argument('--filter_type', type=str, default='none', help="Filter type for the input data ('none', 'strain', 'phage').")
    fs_modeling_group.add_argument('--method', type=str, default='rfe', choices=['rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'],
                                   help="Feature selection method ('rfe', 'shap_rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'; default: rfe).")
    fs_modeling_group.add_argument('--num_features', type=int, default=100, help='Number of features to select (default: 100).')
    fs_modeling_group.add_argument('--num_runs_fs', type=int, default=10, help='Number of feature selection iterations to run (default: 10).')
    fs_modeling_group.add_argument('--num_runs_modeling', type=int, default=10, help='Number of runs per feature table for modeling (default: 10).')
    fs_modeling_group.add_argument('--task_type', type=str, default='classification', choices=['classification', 'regression'], help="Task type for modeling ('classification' or 'regression').")
    fs_modeling_group.add_argument('--max_features', default='none', help='Maximum number of features to include in the feature tables.')

    # General parameters
    general_group = parser.add_argument_group('General')
    general_group.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4).')
    general_group.add_argument('--max_ram', type=float, default=8, help='Maximum RAM usage in GB for feature selection (default: 8).')
    
    args = parser.parse_args()

    # Run the full workflow
    run_protein_family_workflow(
        input_path_strain=args.input_strain,
        input_path_phage=args.input_phage,  # Optional; may be None if not provided
        phenotype_matrix=args.phenotype_matrix,
        output_dir=args.output,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        phage_list=args.phage_list,
        strain_column=args.strain_column,
        phage_column=args.phage_column,
        compare=args.compare,
        source_strain=args.source_strain,
        source_phage=args.source_phage,
        num_features=args.num_features,
        filter_type=args.filter_type,
        num_runs_fs=args.num_runs_fs,
        num_runs_modeling=args.num_runs_modeling,
        sample_column=args.sample_column,
        phenotype_column=args.phenotype_column,
        method=args.method,
        annotation_table_path=args.annotation_table_path,  # Optional
        protein_id_col=args.protein_id_col,
        task_type=args.task_type,
        max_features=args.max_features,
        max_ram=args.max_ram
    )

if __name__ == "__main__":
    main()
