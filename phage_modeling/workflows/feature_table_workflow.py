import os
import argparse
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables

def run_full_feature_workflow(input_path_host, input_path_phage, interaction_matrix, output_dir, tmp_dir="tmp", min_seq_id=0.6, coverage=0.8, sensitivity=7.5, suffix='faa', threads=4, strain_list=None, strain_column='strain', compare=False, source_host='host', source_phage='phage'):
    """
    Combines MMseqs2 clustering, feature assignment for both host and phage genomes, and merges feature tables.

    Args:
        input_path_host (str): Path to the input directory or file for host clustering.
        input_path_phage (str): Path to the input directory or file for phage clustering.
        interaction_matrix (str): Path to the interaction matrix.
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.
        min_seq_id (float): Minimum sequence identity for clustering.
        coverage (float): Minimum coverage for clustering.
        sensitivity (float): Sensitivity for clustering.
        suffix (str): Suffix for input FASTA files.
        threads (int): Number of threads to use.
        strain_list (str or None): Path to a strain list file, or None for no filtering.
        strain_column (str): Column in the strain list file containing strain names.
        compare (bool): Whether to compare original clusters with assigned clusters.
        source_host (str): Prefix for naming selected features for host in the assignment step.
        source_phage (str): Prefix for naming selected features for phage in the assignment step.
    """
    # Run clustering and feature assignment for host
    print("Running clustering workflow for host genomes...")
    host_output_dir = os.path.join(output_dir, "host")
    run_clustering_workflow(input_path_host, host_output_dir, tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list, strain_column, compare)
    presence_absence_host = os.path.join(host_output_dir, "presence_absence_matrix.csv")
    feature_output_dir_host = os.path.join(host_output_dir, "features")
    print("Running feature assignment workflow for host genomes...")
    run_feature_assignment(presence_absence_host, feature_output_dir_host, source=source_host)

    # Run clustering and feature assignment for phage
    print("Running clustering workflow for phage genomes...")
    phage_output_dir = os.path.join(output_dir, "phage")
    run_clustering_workflow(input_path_phage, phage_output_dir, tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list, strain_column, compare)
    presence_absence_phage = os.path.join(phage_output_dir, "presence_absence_matrix.csv")
    feature_output_dir_phage = os.path.join(phage_output_dir, "features")
    print("Running feature assignment workflow for phage genomes...")
    run_feature_assignment(presence_absence_phage, feature_output_dir_phage, source=source_phage)

    # Merge host and phage feature tables
    print("Merging feature tables for host and phage genomes...")
    host_features = os.path.join(feature_output_dir_host, "feature_table.csv")
    phage_features = os.path.join(feature_output_dir_phage, "feature_table.csv")
    merged_output_dir = os.path.join(output_dir, "merged")
    os.makedirs(merged_output_dir, exist_ok=True)
    merge_feature_tables(host_features, phage_features, interaction_matrix, merged_output_dir, remove_suffix=False)

    print(f"Merged feature table saved in: {merged_output_dir}")

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run full feature table generation and merging workflow.')
    parser.add_argument('-ih', '--input_host', type=str, required=True, help='Input path for host clustering (directory or file).')
    parser.add_argument('-ip', '--input_phage', type=str, required=True, help='Input path for phage clustering (directory or file).')
    parser.add_argument('-im', '--interaction_matrix', type=str, required=True, help='Path to the interaction matrix.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory to save results.')
    parser.add_argument('--tmp', type=str, default="tmp", help='Temporary directory for intermediate files.')
    parser.add_argument('--min_seq_id', type=float, default=0.6, help='Minimum sequence identity for clustering.')
    parser.add_argument('--coverage', type=float, default=0.8, help='Minimum coverage for clustering.')
    parser.add_argument('--sensitivity', type=float, default=7.5, help='Sensitivity for clustering.')
    parser.add_argument('--suffix', type=str, default='faa', help='Suffix for input FASTA files.')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use.')
    parser.add_argument('--strain_list', type=str, help='Path to a strain list file for filtering.')
    parser.add_argument('--strain_column', type=str, default='strain', help='Column in the strain list containing strain names.')
    parser.add_argument('--compare', action='store_true', help='Compare original clusters with assigned clusters.')
    parser.add_argument('--source_host', type=str, default='host', help='Prefix for naming selected features for host in the assignment step.')
    parser.add_argument('--source_phage', type=str, default='phage', help='Prefix for naming selected features for phage in the assignment step.')

    args = parser.parse_args()

    # Run the full feature workflow
    run_full_feature_workflow(
        input_path_host=args.input_host,
        input_path_phage=args.input_phage,
        interaction_matrix=args.interaction_matrix,
        output_dir=args.output,
        tmp_dir=args.tmp,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        sensitivity=args.sensitivity,
        suffix=args.suffix,
        threads=args.threads,
        strain_list=args.strain_list,
        strain_column=args.strain_column,
        compare=args.compare,
        source_host=args.source_host,
        source_phage=args.source_phage
    )

if __name__ == "__main__":
    main()
