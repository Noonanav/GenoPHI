import os
import argparse
from phage_modeling.mmseqs2_clustering import run_clustering_workflow, run_feature_assignment, merge_feature_tables

def run_full_feature_workflow(input_path_strain, output_dir, phenotype_matrix, tmp_dir="tmp", 
                              input_path_phage=None, min_seq_id=0.6, coverage=0.8, 
                              sensitivity=7.5, suffix='faa', threads=4, strain_list=None, 
                              strain_column='strain', phage_list=None, phage_column='phage', 
                              compare=False, source_strain='strain', source_phage='phage'):
    """
    Combines MMseqs2 clustering, feature assignment for strain (and optionally phage) genomes, 
    and merges feature tables with the phenotype matrix.
    """
    
    # Define separate temporary directories for strain and phage to avoid conflicts
    strain_tmp_dir = os.path.join(tmp_dir, "strain")
    
    # Run clustering and feature assignment for strain
    print("Running clustering workflow for strain genomes...")
    strain_output_dir = os.path.join(output_dir, "strain")
    run_clustering_workflow(input_path_strain, strain_output_dir, strain_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, strain_list or 'none', strain_column, compare)
    
    presence_absence_strain = os.path.join(strain_output_dir, "presence_absence_matrix.csv")
    feature_output_dir_strain = os.path.join(strain_output_dir, "features")
    
    print("Running feature assignment workflow for strain genomes...")
    strain_input_type = 'directory' if os.path.isdir(input_path_strain) else 'file'
    run_feature_assignment(presence_absence_strain, feature_output_dir_strain, source=source_strain, select=strain_list or 'none', select_column=strain_column, input_type=strain_input_type)

    if input_path_phage:
        # Run clustering and feature assignment for phage if provided
        print("Running clustering workflow for phage genomes...")
        phage_output_dir = os.path.join(output_dir, "phage")
        phage_tmp_dir = os.path.join(tmp_dir, "phage")
        
        run_clustering_workflow(input_path_phage, phage_output_dir, phage_tmp_dir, min_seq_id, coverage, sensitivity, suffix, threads, phage_list or 'none', phage_column, compare)
        
        presence_absence_phage = os.path.join(phage_output_dir, "presence_absence_matrix.csv")
        feature_output_dir_phage = os.path.join(phage_output_dir, "features")
        
        print("Running feature assignment workflow for phage genomes...")
        phage_input_type = 'directory' if os.path.isdir(input_path_phage) else 'file'
        run_feature_assignment(presence_absence_phage, feature_output_dir_phage, source=source_phage, select=phage_list or 'none', select_column=phage_column, input_type=phage_input_type)
        
        if phenotype_matrix:
            # Merge strain and phage feature tables
            print("Merging feature tables for strain and phage genomes...")
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
                remove_suffix=False
            )
            print(f"Merged feature table saved in: {merged_output_dir}")
        else:
            print("No phenotype matrix provided. Skipping merging step.")

    else:
        if phenotype_matrix:     
            # Only strain data: merge with phenotype_matrix
            print("Merging strain features with the phenotype matrix...")
            strain_features = os.path.join(feature_output_dir_strain, "feature_table.csv")   

            merged_output_dir = os.path.join(output_dir, "merged")
            os.makedirs(merged_output_dir, exist_ok=True)

            merge_feature_tables(
                strain_features=strain_features,
                phenotype_matrix=phenotype_matrix,
                output_dir=merged_output_dir,
                sample_column=strain_column,
                remove_suffix=False
            )
            print(f"Strain feature table merged with phenotype matrix and saved at: {merged_output_dir}")
        else:
            print("No phenotype matrix provided. Skipping merging step.")

# Main function for CLI
def main():
    parser = argparse.ArgumentParser(description='Run full feature table generation and merging workflow.')
    parser.add_argument('-ih', '--input_strain', type=str, required=True, help='Input path for strain clustering (directory or file).')
    parser.add_argument('-ip', '--input_phage', type=str, help='Input path for phage clustering (directory or file). Optional; if not provided, only strain data will be used.')
    parser.add_argument('-pm', '--phenotype_matrix', type=str, help='Path to the phenotype matrix.')
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
        source_phage=args.source_phage
    )

if __name__ == "__main__":
    main()
