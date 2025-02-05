import os
import logging
from argparse import ArgumentParser
import pandas as pd
from Bio import SeqIO
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

def create_contig_to_genome_dict(fasta_files, input_type):
    """
    Creates a mapping from contigs to genomes based on input FASTA files.

    Args:
        fasta_files (list): List of FASTA file paths.
        input_type (str): Type of input, either 'directory' or 'file'.

    Returns:
        dict: Dictionary mapping contig IDs to genome names.
    """
    contig_to_genome = {}
    logging.info("Creating contig to genome dictionary...")

    if input_type == 'directory':
        for fasta in fasta_files:
            genome_name = os.path.basename(fasta).split(".")[0]
            for record in SeqIO.parse(fasta, "fasta"):
                contig_to_genome[record.id] = genome_name
    else:
        for record in SeqIO.parse(fasta_files[0], "fasta"):
            protein_id = record.id
            genome_name = '_'.join(protein_id.split(' # ')[0].split('_')[:-1])
            contig_to_genome[protein_id] = genome_name

    logging.info(f"Created contig to genome dictionary with {len(contig_to_genome)} entries.")
    return contig_to_genome

def map_features_with_threshold(best_hits_tsv, feature_map, genome_contig_mapping, genome_type, threshold, select_protein_families=None):
    """
    Maps features for all genomes based on a dynamic threshold for matching protein families.

    Args:
        best_hits_tsv (str): Path to the best hits TSV file.
        feature_map (str): Path to the feature mapping CSV file.
        genome_contig_mapping (dict): Dictionary mapping contigs to genomes.
        genome_type (str): Type of genome ('strain' or 'phage').
        threshold (float): Minimum percentage of protein families (clusters) per feature that need to have at least one protein match.

    Returns:
        pd.DataFrame: Combined feature presence table for all genomes.
    """
    if not os.path.exists(best_hits_tsv):
        logging.error(f"Best hits TSV file does not exist: {best_hits_tsv}")
        return None

    try:
        best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster'])
        feature_mapping = pd.read_csv(feature_map)
        # logging.info("feature_mapping:")
        # logging.info(feature_mapping.head())

        best_hits_df['Cluster'] = best_hits_df['Cluster'].astype(str)
        feature_mapping['Cluster_Label'] = feature_mapping['Cluster_Label'].astype(str)

        logging.info(f"Mapping features for {genome_type}s with threshold: {threshold}...")
    except Exception as e:
        logging.error(f"Error reading input files: {e}")
        return None

    genome_contig_mapping_df = pd.DataFrame(list(genome_contig_mapping.items()), columns=['contig_id', genome_type])
    # logging.info('Genome contig mapping:')
    # logging.info(genome_contig_mapping_df.head())

    merged_df = best_hits_df.merge(feature_mapping, left_on='Cluster', right_on='Cluster_Label')
    merged_df = merged_df.merge(genome_contig_mapping_df, left_on='Query', right_on='contig_id')

    if select_protein_families != 'None':
        selected_protein_families_df = pd.read_csv(select_protein_families)
        selected_protein_families = selected_protein_families_df['protein_family'].tolist()
        
        # First, create a new DataFrame with rows where 'Cluster' is in selected_protein_families
        # This will inherently leave single-cluster features unchanged since they'll either be in or out
        merged_df_mod = pd.DataFrame(columns=merged_df.columns)
        
        for feature, group in merged_df.groupby('Feature'):
            unique_protein_families = group['Cluster'].nunique()
            if unique_protein_families > 1:
                # For features with more than one protein family, filter based on selection
                filtered_group = group[group['Cluster'].isin(selected_protein_families)]
                merged_df_mod = pd.concat([merged_df_mod, filtered_group])
            else:
                # For features with only one protein family, keep all rows
                merged_df_mod = pd.concat([merged_df_mod, group])
        merged_df = merged_df_mod
        merged_df_count = merged_df.groupby('Feature')['Cluster'].nunique().reset_index(name='count')
        merged_df_count = merged_df_count.sort_values('count', ascending=False)
        logging.info(merged_df_count.head())

    # Count unique clusters per feature per genome
    feature_cluster_count = merged_df.groupby([genome_type, 'Feature'])['Cluster'].nunique().reset_index(name='Cluster_Count')

    # Count total unique clusters per feature
    total_clusters_per_feature = merged_df.groupby('Feature')['Cluster'].nunique().reset_index(name='Total_Clusters')
    
    # Merge to get counts in context of each feature per genome
    feature_cluster_count = feature_cluster_count.merge(total_clusters_per_feature, on='Feature')

    # Apply the feature-specific threshold
    feature_cluster_count['Meets_Threshold'] = feature_cluster_count.apply(
        lambda row: (row['Cluster_Count'] / row['Total_Clusters']) >= threshold if threshold <= 1 else row['Cluster_Count'] >= threshold, 
        axis=1
    )

    logging.info('Feature cluster count with threshold:')
    logging.info(feature_cluster_count.head())

    # Pivot to get the feature presence table
    feature_presence = feature_cluster_count.pivot(index=genome_type, columns='Feature', values='Meets_Threshold').fillna(0).astype(int)
    
    # Ensure all features are present in the final table
    all_features = feature_mapping['Feature'].unique()
    for feature in all_features:
        if feature not in feature_presence.columns:
            feature_presence[feature] = 0

    feature_presence = feature_presence.reindex(columns=all_features, fill_value=0).reset_index()

    return feature_presence

def run_assign_predict_with_thresholds(fasta_files, feature_map, best_hits_tsv, phage_feature_table_path, output_dir, model_dir, thresholds, genome_type='strain', select_protein_families=None):
    """
    Runs assignment and prediction workflows for varying feature assignment thresholds.

    Args:
        fasta_files (list): List of FASTA file paths.
        feature_map (str): Path to the feature mapping CSV file.
        best_hits_tsv (str): Path to the best hits TSV file.
        output_dir (str): Directory to save results.
        model_dir (str): Directory containing trained models.
        thresholds (list): List of thresholds to evaluate.
        genome_type (str): Type of genome being processed ('strain' or 'phage').
    """
    os.makedirs(output_dir, exist_ok=True)

    # Generate contig-to-genome mapping
    genome_contig_mapping = create_contig_to_genome_dict(fasta_files, input_type='directory')

    for threshold in thresholds:
        logging.info(f"Running workflow for threshold: {threshold}")
        threshold_dir = os.path.join(output_dir, f"threshold_{threshold:.2f}")
        os.makedirs(threshold_dir, exist_ok=True)

        # Generate feature table for the threshold
        feature_table = map_features_with_threshold(best_hits_tsv, feature_map, genome_contig_mapping, genome_type, threshold, select_protein_families)
        if feature_table is None:
            logging.error(f"Feature table generation failed for threshold {threshold}. Skipping.")
            continue

        # Save the feature table
        feature_table_path = os.path.join(threshold_dir, f"{genome_type}_feature_table.csv")
        feature_table.to_csv(feature_table_path, index=False)
        logging.info(f"Feature table saved for threshold {threshold} at {feature_table_path}")

        # Directory for prediction results
        predict_output_dir = os.path.join(threshold_dir, "predict_results")
        os.makedirs(predict_output_dir, exist_ok=True)

        # Run prediction workflow
        run_prediction_workflow(
            input_dir=threshold_dir,
            phage_feature_table_path=phage_feature_table_path,
            model_dir=model_dir,
            output_dir=predict_output_dir
        )

        logging.info(f"Workflow completed for threshold {threshold}")

def main():
    parser = ArgumentParser(description="Run assignment and prediction workflows with varying feature assignment thresholds.")
    parser.add_argument('--fasta_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--feature_map', type=str, required=True, help="Path to the feature map CSV file.")
    parser.add_argument('--best_hits_tsv', type=str, required=True, help="Path to the best hits TSV file.")
    parser.add_argument('--phage_feature_table_path', type=str, help="Path to the phage feature table for prediction.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory containing trained models.")
    parser.add_argument('--thresholds', type=str, default="0.25,0.5,0.75,1.0", help="Comma-separated thresholds (e.g., 0.25,0.5).")
    parser.add_argument('--select_protein_families', type=str, default=None, help="Path to a file with selected protein families.")

    args = parser.parse_args()
    thresholds = [float(t) for t in args.thresholds.split(',')]

    # Collect FASTA files
    fasta_files = [os.path.join(args.fasta_dir, f) for f in os.listdir(args.fasta_dir) if f.endswith('.faa')]

    run_assign_predict_with_thresholds(
        fasta_files=fasta_files,
        feature_map=args.feature_map,
        best_hits_tsv=args.best_hits_tsv,
        phage_feature_table_path=args.phage_feature_table_path,
        output_dir=args.output_dir,
        model_dir=args.model_dir,
        thresholds=thresholds,
        select_protein_families=args.select_protein_families
    )

if __name__ == "__main__":
    main()
