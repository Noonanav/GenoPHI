import os
import pandas as pd
import logging
import subprocess
from Bio import SeqIO
from argparse import ArgumentParser
from phage_modeling.mmseqs2_clustering import create_mmseqs_database, load_strains, create_contig_to_genome_dict, select_best_hits
from phage_modeling.workflows.prediction_workflow import generate_full_feature_table, predict_interactions, calculate_mean_predictions, load_model

def load_sequences(fasta_file, protein_ids):
    """Load sequences from a FASTA file for given protein IDs."""
    sequences = {}
    if not os.path.exists(fasta_file):
        logging.warning(f"FASTA file not found: {fasta_file}")
        return sequences

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in protein_ids:
            sequences[record.id] = str(record.seq)
    return sequences

def map_features(best_hits_tsv, kmer_map_path, fasta_file, genome_contig_mapping, genome, genome_type):
    """
    Maps features to genomes based on k-mer matches using protein family as an intermediate step.
    
    Args:
        best_hits_tsv (str): Path to the best hits TSV file.
        kmer_map_path (str): Path to the filtered k-mers CSV file.
        fasta_file (str): Path to the FASTA file for the genome.
        genome_contig_mapping (dict): Dictionary mapping contigs to genomes.
        genome (str): Name of the genome being processed.
        genome_type (str): Type of genome ('strain' or 'phage').
    
    Returns:
        pd.DataFrame: A DataFrame with k-mer feature presence for the genome.
    """
    if not os.path.exists(best_hits_tsv):
        logging.error(f"Best hits TSV file does not exist: {best_hits_tsv}")
        return None

    try:
        # Load best hits data
        best_hits_df = pd.read_csv(best_hits_tsv, sep='\t', header=None, names=['Query', 'Cluster'])
        
        # Load k-mer data
        kmer_data = pd.read_csv(kmer_map_path)
        
        # Ensure consistent data types
        best_hits_df['Cluster'] = best_hits_df['Cluster'].astype(str)
        kmer_data['cluster'] = kmer_data['cluster'].astype(str)

    except Exception as e:
        logging.error(f"Error reading input files: {e}")
        return None

    logging.info(f"Mapping features for {genome_type} '{genome}'")

    # Merge with genome mapping to get genome information
    genome_contig_mapping_df = pd.DataFrame.from_dict(genome_contig_mapping, orient='index', columns=['contig_id']).reset_index().rename(columns={'index': genome_type})
    merged_df = best_hits_df.merge(genome_contig_mapping_df, left_on='Query', right_on='contig_id', how='left')

    # Load sequences for proteins in this genome
    protein_ids = merged_df['Query'].unique().tolist()
    sequences = load_sequences(fasta_file, protein_ids)

    # Get all genomes from the mapping
    all_genomes = list(genome_contig_mapping.keys())

    # Assign k-mer features
    kmer_features = []
    for _, row in merged_df.iterrows():
        sequence = sequences.get(row['Query'], "")
        for _, kmer_row in kmer_data[kmer_data['cluster'] == row['Cluster']].iterrows():
            kmer_features.append({
            genome_type: row[genome_type],
            kmer_row['Feature']: int(kmer_row['kmer'] in sequence)
        })

    kmer_features_df = pd.DataFrame(kmer_features)
    logging.info(f"Columns in kmer_features_df: {kmer_features_df.columns}")
    if genome_type not in kmer_features_df.columns:
        logging.error(f"Column '{genome_type}' not found in kmer_features_df. Available columns are: {kmer_features_df.columns.tolist()}")

    # Create base DataFrame with all genomes and all k-mer features set to 0
    all_kmer_features = kmer_data['Feature'].unique()
    base_df = pd.DataFrame({
        genome_type: all_genomes
    })

    for feature in all_kmer_features:
        base_df[feature] = 0

    # Update base_df with actual matches
    feature_presence = base_df.merge(
        kmer_features_df.groupby(genome_type).max().reset_index(),
        on=genome_type,
        how='left'
    ).fillna(0).astype(int)

    return feature_presence

def assign_sequences_to_existing_clusters(query_db, target_db, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv):
    """
    Assigns sequences from the query database to clusters in the existing target database.

    Args:
        query_db (str): Path to the query MMseqs2 database (e.g., validation strains).
        target_db (str): Path to the target MMseqs2 database (e.g., strain feature database).
        output_dir (str): Directory to save results.
        tmp_dir (str): Temporary directory for intermediate files.
        coverage (float): Minimum coverage for assignment.
        min_seq_id (float): Minimum sequence identity for assignment.
        sensitivity (float): Sensitivity for assignment.
        threads (int): Number of threads to use.
        clusters_tsv (str): Path to the clusters TSV file.
    """
    logging.info("Assigning query sequences to existing target clusters...")

    result_db = os.path.join(tmp_dir, "result_db")
    # Run mmseqs2 search from query_db to target_db
    search_command = (
        f"mmseqs search {query_db} {target_db} {result_db} {tmp_dir} "
        f"-c {coverage} --min-seq-id {min_seq_id} -s {sensitivity} "
        f"--threads {threads} -v 3"
    )
    try:
        subprocess.run(search_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during MMseqs search: {e}")
        return None  # Early exit on error

    logging.info("Sequence assignment to existing clusters completed successfully.")

    assigned_tsv = os.path.join(output_dir, "assigned_clusters.tsv")
    # Create a TSV from query_db to target_db results
    createtsv_command = f"mmseqs createtsv {query_db} {target_db} {result_db} {assigned_tsv} --threads {threads} -v 3"
    
    try:
        subprocess.run(createtsv_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error during MMseqs createtsv: {e}")
        return None  # Early exit on error
        
    logging.info(f"Assigned clusters TSV saved to {assigned_tsv}")

    best_hits_tsv = os.path.join(output_dir, "best_hits.tsv")
    select_best_hits(assigned_tsv, best_hits_tsv, clusters_tsv)

    return best_hits_tsv
    
def process_new_genomes(input_dir, mmseqs_db, suffix, tmp_dir, output_dir, clusters_tsv, genome_type, genomes=None, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4, kmer_map_path=None):
    os.makedirs(output_dir, exist_ok=True)
    feature_table_dir = os.path.join(output_dir, "feature_tables")
    os.makedirs(feature_table_dir, exist_ok=True)

    all_genomes_feature_tables = []

    # Create a single query database for all genomes
    all_genomes_query_db = os.path.join(tmp_dir, 'all_genomes_query_db')
    if genomes is None:
        genomes = ['.'.join(f.split('.')[:-1]) for f in os.listdir(input_dir) if f.endswith(suffix)]
    
    fasta_files = create_mmseqs_database(input_dir, all_genomes_query_db, suffix, 'directory', genomes, threads)

    # Single search for all genomes
    best_hits_tsv = assign_sequences_to_existing_clusters(all_genomes_query_db, mmseqs_db, output_dir, tmp_dir, coverage, min_seq_id, sensitivity, threads, clusters_tsv)

    for genome in genomes:
        genome_tmp_dir = os.path.join(tmp_dir, genome)
        os.makedirs(genome_tmp_dir, exist_ok=True)

        logging.info(f"Processing {genome_type} {genome}...")

        try:
            # Here you would process the results from best_hits_tsv for each genome
            genome_contig_mapping, _ = create_contig_to_genome_dict([os.path.join(input_dir, f"{genome}.{suffix}")], 'file')
            genome_fasta_file = os.path.join(input_dir, f"{genome}.{suffix}")
            feature_presence = map_features(best_hits_tsv, kmer_map_path, genome_fasta_file, genome_contig_mapping, genome, genome_type)

            if feature_presence is not None:
                output_path = os.path.join(feature_table_dir, f'{genome}_feature_table.csv')
                feature_presence.to_csv(output_path, index=False)
                logging.info(f"Feature table for {genome_type} '{genome}' saved to {output_path}")
                all_genomes_feature_tables.append(feature_presence)

        except FileNotFoundError as e:
            logging.error(f"Error processing {genome_type} {genome}: {e}. Skipping...")
            continue  
        except Exception as e:
            logging.error(f"Unexpected error while processing {genome_type} {genome}: {str(e)}", exc_info=True)
            continue

    logging.info("Genome processing complete.")

    if all_genomes_feature_tables:
        combined_feature_table = pd.concat(all_genomes_feature_tables)
        combined_output_path = os.path.join(output_dir, 'combined_feature_table.csv')
        combined_feature_table.to_csv(combined_output_path, index=False)
        logging.info(f"Combined feature table saved to {combined_output_path}")

def run_assign_and_predict_workflow(input_dir, mmseqs_db, clusters_tsv, tmp_dir, suffix, genome_list, genome_type, genome_column, model_dir, phage_feature_table, output_dir, sensitivity=7.5, coverage=0.8, min_seq_id=0.6, threads=4, kmer_map_path=None):
    """
    Full workflow to assign features to genomes and predict interactions.
    """
    # Make output and tmp dirs
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    if genome_column is None:
        genome_column = genome_type

    if genome_list and os.path.exists(genome_list):
        genomes = load_strains(genome_list, genome_column)
    else:
        genomes = None

    process_new_genomes(
        input_dir=input_dir,
        genomes=genomes,
        mmseqs_db=mmseqs_db,
        suffix=suffix,
        tmp_dir=tmp_dir,
        output_dir=output_dir,
        clusters_tsv=clusters_tsv,
        genome_type=genome_type,
        sensitivity=sensitivity,
        coverage=coverage,
        min_seq_id=min_seq_id,
        threads=threads,
        kmer_map_path=kmer_map_path
    )

    # Prediction phase
    full_predictions_df = pd.DataFrame()
    phage_feature_table = pd.read_csv(phage_feature_table)

    strain_files = [f for f in os.listdir(output_dir) if f.endswith('_feature_table.csv')]
    for strain_file in strain_files:
        strain_table = pd.read_csv(os.path.join(output_dir, strain_file))
        prediction_feature_table = generate_full_feature_table(strain_table, phage_feature_table)

        all_predictions_df = predict_interactions(model_dir, prediction_feature_table)
        full_predictions_df = pd.concat([full_predictions_df, all_predictions_df])

    # Save predictions
    full_predictions_df.to_csv(os.path.join(output_dir, 'all_predictions.csv'), index=False)
    mean_predictions_df = calculate_mean_predictions(full_predictions_df)
    mean_predictions_df.to_csv(os.path.join(output_dir, 'mean_predictions.csv'), index=False)


# CLI main function
def main():
    parser = ArgumentParser(description="Assign genomic features and predict interactions.")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory containing genome FASTA files.")
    parser.add_argument('--genome_list', type=str, help="CSV file with genome names.")
    parser.add_argument('--genome_type', type=str, choices=['strain', 'phage'], default='strain', help="Type of genome to process.")
    parser.add_argument('--genome_column', type=str, default='strain', help="Column name for genome identifiers in genome_list.")
    parser.add_argument('--mmseqs_db', type=str, required=True, help="Path to the MMseqs2 database.")
    parser.add_argument('--clusters_tsv', type=str, required=True, help="Path to the clusters TSV file.")
    parser.add_argument('--suffix', type=str, default="faa", help="Suffix for FASTA files.")
    parser.add_argument('--model_dir', type=str, required=True, help="Directory with trained models.")
    parser.add_argument('--phage_feature_table', type=str, required=True, help="Path to the phage feature table.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to save results.")
    parser.add_argument('--tmp_dir', type=str, required=True, help="Temporary directory for intermediate files.")
    parser.add_argument('--sensitivity', type=float, default=7.5, help="Sensitivity for MMseqs2 search.")
    parser.add_argument('--coverage', type=float, default=0.8, help="Minimum coverage for assignment.")
    parser.add_argument('--min_seq_id', type=float, default=0.4, help="Minimum sequence identity for assignment.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for MMseqs2.")
    parser.add_argument('--kmer_map', type=str, default=None, help="Path to the filtered k-mers CSV file.")

    args = parser.parse_args()

    run_assign_and_predict_workflow(
        input_dir=args.input_dir,
        genome_list=args.genome_list,
        genome_type=args.genome_type,
        genome_column=args.genome_column,
        mmseqs_db=args.mmseqs_db,
        clusters_tsv=args.clusters_tsv,
        tmp_dir=args.tmp_dir,
        suffix=args.suffix,
        model_dir=args.model_dir,
        phage_feature_table=args.phage_feature_table,
        output_dir=args.output_dir,
        sensitivity=args.sensitivity,
        coverage=args.coverage,
        min_seq_id=args.min_seq_id,
        threads=args.threads,
        kmer_map_path=args.kmer_map
    )


if __name__ == "__main__":
    main()
