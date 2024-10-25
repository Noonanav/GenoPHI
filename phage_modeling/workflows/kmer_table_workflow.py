import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
from phage_modeling.mmseqs2_clustering import merge_feature_tables
from phage_modeling.workflows.select_and_model_workflow import run_modeling_workflow_from_feature_table

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

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

def construct_feature_table(fasta_file, protein_csv, k, id_col, one_gene, output_dir, output_name, k_range=False):
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

    # Updated k-mer filtering for unique counts
    logging.info("Aggregating k-mer data by feature and cluster.")
    full_kmer_df_counts = full_kmer_df_merged.copy()
    full_kmer_df_counts['counts'] = full_kmer_df_counts.groupby(['kmer', 'Feature', 'cluster'])['protein_ID'].transform('nunique')
    full_kmer_df_counts['k_length'] = [len(k) for k in full_kmer_df_counts['kmer']]
    
    # Filter unique k-mers based on 'one_gene' parameter
    if one_gene:
        logging.info("Including features with 1 gene.")
        full_kmer_df_counts = full_kmer_df_counts[full_kmer_df_counts['counts'] >= 1]
    else:
        logging.info("Filtering out features with only 1 gene.")
        full_kmer_df_counts = full_kmer_df_counts[full_kmer_df_counts['counts'] > 1]

    # Create kmer_id by concatenating cluster and kmer, and filter by k length
    logging.info("Creating kmer IDs and filtering by k length.")
    full_kmer_df_counts['kmer_id'] = full_kmer_df_counts['cluster'] + '_' + full_kmer_df_counts['kmer']
    full_kmer_df_counts = full_kmer_df_counts[full_kmer_df_counts['k_length'] == k].drop_duplicates()

    # Construct presence-absence matrix
    logging.info("Constructing final presence-absence matrix.")
    feature_table = full_kmer_df_counts[[id_col, 'kmer_id']].copy()
    feature_table['presence'] = 1
    feature_table = feature_table.pivot_table(index=id_col, columns='kmer_id', values='presence', fill_value=0).reset_index()

    # Save the presence-absence matrix as a CSV
    feature_table_output = os.path.join(output_dir, f"{output_name}_feature_table.csv")
    feature_table.to_csv(feature_table_output, index=False)
    logging.info(f"Presence-absence matrix saved to {feature_table_output}")

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

def feature_selection_optimized(presence_absence, source, genome_column_name, output_dir, prefix=None):
    """Optimizes feature selection by identifying perfect co-occurrence of features."""
    logging.info("Optimizing feature selection...")
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
    
    selected_features = pd.DataFrame(data, columns=["Feature", "Cluster_Label"])
    selected_features['protein_family'] = selected_features['Cluster_Label'].apply(lambda x: '_'.join(x.split('_')[0:-1]))

    if prefix:
        selected_features_output = os.path.join(output_dir, f"{prefix}_selected_features.csv")
    else:
        selected_features_output = os.path.join(output_dir, "selected_features.csv")
    selected_features.to_csv(selected_features_output, index=False)
    logging.info(f"Selected features saved to {selected_features_output}")

    return selected_features

def feature_assignment(genome_assignments, selected_features, genome_column_name, output_dir, prefix=None):
    """Assigns features to genomes based on the selected features."""
    logging.info("Assigning features to genomes...")
    assignment_df = genome_assignments.merge(selected_features, on="Cluster_Label", how="inner")
    assignment_df = assignment_df.drop(columns=["Cluster_Label"]).drop_duplicates()

    # Create feature table in wide format
    feature_table = assignment_df.pivot_table(index=genome_column_name, columns="Feature", aggfunc="size", fill_value=0).reset_index()

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

def run_kmer_table_workflow(strain_fasta, protein_csv, k, id_col, one_gene, output_dir, k_range=False,
                      phenotype_matrix=None, phage_fasta=None, protein_csv_phage=None, remove_suffix=False, 
                      sample_column='strain', phenotype_column='interaction', modeling=False, threads=4):
    """
    Runs the full workflow for generating k-mer feature tables, filtering, optimizing, 
    and optionally merging with phenotype and phage feature tables.
    """
    logging.info("Running the full k-mer feature table workflow...")

    feature_output_dir = os.path.join(output_dir, "feature_tables")
    if not os.path.exists(feature_output_dir):
        os.makedirs(feature_output_dir)

    # Step 1: Construct strain feature table
    strain_feature_table_path = construct_feature_table(strain_fasta, protein_csv, k, id_col, one_gene, feature_output_dir, "strain", k_range)
    strain_feature_table = pd.read_csv(strain_feature_table_path)

    # Step 2: Get genome assignments for strain
    genome_assignments = get_genome_assignments_tables(strain_feature_table, id_col, feature_output_dir)

    # Step 3: Optimize feature selection for strain
    selected_features = feature_selection_optimized(strain_feature_table, "selected", id_col, feature_output_dir)

    # Step 4: Assign features to genomes and generate final feature table for strain
    assignment_df, final_feature_table, final_feature_table_output = feature_assignment(genome_assignments, selected_features, id_col, feature_output_dir)

    # Optionally, construct phage feature table if phage_fasta is provided
    phage_feature_table_path = None
    if phage_fasta:
        # Use protein_csv_phage if provided; otherwise, fall back to protein_csv
        phage_protein_csv = protein_csv_phage if protein_csv_phage else protein_csv
        phage_feature_table_path = construct_feature_table(phage_fasta, phage_protein_csv, k, 'phage', one_gene, feature_output_dir, "phage", k_range)
        phage_feature_table = pd.read_csv(phage_feature_table_path)

        # Phage genome assignments
        phage_genome_assignments = get_genome_assignments_tables(phage_feature_table, 'phage', feature_output_dir, prefix='phage')

        # Optimize feature selection for phage
        phage_selected_features = feature_selection_optimized(phage_feature_table, "phage_selected", 'phage', feature_output_dir, prefix='phage')

        # Phage feature assignment
        phage_assignment_df, phage_final_feature_table, phage_final_feature_table_output = feature_assignment(phage_genome_assignments, phage_selected_features, 'phage', feature_output_dir, prefix='phage')

    # Step 5: Merge feature tables if phenotype_matrix is provided
    if phenotype_matrix:
        logging.info("Merging with phenotype and optional phage feature tables.")
        merged_table_path = merge_feature_tables(
            strain_features=final_feature_table_output,
            phenotype_matrix=phenotype_matrix,
            output_dir=output_dir,
            sample_column=sample_column,
            phage_features=phage_final_feature_table_output,
            remove_suffix=remove_suffix,
            output_file="merged_features"
        )
        logging.info(f"Merged feature table saved to {merged_table_path}")
    else:
        logging.info("Skipping merge as no phenotype matrix provided.")

    logging.info("Feature table construction workflow completed successfully.")

    if modeling and merged_table_path:
        logging.info("Running modeling workflow...")
        modeling_output_dir = os.path.join(output_dir, "modeling")
        if not os.path.exists(modeling_output_dir):
            os.makedirs(modeling_output_dir)
        # Run the modeling workflow using the merged feature tabl
        # Assuming `final_feature_table_output` is the path to the feature table from the k-mer workflow
        run_modeling_workflow_from_feature_table(
            full_feature_table=merged_table_path,
            output_dir=modeling_output_dir,
            threads=threads,
            num_features=100,
            filter_type='none',
            num_runs_fs=10,
            num_runs_modeling=20,
            sample_column=sample_column,
            phenotype_column=phenotype_column,
            method='rfe'
        )

    else:
        logging.warning("Modeling step skipped due to missing merged table or phenotype matrix.")

# Command-line interface

def main():
    """
    Command-line interface for generating k-mer feature tables and running the full workflow.
    """
    parser = argparse.ArgumentParser(prog='AA_seq_to_kmer', description='Generates k-mer feature tables and final presence-absence matrix from AA sequences.')
    parser.add_argument('-i', '--strain_fasta', required=True, help="The path to the FASTA file containing strain AA sequences.")
    parser.add_argument('-ip', '--phage_fasta', default=None, help="The path to the FASTA file containing phage AA sequences.")
    parser.add_argument('-p', '--protein_csv', required=True, help="The path to the CSV file containing protein feature data for strain.")
    parser.add_argument('--protein_csv_phage', default=None, help="The path to the CSV file containing protein feature data for phage.")
    parser.add_argument('--k', type=int, required=True, help="k-mer length.")
    parser.add_argument('--id_col', default="strain", help="The column name for genome ID column.")
    parser.add_argument('--one_gene', action='store_true', help="Include features with 1 gene.")
    parser.add_argument('--k_range', action='store_true', help="Generate range of k-mer lengths from 3 to k.")
    parser.add_argument('-o', '--output_dir', required=True, help="Directory to save all output files.")
    parser.add_argument('--phenotype_matrix', default=None, help="Path to the phenotype matrix CSV for merging.")
    parser.add_argument('--remove_suffix', action='store_true', help="Remove suffix from genome names when merging.")
    parser.add_argument('--sample_column', default='strain', help="Sample identifier column name (default is 'strain').")
    parser.add_argument('--modeling', action='store_true', help="Run modeling workflow after feature table generation.")
    parser.add_argument('--phenotype_column', default='interaction', help="Phenotype column name in the phenotype matrix.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads to use (default: 4).")

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
        threads=args.threads
    )

if __name__ == "__main__":
    main()