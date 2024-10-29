import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from plotnine import *
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Load AA sequences from file
def load_aa_sequences(aa_sequence_file):
    """
    Loads amino acid sequences from a FASTA file into a DataFrame.

    Parameters:
    aa_sequence_file (str): Path to the amino acid sequence file in FASTA format.

    Returns:
    DataFrame: A DataFrame with 'protein_ID' and 'sequence' columns.
    """
    logging.info(f"Loading amino acid sequences from {aa_sequence_file}")
    aa_records = SeqIO.parse(aa_sequence_file, 'fasta')
    aa_sequences_df = pd.DataFrame({
        'protein_ID': [record.id for record in aa_records],
        'sequence': [str(record.seq) for record in SeqIO.parse(aa_sequence_file, 'fasta')]
    })
    logging.info(f"Loaded {len(aa_sequences_df)} sequences.")
    print(aa_sequences_df.head())
    return aa_sequences_df

# Get predictive features based on host or phage
def get_predictive_kmers(feature_file_path, feature2cluster_path, feature_type):
    """
    Filters predictive features based on the specified feature type.

    Parameters:
    feature_file_path (str): Path to the feature CSV file.
    feature2cluster_path (str): Path to the feature-to-cluster mapping CSV file.
    feature_type (str): Either 'host' or 'phage' to indicate feature type.

    Returns:
    DataFrame: A DataFrame containing filtered k-mers with 'kmer' and 'protein_ID' columns.
    """
    logging.info(f"Extracting predictive features of type '{feature_type}' from {feature_file_path}")
    feature_df = pd.read_csv(feature_file_path)
    feature_prefix = feature_type[0] + 'c_'
    select_features = [col for col in feature_df.columns if feature_prefix in col]

    feature2cluster_df = pd.read_csv(feature2cluster_path)
    feature2cluster_df.rename(columns={'Cluster_Label': 'cluster'}, inplace=True)
    filtered_kmers = feature2cluster_df[feature2cluster_df['Feature'].isin(select_features)]
    filtered_kmers['kmer'] = filtered_kmers['cluster'].str.split('_').str[-1]
    filtered_kmers['protein_ID'] = ['_'.join(x.split('_')[:-1]) for x in filtered_kmers['cluster']]
    print(filtered_kmers.head())
    
    logging.info(f"Filtered down to {len(filtered_kmers)} predictive k-mers.")
    return filtered_kmers

# Merge kmers with protein families
def merge_kmers_with_families(kmer_df, protein_families_file, aa_sequences_df, feature_type='strain'):
    """
    Merges k-mer data with protein family information.

    Parameters:
    kmer_df (DataFrame): DataFrame containing k-mers and corresponding proteins.
    protein_families_file (str): Path to the protein families file in CSV format.
    aa_sequences_df (DataFrame): DataFrame containing amino acid sequences.
    feature_type (str): Type of feature, either 'strain' or other.

    Returns:
    DataFrame: A merged DataFrame with k-mers and protein family information.
    """
    logging.info(f"Merging k-mer data with protein families from {protein_families_file}")
    protein_families_df = pd.read_csv(protein_families_file)
    protein_families_df = protein_families_df[[feature_type, 'cluster', 'protein_ID']]
    protein_families_df.rename(columns={'cluster': 'protein_family'}, inplace=True)
    merged_df = protein_families_df.merge(aa_sequences_df, on='protein_ID', how='inner')
    print(merged_df.head())
    logging.info(f"Merged k-mer data with {len(merged_df)} protein family entries.")
    return merged_df

# Construct kmer ID DataFrame for alignment
def construct_kmer_id_df(protein_families_df, kmer_df):
    """
    Constructs a DataFrame of k-mers for each protein family.

    Parameters:
    protein_families_df (DataFrame): DataFrame of protein families.
    kmer_df (DataFrame): DataFrame of k-mers.

    Returns:
    DataFrame: A DataFrame of k-mers with associated protein families.
    """
    logging.info("Constructing k-mer ID DataFrame for alignment.")
    kmer_id_df = pd.DataFrame()
    for protein in kmer_df['protein_ID'].unique():
        family_id = protein_families_df.loc[protein_families_df['protein_ID'] == protein, 'protein_family'].values[0]
        family_df = protein_families_df[protein_families_df['protein_family'] == family_id]
        family_df['kmer_cluster'] = protein
        kmer_id_df = pd.concat([kmer_id_df, family_df])
    print(kmer_id_df.head())
    logging.info(f"Constructed k-mer ID DataFrame with {len(kmer_id_df)} entries.")
    return kmer_id_df

# Perform MSA and extract indices
def align_sequences(sequences):
    """
    Performs multiple sequence alignment and extracts alignment start indices.

    Parameters:
    sequences (list of tuples): List of (header, sequence) tuples.

    Returns:
    DataFrame: DataFrame with 'protein_ID', 'aln_sequence', and 'start_index'.
    """
    logging.info("Performing multiple sequence alignment.")
    with open("temp_sequences.fasta", "w") as f:
        for header, seq in sequences:
            f.write(f">{header}\n{seq}\n")
    clustalw_cline = ClustalwCommandline("clustalw2", infile="temp_sequences.fasta")
    clustalw_cline()
    alignment = AlignIO.read("temp_sequences.aln", "clustal")

    result = []
    for record in alignment:
        seq_str = str(record.seq)
        start_index = seq_str.find(seq_str.lstrip('-'))
        result.append((record.id, seq_str, start_index))
    logging.info("Alignment and index extraction completed.")
    return pd.DataFrame(result, columns=['protein_ID', 'aln_sequence', 'start_index'])

# Find kmer indices within aligned sequences
def find_kmer_indices(row):
    """
    Finds indices of k-mer within an aligned sequence.

    Parameters:
    row (Series): Row containing 'kmer' and 'aln_sequence'.

    Returns:
    Series: Start and stop indices of k-mer in aligned sequence.
    """
    kmer_pattern = '-*'.join(row['kmer'])
    seq = row['aln_sequence']
    matches = [match.start() for match in re.finditer(f'(?={kmer_pattern})', seq)]
    logging.info(f"Found {len(matches)} matches for k-mer pattern in sequence.")
    return pd.Series([matches, [m + len(row['kmer']) for m in matches]])

# Coverage calculation
def calculate_coverage(df):
    """
    Calculates binary coverage for amino acids in aligned sequences.

    Parameters:
    df (DataFrame): DataFrame with aligned sequences and k-mer positions.

    Returns:
    DataFrame: DataFrame with coverage information for each amino acid position.
    """
    logging.info("Calculating coverage for aligned sequences.")
    coverage_data = []
    for _, row in df.iterrows():
        coverage = np.full(len(row['aln_sequence']), np.nan)
        if not pd.isna(row['start_index']) and not pd.isna(row['stop_index']):
            coverage[int(row['start_index']):int(row['stop_index']) + 1] = 1
        for idx, residue in enumerate(row['aln_sequence']):
            if residue != '-':
                coverage_data.append({'protein_family': row['protein_family'], 'protein_ID': row['protein_ID'],
                                      'AA_index': idx, 'Residue': residue, 'coverage': coverage[idx] if not pd.isna(coverage[idx]) else 0})
    logging.info(f"Calculated coverage for {len(coverage_data)} amino acid positions.")
    print(pd.DataFrame(coverage_data).head())
    return pd.DataFrame(coverage_data)

# Identify segments from binary coverage
def identify_segments(df):
    """
    Identifies contiguous coverage segments in amino acid sequences.

    Parameters:
    df (DataFrame): DataFrame with binary coverage information.

    Returns:
    DataFrame: DataFrame with segment start and stop positions.
    """
    logging.info("Identifying coverage segments.")
    df = df.sort_values(by=['protein_family', 'protein_ID', 'AA_index'])
    df['segment_change'] = (df['coverage'] != df['coverage'].shift(1)).cumsum()
    segments_df = df.groupby(['protein_family', 'protein_ID', 'coverage', 'segment_change']).agg(
        start=('AA_index', 'min'), stop=('AA_index', 'max')).reset_index()
    segments_df['stop'] += 1  # Inclusive stop
    logging.info(f"Identified {len(segments_df)} segments.")
    print(segments_df.head())
    return segments_df.drop_duplicates()

# Plot segments
def plot_segments(segment_summary_df, output_dir):
    """
    Plots segments by protein family and saves the plots.

    Parameters:
    segment_summary_df (DataFrame): DataFrame with segment summary data.
    output_dir (str): Directory to save output plots.

    Returns:
    None
    """
    logging.info(f"Plotting segments to {output_dir}.")
    os.makedirs(output_dir, exist_ok=True)
    for family, group in segment_summary_df.groupby('protein_family'):
        plot = (
            ggplot() +
            geom_segment(data=group[group['coverage'] == 1], mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID'), color='green', size=5) +
            geom_segment(data=group[group['coverage'] == 0], mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID'), color='grey', size=5) +
            labs(title=f'Protein Family: {family}', x='AA Index', y='protein_ID') +
            theme(axis_text_x=element_text(rotation=90), panel_background=element_rect(fill='white'))
        )
        plot.save(f"{output_dir}/{family}_coverage_plot.png")
        logging.info(f"Saved plot for protein family {family} at {output_dir}")
