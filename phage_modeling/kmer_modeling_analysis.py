import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from plotnine import *
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Load AA sequences from file
def load_aa_sequences(aa_sequence_file, feature_type='strain'):
    """
    Loads amino acid sequences from a FASTA file into a DataFrame.

    Parameters:
    aa_sequence_file (str): Path to the amino acid sequence file in FASTA format.

    Returns:
    DataFrame: A DataFrame with 'protein_ID' and 'sequence' columns.
    """
    logging.info(f"Loading amino acid sequences from {aa_sequence_file}")
    
    # Convert the iterator to a list to prevent re-consumption
    aa_records = list(SeqIO.parse(aa_sequence_file, 'fasta'))
    
    # Create the DataFrame with consistent-length lists
    aa_sequences_df = pd.DataFrame({
        'protein_ID': [record.id for record in aa_records],
        'sequence': [str(record.seq) for record in aa_records]
    })
    
    logging.info(f"Loaded {len(aa_sequences_df)} sequences.")
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
    filtered_kmers['protein_family'] = filtered_kmers['cluster'].str.split('_').str[:-1].str.join('_')
    # filtered_kmers['protein_ID'] = ['_'.join(x.split('_')[:-1]) for x in filtered_kmers['cluster']]
    # print(filtered_kmers.head())
    
    logging.info(f"Filtered down to {len(filtered_kmers)} predictive k-mers.")
    return filtered_kmers

# Merge kmers with protein families
def merge_kmers_with_families(protein_families_file, aa_sequences_df, feature_type='strain'):
    """
    Merges k-mer data with protein family information.

    Parameters:
    protein_families_file (str): Path to the protein families file in CSV format.
    aa_sequences_df (DataFrame): DataFrame containing amino acid sequences.
    feature_type (str): Type of feature, either 'strain' or other.

    Returns:
    DataFrame: A merged DataFrame with k-mers and protein family information.
    """
    logging.info(f"Merging k-mer data with protein families from {protein_families_file}")
    protein_families_df = pd.read_csv(protein_families_file)
    protein_families_df = protein_families_df[[feature_type, 'cluster', 'protein_ID']].drop_duplicates()
    protein_families_df.rename(columns={'cluster': 'protein_family'}, inplace=True)
    merged_df = protein_families_df.merge(aa_sequences_df, on='protein_ID', how='inner')
    # print(merged_df.head())
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
    # print(kmer_id_df.head())
    logging.info(f"Constructed k-mer ID DataFrame with {len(kmer_id_df)} entries.")
    return kmer_id_df

# Perform MSA and extract indices
def align_sequences(sequences, output_dir, family_name):
    """
    Aligns sequences within a protein family and removes excessive leading gaps.

    Parameters:
    sequences (list of tuples): List of (header, sequence) tuples.
    output_dir (str): Directory to save temporary files for alignment.
    family_name (str): Name of the protein family for unique file handling.

    Returns:
    DataFrame: DataFrame with 'protein_ID', 'aln_sequence', and 'start_index'.
    """
    logging.info(f"Aligning sequences for protein family: {family_name}")
    
    # Paths for temporary files
    temp_fasta_path = os.path.join(output_dir, f"{family_name}_temp_sequences.fasta")
    temp_aln_path = os.path.join(output_dir, f"{family_name}_temp_sequences.aln")

    # Create temporary IDs and map to original IDs
    seq_records = []
    temp_id_map = {}
    for i, (header, seq) in enumerate(sequences):
        temp_id = f"seq_{i}"
        record = SeqRecord(Seq(seq), id=temp_id, description="")
        seq_records.append(record)
        temp_id_map[temp_id] = header

    # Write sequences to the temporary FASTA file
    with open(temp_fasta_path, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")

    # Run ClustalW for alignment
    clustalw_cline = ClustalwCommandline("clustalw2", infile=temp_fasta_path)
    clustalw_cline()

    # Read alignment results
    alignment = AlignIO.read(temp_aln_path, "clustal")

    # Remove leading gaps from alignment
    aln_length = alignment.get_alignment_length()
    non_gap_positions = [any(record.seq[i] != '-' for record in alignment) for i in range(aln_length)]
    first_non_gap = non_gap_positions.index(True)
    
    # Construct DataFrame with trimmed alignment
    start_positions = {}
    aligned_sequences = {}

    for record in alignment:
        trimmed_seq = str(record.seq[first_non_gap:])  # Trim leading gaps
        original_id = temp_id_map[record.id]
        start_pos = trimmed_seq.find(trimmed_seq.lstrip('-'))
        start_positions[original_id] = start_pos
        aligned_sequences[original_id] = trimmed_seq

    # Construct the DataFrame for this protein family
    result_df = pd.DataFrame({
        'protein_ID': list(aligned_sequences.keys()),
        'aln_sequence': list(aligned_sequences.values()),
        'start_index': list(start_positions.values())
    })

    # Clean up temporary files
    os.remove(temp_fasta_path)
    os.remove(temp_aln_path)
    os.remove(os.path.join(output_dir, f"{family_name}_temp_sequences.dnd"))

    logging.info(f"Alignment completed for protein family: {family_name}")
    return result_df

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
    start_indices = matches
    stop_indices = [m + len(row['kmer']) - 1 for m in matches]  # End of each kmer match
    # logging.info(f"Found {len(matches)} matches for k-mer pattern in sequence.")
    return pd.Series([start_indices, stop_indices], index=['start_indices', 'stop_indices'])

# Coverage calculation
def calculate_coverage(df):
    """
    Calculates binary coverage for amino acids in aligned sequences with multiple k-mer matches.

    Parameters:
    df (DataFrame): DataFrame with aligned sequences and lists of k-mer start and stop positions.

    Returns:
    DataFrame: DataFrame with coverage information for each amino acid position.
    """
    logging.info("Calculating coverage for aligned sequences.")
    coverage_data = []
    for _, row in df.iterrows():
        coverage = np.full(len(row['aln_sequence']), 0)  # Initialize with 0 (absence)
        
        # Process each k-mer match for the sequence
        for start, stop in zip(row['start_indices'], row['stop_indices']):
            coverage[start:stop + 1] = 1  # Mark presence from start to stop index

        for idx, residue in enumerate(row['aln_sequence']):
            if residue != '-':  # Skip gap positions
                coverage_data.append({
                    'Feature': row['Feature'],
                    'protein_family': row['protein_family'],
                    'protein_ID': row['protein_ID'],
                    'AA_index': idx,
                    'Residue': residue,
                    'coverage': coverage[idx]
                })
    
    coverage_data = pd.DataFrame(coverage_data)
    coverage_data = coverage_data.groupby(['Feature', 'protein_family', 'protein_ID', 'AA_index', 'Residue'])['coverage'].max().reset_index()

    logging.info(f"Calculated coverage for {len(coverage_data)} amino acid positions.")
    return coverage_data

# Identify segments from binary coverage
def identify_segments(df):
    """
    Identifies contiguous coverage segments in amino acid sequences for each unique combination
    of Feature, protein family, and protein ID.

    Parameters:
    df (DataFrame): Input DataFrame containing the following columns:
        - Feature: Identifier for the feature category.
        - protein_family: Identifier for the protein family.
        - protein_ID: Identifier for the specific protein.
        - AA_index: Amino acid index within the protein sequence.
        - coverage: Binary indicator of coverage (1 for covered, 0 for uncovered).
    """
    logging.info("Identifying coverage segments.")
    
    print('Step 0')
    print(df.head())
    # Sort values to ensure segment detection aligns with amino acid sequence order
    df = df.sort_values(by=['Feature', 'protein_family', 'protein_ID', 'AA_index'])
    print('Step 1')
    print(df.head())
    
    # Calculate changes in coverage, indicating the start of new segments
    df['prev_coverage'] = df.groupby(['Feature', 'protein_family', 'protein_ID'])['coverage'].shift(1)
    print('Step 2')
    print(df.head())
    df['segment_change'] = (df['coverage'] != df['prev_coverage']) | (df['prev_coverage'].isna())
    print('Step 3')
    print(df.head())
    df['segment_id'] = df.groupby(['Feature', 'protein_family', 'protein_ID'])['segment_change'].cumsum()
    print('Step 4')
    print(df.head())
    
    # Aggregate to find start and stop indices of each segment
    segments_df = df.groupby(['Feature', 'protein_family', 'protein_ID', 'coverage', 'segment_id']).agg(
        start=('AA_index', 'min'),
        stop=('AA_index', 'max')
    ).reset_index()
    print('Step 5')
    print(segments_df.head())
    
    # Adjust stop index to be inclusive
    segments_df['stop'] += 1
    print('Step 6')
    print(segments_df.head())
    
    # Clean up by removing temporary columns and handling missing values
    segments_df = segments_df[['Feature', 'protein_family', 'protein_ID', 'coverage', 'segment_id', 'start', 'stop']]
    segments_df = segments_df.dropna()
    
    logging.info(f"Identified {len(segments_df)} segments.")
    return segments_df

def merge_no_coverage_proteins(coverage_segments_df, aligned_df):
    """
    Merges proteins with no coverage into the final segment summary.

    Parameters:
    coverage_segments_df (DataFrame): DataFrame with identified segments and coverage.
    aligned_df (DataFrame): DataFrame with aligned sequences and k-mer indices.

    Returns:
    DataFrame: The final coverage summary DataFrame, including proteins with no segments.
    """
    logging.info("Adding proteins with no coverage segments to the final summary.")
    
    # Identify proteins without any coverage segments
    no_coverage_proteins_df = aligned_df[~aligned_df['protein_ID'].isin(coverage_segments_df['protein_ID'].unique())]
    no_coverage_proteins_df = no_coverage_proteins_df[['protein_family', 'protein_ID', 'start_index', 'aln_sequence']].drop_duplicates()
    
    # Clean up and prepare no-coverage proteins for merging
    no_coverage_proteins_df['aln_sequence'] = no_coverage_proteins_df['aln_sequence'].str.strip('-')
    no_coverage_proteins_df['stop'] = no_coverage_proteins_df['aln_sequence'].str.len()
    no_coverage_proteins_df = no_coverage_proteins_df.rename(columns={'start_index': 'start'})
    no_coverage_proteins_df = no_coverage_proteins_df.drop(columns=['aln_sequence'])
    no_coverage_proteins_df['coverage'] = 0  # Set coverage to 0 for proteins with no segments

    # Concatenate the no-coverage proteins with the coverage segments
    final_coverage_summary_df = pd.concat([coverage_segments_df, no_coverage_proteins_df], ignore_index=True)
    logging.info(f"Final coverage summary includes {len(final_coverage_summary_df)} proteins.")

    return final_coverage_summary_df

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
            geom_segment(data=group[group['coverage'] == 0], mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID'), color='grey', size=5) +
            geom_segment(data=group[group['coverage'] == 1], mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID', color='Feature'), size=5) +
            labs(title=f'Protein Family: {family}', x='AA Index', y='protein_ID') +
            theme(axis_text_x=element_text(rotation=90), panel_background=element_rect(fill='white'))
        )
        plot.save(f"{output_dir}/{family}_coverage_plot.png")
        logging.info(f"Saved plot for protein family {family} at {output_dir}")
