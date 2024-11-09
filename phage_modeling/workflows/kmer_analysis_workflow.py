import argparse
import pandas as pd
import logging
import os
from phage_modeling.kmer_modeling_analysis import (
    load_aa_sequences,
    get_predictive_kmers,
    merge_kmers_with_families,
    construct_kmer_id_df,
    align_sequences,
    find_kmer_indices,
    calculate_coverage,
    identify_segments,
    merge_no_coverage_proteins,
    plot_segments,
)

def kmer_analysis_workflow(
    aa_sequence_file,
    feature_file_path,
    feature2cluster_path,
    protein_families_file,
    output_dir,
    feature_type='strain',
    annotation_file=None,
):
    """
    Complete workflow function to load amino acid sequences, extract and align predictive kmers,
    calculate coverage, identify segments, and plot segment mappings on sequences.

    Args:
        aa_sequence_file (str): Path to the amino acid sequences file in FASTA format.
        feature_file_path (str): Path to the feature selection output file (CSV).
        feature2cluster_path (str): Path to the feature-to-cluster mapping file (CSV).
        protein_families_file (str): Path to the protein families file (CSV).
        output_dir (str): Directory where output files (plots, coverage summaries) will be saved.
        feature_type (str): Type of features to extract ('host' or 'phage').
        annotation_file (str, optional): Optional path to annotations file for enhanced plotting.

    Outputs:
        Segment plots and coverage summaries in the specified output directory.
    """
    # Step 1: Load amino acid sequences
    aa_sequences_df = load_aa_sequences(aa_sequence_file, feature_type)
    aa_sequences_df.to_csv(os.path.join(output_dir, 'aa_sequences_df.csv'), index=False)
    print('Printing aa_sequences_df')
    print(aa_sequences_df.head())
    
    # Step 2: Get predictive kmers
    filtered_kmers = get_predictive_kmers(feature_file_path, feature2cluster_path, feature_type)
    filtered_kmers.to_csv(os.path.join(output_dir, 'filtered_kmers.csv'), index=False)
    print('Printing filtered_kmers')
    print(filtered_kmers.head())
    
    # Step 3: Merge with protein families
    protein_families_df = merge_kmers_with_families(protein_families_file, aa_sequences_df, feature_type)
    protein_families_df.to_csv(os.path.join(output_dir, 'protein_families_df.csv'), index=False)
    print('Printing protein_families_df')
    print(protein_families_df.head())


    kmer_full_df = filtered_kmers.merge(protein_families_df, on='protein_family', how='inner')
    # kmer_full_df.to_csv(os.path.join(output_dir, 'kmer_full_df.csv'), index=False)
    print('Printing kmer_full_df')
    print(kmer_full_df.head())
    
    # Step 4: Construct kmer ID DataFrame
    kmer_id_df = construct_kmer_id_df(protein_families_df, kmer_full_df)
    print('Printing kmer_id_df')
    print(kmer_id_df.head())

    # Step 5: Align sequences by protein family
    aligned_df = pd.DataFrame()
    for family_name, family_group in protein_families_df.groupby("protein_family"):
        seqs_for_family = [(row["protein_ID"], row["sequence"]) for _, row in family_group.iterrows()]
        
        # Align and process each family separately
        aligned_family_df = align_sequences(seqs_for_family, output_dir, family_name)
        
        # Merge alignment results with the main kmer data
        aligned_family_df['protein_family'] = family_name
        # aligned_family_df = aligned_family_df.merge(kmer_full_df[['Feature', 'cluster', 'protein_ID', 'kmer']],
        #                                             on='protein_ID', how='inner')
        aligned_df = pd.concat([aligned_df, aligned_family_df], ignore_index=True)

    print('Printing aligned_df')
    print(aligned_df.head())
    aligned_df = aligned_df.merge(kmer_full_df[['Feature', 'cluster', 'protein_ID', 'kmer']], on='protein_ID', how='inner')
    aligned_df[['start_indices', 'stop_indices']] = aligned_df.apply(find_kmer_indices, axis=1)
    print('Printing aligned_df')
    print(aligned_df.head())
    
    # Step 6: Calculate coverage and identify segments
    coverage_summary = calculate_coverage(aligned_df)
    coverage_summary = coverage_summary.drop_duplicates()
    segments_df = identify_segments(coverage_summary)
    print('Printing segments_df')
    print(segments_df.head())
    segments_df.to_csv(os.path.join(output_dir, 'segments_df.csv'), index=False)

    # Step 6.1: Merge proteins without coverage segments
    # Adding proteins with no coverage segments to ensure complete representation
    final_segments_df = merge_no_coverage_proteins(segments_df, aligned_df)
    
    # Step 7: Plot segments with optional annotations
    if annotation_file:
        logging.info(f"Merging protein_families_df with annotation data from {annotation_file}")
        annotation_df = pd.read_csv(annotation_file)
        protein_families_annotated_df = protein_families_df.merge(annotation_df, on='protein_ID', how='left')
        protein_families_annotated_df.to_csv(os.path.join(output_dir, 'protein_families_annotated.csv'), index=False)
    
    logging.info(f"Plotting segments in {output_dir}")
    plot_segments(final_segments_df, output_dir)

def main():
    parser = argparse.ArgumentParser(
        description="Workflow script for analyzing predictive kmers in amino acid sequences."
    )
    parser.add_argument("--aa_sequence_file", required=True, help="Path to the AA sequence FASTA file.")
    parser.add_argument("--feature_file_path", required=True, help="Path to the feature selection output CSV file.")
    parser.add_argument("--feature2cluster_path", required=True, help="Path to the feature-to-cluster mapping CSV file.")
    parser.add_argument("--protein_families_file", required=True, help="Path to the protein families CSV file.")
    parser.add_argument("--output_dir", required=True, help="Directory for output files.")
    parser.add_argument("--feature_type", default="strain", help="Feature type to analyze.")
    parser.add_argument("--annotation_file", help="Path to optional annotation CSV file for plots.")

    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Run workflow
    kmer_analysis_workflow(
        aa_sequence_file=args.aa_sequence_file,
        feature_file_path=args.feature_file_path,
        feature2cluster_path=args.feature2cluster_path,
        protein_families_file=args.protein_families_file,
        output_dir=args.output_dir,
        feature_type=args.feature_type,
        annotation_file=args.annotation_file,
    )

if __name__ == "__main__":
    main()
