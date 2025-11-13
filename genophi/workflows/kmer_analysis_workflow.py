def kmer_analysis_workflow(
    aa_sequence_file,
    feature_file_path,
    feature2cluster_path,
    protein_families_file,
    output_dir,
    feature_type='strain',
    annotation_file=None,
    model_output_dir=None,
    quick_run=False,
    ignore_families=False  # NEW PARAMETER - MUST match setting from k-mer table generation
):
    """
    Complete workflow function to load amino acid sequences, extract and align predictive kmers,
    calculate coverage, identify segments, and plot segment mappings on sequences.

    Optionally aggregates SHAP importance values if `model_output_dir` is provided.

    IMPORTANT: The ignore_families parameter must match the setting used during k-mer table generation.
    If k-mer tables were generated with ignore_families=True, this must also be True here.

    Args:
        aa_sequence_file (str): Path to the amino acid sequences file in FASTA format.
        feature_file_path (str): Path to the feature selection output file (CSV).
        feature2cluster_path (str): Path to the feature-to-cluster mapping file (CSV).
        protein_families_file (str): Path to the protein families file (CSV).
        output_dir (str): Directory where output files (plots, coverage summaries) will be saved.
        feature_type (str): Type of features to extract ('host', 'strain', or 'phage').
        annotation_file (str, optional): Optional path to annotations file for enhanced plotting.
        model_output_dir (str, optional): Directory containing SHAP importance files from model runs.
        quick_run (bool): If True, skip alignment and coverage calculation.
        ignore_families (bool): CRITICAL - Must match the setting used during k-mer table generation.
                               If True, treats each k-mer independently without protein family context.
                               When True, alignments will likely fail since proteins aren't biologically related.

    Outputs:
        Segment plots and coverage summaries in the specified output directory.
        Aggregated SHAP values saved as `full_SHAP_values.csv` if `model_output_dir` is provided.
        
    Returns:
        None
    """
    logging.info("Starting k-mer analysis workflow...")
    logging.info(f"Parameters: feature_type={feature_type}, quick_run={quick_run}, ignore_families={ignore_families}")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    type_output_dir = os.path.join(output_dir, feature_type)
    if not os.path.exists(type_output_dir):
        os.makedirs(type_output_dir)

    # Step 1: Load amino acid sequences
    logging.info("Step 1: Loading amino acid sequences...")
    aa_sequences_df = load_aa_sequences(aa_sequence_file)
    aa_sequences_df.to_csv(os.path.join(type_output_dir, 'aa_sequences_df.csv'), index=False)
    logging.info(f"Amino acid sequences loaded and saved: {len(aa_sequences_df)} sequences")

    # Step 2: Get predictive kmers with ignore_families parameter
    logging.info("Step 2: Extracting predictive kmers...")
    filtered_kmers = get_predictive_kmers(
        feature_file_path, 
        feature2cluster_path, 
        feature_type,
        ignore_families=ignore_families
    )
    
    # Check if any predictive kmers were found
    if filtered_kmers.empty:
        logging.warning(f"No predictive {feature_type} kmers found. Saving empty kmers file and exiting workflow.")
        filtered_kmers.to_csv(os.path.join(type_output_dir, 'filtered_kmers.csv'), index=False)
        return
    
    filtered_kmers.to_csv(os.path.join(type_output_dir, 'filtered_kmers.csv'), index=False)
    logging.info(f"Predictive kmers extracted and saved: {len(filtered_kmers)} k-mers")

    # Step 3: Merge with protein families
    logging.info("Step 3: Merging kmers with protein families...")
    protein_families_df = merge_kmers_with_families(protein_families_file, aa_sequences_df, feature_type)
    
    # Check if any protein families were found
    if protein_families_df.empty:
        logging.warning(f"No protein families found for {feature_type} kmers. Exiting workflow.")
        protein_families_df.to_csv(os.path.join(type_output_dir, 'protein_families_df.csv'), index=False)
        return
    
    protein_families_df.to_csv(os.path.join(type_output_dir, 'protein_families_df.csv'), index=False)
    logging.info("Protein families merged and saved.")

    kmer_full_df = filtered_kmers.merge(protein_families_df, on='protein_family', how='inner')
    
    # Check if kmer_full_df is empty after merging
    if kmer_full_df.empty:
        logging.warning(f"No matching kmers found after merging with protein families. Exiting workflow.")
        with open(os.path.join(type_output_dir, f'{feature_type}_protein_sequences.faa'), 'w') as f:
            f.write('')
        return

    logging.info(f"Step 4: Saving protein sequences as FASTA file...")
    protein_seqs_df = kmer_full_df[['protein_ID', 'sequence']].drop_duplicates()
    seqrecords = []
    for _, row in protein_seqs_df.iterrows():
        seqrecords.append(SeqRecord(Seq(row['sequence']), id=row['protein_ID'], description=''))
    SeqIO.write(seqrecords, os.path.join(type_output_dir, f'{feature_type}_protein_sequences.faa'), 'fasta')
    logging.info(f"Saved {len(protein_seqs_df)} protein sequences to FASTA file")

    if quick_run:
        logging.info("Quick run mode enabled. Skipping sequence alignment and coverage calculation.")
        return

    if ignore_families:
        logging.warning("=" * 80)
        logging.warning("WARNING: ignore_families=True was used during k-mer generation.")
        logging.warning("Alignment may fail or produce meaningless results because proteins")
        logging.warning("sharing the same k-mer may not be biologically related.")
        logging.warning("Consider using --quick_run to skip alignment for ignore_families=True.")
        logging.warning("=" * 80)

    # Step 5: Construct kmer ID DataFrame
    logging.info("Step 5: Constructing kmer ID DataFrame...")
    kmer_id_df = construct_kmer_id_df(protein_families_df, kmer_full_df)
    logging.info("kmer ID DataFrame constructed.")

    # Step 6: Align sequences by protein family
    logging.info("Step 6: Aligning sequences by protein family...")
    aligned_df = pd.DataFrame()
    
    # Track families that need alignment
    families_by_size = protein_families_df.groupby("protein_family").size().reset_index(name='count')
    single_protein_families = len(families_by_size[families_by_size['count'] == 1])
    multi_protein_families = len(families_by_size[families_by_size['count'] > 1])
    
    logging.info(f"Found {len(families_by_size)} unique protein families:")
    logging.info(f"  - {single_protein_families} families with 1 protein (will skip)")
    logging.info(f"  - {multi_protein_families} families with 2+ proteins (will attempt alignment)")
    
    successful_alignments = 0
    failed_alignments = 0
    
    for family_name, family_group in protein_families_df.groupby("protein_family"):
        seqs_for_family = [(row["protein_ID"], row["sequence"]) for _, row in family_group.iterrows()]
        
        if not seqs_for_family or len(seqs_for_family) < 2:
            logging.debug(f"Skipping protein family {family_name}: Only {len(seqs_for_family)} sequence(s)")
            continue

        aligned_family_df = align_sequences(seqs_for_family, type_output_dir, family_name)
        if aligned_family_df.empty:
            failed_alignments += 1
            logging.debug(f"Alignment failed for protein family {family_name}")
            continue

        aligned_family_df['protein_family'] = family_name
        aligned_df = pd.concat([aligned_df, aligned_family_df], ignore_index=True)
        successful_alignments += 1

    logging.info(f"Alignment results: {successful_alignments} successful, {failed_alignments} failed")

    # Check if we have any successful alignments
    if aligned_df.empty:
        logging.error("=" * 80)
        logging.error("ERROR: No successful alignments performed!")
        logging.error("Possible reasons:")
        logging.error("  1. ignore_families=True was used (k-mers not grouped by protein families)")
        logging.error("  2. All protein families have only 1 sequence")
        logging.error("  3. All alignment attempts failed")
        logging.error("")
        logging.error("If using ignore_families=True, use --quick_run to skip alignment.")
        logging.error("=" * 80)
        return

    # Step 7: Find k-mer positions in aligned sequences
    logging.info("Step 7: Finding k-mer positions in aligned sequences...")
    aligned_df = aligned_df.merge(
        kmer_full_df[['Feature', 'cluster', 'protein_ID', 'kmer']], 
        on='protein_ID', 
        how='inner'
    )
    
    if aligned_df.empty:
        logging.error("No k-mers found in aligned sequences after merge. Exiting.")
        return
    
    aligned_df[['start_indices', 'stop_indices']] = aligned_df.apply(find_kmer_indices, axis=1)
    logging.info(f"Sequence alignment completed for {len(aligned_df)} proteins")

    # Step 8: Calculate coverage and identify segments
    logging.info("Step 8: Calculating coverage and identifying segments...")
    coverage_summary = calculate_coverage(aligned_df)
    coverage_summary = coverage_summary.drop_duplicates()
    segments_df = identify_segments(coverage_summary)
    segments_df.to_csv(os.path.join(type_output_dir, 'segments_df.csv'), index=False)
    logging.info(f"Identified {len(segments_df)} segments")

    # Step 9: Merge proteins without coverage segments
    logging.info("Step 9: Merging proteins without coverage segments...")
    final_segments_df = merge_no_coverage_proteins(segments_df, aligned_df)
    logging.info(f"Final segments DataFrame has {len(final_segments_df)} entries")

    # Step 10: Plot segments
    logging.info("Step 10: Plotting segments...")
    plotting_output_dir = os.path.join(type_output_dir, 'alignment_plots')
    plot_segments(final_segments_df, plotting_output_dir)
    logging.info("Segment plotting completed.")

    # Step 11: SHAP aggregation (optional)
    if model_output_dir:
        logging.info("Step 11: Aggregating SHAP values...")
        top_models_shap_df = aggregate_shap_values(model_output_dir)
        if not top_models_shap_df.empty:
            shap_output_path = os.path.join(output_dir, 'full_SHAP_values.csv')
            top_models_shap_df.to_csv(shap_output_path, index=False)
            logging.info(f"Aggregated SHAP values saved to {shap_output_path}")
        else:
            logging.warning("No SHAP importance data found in the specified directory.")
    
    logging.info("=" * 80)
    logging.info("K-mer analysis workflow completed successfully!")
    logging.info("=" * 80)import argparse
import pandas as pd
import logging
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings

warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

from genophi.kmer_modeling_analysis import (
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

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

def kmer_analysis_workflow(
    aa_sequence_file,
    feature_file_path,
    feature2cluster_path,
    protein_families_file,
    output_dir,
    feature_type='strain',
    annotation_file=None,
    model_output_dir=None,
    quick_run=False,
    ignore_families=False  # NEW PARAMETER - MUST match setting from k-mer table generation
):
    """
    Complete workflow function to load amino acid sequences, extract and align predictive kmers,
    calculate coverage, identify segments, and plot segment mappings on sequences.

    Optionally aggregates SHAP importance values if `model_output_dir` is provided.

    Args:
        aa_sequence_file (str): Path to the amino acid sequences file in FASTA format.
        feature_file_path (str): Path to the feature selection output file (CSV).
        feature2cluster_path (str): Path to the feature-to-cluster mapping file (CSV).
        protein_families_file (str): Path to the protein families file (CSV).
        output_dir (str): Directory where output files (plots, coverage summaries) will be saved.
        feature_type (str): Type of features to extract ('host', 'strain', or 'phage').
        annotation_file (str, optional): Optional path to annotations file for enhanced plotting.
        model_output_dir (str, optional): Directory containing SHAP importance files from model runs.
        quick_run (bool): If True, skip alignment and coverage calculation.
        ignore_families (bool): CRITICAL - Must match the setting used during k-mer table generation.
                               If True, treats each k-mer independently without protein family context.

    Outputs:
        Segment plots and coverage summaries in the specified output directory.
        Aggregated SHAP values saved as `full_SHAP_values.csv` if `model_output_dir` is provided.
    """
    logging.info("Starting k-mer analysis workflow...")
    logging.info(f"ignore_families setting: {ignore_families}")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    type_output_dir = os.path.join(output_dir, feature_type)
    if not os.path.exists(type_output_dir):
        os.makedirs(type_output_dir)

    # Step 1: Load amino acid sequences
    logging.info("Loading amino acid sequences...")
    aa_sequences_df = load_aa_sequences(aa_sequence_file)
    aa_sequences_df.to_csv(os.path.join(type_output_dir, 'aa_sequences_df.csv'), index=False)
    logging.info("Amino acid sequences loaded and saved.")

    # Step 2: Get predictive kmers with ignore_families parameter
    logging.info("Extracting predictive kmers...")
    filtered_kmers = get_predictive_kmers(
        feature_file_path, 
        feature2cluster_path, 
        feature_type,
        ignore_families=ignore_families
    )
    
    # Check if any predictive kmers were found
    if filtered_kmers.empty:
        logging.warning(f"No predictive {feature_type} kmers found. Saving empty kmers file and exiting workflow.")
        filtered_kmers.to_csv(os.path.join(type_output_dir, 'filtered_kmers.csv'), index=False)
        return
    
    filtered_kmers.to_csv(os.path.join(type_output_dir, 'filtered_kmers.csv'), index=False)
    logging.info("Predictive kmers extracted and saved.")

    # Step 3: Merge with protein families
    logging.info("Merging kmers with protein families...")
    protein_families_df = merge_kmers_with_families(protein_families_file, aa_sequences_df, feature_type)
    
    # Check if any protein families were found
    if protein_families_df.empty:
        logging.warning(f"No protein families found for {feature_type} kmers. Exiting workflow.")
        protein_families_df.to_csv(os.path.join(type_output_dir, 'protein_families_df.csv'), index=False)
        return
    
    protein_families_df.to_csv(os.path.join(type_output_dir, 'protein_families_df.csv'), index=False)
    logging.info("Protein families merged and saved.")

    kmer_full_df = filtered_kmers.merge(protein_families_df, on='protein_family', how='inner')
    
    # Check if kmer_full_df is empty after merging
    if kmer_full_df.empty:
        logging.warning(f"No matching kmers found after merging with protein families. Exiting workflow.")
        with open(os.path.join(type_output_dir, f'{feature_type}_protein_sequences.faa'), 'w') as f:
            f.write('')
        return

    logging.info("Saving proteins sequences as fasta file...")
    protein_seqs_df = kmer_full_df[['protein_ID', 'sequence']].drop_duplicates()
    seqrecords = []
    for _, row in protein_seqs_df.iterrows():
        seqrecords.append(SeqRecord(Seq(row['sequence']), id=row['protein_ID'], description=''))
    SeqIO.write(seqrecords, os.path.join(type_output_dir, f'{feature_type}_protein_sequences.faa'), 'fasta')

    if quick_run:
        logging.info("Quick run mode enabled. Skipping sequence alignment and coverage calculation.")
        return

    # Step 4: Construct kmer ID DataFrame
    logging.info("Constructing kmer ID DataFrame...")
    kmer_id_df = construct_kmer_id_df(protein_families_df, kmer_full_df)
    logging.info("kmer ID DataFrame constructed.")

    # Step 5: Align sequences by protein family
    logging.info("Aligning sequences by protein family...")
    aligned_df = pd.DataFrame()
    
    # Track families that need alignment
    families_needing_alignment = protein_families_df.groupby("protein_family").size()
    logging.info(f"Found {len(families_needing_alignment)} unique protein families")
    
    for family_name, family_group in protein_families_df.groupby("protein_family"):
        seqs_for_family = [(row["protein_ID"], row["sequence"]) for _, row in family_group.iterrows()]
        
        if not seqs_for_family or len(seqs_for_family) < 2:
            logging.warning(f"Skipping protein family {family_name}: Insufficient sequences for alignment ({len(seqs_for_family)} sequences).")
            continue

        aligned_family_df = align_sequences(seqs_for_family, type_output_dir, family_name)
        if aligned_family_df.empty:
            logging.warning(f"Skipping protein family {family_name}: Alignment failed or not performed.")
            continue

        aligned_family_df['protein_family'] = family_name
        aligned_df = pd.concat([aligned_df, aligned_family_df], ignore_index=True)

    # Check if we have any successful alignments
    if aligned_df.empty:
        logging.warning("No successful alignments performed. This may indicate:")
        logging.warning("  1. ignore_families=True was used (each k-mer is independent)")
        logging.warning("  2. All protein families have only 1 sequence")
        logging.warning("  3. Alignment failures occurred")
        logging.warning("Exiting workflow - cannot perform coverage analysis without alignments.")
        return

    aligned_df = aligned_df.merge(kmer_full_df[['Feature', 'cluster', 'protein_ID', 'kmer']], on='protein_ID', how='inner')
    aligned_df[['start_indices', 'stop_indices']] = aligned_df.apply(find_kmer_indices, axis=1)
    logging.info("Sequence alignment completed.")

    # Step 6: Calculate coverage and identify segments
    logging.info("Calculating coverage and identifying segments...")
    coverage_summary = calculate_coverage(aligned_df)
    coverage_summary = coverage_summary.drop_duplicates()
    segments_df = identify_segments(coverage_summary)
    segments_df.to_csv(os.path.join(type_output_dir, 'segments_df.csv'), index=False)
    logging.info("Coverage calculation and segment identification completed.")

    # Step 6.1: Merge proteins without coverage segments
    logging.info("Merging proteins without coverage segments...")
    final_segments_df = merge_no_coverage_proteins(segments_df, aligned_df)

    # Step 7: Plot segments with optional annotations
    logging.info("Plotting segments...")
    plotting_output_dir = os.path.join(type_output_dir, 'alignment_plots')
    plot_segments(final_segments_df, plotting_output_dir)
    logging.info("Segment plotting completed.")

    # Step 8: SHAP aggregation (optional)
    if model_output_dir:
        logging.info("Starting SHAP value aggregation...")
        top_models_shap_df = aggregate_shap_values(model_output_dir)
        if not top_models_shap_df.empty:
            shap_output_path = os.path.join(output_dir, 'full_SHAP_values.csv')
            top_models_shap_df.to_csv(shap_output_path, index=False)
            logging.info(f"Aggregated SHAP values saved to {shap_output_path}")
        else:
            logging.warning("No SHAP importance data found in the specified directory.")


def aggregate_shap_values(model_output_dir):
    """
    Aggregates SHAP importance values from multiple model runs.

    Args:
        model_output_dir (str): Directory containing model runs (e.g., 'run_0', 'run_1', etc.).

    Returns:
        DataFrame: Aggregated SHAP importance values.
    """
    runs = [x for x in os.listdir(model_output_dir) if 'run' in x]

    top_models_shap_df = pd.DataFrame()
    for run in runs:
        shap_values_csv_path = os.path.join(model_output_dir, run, 'shap_importances.csv')
        if os.path.exists(shap_values_csv_path):
            shap_values_temp = pd.read_csv(shap_values_csv_path)
            shap_values_temp = shap_values_temp.groupby(['feature', 'value']).agg({'shap_value': 'median'}).reset_index()
            top_models_shap_df = pd.concat([top_models_shap_df, shap_values_temp], ignore_index=True)
        else:
            logging.warning(f"SHAP importance file not found in {run}")

    return top_models_shap_df

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Workflow script for analyzing predictive kmers in amino acid sequences and optionally aggregating SHAP values.")
    parser.add_argument("--aa_sequence_file", required=True, help="Path to the AA sequence FASTA file.")
    parser.add_argument("--feature_file_path", required=True, help="Path to the feature selection output CSV file.")
    parser.add_argument("--feature2cluster_path", required=True, help="Path to the feature-to-cluster mapping CSV file.")
    parser.add_argument("--protein_families_file", required=True, help="Path to the protein families CSV file.")
    parser.add_argument("--output_dir", required=True, help="Directory for output files.")
    parser.add_argument("--feature_type", default="strain", help="Feature type to analyze (strain/phage).")
    parser.add_argument("--annotation_file", help="Path to optional annotation CSV file for plots.")
    parser.add_argument("--model_output_dir", help="Optional path to the directory containing SHAP importance files.")
    parser.add_argument("--quick_run", action="store_true", help="Run a quick version of the workflow without alignment and coverage calculation.")
    parser.add_argument("--ignore_families", action="store_true", help="CRITICAL: Set to True if k-mer table was generated with ignore_families=True.")

    args = parser.parse_args()

    kmer_analysis_workflow(
        aa_sequence_file=args.aa_sequence_file,
        feature_file_path=args.feature_file_path,
        feature2cluster_path=args.feature2cluster_path,
        protein_families_file=args.protein_families_file,
        output_dir=args.output_dir,
        feature_type=args.feature_type,
        annotation_file=args.annotation_file,
        model_output_dir=args.model_output_dir,
        quick_run=args.quick_run,
        ignore_families=args.ignore_families
    )

if __name__ == "__main__":
    main()