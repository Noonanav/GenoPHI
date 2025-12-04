import os
import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from genophi.kmer_modeling_analysis import (
    load_aa_sequences,
    get_predictive_kmers,
    merge_kmers_with_families,
    cluster_sequences_mmseqs,
    align_sequences,
    find_kmer_indices,
    calculate_coverage,
    identify_segments,
    plot_segments,
    aggregate_shap_values
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
    ignore_families=False
):
    """
    Workflow that supports both standard protein family alignment 
    AND 'ignore_families' mode which uses MMseqs2 clustering.
    """
    logging.info("Starting k-mer analysis workflow...")
    
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    type_output_dir = os.path.join(output_dir, feature_type)
    if not os.path.exists(type_output_dir): os.makedirs(type_output_dir)

    # 1. Load Sequences
    aa_sequences_df = load_aa_sequences(aa_sequence_file)
    # SAVE OUTPUT
    aa_sequences_df.to_csv(os.path.join(type_output_dir, 'aa_sequences_df.csv'), index=False)

    # 2. Extract K-mers
    filtered_kmers = get_predictive_kmers(feature_file_path, feature2cluster_path, feature_type, ignore_families)
    if filtered_kmers.empty: return
    # SAVE OUTPUT
    filtered_kmers.to_csv(os.path.join(type_output_dir, 'filtered_kmers.csv'), index=False)

    # 3. Merge / Search
    full_df = merge_kmers_with_families(
        protein_families_file, 
        aa_sequences_df, 
        feature_type, 
        ignore_families=ignore_families, 
        filtered_kmers=filtered_kmers
    )
    if full_df.empty: return

    # 3.5. Ensure Feature/Kmer Data is present in Standard Mode
    if not ignore_families:
        logging.info("Standard Mode: Merging predictive K-mers into protein dataframe...")
        full_df = full_df.merge(
            filtered_kmers[['protein_family', 'kmer', 'Feature']], 
            on='protein_family', 
            how='inner'
        )
        if full_df.empty:
            logging.warning("No intersection between protein families and predictive k-mers.")
            return
            
    # SAVE OUTPUT (This corresponds to protein_families_df.csv)
    full_df.to_csv(os.path.join(type_output_dir, 'protein_families_df.csv'), index=False)

    # Save protein sequences as FASTA (Useful for external tools)
    logging.info("Saving protein sequences to FASTA...")
    protein_seqs_df = full_df[['protein_ID', 'sequence']].drop_duplicates()
    seqrecords = [SeqRecord(Seq(row['sequence']), id=row['protein_ID'], description='') for _, row in protein_seqs_df.iterrows()]
    SeqIO.write(seqrecords, os.path.join(type_output_dir, f'{feature_type}_protein_sequences.faa'), 'fasta')

    # 4. Define Groups (Clustering vs Family)
    if ignore_families:
        # CLUSTERING MODE
        logging.info("ignore_families=True: Running MMseqs2 clustering...")
        clusters_df = cluster_sequences_mmseqs(full_df, type_output_dir)
        
        if clusters_df.empty:
            logging.error("Clustering failed.")
            return
        
        # SAVE OUTPUT
        clusters_df.to_csv(os.path.join(type_output_dir, 'mmseqs_clusters.csv'), index=False)
            
        # Merge cluster IDs back
        full_df = full_df.merge(clusters_df, on='protein_ID', how='inner')
        grouping_col = 'cluster_id'
        
    else:
        # STANDARD MODE
        logging.info("Standard Mode: Grouping by original Protein Family.")
        grouping_col = 'protein_family'
        full_df['cluster_id'] = full_df['protein_family']

    if quick_run: 
        logging.info("Quick run selected. Skipping alignment and coverage calculation.")
        return

    # 5. Alignment
    logging.info(f"Aligning sequences grouped by {grouping_col}...")
    aligned_dfs = []
    
    for group_id, group_data in full_df.groupby(grouping_col):
        seqs_for_group = group_data[['protein_ID', 'sequence']].drop_duplicates()
        seq_tuples = [(row['protein_ID'], row['sequence']) for _, row in seqs_for_group.iterrows()]
        
        if len(seq_tuples) < 2:
            continue
            
        aln_df = align_sequences(seq_tuples, type_output_dir, f"group_{group_id}")
        
        if not aln_df.empty:
            merged_aln = aln_df.merge(
                group_data[['protein_ID', 'cluster_id', 'kmer', 'Feature']].drop_duplicates(),
                on='protein_ID',
                how='inner'
            )
            aligned_dfs.append(merged_aln)

    if not aligned_dfs:
        logging.warning("No successful alignments.")
        return

    final_aligned_df = pd.concat(aligned_dfs, ignore_index=True)
    
    # 6. Indices & Coverage
    logging.info("Finding k-mer indices in aligned sequences...")
    final_aligned_df[['start_indices', 'stop_indices']] = final_aligned_df.apply(find_kmer_indices, axis=1)
    
    # SAVE OUTPUT
    final_aligned_df.to_csv(os.path.join(type_output_dir, 'aligned_df.csv'), index=False)
    
    coverage_df = calculate_coverage(final_aligned_df)
    
    # 7. Segments & Plotting
    segments_df = identify_segments(coverage_df)
    
    # SAVE OUTPUT
    segments_df.to_csv(os.path.join(type_output_dir, 'segments_df.csv'), index=False)
    
    plot_dir = os.path.join(type_output_dir, 'plots')
    plot_segments(segments_df, plot_dir)
    
    # 8. Optional SHAP
    if model_output_dir:
        shap_df = aggregate_shap_values(model_output_dir)
        if not shap_df.empty:
            shap_df.to_csv(os.path.join(output_dir, 'full_SHAP_values.csv'), index=False)

    logging.info("Workflow complete.")