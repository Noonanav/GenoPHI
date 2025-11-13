#!/usr/bin/env python3
"""
Standalone k-mer analysis workflow for ignore_families=True k-mer modeling.

This script:
1. Identifies predictive k-mers from feature selection results
2. Finds proteins containing those k-mers
3. Optionally clusters proteins with MMseqs2
4. Aligns proteins and calculates k-mer coverage
5. Identifies segments and outputs results

No protein families or predictive_proteins directory required.
"""

import os
import re
import logging
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from plotnine import *
import sys

# Optional MMseqs2 import
try:
    from phage_modeling.mmseqs2_clustering import create_mmseqs_database, run_mmseqs_cluster
    MMSEQS_AVAILABLE = True
except ImportError:
    MMSEQS_AVAILABLE = False
    logging.warning("MMseqs2 functions not available. Clustering will be disabled.")

logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)

def get_best_cutoff(modeling_dir):
    """Get the best performing cutoff from model performance metrics."""
    metrics_file = os.path.join(modeling_dir, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
    
    if not os.path.exists(metrics_file):
        raise FileNotFoundError(f"Model performance metrics not found: {metrics_file}")
    
    df = pd.read_csv(metrics_file)
    best_cutoff = df.iloc[0]['cut_off'].split('_')[-1]
    logging.info(f"Best performing cutoff: {best_cutoff}")
    return best_cutoff


def load_predictive_kmers(feature_table_path, feature2cluster_path, feature_type='strain'):
    """
    Load predictive k-mers from feature selection results.
    Diagnostic version with aggressive logging.
    """
    import time
    start_time = time.time()
    
    print(f"[DEBUG] Entering load_predictive_kmers at {time.time()}", flush=True)
    logging.info(f"Loading predictive {feature_type} k-mers...")
    logging.info(f"Feature table path: {feature_table_path}")
    logging.info(f"Feature2cluster path: {feature2cluster_path}")
    
    # Check file exists and size
    import os
    if os.path.exists(feature_table_path):
        size_mb = os.path.getsize(feature_table_path) / 1024**2
        logging.info(f"Feature table exists: {size_mb:.2f} MB")
        print(f"[DEBUG] Feature table size: {size_mb:.2f} MB", flush=True)
    else:
        logging.error(f"Feature table NOT FOUND: {feature_table_path}")
        raise FileNotFoundError(f"Feature table not found: {feature_table_path}")
    
    # Load feature table
    logging.info("Loading feature table...")
    print("[DEBUG] About to call pd.read_csv on feature table...", flush=True)
    
    t0 = time.time()
    feature_df = pd.read_csv(feature_table_path)
    t1 = time.time()
    
    print(f"[DEBUG] pd.read_csv completed in {t1-t0:.2f} seconds", flush=True)
    logging.info(f"Feature table loaded in {t1-t0:.2f} seconds")
    logging.info(f"Feature table shape: {feature_df.shape}")
    
    feature_prefix = feature_type[0] + 'c_'
    predictive_features = [col for col in feature_df.columns if col.startswith(feature_prefix)]
    logging.info(f"Found {len(predictive_features)} predictive {feature_type} features")
    print(f"[DEBUG] Found {len(predictive_features)} predictive features", flush=True)
    
    # Load feature to cluster mapping
    logging.info("Loading feature2cluster mapping...")
    print("[DEBUG] Loading feature2cluster...", flush=True)
    
    feature2cluster = pd.read_csv(feature2cluster_path)
    print(f"[DEBUG] Feature2cluster loaded: {feature2cluster.shape}", flush=True)
    
    feature2cluster.rename(columns={'Cluster_Label': 'cluster'}, inplace=True)
    predictive_kmers = feature2cluster[feature2cluster['Feature'].isin(predictive_features)].copy()
    predictive_kmers['kmer'] = predictive_kmers['cluster']
    
    logging.info(f"Loaded {len(predictive_kmers)} predictive k-mer mappings")
    print(f"[DEBUG] Predictive k-mers: {len(predictive_kmers)}", flush=True)
    
    # Build feature presence WITHOUT stack()
    logging.info("Building feature presence (optimized method without stack)...")
    print("[DEBUG] Starting feature presence building...", flush=True)
    
    genome_col = feature_type
    results = []
    
    t0 = time.time()
    for i, feature in enumerate(predictive_features):
        if i % 100 == 0:
            elapsed = time.time() - t0
            rate = i / elapsed if elapsed > 0 else 0
            eta = (len(predictive_features) - i) / rate if rate > 0 else 0
            print(f"[DEBUG] Processing feature {i}/{len(predictive_features)} - ETA: {eta:.1f}s", flush=True)
            logging.info(f"  Feature {i}/{len(predictive_features)}")
        
        # Get genomes where feature is present
        present_mask = feature_df[feature] == 1
        present_genomes = feature_df.loc[present_mask, genome_col].unique()
        
        for genome in present_genomes:
            results.append({
                genome_col: genome,
                'Feature': feature
            })
    
    print(f"[DEBUG] Feature presence complete: {len(results)} pairs", flush=True)
    logging.info(f"Created {len(results):,} feature-genome pairs")
    
    feature_presence = pd.DataFrame(results)
    
    # Merge with k-mer info
    logging.info("Merging with k-mer information...")
    print("[DEBUG] Merging...", flush=True)
    
    predictive_kmers_with_genomes = feature_presence.merge(
        predictive_kmers[['Feature', 'cluster', 'kmer']], 
        on='Feature', 
        how='inner'
    )
    
    elapsed_total = time.time() - start_time
    logging.info(f"Completed in {elapsed_total:.2f} seconds ({elapsed_total/60:.2f} minutes)")
    print(f"[DEBUG] load_predictive_kmers complete in {elapsed_total:.2f}s", flush=True)
    
    return predictive_kmers_with_genomes

def load_aa_sequences(aa_sequence_file):
    """Load amino acid sequences from FASTA file."""
    import time
    
    print(f"[DEBUG] Starting load_aa_sequences at {time.time()}", flush=True)
    logging.info(f"Loading sequences from {aa_sequence_file}")
    
    # Check file size first
    import os
    if os.path.exists(aa_sequence_file):
        size_mb = os.path.getsize(aa_sequence_file) / 1024**2
        print(f"[DEBUG] FASTA file size: {size_mb:.2f} MB", flush=True)
        logging.info(f"FASTA file size: {size_mb:.2f} MB")
    
    print("[DEBUG] About to parse FASTA file with SeqIO.parse()...", flush=True)
    t0 = time.time()
    
    aa_records = list(SeqIO.parse(aa_sequence_file, 'fasta'))
    
    t1 = time.time()
    print(f"[DEBUG] FASTA parsing complete in {t1-t0:.2f} seconds", flush=True)
    print(f"[DEBUG] Parsed {len(aa_records)} records", flush=True)
    
    print("[DEBUG] Converting to DataFrame...", flush=True)
    t2 = time.time()
    
    aa_sequences_df = pd.DataFrame({
        'protein_ID': [record.id for record in aa_records],
        'sequence': [str(record.seq) for record in aa_records]
    })
    
    t3 = time.time()
    print(f"[DEBUG] DataFrame conversion complete in {t3-t2:.2f} seconds", flush=True)
    
    logging.info(f"Loaded {len(aa_sequences_df)} protein sequences")
    print(f"[DEBUG] load_aa_sequences complete", flush=True)
    
    return aa_sequences_df


def map_kmers_to_proteins(predictive_kmers_df, protein_mapping_csv, aa_sequences_df, feature_type='strain'):
    """
    Map k-mers to proteins using k-mer index for ultra-fast lookup.
    
    Returns DataFrame with: feature_type, Feature, cluster, kmer, protein_ID, sequence
    """
    import time
    from collections import defaultdict
    
    start_time = time.time()
    
    print("[DEBUG] Entering map_kmers_to_proteins", flush=True)  # ← ADD THIS
    
    logging.info("="*60)
    logging.info("Starting map_kmers_to_proteins (INDEXED VERSION)")
    logging.info("="*60)
    logging.info(f"Input predictive_kmers_df: {predictive_kmers_df.shape}")
    logging.info(f"Input aa_sequences_df: {aa_sequences_df.shape}")
    
    print("[DEBUG] Creating unique_kmers set...", flush=True)  # ← ADD THIS
    
    # Step 1: Get unique predictive k-mers and their lengths
    unique_kmers = set(predictive_kmers_df['kmer'].unique())
    kmer_lengths = set(len(k) for k in unique_kmers)
    
    print(f"[DEBUG] Found {len(unique_kmers):,} unique k-mers", flush=True)  # ← ADD THIS
    print(f"[DEBUG] K-mer lengths: {sorted(kmer_lengths)}", flush=True)  # ← ADD THIS
    
    logging.info(f"Found {len(unique_kmers):,} unique predictive k-mers")
    logging.info(f"K-mer lengths: {sorted(kmer_lengths)}")
    
    # Step 2: Build k-mer index (hash table: kmer -> set of protein_IDs)
    logging.info("Building k-mer index (this is a one-time cost)...")
    print("[DEBUG] Starting k-mer indexing...", flush=True)  # ← ADD THIS
    
    kmer_to_proteins = defaultdict(set)
    
    for idx, (protein_id, seq) in enumerate(zip(aa_sequences_df['protein_ID'], 
                                                  aa_sequences_df['sequence']), 1):
        if idx % 50000 == 0 or idx == 1:
            elapsed = time.time() - start_time
            rate = idx / elapsed if elapsed > 0 else 0
            eta = (len(aa_sequences_df) - idx) / rate if rate > 0 else 0
            
            # ← ADD THIS PRINT:
            print(f"[DEBUG] Indexing: {idx:,}/{len(aa_sequences_df):,} "
                  f"({100*idx/len(aa_sequences_df):.1f}%) | "
                  f"Rate: {rate:.1f} prot/s | ETA: {eta/60:.1f} min", flush=True)
            
            logging.info(f"  Indexing: {idx:,}/{len(aa_sequences_df):,} proteins "
                        f"({100*idx/len(aa_sequences_df):.1f}%) | "
                        f"ETA: {eta:.1f}s")
        
        # Extract all k-mers of relevant lengths from this protein
        for k in kmer_lengths:
            if len(seq) < k:
                continue
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                # Only index if this k-mer is predictive
                if kmer in unique_kmers:
                    kmer_to_proteins[kmer].add(protein_id)
    
    index_time = time.time() - start_time
    print(f"[DEBUG] Indexing COMPLETE in {index_time/60:.2f} minutes", flush=True)  # ADD THIS
    logging.info(f"Indexing complete in {index_time:.2f}s ({index_time/60:.2f} min)")
    logging.info(f"Index contains {len(kmer_to_proteins):,} k-mers that appear in proteins")

    # Step 3: Load protein mapping AND BUILD LOOKUP TABLE
    print("[DEBUG] Loading protein mapping CSV...", flush=True)
    logging.info(f"Loading protein mapping from: {protein_mapping_csv}")
    protein_mapping = pd.read_csv(protein_mapping_csv)
    print(f"[DEBUG] Protein mapping loaded: {protein_mapping.shape}", flush=True)
    logging.info(f"Protein mapping shape: {protein_mapping.shape}")

    # BUILD GENOME -> PROTEINS LOOKUP (this is the key optimization!)
    print("[DEBUG] Building genome->proteins lookup...", flush=True)
    genome_to_proteins = {}
    for genome, group in protein_mapping.groupby(feature_type):
        genome_to_proteins[genome] = set(group['protein_ID'])
    print(f"[DEBUG] Genome lookup built for {len(genome_to_proteins)} genomes", flush=True)

    # Step 4: Fast lookup - use BOTH indexes
    print("[DEBUG] Starting lookup phase...", flush=True)
    logging.info("Mapping k-mers to proteins using index (ultra-fast lookup)...")

    results = []
    total_kmers = len(predictive_kmers_df)
    print(f"[DEBUG] Total k-mer-genome pairs to process: {total_kmers:,}", flush=True)

    for idx, (_, row) in enumerate(predictive_kmers_df.iterrows(), 1):
        if idx % 100000 == 0 or idx == 1:  # Changed to 100K for less frequent updates
            elapsed = time.time() - start_time
            rate = idx / elapsed if elapsed > 0 else 0
            eta = (total_kmers - idx) / rate if rate > 0 else 0
            
            print(f"[DEBUG] Lookup: {idx:,}/{total_kmers:,} ({100*idx/total_kmers:.1f}%) | "
                f"Found: {len(results):,} mappings | ETA: {eta/60:.1f} min", flush=True)
            
            logging.info(f"  Progress: {idx:,}/{total_kmers:,} k-mer-genome pairs "
                        f"({100*idx/total_kmers:.1f}%) | "
                        f"ETA: {eta:.1f}s | Found: {len(results):,} mappings")
        
        kmer = row['kmer']
        genome = row[feature_type]
        feature = row['Feature']
        cluster = row['cluster']
        
        # Fast O(1) lookup: which proteins have this k-mer?
        candidate_proteins = kmer_to_proteins.get(kmer, set())
        
        if not candidate_proteins:
            continue
        
        # Fast O(1) lookup: which proteins are in this genome?
        genome_proteins_set = genome_to_proteins.get(genome, set())
        
        # Set intersection to find proteins that have the k-mer AND are in this genome
        matching_proteins = candidate_proteins & genome_proteins_set
        
        # Add all matching proteins
        for protein_id in matching_proteins:
            results.append({
                feature_type: genome,
                'Feature': feature,
                'cluster': cluster,
                'kmer': kmer,
                'protein_ID': protein_id
            })
    
    # Convert to DataFrame
    kmer_protein = pd.DataFrame(results)
    
    if kmer_protein.empty:
        logging.error("No k-mer-protein mappings found!")
        return pd.DataFrame()
    
    # Add sequences
    logging.info("Adding sequences to mappings...")
    kmer_protein = kmer_protein.merge(aa_sequences_df, on='protein_ID', how='left')
    
    elapsed_total = time.time() - start_time
    logging.info("="*60)
    logging.info("MAPPING COMPLETE")
    logging.info("="*60)
    logging.info(f"Total time: {elapsed_total:.2f}s ({elapsed_total/60:.2f} min)")
    logging.info(f"  - Indexing time: {index_time:.2f}s ({100*index_time/elapsed_total:.1f}%)")
    logging.info(f"  - Lookup time: {elapsed_total-index_time:.2f}s ({100*(elapsed_total-index_time)/elapsed_total:.1f}%)")
    logging.info(f"Final dataset:")
    logging.info(f"  - Rows: {len(kmer_protein):,}")
    logging.info(f"  - Proteins: {kmer_protein['protein_ID'].nunique():,}")
    logging.info(f"  - Genomes: {kmer_protein[feature_type].nunique():,}")
    logging.info(f"  - Features: {kmer_protein['Feature'].nunique():,}")
    logging.info(f"  - Unique k-mers: {kmer_protein['kmer'].nunique():,}")
    logging.info(f"  - Memory: {kmer_protein.memory_usage(deep=True).sum() / 1024**2:.2f} MB")
    logging.info("="*60)
    
    return kmer_protein

def find_kmer_positions_in_sequence(kmer, sequence):
    """
    Find all positions where a k-mer occurs in a sequence.
    
    Returns list of (start, end) tuples (end is exclusive, like Python slicing).
    """
    positions = []
    start = 0
    
    while True:
        pos = sequence.find(kmer, start)
        if pos == -1:
            break
        positions.append((pos, pos + len(kmer)))
        start = pos + 1  # Allow overlapping matches
    
    return positions


def extract_segments_from_raw_sequences(kmer_protein_df, feature_type='strain'):
    """
    Extract k-mer segment positions directly from raw (unaligned) sequences.
    No alignment required!
    """
    import time
    start_time = time.time()
    
    print("[DEBUG] Entering extract_segments_from_raw_sequences", flush=True)
    
    logging.info("=" * 60)
    logging.info("Extracting k-mer segments from raw sequences (no alignment)")
    logging.info("=" * 60)
    
    segments = []
    
    # DEDUPLICATE first to avoid processing same k-mer-protein pair multiple times
    print("[DEBUG] Deduplicating k-mer-protein pairs...", flush=True)
    kmer_protein_unique = kmer_protein_df[[
        feature_type, 'protein_ID', 'Feature', 'cluster', 'kmer', 'sequence'
    ]].drop_duplicates()
    
    total_rows = len(kmer_protein_unique)
    print(f"[DEBUG] Processing {total_rows:,} unique k-mer-protein pairs", flush=True)
    logging.info(f"Processing {total_rows:,} unique k-mer-protein pairs...")
    
    for idx, (_, data) in enumerate(kmer_protein_unique.iterrows(), 1):
        if idx % 10000 == 0 or idx == 1:
            elapsed = time.time() - start_time
            rate = idx / elapsed if elapsed > 0 else 0
            eta = (total_rows - idx) / rate if rate > 0 else 0
            
            print(f"[DEBUG] Segment extraction: {idx:,}/{total_rows:,} ({100*idx/total_rows:.1f}%) | "
                  f"Rate: {rate:.1f}/sec | ETA: {eta/60:.1f} min | Segments: {len(segments):,}", flush=True)
            
            logging.info(f"  Progress: {idx:,}/{total_rows:,} ({100*idx/total_rows:.1f}%) | "
                        f"Rate: {rate:.1f} pairs/sec | ETA: {eta:.1f} sec")
        
        kmer = data['kmer']
        sequence = data['sequence']
        protein_id = data['protein_ID']
        feature = data['Feature']
        genome = data[feature_type]
        cluster = data['cluster']
        
        # Find all positions of this k-mer in this protein
        positions = find_kmer_positions_in_sequence(kmer, sequence)
        
        if not positions:
            logging.warning(f"K-mer '{kmer}' not found in protein {protein_id} (unexpected)")
            continue
        
        # Create a segment for each occurrence
        for start, end in positions:
            segments.append({
                feature_type: genome,
                'protein_ID': protein_id,
                'Feature': feature,
                'cluster': cluster,
                'kmer': kmer,
                'start': start,
                'end': end,
                'segment_sequence': sequence[start:end],
                'segment_length': end - start,
                'full_protein_sequence': sequence  # Keep for merged segments
            })
    
    elapsed_total = time.time() - start_time
    print(f"[DEBUG] Segment extraction COMPLETE in {elapsed_total/60:.2f} minutes", flush=True)
    logging.info(f"Completed in {elapsed_total:.2f} seconds")
    
    if not segments:
        logging.warning("No segments found!")
        return pd.DataFrame()
    
    print(f"[DEBUG] Creating DataFrame from {len(segments):,} segments...", flush=True)
    segments_df = pd.DataFrame(segments)
    print(f"[DEBUG] DataFrame created", flush=True)
    
    logging.info("=" * 60)
    logging.info("SEGMENT EXTRACTION SUMMARY")
    logging.info("=" * 60)
    logging.info(f"Total segments found: {len(segments_df):,}")
    logging.info(f"Unique proteins: {segments_df['protein_ID'].nunique():,}")
    logging.info(f"Unique k-mers: {segments_df['kmer'].nunique():,}")
    logging.info(f"Unique features: {segments_df['Feature'].nunique():,}")
    logging.info(f"Average segments per protein: {len(segments_df)/segments_df['protein_ID'].nunique():.1f}")
    
    # Show k-mer length distribution
    kmer_lengths = segments_df['kmer'].str.len().value_counts().sort_index()
    logging.info(f"K-mer length distribution:")
    for length, count in kmer_lengths.items():
        logging.info(f"  {length}-mers: {count:,} segments")
    
    logging.info("=" * 60)
    
    return segments_df


def merge_overlapping_segments(segments_df, feature_type='strain'):
    """
    Merge overlapping k-mer segments within each protein.
    This creates "covered regions" where one or more k-mers are present.
    """
    logging.info("Merging overlapping segments...")
    
    if segments_df.empty:
        logging.warning("No segments to merge")
        return pd.DataFrame()
    
    merged_segments = []
    
    # Group by protein
    for (genome, protein_id), protein_group in segments_df.groupby([feature_type, 'protein_ID']):
        # Get the full protein sequence (should be same for all rows)
        full_sequence = protein_group['full_protein_sequence'].iloc[0]
        
        # Sort by start position
        protein_group = protein_group.sort_values('start')
        
        # Track current merged segment
        current_start = None
        current_end = None
        current_features = set()
        current_kmers = set()
        
        for _, row in protein_group.iterrows():
            if current_start is None:
                # Start new segment
                current_start = row['start']
                current_end = row['end']
                current_features.add(row['Feature'])
                current_kmers.add(row['kmer'])
            elif row['start'] <= current_end:
                # Overlaps with current segment - extend it
                current_end = max(current_end, row['end'])
                current_features.add(row['Feature'])
                current_kmers.add(row['kmer'])
            else:
                # No overlap - save current and start new
                merged_segments.append({
                    feature_type: genome,
                    'protein_ID': protein_id,
                    'start': current_start,
                    'end': current_end,
                    'length': current_end - current_start,
                    'segment_sequence': full_sequence[current_start:current_end],
                    'n_features': len(current_features),
                    'n_kmers': len(current_kmers),
                    'features': ','.join(sorted(current_features)),
                    'kmers': ','.join(sorted(current_kmers))
                })
                
                current_start = row['start']
                current_end = row['end']
                current_features = {row['Feature']}
                current_kmers = {row['kmer']}
        
        # Don't forget the last segment
        if current_start is not None:
            merged_segments.append({
                feature_type: genome,
                'protein_ID': protein_id,
                'start': current_start,
                'end': current_end,
                'length': current_end - current_start,
                'segment_sequence': full_sequence[current_start:current_end],
                'n_features': len(current_features),
                'n_kmers': len(current_kmers),
                'features': ','.join(sorted(current_features)),
                'kmers': ','.join(sorted(current_kmers))
            })
    
    if not merged_segments:
        logging.warning("No merged segments created")
        return pd.DataFrame()
    
    merged_df = pd.DataFrame(merged_segments)
    
    logging.info(f"Created {len(merged_df):,} merged segments from {len(segments_df):,} original segments")
    logging.info(f"Reduction: {100*(1 - len(merged_df)/len(segments_df)):.1f}%")
    
    return merged_df

def cluster_proteins_mmseqs(protein_df, output_dir, feature_type, threads=4, 
                            min_seq_id=0.5, coverage=0.8, sensitivity=7.5):
    """Cluster proteins using MMseqs2."""
    if not MMSEQS_AVAILABLE:
        logging.warning("MMseqs2 not available. Assigning all proteins to cluster 0.")
        protein_df['cluster_id'] = 0
        return protein_df
    
    logging.info("Clustering proteins with MMseqs2...")
    
    cluster_dir = os.path.join(output_dir, "clustering")
    os.makedirs(cluster_dir, exist_ok=True)
    
    # Get unique proteins
    unique_proteins = protein_df[['protein_ID', 'sequence']].drop_duplicates()
    
    # Write temp FASTA
    temp_fasta = os.path.join(cluster_dir, f"{feature_type}_temp.fasta")
    records = [
        SeqRecord(Seq(row['sequence']), id=row['protein_ID'], description='')
        for _, row in unique_proteins.iterrows()
    ]
    SeqIO.write(records, temp_fasta, 'fasta')
    
    # Create MMseqs2 database
    db_name = os.path.join(cluster_dir, f"{feature_type}_db")
    create_mmseqs_database(temp_fasta, db_name, suffix="fasta", 
                          input_type="file", strains=None, threads=threads)
    
    # Run clustering
    tmp_dir = os.path.join(cluster_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    
    clusters_tsv = run_mmseqs_cluster(db_name, cluster_dir, tmp_dir, 
                                     coverage, min_seq_id, sensitivity, threads)
    
    # Parse clusters
    clusters_df = pd.read_csv(clusters_tsv, sep="\t", header=None, 
                             names=["cluster_rep", "protein_ID"])
    clusters_df['cluster_id'] = clusters_df.groupby("cluster_rep").ngroup()
    
    # Merge back
    protein_df = protein_df.merge(
        clusters_df[['protein_ID', 'cluster_id']], 
        on='protein_ID', 
        how='left'
    )
    
    # Fill any missing cluster IDs
    protein_df['cluster_id'] = protein_df['cluster_id'].fillna(0).astype(int)
    
    logging.info(f"Clustered into {protein_df['cluster_id'].nunique()} clusters")
    
    # Cleanup
    os.remove(temp_fasta)
    
    return protein_df


def align_sequences(sequences, output_dir, cluster_name):
    """Align sequences within a cluster using MAFFT."""
    logging.info(f"Aligning {len(sequences)} sequences for {cluster_name}")
    
    if len(sequences) < 2:
        logging.warning(f"Skipping {cluster_name}: Need at least 2 sequences")
        return pd.DataFrame()
    
    # Setup temp files - sanitize all special characters
    safe_name = (cluster_name.replace('/', '_').replace('|', '_')
                           .replace(':', '_').replace('.', '_')
                           .replace(' ', '_'))
    temp_fasta = os.path.join(output_dir, f"{safe_name}_temp.fasta")
    temp_aln = os.path.join(output_dir, f"{safe_name}_temp.aln")
    
    # Create sequence records with temp IDs
    seq_records = []
    temp_id_map = {}
    for i, (protein_id, seq) in enumerate(sequences):
        if not seq or len(seq.strip('-')) == 0:
            continue
        temp_id = f"seq_{i}"
        seq_records.append(SeqRecord(Seq(seq), id=temp_id, description=""))
        temp_id_map[temp_id] = protein_id
    
    if len(seq_records) < 2:
        logging.warning(f"Skipping {cluster_name}: Insufficient valid sequences")
        return pd.DataFrame()
    
    # Write FASTA
    SeqIO.write(seq_records, temp_fasta, "fasta")
    
    # Run MAFFT
    try:
        mafft_cline = MafftCommandline(input=temp_fasta)
        stdout, stderr = mafft_cline()
        
        with open(temp_aln, "w") as handle:
            handle.write(stdout)
    except Exception as e:
        logging.error(f"MAFFT failed for {cluster_name}: {e}")
        return pd.DataFrame()
    
    # Parse alignment
    try:
        alignment = AlignIO.read(temp_aln, "fasta")
    except Exception as e:
        logging.error(f"Failed to parse alignment for {cluster_name}: {e}")
        return pd.DataFrame()
    
    # Remove leading gaps
    aln_length = alignment.get_alignment_length()
    non_gap_positions = [any(record.seq[i] != '-' for record in alignment) 
                        for i in range(aln_length)]
    first_non_gap = non_gap_positions.index(True) if True in non_gap_positions else 0
    
    # Extract aligned sequences
    aligned_sequences = []
    for record in alignment:
        trimmed_seq = str(record.seq[first_non_gap:])
        original_id = temp_id_map[record.id]
        start_pos = len(trimmed_seq) - len(trimmed_seq.lstrip('-'))
        aligned_sequences.append({
            'protein_ID': original_id,
            'aln_sequence': trimmed_seq,
            'start_index': start_pos
        })
    
    # Cleanup
    for f in [temp_fasta, temp_aln, temp_aln.replace('.aln', '.dnd')]:
        if os.path.exists(f):
            os.remove(f)
    
    return pd.DataFrame(aligned_sequences)


def find_kmer_indices(row):
    """Find k-mer positions in aligned sequence."""
    kmer = row['kmer']
    seq = row['aln_sequence']
    
    # Create pattern with optional gaps
    pattern = "-*".join(kmer)
    
    start_indices = []
    stop_indices = []
    
    for match in re.finditer(pattern, seq):
        start = match.start()
        stop = match.end() - 1
        
        # Verify correct number of AAs
        segment = seq[start:stop+1]
        aa_count = sum(1 for c in segment if c != '-')
        
        if aa_count == len(kmer):
            start_indices.append(start)
            stop_indices.append(stop)
    
    return pd.Series([start_indices, stop_indices], 
                    index=['start_indices', 'stop_indices'])


def calculate_coverage(aligned_df):
    """Calculate k-mer coverage for aligned sequences."""
    logging.info("Calculating k-mer coverage...")
    
    coverage_data = []
    
    for _, row in aligned_df.iterrows():
        if len(row['start_indices']) == 0:
            continue
        
        coverage = np.zeros(len(row['aln_sequence']), dtype=int)
        
        # Mark coverage
        for start, stop in zip(row['start_indices'], row['stop_indices']):
            coverage[start:stop+1] = 1
        
        # Create records for each position
        for idx, residue in enumerate(row['aln_sequence']):
            coverage_data.append({
                'Feature': row['Feature'],
                'protein_ID': row['protein_ID'],
                'AA_index': idx,
                'Residue': residue,
                'coverage': coverage[idx],
                'is_gap': 1 if residue == '-' else 0
            })
    
    if not coverage_data:
        logging.warning("No coverage data generated")
        return pd.DataFrame()
    
    coverage_df = pd.DataFrame(coverage_data)
    coverage_df = coverage_df.groupby(
        ['Feature', 'protein_ID', 'AA_index', 'Residue']
    ).agg({'coverage': 'max', 'is_gap': 'max'}).reset_index()
    
    logging.info(f"Calculated coverage for {len(coverage_df)} positions")
    return coverage_df


def identify_segments(coverage_df):
    """Identify contiguous segments from coverage data."""
    logging.info("Identifying coverage segments...")
    
    # Sort by position
    df = coverage_df.sort_values(by=['Feature', 'protein_ID', 'AA_index'])
    
    # Detect segment changes
    df['prev_coverage'] = df.groupby(['Feature', 'protein_ID'])['coverage'].shift(1).fillna(-1)
    df['segment_change'] = (df['coverage'] != df['prev_coverage'])
    df['segment_id'] = df.groupby(['Feature', 'protein_ID'])['segment_change'].cumsum()
    
    # Aggregate segments
    segments_df = df.groupby(
        ['Feature', 'protein_ID', 'coverage', 'segment_id']
    ).agg(
        start=('AA_index', 'min'),
        stop=('AA_index', 'max'),
        is_gap=('is_gap', 'max')
    ).reset_index()
    
    # Adjust stop to be inclusive
    segments_df['stop'] += 1
    
    logging.info(f"Identified {len(segments_df)} segments")
    return segments_df


def add_genome_info(segments_df, kmer_protein_df, feature_type='strain'):
    """Add genome (strain/phage) information to segments."""
    # Get protein -> genome mapping
    protein_genome = kmer_protein_df[[feature_type, 'protein_ID']].drop_duplicates()
    
    # Merge
    segments_df = segments_df.merge(protein_genome, on='protein_ID', how='left')
    
    return segments_df


def extract_segment_sequences(segments_df, aligned_df):
    """Extract actual sequences for segments."""
    logging.info("Extracting segment sequences...")
    
    # Only process covered segments
    covered_segments = segments_df[segments_df['coverage'] > 0].copy()
    
    # Merge with aligned sequences
    covered_segments = covered_segments.merge(
        aligned_df[['protein_ID', 'aln_sequence']].drop_duplicates(),
        on='protein_ID',
        how='left'
    )
    
    # Extract sequences
    covered_segments['segment_sequence'] = covered_segments.apply(
        lambda x: x['aln_sequence'][int(x['start']):int(x['stop'])] 
        if pd.notna(x['aln_sequence']) else '',
        axis=1
    )
    
    # Remove gaps
    covered_segments['segment_sequence_nogaps'] = covered_segments['segment_sequence'].str.replace('-', '')
    
    # Calculate metrics
    covered_segments['segment_length'] = covered_segments['stop'] - covered_segments['start']
    covered_segments['aa_count'] = covered_segments['segment_sequence_nogaps'].str.len()
    
    logging.info(f"Extracted sequences for {len(covered_segments)} segments")
    
    return covered_segments


def plot_segments(segments_df, aligned_df, output_dir, feature_type):
    """Create segment visualization plots."""
    logging.info(f"Creating segment plots in {output_dir}...")
    os.makedirs(output_dir, exist_ok=True)
    
    # Merge cluster info if available
    if 'cluster_id' in aligned_df.columns:
        segments_df = segments_df.merge(
            aligned_df[['protein_ID', 'cluster_id']].drop_duplicates(),
            on='protein_ID',
            how='left'
        )
    else:
        segments_df['cluster_id'] = 0
    
    # Separate covered and uncovered
    uncovered = segments_df[segments_df['coverage'] == 0].copy()
    covered = segments_df[segments_df['coverage'] == 1].copy()
    
    # Order proteins by feature presence
    if not covered.empty:
        feature_matrix = covered.pivot_table(
            index='protein_ID',
            columns='Feature',
            values='coverage',
            aggfunc='sum',
            fill_value=0
        )
        
        if len(feature_matrix) > 1:
            try:
                distances = pdist(feature_matrix, metric='euclidean')
                linkage_matrix = linkage(distances, method='ward')
                ordered_indices = leaves_list(linkage_matrix)
                ordered_proteins = feature_matrix.index[ordered_indices].tolist()
            except:
                ordered_proteins = feature_matrix.index.tolist()
        else:
            ordered_proteins = feature_matrix.index.tolist()
    else:
        ordered_proteins = segments_df['protein_ID'].unique().tolist()
    
    segments_df['protein_ID'] = pd.Categorical(
        segments_df['protein_ID'],
        categories=ordered_proteins,
        ordered=True
    )
    
    uncovered = segments_df[segments_df['coverage'] == 0]
    covered = segments_df[segments_df['coverage'] == 1]
    
    # Calculate plot height
    n_proteins = segments_df['protein_ID'].nunique()
    fig_height = min(max(4, n_proteins * 0.3), 20)
    
    # Create plot
    plot = (
        ggplot() +
        (geom_segment(
            data=uncovered,
            mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID'),
            color='grey',
            size=3
        ) if not uncovered.empty else geom_blank()) +
        (geom_segment(
            data=covered,
            mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID', 
                       color='Feature'),
            size=3
        ) if not covered.empty else geom_blank()) +
        labs(
            title=f'{feature_type.capitalize()} K-mer Coverage',
            x='Aligned Position',
            y='Protein ID'
        ) +
        theme(
            axis_text_x=element_text(rotation=90),
            panel_background=element_rect(fill='white'),
            figure_size=(12, fig_height),
            axis_text_y=element_text(size=6)
        )
    )
    
    # Add faceting if clusters exist
    if 'cluster_id' in segments_df.columns and segments_df['cluster_id'].nunique() > 1:
        plot = plot + facet_wrap('~cluster_id', dir='v', scales='free_y', ncol=1)
    
    try:
        plot_file = os.path.join(output_dir, f'{feature_type}_kmer_coverage.png')
        plot.save(plot_file)
        logging.info(f"Saved plot: {plot_file}")
    except Exception as e:
        logging.error(f"Failed to save plot: {e}")

def save_proteins_with_kmers_fasta(kmer_protein_df, aa_sequences_df, output_path, feature_type):
    """
    Save FASTA file containing only proteins with predictive k-mers.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import time
    
    start_time = time.time()
    print("[DEBUG] Entering save_proteins_with_kmers_fasta", flush=True)
    
    logging.info(f"Creating FASTA file of proteins with predictive k-mers...")
    
    # Get unique proteins with k-mers
    print("[DEBUG] Getting unique protein IDs...", flush=True)
    proteins_with_kmers = set(kmer_protein_df['protein_ID'].unique())
    print(f"[DEBUG] Found {len(proteins_with_kmers):,} unique proteins", flush=True)
    
    # Filter sequences
    print("[DEBUG] Filtering and writing sequences...", flush=True)
    filtered_seqs = aa_sequences_df[aa_sequences_df['protein_ID'].isin(proteins_with_kmers)]
    
    # Create SeqRecord objects - simple, no metadata needed
    records = []
    for _, row in filtered_seqs.iterrows():
        record = SeqRecord(
            Seq(row['sequence']),
            id=row['protein_ID'],
            description=""  # Keep it simple
        )
        records.append(record)
    
    # Write FASTA
    print(f"[DEBUG] Writing {len(records):,} records to FASTA...", flush=True)
    SeqIO.write(records, output_path, 'fasta')
    print(f"[DEBUG] FASTA written successfully", flush=True)
    
    elapsed = time.time() - start_time
    logging.info(f"Saved {len(records):,} proteins to: {output_path}")
    logging.info(f"Reduction: {100*(1 - len(records)/len(aa_sequences_df)):.1f}% of original proteins")
    print(f"[DEBUG] save_proteins_with_kmers_fasta complete in {elapsed:.2f}s", flush=True)

def main():
    parser = argparse.ArgumentParser(
        description='Standalone k-mer analysis for ignore_families=True workflows'
    )
    
    # Required arguments
    parser.add_argument(
        '--modeling_dir',
        required=True,
        help='Path to modeling directory (contains modeling_results/)'
    )
    parser.add_argument(
        '--feature_selection_dir',
        required=True,
        help='Path to feature_selection directory (contains filtered_feature_tables/)'
    )
    parser.add_argument(
        '--feature2cluster_path',
        required=True,
        help='Path to selected_features.csv (or strain_selected_features.csv)'
    )
    parser.add_argument(
        '--protein_mapping_csv',
        required=True,
        help='Path to protein mapping CSV (strain_proteins.csv or phage_proteins.csv)'
    )
    parser.add_argument(
        '--aa_sequence_file',
        required=True,
        help='Path to combined FASTA file (strain_combined.faa or phage_combined.faa)'
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='Output directory for results'
    )
    
    # Feature type
    parser.add_argument(
        '--feature_type',
        required=True,
        choices=['strain', 'phage'],
        help='Type of features to analyze'
    )
    
    # Optional arguments
    parser.add_argument(
        '--cutoff',
        help='Specific cutoff to use (default: auto-detect best)'
    )
    parser.add_argument(
        '--use_clustering',
        action='store_true',
        help='Use MMseqs2 to cluster proteins before alignment'
    )
    parser.add_argument(
        '--min_seq_id',
        type=float,
        default=0.5,
        help='MMseqs2 minimum sequence identity (default: 0.5)'
    )
    parser.add_argument(
        '--coverage',
        type=float,
        default=0.8,
        help='MMseqs2 minimum coverage (default: 0.8)'
    )
    parser.add_argument(
        '--sensitivity',
        type=float,
        default=7.5,
        help='MMseqs2 sensitivity (default: 7.5)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads (default: 4)'
    )
    parser.add_argument(
        '--create_plots',
        action='store_true',
        help='Create visualization plots (requires alignment)'
    )
    parser.add_argument(
        '--skip_alignment',
        action='store_true',
        help='Skip alignment and segment analysis (faster, only produces k-mer-protein mappings)'
    )
    parser.add_argument(
        '--max_proteins_per_cluster',
        type=int,
        default=200,
        help='Maximum proteins per cluster to align (default: 200, prevents memory issues)'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.create_plots and args.skip_alignment:
        logging.error("Cannot create plots with --skip_alignment. Remove one of these flags.")
        return
    
    if args.create_plots and not args.use_clustering:
        logging.warning("--create_plots without --use_clustering may produce poor results.")
        logging.warning("Consider using --use_clustering to group related proteins before alignment.")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Get best cutoff
    if args.cutoff:
        cutoff = args.cutoff
    else:
        cutoff = get_best_cutoff(args.modeling_dir)
    
    feature_table_path = os.path.join(
        args.feature_selection_dir,
        'filtered_feature_tables',
        f'select_feature_table_cutoff_{cutoff}.csv'
    )
    
    logging.info(f"Using cutoff: {cutoff}")
    logging.info(f"Feature table: {feature_table_path}")
    
    # Step 1: Load predictive k-mers
    predictive_kmers = load_predictive_kmers(
        feature_table_path,
        args.feature2cluster_path,
        args.feature_type
    )
    
    # Step 2: Load protein sequences
    aa_sequences = load_aa_sequences(args.aa_sequence_file)
    
    # Step 3: Map k-mers to proteins
    kmer_protein_df = map_kmers_to_proteins(
        predictive_kmers,
        args.protein_mapping_csv,
        aa_sequences,
        args.feature_type
    )
    
    if kmer_protein_df.empty:
        logging.error("No k-mer-protein mappings found. Exiting.")
        return
    
    # Save k-mer-protein mappings
    kmer_protein_output = os.path.join(args.output_dir, f'{args.feature_type}_kmer_proteins.csv')
    kmer_protein_df.to_csv(kmer_protein_output, index=False)
    logging.info(f"Saved k-mer-protein mappings to: {kmer_protein_output}")

    # Save FASTA of proteins with k-mers
    fasta_output = os.path.join(args.output_dir, f'{args.feature_type}_proteins_with_kmers.faa')
    save_proteins_with_kmers_fasta(kmer_protein_df, aa_sequences, fasta_output, args.feature_type)
    logging.info(f"Saved FASTA of proteins with k-mers to: {fasta_output}")
    
    # Step 3b: Extract segments from raw sequences (NO ALIGNMENT NEEDED)
    logging.info("Extracting k-mer segments from raw sequences...")
    segments_raw = extract_segments_from_raw_sequences(kmer_protein_df, args.feature_type)
    
    if not segments_raw.empty:
        # Save individual segment positions (without full_protein_sequence column)
        segments_raw_output_df = segments_raw.drop(columns=['full_protein_sequence'])
        segments_raw_output = os.path.join(args.output_dir, f'{args.feature_type}_kmer_segments_raw.csv')
        segments_raw_output_df.to_csv(segments_raw_output, index=False)
        
        # Create merged/overlapping segments
        merged_segments = merge_overlapping_segments(segments_raw, args.feature_type)
        merged_output = os.path.join(args.output_dir, f'{args.feature_type}_kmer_segments_merged.csv')
        merged_segments.to_csv(merged_output, index=False)
        logging.info(f"Saved merged segments to: {merged_output}")
    else:
        logging.warning("No segments extracted from sequences")
    
    # Step 4: Optional clustering
    if args.use_clustering:
        kmer_protein_df = cluster_proteins_mmseqs(
            kmer_protein_df,
            args.output_dir,
            args.feature_type,
            args.threads,
            args.min_seq_id,
            args.coverage,
            args.sensitivity
        )
    else:
        logging.info("Clustering disabled. All proteins assigned to cluster 0.")
        kmer_protein_df['cluster_id'] = 0
    
    # Decide whether to run alignment
    run_alignment = not args.skip_alignment
    
    if run_alignment:
        # Check if alignment makes sense
        if not args.use_clustering:
            total_proteins = kmer_protein_df['protein_ID'].nunique()
            if total_proteins > args.max_proteins_per_cluster:
                logging.warning("=" * 60)
                logging.warning(f"WARNING: {total_proteins} proteins in single cluster!")
                logging.warning("Aligning unrelated proteins is not recommended.")
                logging.warning("Options:")
                logging.warning("  1. Use --use_clustering to group related proteins")
                logging.warning("  2. Use --skip_alignment to skip this step")
                logging.warning(f"  3. Increase --max_proteins_per_cluster (currently {args.max_proteins_per_cluster})")
                logging.warning("=" * 60)
                logging.warning("Skipping alignment due to cluster size limit.")
                run_alignment = False
    
    if run_alignment:
        # Step 5: Align sequences
        logging.info("Aligning sequences...")
        aligned_dfs = []
        
        total_clusters = kmer_protein_df['cluster_id'].nunique()
        successful_alignments = 0
        skipped_large = 0
        skipped_small = 0
        
        for i, (cluster_id, cluster_group) in enumerate(kmer_protein_df.groupby('cluster_id'), 1):
            unique_proteins = cluster_group[['protein_ID', 'sequence']].drop_duplicates()
            n_proteins = len(unique_proteins)
            
            if n_proteins < 2:
                logging.info(f"Cluster {i}/{total_clusters} (id={cluster_id}): Skipping - only {n_proteins} protein(s)")
                skipped_small += 1
                continue
            
            if n_proteins > args.max_proteins_per_cluster:
                logging.warning(f"Cluster {i}/{total_clusters} (id={cluster_id}): Skipping - {n_proteins} proteins exceeds max ({args.max_proteins_per_cluster})")
                skipped_large += 1
                continue
            
            logging.info(f"Cluster {i}/{total_clusters} (id={cluster_id}): Aligning {n_proteins} proteins")
            
            seqs = [(row['protein_ID'], row['sequence']) 
                    for _, row in unique_proteins.iterrows()]
            
            aligned_df = align_sequences(
                seqs,
                args.output_dir,
                f"{args.feature_type}_cluster{cluster_id}"
            )
            
            if not aligned_df.empty:
                aligned_df['cluster_id'] = cluster_id
                aligned_dfs.append(aligned_df)
                successful_alignments += 1
        
        logging.info("=" * 60)
        logging.info(f"Alignment summary:")
        logging.info(f"  Successful: {successful_alignments}")
        logging.info(f"  Skipped (too small): {skipped_small}")
        logging.info(f"  Skipped (too large): {skipped_large}")
        logging.info("=" * 60)
        
        if not aligned_dfs:
            logging.error("No successful alignments. Cannot proceed with segment analysis.")
            if args.create_plots:
                logging.error("Cannot create plots without alignments.")
            logging.info("Main output (k-mer-protein mappings) has been saved successfully.")
            return
        
        aligned_combined = pd.concat(aligned_dfs, ignore_index=True)
        
        # Step 6: Find k-mer positions
        logging.info("Finding k-mer positions in aligned sequences...")
        aligned_with_kmers = aligned_combined.merge(
            kmer_protein_df[['protein_ID', 'Feature', 'cluster', 'kmer', 'cluster_id']],
            on=['protein_ID', 'cluster_id'],
            how='inner'
        )
        
        aligned_with_kmers[['start_indices', 'stop_indices']] = aligned_with_kmers.apply(
            find_kmer_indices, axis=1
        )
        
        # Save aligned data
        aligned_output = os.path.join(args.output_dir, f'{args.feature_type}_aligned_kmers.csv')
        aligned_with_kmers.to_csv(aligned_output, index=False)
        logging.info(f"Saved aligned k-mers to: {aligned_output}")
        
        # Step 7: Calculate coverage
        coverage_df = calculate_coverage(aligned_with_kmers)
        
        if coverage_df.empty:
            logging.error("No coverage data generated.")
            logging.info("Main output (k-mer-protein mappings) has been saved successfully.")
            return
        
        # Step 8: Identify segments
        segments_df = identify_segments(coverage_df)
        
        # Step 9: Add genome information
        segments_df = add_genome_info(segments_df, kmer_protein_df, args.feature_type)
        
        # Step 10: Extract segment sequences
        segments_with_seqs = extract_segment_sequences(segments_df, aligned_combined)
        
        # Save segments
        segments_output = os.path.join(args.output_dir, f'{args.feature_type}_segments.csv')
        segments_with_seqs.to_csv(segments_output, index=False)
        logging.info(f"Saved segments to: {segments_output}")
        
        # Step 11: Optional plotting
        if args.create_plots:
            plot_dir = os.path.join(args.output_dir, 'plots')
            plot_segments(segments_df, aligned_combined, plot_dir, args.feature_type)
        
        # Summary statistics with alignment
        logging.info("=" * 60)
        logging.info("ANALYSIS SUMMARY (with alignment)")
        logging.info("=" * 60)
        logging.info(f"Feature type: {args.feature_type}")
        logging.info(f"Predictive features: {kmer_protein_df['Feature'].nunique()}")
        logging.info(f"Unique k-mers: {kmer_protein_df['kmer'].nunique()}")
        logging.info(f"Proteins analyzed: {kmer_protein_df['protein_ID'].nunique()}")
        logging.info(f"Genomes represented: {kmer_protein_df[args.feature_type].nunique()}")
        logging.info(f"Raw segments found: {len(segments_raw) if not segments_raw.empty else 0}")
        logging.info(f"Merged segments: {len(merged_segments) if not segments_raw.empty else 0}")
        logging.info(f"Proteins aligned: {aligned_combined['protein_ID'].nunique()}")
        logging.info(f"Alignment-based segments: {len(segments_df)}")
        logging.info(f"Covered segments: {len(segments_df[segments_df['coverage'] > 0])}")
        logging.info("=" * 60)
    else:
        # No alignment performed
        logging.info("=" * 60)
        logging.info("ANALYSIS SUMMARY (k-mer-protein mappings only)")
        logging.info("=" * 60)
        logging.info(f"Feature type: {args.feature_type}")
        logging.info(f"Predictive features: {kmer_protein_df['Feature'].nunique()}")
        logging.info(f"Unique k-mers: {kmer_protein_df['kmer'].nunique()}")
        logging.info(f"Proteins with k-mers: {kmer_protein_df['protein_ID'].nunique()}")
        logging.info(f"Genomes represented: {kmer_protein_df[args.feature_type].nunique()}")
        logging.info(f"Raw segments found: {len(segments_raw) if not segments_raw.empty else 0}")
        logging.info(f"Merged segments: {len(merged_segments) if not segments_raw.empty else 0}")
        logging.info("")
        logging.info("Alignment and alignment-based segment analysis skipped.")
        logging.info("Raw segment coordinates (in original sequences) have been saved.")
        logging.info("To perform alignment:")
        logging.info("  - Remove --skip_alignment flag")
        logging.info("  - Use --use_clustering to group related proteins")
        logging.info("=" * 60)


if __name__ == "__main__":
    main()
