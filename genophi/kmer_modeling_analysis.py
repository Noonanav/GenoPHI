import os
import re
import shutil
import subprocess
import logging
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plotnine import *
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ==========================
# 1. DATA LOADING
# ==========================

def load_aa_sequences(aa_sequence_file):
    logging.info(f"Loading amino acid sequences from {aa_sequence_file}")
    aa_records = list(SeqIO.parse(aa_sequence_file, 'fasta'))
    aa_sequences_df = pd.DataFrame({
        'protein_ID': [record.id for record in aa_records],
        'sequence': [str(record.seq) for record in aa_records]
    })
    logging.info(f"Loaded {len(aa_sequences_df)} sequences.")
    return aa_sequences_df

def aggregate_shap_values(model_output_dir):
    """Aggregates SHAP importance values from multiple model runs."""
    runs = [x for x in os.listdir(model_output_dir) if 'run' in x]
    top_models_shap_df = pd.DataFrame()
    for run in runs:
        shap_values_csv_path = os.path.join(model_output_dir, run, 'shap_importances.csv')
        if os.path.exists(shap_values_csv_path):
            shap_values_temp = pd.read_csv(shap_values_csv_path)
            shap_values_temp = shap_values_temp.groupby(['feature', 'value']).agg({'shap_value': 'median'}).reset_index()
            top_models_shap_df = pd.concat([top_models_shap_df, shap_values_temp], ignore_index=True)
    return top_models_shap_df

# ==========================
# 2. FEATURE EXTRACTION & MAPPING
# ==========================

def get_predictive_kmers(feature_file_path, feature2cluster_path, feature_type, ignore_families=False):
    logging.info(f"Extracting predictive features of type '{feature_type}'")
    feature_df = pd.read_csv(feature_file_path)
    
    feature_prefix = feature_type[0] + 'c_' 
    select_features = [col for col in feature_df.columns if feature_prefix in col]

    if not select_features:
        logging.warning(f"No predictive {feature_type} features found.")
        return pd.DataFrame(columns=['cluster', 'kmer', 'protein_family', 'Feature'])

    feature2cluster_df = pd.read_csv(feature2cluster_path)
    feature2cluster_df.rename(columns={'Cluster_Label': 'cluster'}, inplace=True)
    
    filtered_kmers = feature2cluster_df[feature2cluster_df['Feature'].isin(select_features)].copy()
    
    if ignore_families:
        # In ignore_families mode, the 'cluster' column IS the k-mer sequence
        # We preserve the 'Feature' column (sc_*) as is
        filtered_kmers['kmer'] = filtered_kmers['cluster']
        filtered_kmers['protein_family'] = 'Unclustered' 
    else:
        # Standard mode: cluster = family_kmer
        filtered_kmers['kmer'] = filtered_kmers['cluster'].str.split('_').str[-1]
        filtered_kmers['protein_family'] = filtered_kmers['cluster'].str.split('_').str[:-1].str.join('_')
    
    logging.info(f"Filtered down to {len(filtered_kmers)} predictive k-mers.")
    return filtered_kmers

def merge_kmers_with_families(protein_families_file, aa_sequences_df, feature_type='strain', ignore_families=False, filtered_kmers=None):
    """
    Merges k-mer data with proteins.
    In Standard Mode: Joins using the protein_families_file.
    In Ignore Mode: Scans sequences for k-mers and maps back to Feature ID.
    """
    
    # --- PATH A: STANDARD (Use Lookup File) ---
    if not ignore_families:
        logging.info(f"Merging k-mer data via Protein Families file...")
        if not protein_families_file:
            logging.error("Protein families file is missing!")
            return pd.DataFrame()

        protein_families_df = pd.read_csv(protein_families_file)
        protein_families_df = protein_families_df[[feature_type, 'cluster', 'protein_ID']].drop_duplicates()
        protein_families_df.rename(columns={'cluster': 'protein_family'}, inplace=True)
        
        # Return Families + Sequences (Features merged later in workflow)
        merged_df = protein_families_df.merge(aa_sequences_df, on='protein_ID', how='inner')
        return merged_df

    # --- PATH B: SEARCH MODE (No Lookup File) ---
    else:
        logging.info("ignore_families=True: Searching AA sequences for k-mers...")
        if filtered_kmers is None or filtered_kmers.empty:
            return pd.DataFrame()

        # Unique K-mers to search for
        target_kmers = filtered_kmers['kmer'].unique()
        hits = []
        
        for _, row in aa_sequences_df.iterrows():
            seq = row['sequence']
            pid = row['protein_ID']
            
            for kmer in target_kmers:
                if kmer in seq:
                    # Get the Feature ID (e.g. sc_123) corresponding to this kmer
                    # Handle case where one kmer might map to multiple Features (rare but possible)
                    matching_features = filtered_kmers[filtered_kmers['kmer'] == kmer]['Feature'].values
                    
                    for feat in matching_features:
                        hits.append({
                            'protein_family': 'Unclustered', 
                            'protein_ID': pid,
                            'sequence': seq,
                            'kmer': kmer,
                            'Feature': feat 
                        })
        
        if not hits:
            logging.warning("No proteins found containing predictive k-mers.")
            return pd.DataFrame()

        merged_df = pd.DataFrame(hits)
        logging.info(f"Found {len(merged_df)} hits across {merged_df['protein_ID'].nunique()} proteins.")
        return merged_df

# ==========================
# 3. CLUSTERING & ALIGNMENT
# ==========================

def cluster_sequences_mmseqs(seq_df, output_dir, threads=4, min_seq_id=0.4, coverage=0.8):
    """Clusters sequences using MMseqs2."""
    logging.info("Clustering sequences with MMseqs2...")
    cluster_dir = os.path.join(output_dir, "mmseqs_clustering")
    if os.path.exists(cluster_dir):
        shutil.rmtree(cluster_dir)
    os.makedirs(cluster_dir, exist_ok=True)

    input_fasta = os.path.join(cluster_dir, "input.fasta")
    unique_seqs = seq_df[['protein_ID', 'sequence']].drop_duplicates()
    
    with open(input_fasta, 'w') as f:
        for _, row in unique_seqs.iterrows():
            f.write(f">{row['protein_ID']}\n{row['sequence']}\n")

    db_path = os.path.join(cluster_dir, "DB")
    cluster_db = os.path.join(cluster_dir, "DB_clu")
    tsv_path = os.path.join(cluster_dir, "clusters.tsv")
    tmp_dir = os.path.join(cluster_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    try:
        subprocess.run(['mmseqs', 'createdb', input_fasta, db_path], check=True, stdout=subprocess.DEVNULL)
        subprocess.run(['mmseqs', 'cluster', db_path, cluster_db, tmp_dir, 
                        '--min-seq-id', str(min_seq_id), '-c', str(coverage), '--threads', str(threads)], 
                        check=True, stdout=subprocess.DEVNULL)
        subprocess.run(['mmseqs', 'createtsv', db_path, db_path, cluster_db, tsv_path], 
                        check=True, stdout=subprocess.DEVNULL)
        
        clusters_df = pd.read_csv(tsv_path, sep="\t", header=None, names=["cluster_rep", "protein_ID"])
        clusters_df["cluster_id"] = clusters_df.groupby("cluster_rep").ngroup()
        
        logging.info(f"MMseqs2 created {clusters_df['cluster_id'].nunique()} clusters.")
        return clusters_df[["protein_ID", "cluster_id"]]

    except Exception as e:
        logging.error(f"MMseqs2 clustering failed: {e}")
        return pd.DataFrame({'protein_ID': unique_seqs['protein_ID'], 'cluster_id': 0})

def align_sequences(sequences, output_dir, family_name):
    """Aligns sequences using MAFFT via subprocess."""
    logging.info(f"Aligning {len(sequences)} sequences for group: {family_name}")
    
    if len(sequences) < 2:
        return pd.DataFrame()

    safe_name = str(family_name).replace('|', '_').replace('/', '_')
    temp_fasta_path = os.path.join(output_dir, f"{safe_name}_temp.fasta")
    temp_aln_path = os.path.join(output_dir, f"{safe_name}_temp.aln")

    with open(temp_fasta_path, "w") as f:
        for pid, seq in sequences:
            f.write(f">{pid}\n{seq}\n")

    try:
        cmd = ["mafft", "--auto", "--quiet", temp_fasta_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        with open(temp_aln_path, "w") as f:
            f.write(result.stdout)
            
        alignment = AlignIO.read(temp_aln_path, "fasta")
        
        aln_len = alignment.get_alignment_length()
        non_gap_pos = [any(rec.seq[i] != '-' for rec in alignment) for i in range(aln_len)]
        first_non_gap = non_gap_pos.index(True) if True in non_gap_pos else 0
        
        aligned_data = []
        for record in alignment:
            trimmed = str(record.seq[first_non_gap:])
            start_pos = trimmed.find(trimmed.lstrip('-'))
            aligned_data.append({
                'protein_ID': record.id,
                'aln_sequence': trimmed,
                'start_index': start_pos
            })
            
        if os.path.exists(temp_fasta_path): os.remove(temp_fasta_path)
        if os.path.exists(temp_aln_path): os.remove(temp_aln_path)
        return pd.DataFrame(aligned_data)

    except Exception as e:
        logging.error(f"Alignment failed for {family_name}: {e}")
        return pd.DataFrame()

# ==========================
# 4. ANALYSIS (INDICES & COVERAGE)
# ==========================

def construct_kmer_id_df(protein_families_df, kmer_df):
    """Legacy helper."""
    return kmer_df

def find_kmer_indices(row):
    """Finds k-mer indices allowing for gaps (e.g. A-B-C)."""
    kmer = row['kmer']
    seq = row['aln_sequence']
    
    pattern = ""
    for char in kmer:
        pattern += char + "-*"
    pattern = pattern[:-2]
    
    start_indices = []
    stop_indices = []
    
    for match in re.finditer(pattern, seq):
        start = match.start()
        stop = match.end() - 1
        segment = seq[start:stop+1]
        aa_count = sum(1 for c in segment if c != '-')
        if aa_count == len(kmer):
            start_indices.append(start)
            stop_indices.append(stop)
            
    return pd.Series([start_indices, stop_indices], index=['start_indices', 'stop_indices'])

def calculate_coverage(df):
    """Calculates coverage and tracks gaps."""
    logging.info("Calculating coverage...")
    coverage_data = []
    
    for _, row in df.iterrows():
        seq_len = len(row['aln_sequence'])
        coverage = np.zeros(seq_len, dtype=int)
        
        # Only fill if matches exist
        if len(row['start_indices']) > 0:
            for start, stop in zip(row['start_indices'], row['stop_indices']):
                coverage[start:stop+1] = 1
            
        for idx, residue in enumerate(row['aln_sequence']):
            coverage_data.append({
                'Feature': row.get('Feature', 'Unknown'),
                'protein_ID': row['protein_ID'],
                'cluster_id': row.get('cluster_id', 0),
                'AA_index': idx,
                'coverage': coverage[idx],
                'is_gap': 1 if residue == '-' else 0
            })
            
    if not coverage_data:
        return pd.DataFrame()
        
    cov_df = pd.DataFrame(coverage_data)
    cov_df = cov_df.groupby(['Feature', 'protein_ID', 'cluster_id', 'AA_index']).agg({
        'coverage': 'max', 
        'is_gap': 'max'
    }).reset_index()
    
    return cov_df

def identify_segments(df):
    """Identifies segments (0=Uncovered, 1=Covered, 2=Gap)."""
    logging.info("Identifying segments...")
    if df.empty: return pd.DataFrame()
    
    df = df.sort_values(by=['Feature', 'protein_ID', 'AA_index'])
    
    # Segment type logic
    df['segment_type'] = df['coverage']
    df.loc[df['is_gap'] == 1, 'segment_type'] = 2
    
    df['prev_type'] = df.groupby(['Feature', 'protein_ID'])['segment_type'].shift(1).fillna(-1)
    df['segment_change'] = (df['segment_type'] != df['prev_type'])
    df['segment_id'] = df.groupby(['Feature', 'protein_ID'])['segment_change'].cumsum()
    
    segments = df.groupby(['Feature', 'protein_ID', 'cluster_id', 'segment_type', 'segment_id']).agg(
        start=('AA_index', 'min'),
        stop=('AA_index', 'max')
    ).reset_index()
    
    return segments

# ==========================
# 5. PLOTTING
# ==========================

def plot_segments(segments_df, output_dir):
    """Creates faceted plots based on cluster_id."""
    logging.info(f"Plotting segments to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    
    if segments_df.empty:
        return

    # Separate Types
    uncovered = segments_df[segments_df['segment_type'] == 0]
    covered = segments_df[segments_df['segment_type'] == 1]
    
    n_proteins = segments_df['protein_ID'].nunique()
    fig_height = min(max(5, n_proteins * 0.4), 25)
    
    plot = (
        ggplot() +
        geom_segment(
            data=uncovered,
            mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID'),
            color='grey', size=5
        ) +
        geom_segment(
            data=covered,
            mapping=aes(x='start', xend='stop', y='protein_ID', yend='protein_ID', color='Feature'),
            size=5
        ) +
        facet_wrap('~cluster_id', dir='v', scales='free_y', ncol=1) +
        labs(title="K-mer Coverage by Sequence Cluster", x="AA Index", y="Protein ID") +
        theme(
            axis_text_x=element_text(rotation=90),
            panel_background=element_rect(fill='white'),
            figure_size=(12, fig_height),
            strip_background=element_rect(fill='#f0f0f0')
        )
    )
    
    try:
        plot.save(os.path.join(output_dir, "coverage_plot.png"))
        logging.info("Plot saved successfully.")
    except Exception as e:
        logging.error(f"Plotting failed: {e}")

def merge_no_coverage_proteins(segments, aligned):
    return segments