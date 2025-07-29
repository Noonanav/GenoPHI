from hdbscan import HDBSCAN
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt
# from umap import UMAP
from scipy.cluster.hierarchy import linkage, fcluster
import os
from plotnine import ggplot, aes, geom_point, labs, theme, element_text, geom_col, scale_fill_gradient, lims, facet_grid
import argparse
from IPython.display import display

def plot_phage_umap(feature_table, select_phages, phage_col='phage', 
                   n_neighbors=15, n_components=2, min_dist=0.5, 
                   split_by_strain=False,
                   plot_title="UMAP Projection of Phages by Cluster", 
                   save_plot=False, plot_filename='phage_umap_projection.png'):
    """
    Plots a UMAP projection of phages, colored by their pre-existing cluster assignments based on feature content.

    Args:
        feature_table (DataFrame): DataFrame containing phage features and cluster labels.
        phage_col (str): Column name for phage identifier. Default is 'phage'.
        n_neighbors (int): The size of local neighborhood used for manifold approximation. Default is 15.
        n_components (int): The dimension of the space to embed into. Default is 2.
        min_dist (float): The effective minimum distance between embedded points. Default is 0.1.
        plot_title (str): Title for the plot. Default is "UMAP Projection of Phages by Cluster".
        save_plot (bool): Whether to save the plot. Default is False.
        plot_filename (str): Filename to save the plot if save_plot is True. Default is 'phage_umap_projection.png'.

    Returns:
        None: This function plots the UMAP projection directly or saves it.
    """
    # Extract feature columns (e.g., those starting with 'pc_')
    phage_feature_columns = [col for col in feature_table.columns if 'pc_' in col]
    
    # Ensure 'cluster' column exists
    if 'cluster' not in feature_table.columns:
        raise ValueError("The feature_table must contain a 'cluster' column.")
    
    # Apply UMAP for dimensionality reduction
    umap_reducer = UMAP(n_neighbors=n_neighbors, n_components=n_components, 
                        min_dist=min_dist, random_state=42)
    umap_results = umap_reducer.fit_transform(feature_table[phage_feature_columns])
    print(umap_results.shape)
    
    # Prepare the DataFrame for plotting. Set the index to the phage identifier.
    umap_df = pd.DataFrame({
        'UMAP_1': umap_results[:, 0],
        'UMAP_2': umap_results[:, 1],
        'cluster': feature_table['cluster'],
        phage_col: feature_table[phage_col]
    })
    # display(umap_df.head())
    
    # Avoid modifying select_phages in place:
    select_phages = select_phages.copy()
    select_phages['select'] = True
    # display(select_phages.head())
    
    # Build the plot using Plotnine (ggplot)
    if split_by_strain:
        # Merge on the phage identifier (which is now both the index of umap_df and a column in select_phages)
        umap_df = umap_df.merge(select_phages[[phage_col, 'strain', 'select']], on=phage_col, how='left')
        umap_df['cluster'] = umap_df['cluster'].astype(str)
        # display(umap_df.head())

        for strain in umap_df['strain'].unique():
            umap_df_strain = umap_df.copy()
            umap_df_strain['strain'] = umap_df_strain['strain'].fillna(strain)
            umap_df_strain = umap_df_strain[umap_df_strain['strain'] == strain]
            umap_plot = (
                ggplot() +
                geom_point(umap_df_strain[umap_df_strain['select'] != True], aes(x='UMAP_1', y='UMAP_2', color='cluster')) +
                geom_point(umap_df_strain[umap_df_strain['select'] == True], aes(x='UMAP_1', y='UMAP_2'), color='black', size=3) +
                labs(title=strain, x='UMAP 1', y='UMAP 2') +
                theme(figure_size=(3, 3),
                      legend_position='none'
                      )
            )
    
            print(umap_plot)
    else:
        # Merge on the phage identifier (which is now both the index of umap_df and a column in select_phages)
        umap_df = umap_df.merge(select_phages[[phage_col, 'select']].drop_duplicates(), on=phage_col, how='left')
        umap_df['cluster'] = umap_df['cluster'].astype(str)
        display(umap_df.head())
        
        umap_plot = (
                ggplot() +
                geom_point(umap_df[umap_df['select'] != True], aes(x='UMAP_1', y='UMAP_2', color='cluster')) +
                geom_point(umap_df[umap_df['select'] == True], aes(x='UMAP_1', y='UMAP_2'), color='black', size=3) +
                labs(title=strain, x='UMAP 1', y='UMAP 2') +
                theme(figure_size=(8, 7))
            )
    
        print(umap_plot)
    
    if save_plot:
        umap_plot.save(plot_filename, dpi=300)
    
def plot_phage_clustermap(feature_table, phage_col='phage', plot_title="Phage-Phage Clustermap", save_plot=False, plot_filename='phage_clustermap.png'):
    """
    Plots a clustermap of phages based on Jaccard distances of their feature content.

    Args:
        feature_table (DataFrame): DataFrame containing phage features.
        phage_col (str): Column name for phage identifier. Default is 'phage'.
        plot_title (str): Title for the plot. Default is "Phage-Phage Clustermap".
        save_plot (bool): Whether to save the plot. Default is False.
        plot_filename (str): Filename to save the plot if save_plot is True. Default is 'phage_clustermap.png'.

    Returns:
        None: This function plots the clustermap directly or saves it.
    """
    # Extract phage feature columns
    phage_feature_columns = [col for col in feature_table.columns if 'pc_' in col]
    
    # Create a subset of the feature table for phage clustering
    phage_feature_table = feature_table[[phage_col] + phage_feature_columns].drop_duplicates()
    
    # Convert feature table to a binary matrix
    phage_matrix = phage_feature_table.set_index(phage_col).applymap(lambda x: 1 if x > 0 else 0)
    
    # Compute Jaccard distances
    jaccard_distances = pdist(phage_matrix, metric='jaccard')
    jaccard_distances_square = squareform(jaccard_distances)
    
    # Create a DataFrame from the distance matrix
    jaccard_df = pd.DataFrame(jaccard_distances_square, index=phage_matrix.index, columns=phage_matrix.index)
    
    # Plot the clustermap
    plt.figure(figsize=(6, 6))
    # Use 'jaccard_df.values' instead of 'jaccard_df' to pass the actual distance matrix
    clustermap = sns.clustermap(jaccard_df.values, cmap='viridis', method='ward', metric='euclidean')
    
    # Set the title
    clustermap.fig.suptitle(plot_title, fontsize=16)
    plt.subplots_adjust(top=0.95)  # Adjust the top to make room for the title
    
    # Optionally save the plot
    if save_plot:
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    
    # Show the plot
    plt.show()

def cluster_phages(
    phage_feature_table: pd.DataFrame, 
    method: str = 'HDBSCAN',
    top_n: int = 5,
    verbose=False,
    **kwargs
) -> pd.DataFrame:
    """
    Clusters phages based on provided feature columns using HDBSCAN.
    
    Args:
        phage_feature_table (DataFrame): DataFrame with phage features.
        feature_columns (list): List of columns in phage_feature_table to use for clustering.
        **kwargs: Additional keyword arguments for HDBSCAN.
        
    Returns:
        DataFrame: A copy of phage_feature_table with an added 'cluster' column.
    """
    feature_columns = [x for x in phage_feature_table.columns if 'pc_' in x]
    
    # Add a check to ensure we have feature columns
    if not feature_columns:
        raise ValueError("No feature columns found with 'pc_' in their name. Please ensure your feature table has properly named columns.")
    
    phage_feature_table = phage_feature_table[['phage'] + feature_columns].drop_duplicates()

    if method == 'HDBSCAN':
        if verbose:
            print(f"Clustering phages using HDBSCAN with parameters: {kwargs}")
        clusterer = HDBSCAN(**kwargs)
        cluster_labels = clusterer.fit_predict(phage_feature_table[feature_columns])

        if verbose:
            print(f"Number of clusters: {len(np.unique(cluster_labels))}")
            print(f"Number of noise points: {np.sum(cluster_labels == -1)}")
        
        # Handle noise points (label -1) by assigning unique cluster IDs
        max_cluster_label = cluster_labels.max()
        cluster_labels = list(cluster_labels)  # ensure mutable list
        for i, label in enumerate(cluster_labels):
            if label == -1:
                max_cluster_label += 1
                cluster_labels[i] = max_cluster_label
                
    if method == 'hierarchical':
        if verbose:
            print(f"Clustering phages using hierarchical clustering")

        # Calculate pairwise distances
        phage_distances = pdist(phage_feature_table[feature_columns], metric='euclidean')
        
        # Perform hierarchical clustering
        linkage_matrix = linkage(phage_distances, method='ward')
        
        # Asing top_n clusters
        cluster_n = top_n
        cluster_labels = fcluster(linkage_matrix, cluster_n, criterion='maxclust')
            
    phage_feature_table = phage_feature_table.copy()
    phage_feature_table['cluster'] = cluster_labels
    return phage_feature_table

def cocktail_design(df, feature_table, 
                   strain_col='strain', 
                   phage_col='phage', 
                   interaction_col='interaction', 
                   confidence_col='Confidence', 
                   grouping_col=None, 
                   top_n=5, 
                   clustering=False,
                   clustering_method='HDBSCAN',
                   plot_title="Cocktail Design Analysis", 
                   display_intermediate=False, 
                   generate_plot=True,
                   plot_umap=False,
                   split_by_strain=False,
                   save_plot=False, 
                   plot_filename='cocktail_design_plot.png',
                   umap_filename='phage_umap_projection.png',
                   umap_save_plot=False,
                   default_cocktails=False,
                   balanced_selection=False,
                   lambda_=0.5,
                   verbose=False):
    """
    Designs and analyzes phage cocktails based on interaction data.
    """
    # Define grouping columns
    grouping_cols = [strain_col] + ([grouping_col] if grouping_col else [])
    max_interactions = top_n
    
    # Calculate total interactions per group
    cocktails_filter = df.copy()
    cocktails_filter['total_interactions'] = cocktails_filter.groupby(grouping_cols)[interaction_col].transform('sum')
    
    if verbose:
        display(f'Unique {strain_col}s:' + str(cocktails_filter[strain_col].nunique()))
        display(f'Unique {phage_col}s:' + str(cocktails_filter[phage_col].nunique()))

    if display_intermediate:
        display(cocktails_filter.sort_values(by=grouping_cols, ascending=False).head(30))

    # Filter out strains/phages with no interactions
    cocktails_filter = cocktails_filter[cocktails_filter['total_interactions'] > 0]
    if verbose:
        display(f'Unique {strain_col}s with interactions:' + str(cocktails_filter[strain_col].nunique()))
        display(f'Unique {phage_col}s:' + str(cocktails_filter[phage_col].nunique()))

    # Cap total_interactions at max_interactions
    cocktails_filter['total_interactions'] = cocktails_filter['total_interactions'].apply(lambda x: min(x, max_interactions))

    if default_cocktails:
        phage_percent = feature_table[[strain_col, phage_col, interaction_col]].drop_duplicates()
        phage_percent = phage_percent.groupby(phage_col).agg(
            sum_interaction=(interaction_col, 'sum'),
            total_count=(interaction_col, 'count')
        ).reset_index()
        phage_percent[confidence_col] = phage_percent['sum_interaction'] / phage_percent['total_count']
        phage_percent = phage_percent[[phage_col, confidence_col]].drop_duplicates()
        cocktails_filter = cocktails_filter.drop(columns=[confidence_col])
        cocktails_filter = cocktails_filter.merge(phage_percent, on=phage_col, how='left')

    if balanced_selection:
        phage_percent = feature_table[[strain_col, phage_col, interaction_col]].drop_duplicates()
        phage_percent = phage_percent.groupby(phage_col).agg(
            sum_interaction=(interaction_col, 'sum'),
            total_count=(interaction_col, 'count')
        ).reset_index()
        phage_percent['coverage'] = phage_percent['sum_interaction'] / phage_percent['total_count']
        phage_percent = phage_percent[[phage_col, 'coverage']].drop_duplicates()

        cocktails_filter = cocktails_filter.merge(phage_percent, on=phage_col, how='left')
        cocktails_filter[confidence_col] = lambda_ * cocktails_filter[confidence_col] + (1 - lambda_) * cocktails_filter['coverage']

    if clustering:
        feature_columns = [x for x in feature_table.columns if 'pc_' in x]
        if not feature_columns:
            print("Warning: No feature columns found for clustering. Falling back to non-clustering method.")
            clustering = False
        if clustering_method == 'HDBSCAN':
            phage_clusters = cluster_phages(feature_table, method = 'HDBSCAN', min_cluster_size=2, verbose=verbose)
            cocktails_filter = cocktails_filter.merge(phage_clusters[['phage', 'cluster']], on='phage', how='left')
            cocktails_filter = cocktails_filter.sort_values(by=confidence_col, ascending=False)
            cocktails_filter_select = cocktails_filter.groupby(grouping_cols + ['cluster']).head(1).reset_index(drop=True)
            cocktails_top_n = cocktails_filter_select.groupby(grouping_cols).head(top_n).reset_index(drop=True)

            # Check if any strain has less than top_n phages
            while any(count < top_n for count in cocktails_top_n.groupby(strain_col).size()):
                print('Not enough phages to fill top_n for all strains')
                cocktails_filter_select['select'] = True
                cocktail_filter_remaining = cocktails_filter.merge(
                    cocktails_filter_select[grouping_cols + [phage_col, 'select']], 
                    on=[strain_col, phage_col], 
                    how='left'
                )
                cocktail_filter_remaining = cocktail_filter_remaining[cocktail_filter_remaining['select'] != True]
                display(cocktail_filter_remaining.head())
                cocktail_filter_remaining = cocktail_filter_remaining.groupby(grouping_cols).head(1).reset_index(drop=True)
                print('Selected ', len(cocktail_filter_remaining), 'phages')
                cocktails_top_n = pd.concat([cocktails_top_n, cocktail_filter_remaining], ignore_index=True)
                cocktails_filter_select = pd.concat([cocktails_filter_select, cocktail_filter_remaining], ignore_index=True)

        elif clustering_method == 'hierarchical':
            phage_clusters = cluster_phages(feature_table, method = 'hierarchical', top_n=top_n)
            cocktails_filter = cocktails_filter.merge(phage_clusters[['phage', 'cluster']], on='phage', how='left')
            cocktails_filter = cocktails_filter.sort_values(by=confidence_col, ascending=False)
            cocktails_filter_select = cocktails_filter.groupby(grouping_cols + ['cluster']).head(1).reset_index(drop=True)
            cocktails_top_n = cocktails_filter_select.groupby(grouping_cols).head(top_n).reset_index(drop=True)

    else:
        cocktails_filter = cocktails_filter.sort_values(by=confidence_col, ascending=False)
        cocktails_top_n = cocktails_filter.groupby(grouping_cols).head(top_n).reset_index(drop=True)

    if plot_umap: 
        plot_phage_umap(phage_clusters, cocktails_top_n, split_by_strain=split_by_strain, n_neighbors=10, n_components=2, min_dist=0.4, plot_filename=umap_filename, save_plot=umap_save_plot)

    if display_intermediate:
        # Fix: Only include grouping_col in sort if it's not None
        sort_columns = [strain_col, confidence_col]
        if grouping_col:
            sort_columns = [grouping_col] + sort_columns
        display(cocktails_top_n.sort_values(by=sort_columns, ascending=False).head(30))

    # Aggregate interaction counts
    if grouping_col:
        cocktails_sorted_top_n_tally = cocktails_top_n.groupby([grouping_col, strain_col, 'total_interactions']).agg({interaction_col: 'sum'}).reset_index()
        cocktails_sorted_top_n_tally = cocktails_sorted_top_n_tally.rename(columns={interaction_col: 'interaction_count'})

        # Filter for incomplete interactions if needed
        counts_df_cocktail_filter = cocktails_sorted_top_n_tally[cocktails_sorted_top_n_tally['interaction_count'] != cocktails_sorted_top_n_tally['total_interactions']]
        if display_intermediate:
            display(counts_df_cocktail_filter.sort_values(by=[grouping_col, 'interaction_count'], ascending=False))

        # Calculate distribution of interaction counts
        cocktails_sorted_top_n_tally_size = cocktails_sorted_top_n_tally.groupby([grouping_col, 'interaction_count']).size().reset_index()

    else:
        cocktails_sorted_top_n_tally = cocktails_top_n.groupby([strain_col, 'total_interactions']).agg({interaction_col: 'sum'}).reset_index()
        cocktails_sorted_top_n_tally = cocktails_sorted_top_n_tally.rename(columns={interaction_col: 'interaction_count'})

        # Filter for incomplete interactions if needed
        counts_df_cocktail_filter = cocktails_sorted_top_n_tally[cocktails_sorted_top_n_tally['interaction_count'] != cocktails_sorted_top_n_tally['total_interactions']]
        if display_intermediate:
            display(counts_df_cocktail_filter.sort_values(by=['interaction_count'], ascending=False))

        # Calculate distribution of interaction counts
        cocktails_sorted_top_n_tally_size = cocktails_sorted_top_n_tally.groupby(['interaction_count']).size().reset_index()

    cocktails_sorted_top_n_tally_size = cocktails_sorted_top_n_tally_size.rename(columns={0: 'count'})
    cocktails_sorted_top_n_tally_size['percent'] = cocktails_sorted_top_n_tally_size['count'] / cocktails_sorted_top_n_tally_size['count'].sum()

    # Calculate fraction identified
    counts_df_cocktail = cocktails_sorted_top_n_tally.copy()
    counts_df_cocktail['fraction_identified'] = counts_df_cocktail['interaction_count'] / counts_df_cocktail['total_interactions']

    if display_intermediate:
        sort_columns = ['interaction_count']
        if grouping_col:
            sort_columns = [grouping_col] + sort_columns
        display(counts_df_cocktail.sort_values(by=sort_columns, ascending=False))

    # Prepare for bootstrap metrics
    counts_df_cocktail_summary = counts_df_cocktail.copy()
    counts_df_cocktail_summary['True'] = (counts_df_cocktail_summary['interaction_count'] > 0).astype(int)

    reorder_str = f'reorder({strain_col}, fraction_identified)'
    facet_str = f'{grouping_col}~.'

    # Plotting
    if generate_plot:
        if grouping_col:
            host_plots = (
                ggplot(counts_df_cocktail, aes(y='fraction_identified', x=reorder_str, fill='interaction_count')) +
                geom_col() +
                theme(axis_text_x=element_text(rotation=90),
                    figure_size=(10, 10)) +
                labs(title=plot_title, y="Possible Fraction", x=strain_col) +
                scale_fill_gradient(low="white", high="#4e6156", limits=(0, max_interactions)) +
                lims(y=(0, 1)) +
                facet_grid(facet_str)
            )
            print(host_plots)

        else:
            host_plots = (
                    ggplot(counts_df_cocktail, aes(y = "fraction_identified", x = "reorder(strain, fraction_identified)", fill = 'interaction_count')) +
                    geom_col() +
                    theme(axis_text_x=element_text(rotation=90),
                        # figure_size=(10, 2.5)) +
                        figure_size=(8, 2)) +
                    labs(title=plot_title, y = "Possible Fraction", x = "Strain") +
                    scale_fill_gradient(low="white", high="#4e6156", limits=(0,max_interactions)) +
                    lims(y=(0,1))
                )
            print(host_plots)

        # Save plot if specified
        if save_plot:
            host_plots.save(plot_filename, dpi=300)
    
    else:
        host_plots = None

    # Return the results
    return counts_df_cocktail, counts_df_cocktail_summary, host_plots


def cluster_and_select_strains(strain_features_df, n_clusters=10, random_state=42):
    """
    Clusters strains based on feature content and selects one representative strain from each cluster.
    
    Args:
        strain_features_df (DataFrame): DataFrame containing strain features
        n_clusters (int): Number of clusters to create
        random_state (int): Random seed for reproducibility
        
    Returns:
        tuple: (selected_strains, cluster_assignments)
            - selected_strains: list of selected strain IDs
            - cluster_assignments: DataFrame with strain IDs and their cluster assignments
    """
    # Get feature columns (assuming first column is strain ID)
    feature_columns = strain_features_df.columns[1:]
    
    # Calculate pairwise distances
    strain_distances = pdist(strain_features_df[feature_columns], metric='euclidean')
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(strain_distances, method='ward')
    
    # Assign clusters
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    # Create DataFrame with cluster assignments
    cluster_assignments = pd.DataFrame({
        'strain': strain_features_df.iloc[:, 0],
        'cluster': cluster_labels
    })
    
    # Randomly select one strain from each cluster
    np.random.seed(random_state)
    selected_strains = []
    for cluster in range(1, n_clusters + 1):
        cluster_strains = cluster_assignments[cluster_assignments['cluster'] == cluster]['strain'].values
        if len(cluster_strains) > 0:
            selected_strain = np.random.choice(cluster_strains, 1)[0]
            selected_strains.append(selected_strain)
    
    return selected_strains, cluster_assignments

def evalutate_cocktail_performance(
    modeling_dir,
    full_predictions_df,
    output_path):


    bootstrapping_dir = os.path.join(modeling_dir, 'bootstrapping')

    bootstrap_metrics_df = pd.DataFrame()

    if not full_predictions_df.empty:
        for (iteration, run_label), group in full_predictions_df.groupby(['iteration', 'modeling_method']):
            if run_label != 'inverse_frequency_hierarchical_filter':
                print('Skipping run_label:', run_label)
                continue
            # if iteration != 'iteration_1':
            #     continue
            print('Running cocktail design analysis for:', run_label)
            bootstrap_dir = os.path.join(bootstrapping_dir, run_label)
            performance_path = os.path.join(bootstrap_dir, iteration, 'modeling_results', 'model_performance', 'model_performance_metrics.csv')
            if not os.path.exists(performance_path):
                print('No performance file found for:', iteration, ' - ', run_label)
                continue

            performance_df = pd.read_csv(performance_path)
            # display(performance_df.head())
            max_cutoff = performance_df['cut_off'].iloc[0]
            feature_table_path = os.path.join(bootstrap_dir, iteration, 'feature_selection', 'filtered_feature_tables', f'select_feature_table_{max_cutoff}.csv')

            if not os.path.exists(feature_table_path):
                print('No feature table found for:', iteration, ' - ', run_label)
                continue
            feature_table_df = pd.read_csv(feature_table_path)

            # display(feature_table_df.head())

            # for n in range(1,6):
            for n in [1, 3, 5]:
                print('Running cocktail design analysis for:', iteration)
                cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=False, clustering_method='hierarchical', display_intermediate=False, save_plot=False, top_n=n, generate_plot=False, plot_umap=False)
                cocktail_design_summary['clustering'] = 'None'
                cocktail_design_summary['top_n'] = n
                cocktail_design_summary['run_label'] = run_label
                bootstrap_metrics_df = pd.concat([bootstrap_metrics_df, cocktail_design_summary])

                cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=True, clustering_method='hierarchical', display_intermediate=False, save_plot=False, top_n=n, generate_plot=False, plot_umap=False)
                cocktail_design_summary['clustering'] = 'hierarchical'
                cocktail_design_summary['top_n'] = n
                cocktail_design_summary['run_label'] = run_label
                bootstrap_metrics_df = pd.concat([bootstrap_metrics_df, cocktail_design_summary])

                # if run_label == 'bootstrap_log10_hierarchical_split':
                #     cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=True, clustering_method='HDBSCAN', display_intermediate=False, save_plot=False, top_n=n, generate_plot=True, plot_umap=False, split_by_strain=False)
                # else:
                cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=True, clustering_method='HDBSCAN', display_intermediate=False, save_plot=False, top_n=n, generate_plot=False, plot_umap=False, split_by_strain=False)
                cocktail_design_summary['clustering'] = 'HDBSCAN'
                cocktail_design_summary['top_n'] = n
                cocktail_design_summary['run_label'] = run_label
                bootstrap_metrics_df = pd.concat([bootstrap_metrics_df, cocktail_design_summary])

                
                if run_label == 'inverse_frequency_hierarchical_filter':   
                    # print('Running cocktail design analysis for:', iteration)
                    cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=False, clustering_method='hierarchical', display_intermediate=False, save_plot=False, top_n=n, generate_plot=False, plot_umap=False, default_cocktails=True)
                    cocktail_design_summary['clustering'] = 'None'
                    cocktail_design_summary['top_n'] = n
                    cocktail_design_summary['run_label'] = 'Default cocktail'
                    bootstrap_metrics_df = pd.concat([bootstrap_metrics_df, cocktail_design_summary])

                    cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=True, clustering_method='hierarchical', display_intermediate=False, save_plot=False, top_n=n, generate_plot=False, plot_umap=False, default_cocktails=True)
                    cocktail_design_summary['clustering'] = 'hierarchical'
                    cocktail_design_summary['top_n'] = n
                    cocktail_design_summary['run_label'] = 'Default cocktail'
                    bootstrap_metrics_df = pd.concat([bootstrap_metrics_df, cocktail_design_summary])

                    cocktail_design_df, cocktail_design_summary, cocktail_design_plot = cocktail_design(group, feature_table_df, clustering=True, clustering_method='HDBSCAN', display_intermediate=False, save_plot=False, top_n=n, generate_plot=False, plot_umap=False, split_by_strain=False, default_cocktails=True)
                    cocktail_design_summary['clustering'] = 'HDBSCAN'
                    cocktail_design_summary['top_n'] = n
                    cocktail_design_summary['run_label'] = 'Default cocktail'
                    bootstrap_metrics_df = pd.concat([bootstrap_metrics_df, cocktail_design_summary])

    else:
        print('No predictions parsed')

    bootstrap_metrics_df = bootstrap_metrics_df.groupby(['run_label', 'clustering', 'top_n']).agg({'True': 'sum', 'strain': 'count'}).reset_index()
    bootstrap_metrics_df['accuracy'] = bootstrap_metrics_df['True'] / bootstrap_metrics_df['strain']

    bootstrap_metrics_df.to_csv(output_path, index=False)

    display(bootstrap_metrics_df)

def main():
    parser = argparse.ArgumentParser(description='Evaluate cocktail performance based on phage modeling predictions.')
    parser.add_argument('--modeling_dir', type=str, required=True, help='Path to the modeling directory containing bootstrap predictions.')
    parser.add_argument('--full_predictions_path', type=str, required=True, help='Path to the full predictions CSV file containing phage modeling results.')
    parser.add_argument('--output_path', type=str, required=True, help='Path to save the output metrics CSV file.')
    args = parser.parse_args()

    modeling_dir = args.modeling_dir
    full_predictions_path = args.full_predictions_path
    output_path = args.output_path

    # Load full predictions DataFrame
    if os.path.exists(full_predictions_path):
        full_predictions_df = pd.read_csv(full_predictions_path)
    else:
        print(f"Full predictions file not found: {full_predictions_path}")
        return

    evalutate_cocktail_performance(
        modeling_dir=modeling_dir,
        full_predictions_df=full_predictions_df,
        output_path=output_path
    )

if __name__ == "__main__":
    main()
