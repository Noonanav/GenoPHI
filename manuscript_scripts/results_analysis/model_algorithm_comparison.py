#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, accuracy_score, precision_score, recall_score,
    f1_score, matthews_corrcoef, roc_curve, precision_recall_curve
)
from itertools import combinations
from scipy.stats import mannwhitneyu
import argparse
from plotnine import ggplot, aes, geom_boxplot, geom_jitter, labs, theme, element_text
from IPython.display import display

def calculate_metrics(df):
    """Calculate performance metrics for model predictions."""
    y_true = df['interaction']
    y_pred_prob = df['Confidence']
    y_pred = df['Prediction']
    
    metrics = {
        'AUC': roc_auc_score(y_true, y_pred_prob),
        'Accuracy': accuracy_score(y_true, y_pred),
        'Precision': precision_score(y_true, y_pred),
        'Recall': recall_score(y_true, y_pred),
        'F1': f1_score(y_true, y_pred),
        'MCC': matthews_corrcoef(y_true, y_pred)
    }
    
    fpr, tpr, roc_thresholds = roc_curve(y_true, y_pred_prob)
    roc_df = pd.DataFrame({
        'fpr': fpr, 
        'tpr': tpr, 
        'thresholds': roc_thresholds, 
        'cut_off': df['cut_off'].iloc[0], 
        'model': df['model'].iloc[0]
    })
    
    if 'run' in df.columns:
        roc_df['run'] = df['run'].iloc[0]
    
    precision, recall, pr_thresholds = precision_recall_curve(y_true, y_pred_prob)
    pr_df = pd.DataFrame({
        'precision': precision, 
        'recall': recall, 
        'thresholds': list(pr_thresholds) + [1.0], 
        'cut_off': df['cut_off'].iloc[0], 
        'model': df['model'].iloc[0]
    })
    
    if 'run' in df.columns:
        pr_df['run'] = df['run'].iloc[0]
    
    return metrics, roc_df, pr_df

def perform_pairwise_mannwhitney(df, metric='MCC'):
    """Perform pairwise Mann-Whitney U tests between models."""
    methods = df['model'].unique()
    results = []
    
    # Get all unique pairs of methods
    for m1, m2 in combinations(methods, 2):
        group1 = df[df['model'] == m1][metric]
        group2 = df[df['model'] == m2][metric]
        
        # Perform Mann-Whitney U test
        try:
            u_stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
            
            results.append({
                'model1': m1,
                'model2': m2,
                'u_statistic': u_stat,
                'p_value': p_value,
                'significant': p_value < 0.05  # Using 0.05 as significance threshold
            })
        except Exception as e:
            print(f"Error in Mann-Whitney test for {m1} vs {m2}: {e}")
    
    return pd.DataFrame(results)

def process_dataset(path, output_dir, dataset_name, metric='MCC'):
    """Process a single dataset and save results."""
    if not os.path.exists(path):
        print(f'Model predictions file not found: {path}')
        return None
    
    print(f'Parsing: {path}')
    model_predictions_df_full = pd.read_csv(path)
    model_predictions_df_full['cut_off'] = model_predictions_df_full['cut_off'].astype(str)
    model_predictions_df_full['dataset'] = dataset_name
    
    # Create output directory for this dataset
    dataset_output_dir = os.path.join(output_dir, dataset_name.replace(" ", "_"))
    os.makedirs(dataset_output_dir, exist_ok=True)
    
    # Aggregate predictions
    print('Aggregating predictions')
    model_predictions_df_calcs = model_predictions_df_full.groupby(['model', 'cut_off', 'strain', 'phage', 'interaction']).agg({
        'Confidence': 'median'
    }).reset_index()
    
    model_predictions_df_calcs['Prediction'] = [1 if x > 0.5 else 0 for x in model_predictions_df_calcs['Confidence']]
    model_predictions_df_calcs['dataset'] = dataset_name
    
    # Calculate metrics to identify top cut_off for each model
    print('Identifying top cut_off based on MCC')
    metrics_list = []
    roc_list = []
    pr_list = []
    
    for (cut_off, model), group in model_predictions_df_calcs.groupby(['cut_off', 'model']):
        try:
            metrics, roc_df, pr_df = calculate_metrics(group)
            metrics['cut_off'] = cut_off
            metrics['model'] = model
            metrics['dataset'] = dataset_name
            metrics_list.append(metrics)
            roc_list.append(roc_df)
            pr_list.append(pr_df)
        except Exception as e:
            print(f"Error calculating metrics for cut_off={cut_off}, model={model}: {e}")
    
    if not metrics_list:
        print(f"No valid metrics calculated for {dataset_name}")
        return None
    
    # Convert results to dataframe
    metrics_df = pd.DataFrame(metrics_list)
    
    # Save cut_off selection metrics
    cutoff_metrics_path = os.path.join(dataset_output_dir, "cutoff_selection_metrics.csv")
    metrics_df.to_csv(cutoff_metrics_path, index=False)
    print(f"Saved cut_off selection metrics to: {cutoff_metrics_path}")
    
    # Identify best cut_off for each model based on MCC
    top_model_cutoff = metrics_df.copy()
    top_model_cutoff['max_MCC'] = top_model_cutoff.groupby('model')['MCC'].transform('max')
    top_model_cutoff = top_model_cutoff[top_model_cutoff['MCC'] == top_model_cutoff['max_MCC']]
    
    # Save best cut_offs
    best_cutoffs_path = os.path.join(dataset_output_dir, "best_cutoffs.csv")
    top_model_cutoff[['model', 'cut_off', 'MCC']].to_csv(best_cutoffs_path, index=False)
    print(f"Saved best cut_offs to: {best_cutoffs_path}")
    
    # Filter to top cut_off and include run information
    print('Filtering to top cut_off')
    try:
        model_predictions_df_calcs_runs = model_predictions_df_full.merge(
            top_model_cutoff[['model', 'cut_off']], 
            on=['model', 'cut_off'], 
            how='inner'
        )
        
        model_predictions_df_calcs_runs = model_predictions_df_calcs_runs.groupby(
            ['model', 'run', 'cut_off', 'strain', 'phage', 'interaction']
        ).agg({
            'Confidence': 'median'
        }).reset_index()
        
        model_predictions_df_calcs_runs['Prediction'] = [
            1 if x > 0.5 else 0 for x in model_predictions_df_calcs_runs['Confidence']
        ]
        model_predictions_df_calcs_runs['dataset'] = dataset_name
        
        # Calculate metrics by run
        print('Calculating metrics by run')
        run_metrics_list = []
        run_roc_list = []
        run_pr_list = []
        
        for (cut_off, model, run), group in model_predictions_df_calcs_runs.groupby(['cut_off', 'model', 'run']):
            try:
                metrics, roc_df, pr_df = calculate_metrics(group)
                metrics['cut_off'] = cut_off
                metrics['model'] = model
                metrics['run'] = run
                metrics['dataset'] = dataset_name
                run_metrics_list.append(metrics)
                run_roc_list.append(roc_df)
                run_pr_list.append(pr_df)
            except Exception as e:
                print(f"Error calculating metrics for cut_off={cut_off}, model={model}, run={run}: {e}")
        
        if run_metrics_list:
            # Convert results to dataframe
            run_metrics_df = pd.DataFrame(run_metrics_list)
            
            # Save run-level metrics
            run_metrics_path = os.path.join(dataset_output_dir, "run_level_metrics.csv")
            run_metrics_df.to_csv(run_metrics_path, index=False)
            print(f"Saved run-level metrics to: {run_metrics_path}")
            
            # Create plotting dataframe
            metrics_df_plotting = run_metrics_df.copy()
            metrics_df_plotting[f'median_{metric}'] = metrics_df_plotting.groupby(['cut_off', 'model'])[metric].transform('median')
            
            # Generate boxplot
            try:
                metric_plot = (
                    ggplot(metrics_df_plotting, aes(x=f'reorder(model,median_{metric})', y=metric)) +
                    geom_boxplot() +
                    geom_jitter() +
                    labs(x='Feature Selection Method', y=metric, title=f"{dataset_name} - {metric} by Model") +
                    theme(figure_size=(8, 6), axis_text_x=element_text(rotation=90))
                )
                
                # Save plot
                plot_path = os.path.join(dataset_output_dir, f"{dataset_name.replace(' ', '_')}_{metric.lower()}_boxplot.png")
                metric_plot.save(plot_path, dpi=300)
                print(f"Saved boxplot to: {plot_path}")
            except Exception as e:
                print(f"Error generating boxplot for {dataset_name}: {e}")
            
            # Perform statistical tests
            print(f'Performing Mann-Whitney tests for {metric}')
            try:
                mw_results = perform_pairwise_mannwhitney(run_metrics_df, metric)
                mw_results_sorted = mw_results.sort_values('p_value')
                
                # Save statistical test results
                stats_path = os.path.join(dataset_output_dir, f"mannwhitney_test_results_{metric.lower()}.csv")
                mw_results_sorted.to_csv(stats_path, index=False)
                print(f"Saved statistical test results to: {stats_path}")
            except Exception as e:
                print(f"Error performing Mann-Whitney tests for {dataset_name}: {e}")
            
            return run_metrics_df
        else:
            print(f"No valid run-level metrics calculated for {dataset_name}")
            return None
    except Exception as e:
        print(f"Error filtering to top cut_off for {dataset_name}: {e}")
        return metrics_df  # Return the cutoff-level metrics at least

def main():
    parser = argparse.ArgumentParser(description='Analyze model predictions across multiple datasets.')
    parser.add_argument('--output', '-o', required=True, help='Output directory for tables and figures')
    parser.add_argument('--datasets', '-d', nargs='+', help='Optional: Names for each dataset')
    parser.add_argument('--paths', '-p', nargs='+', help='Optional: Custom paths to prediction files')
    parser.add_argument('--metric', '-m', default='MCC', choices=['AUC', 'Accuracy', 'Precision', 'Recall', 'F1', 'MCC'],
                        help='Metric to use for statistical tests and visualizations (default: MCC)')
    
    args = parser.parse_args()
    output_dir = args.output
    metric = args.metric
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Default paths and dataset names
    default_paths = [
        '/usr2/people/anoonan/BRaVE/resources/genome_data/e_coli/ml_features/Gaborieau/modeling/model_testing/select_features_full_1s_model_predictions.csv',
        '/usr2/people/anoonan/BRaVE/resources/genome_data/klebsiella/beatriz/ml_features/modeling/model_testing/select_features/select_features_model_predictions.csv',
        '/usr2/people/anoonan/BRaVE/resources/genome_data/klebsiella/demi/ml_features/modeling/set1/model_testing/select_features_model_predictions.csv',
        '/usr2/people/anoonan/BRaVE/resources/genome_data/klebsiella/demi/ml_features/modeling/set2/model_testing/select_features_model_predictions.csv',
        '/usr2/people/anoonan/BRaVE/resources/genome_data/pseudomonas/anarita/ml_features/modeling/model_testing/select_features/select_features_model_predictions.csv'
    ]
    
    default_names = [
        'E_coli',
        'Klebsiella_Beatriz',
        'Klebsiella_Dimi1',
        'Klebsiella_Dimi2',
        'Pseudomonas'
    ]
    
    # Use provided paths or defaults
    prediction_paths = args.paths if args.paths else default_paths
    
    # Use provided dataset names or defaults (matching the number of paths)
    if args.datasets and len(args.datasets) == len(prediction_paths):
        dataset_names = args.datasets
    else:
        # If no datasets provided or count doesn't match, use defaults or generate names
        if len(prediction_paths) == len(default_names):
            dataset_names = default_names
        else:
            dataset_names = [f"Dataset_{i+1}" for i in range(len(prediction_paths))]
    
    print(f"Processing {len(prediction_paths)} prediction files")
    for path, name in zip(prediction_paths, dataset_names):
        print(f"Dataset: {name}, Path: {path}")
    print(f"Results will be saved to: {output_dir}")
    print(f"Using MCC for cut-off selection and {metric} for statistical tests and visualizations")
    
    # Process each dataset
    all_metrics = []
    for path, name in zip(prediction_paths, dataset_names):
        metrics_df = process_dataset(path, output_dir, name, metric)
        if metrics_df is not None:
            all_metrics.append(metrics_df)
    
    # Combine results across datasets if we have multiple datasets
    if len(all_metrics) > 1:
        combined_metrics = pd.concat(all_metrics, ignore_index=True)
        
        # Save combined metrics
        combined_path = os.path.join(output_dir, "combined_metrics.csv")
        combined_metrics.to_csv(combined_path, index=False)
        print(f"Saved combined metrics to: {combined_path}")
        
        # Create comparison plots
        try:
            # Boxplot comparing models across datasets
            plt.figure(figsize=(12, 8))
            sns.boxplot(data=combined_metrics, x='model', y=metric, hue='dataset')
            plt.title(f'{metric} by Model and Dataset')
            plt.xlabel('Model')
            plt.ylabel(metric)
            plt.xticks(rotation=90)
            plt.tight_layout()
            
            comparison_path = os.path.join(output_dir, f"dataset_comparison_boxplot_{metric.lower()}.png")
            plt.savefig(comparison_path, dpi=300)
            print(f"Saved dataset comparison boxplot to: {comparison_path}")
            plt.close()
            
            # Create heatmap of median metric values
            pivot_data = combined_metrics.groupby(['dataset', 'model'])[metric].median().reset_index()
            pivot_table = pivot_data.pivot(index='dataset', columns='model', values=metric)
            
            plt.figure(figsize=(12, 8))
            sns.heatmap(pivot_table, annot=True, cmap='viridis', fmt='.3f')
            plt.title(f'Median {metric} by Model and Dataset')
            plt.tight_layout()
            
            heatmap_path = os.path.join(output_dir, f"median_{metric.lower()}_heatmap.png")
            plt.savefig(heatmap_path, dpi=300)
            print(f"Saved median {metric} heatmap to: {heatmap_path}")
            plt.close()
            
            # Create line plot of median metric by dataset
            plt.figure(figsize=(12, 8))
            
            for dataset in combined_metrics['dataset'].unique():
                dataset_data = combined_metrics[combined_metrics['dataset'] == dataset]
                dataset_medians = dataset_data.groupby('model')[metric].median().reset_index()
                dataset_medians = dataset_medians.sort_values(metric)
                
                plt.plot(
                    dataset_medians['model'], 
                    dataset_medians[metric], 
                    marker='o', 
                    label=dataset
                )
            
            plt.xlabel('Model')
            plt.ylabel(f'Median {metric}')
            plt.xticks(rotation=90)
            plt.legend(title='Dataset')
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.tight_layout()
            
            line_plot_path = os.path.join(output_dir, f"dataset_model_comparison_line_{metric.lower()}.png")
            plt.savefig(line_plot_path, dpi=300)
            print(f"Saved dataset comparison line plot to: {line_plot_path}")
            plt.close()
            
        except Exception as e:
            print(f"Error generating comparison plots: {e}")

if __name__ == "__main__":
    main()

# python /usr2/people/anoonan/BRaVE/machine_learning/phage_modeling/manuscript_scripts/results_analysis/model_algorithm_comparison.py --output algorithm_comparison/AUC --metric AUC