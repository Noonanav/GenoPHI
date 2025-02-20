import os
import pandas as pd
import argparse
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

def rerun_predictions(base_output_dir, threads=4):
    """
    Re-runs predictions for all completed iterations using existing feature assignments.
    
    Args:
        base_output_dir (str): Base directory containing iteration results
        threads (int): Number of threads to use for prediction
    """
    # Initialize DataFrame for final predictions
    final_predictions = pd.DataFrame()
    
    # Iterate through all directories
    for dir_name in os.listdir(base_output_dir):
        if not dir_name.startswith('iteration_'):
            continue
            
        iteration_dir = os.path.join(base_output_dir, dir_name)
        if not os.path.isdir(iteration_dir):
            continue
            
        iteration_num = int(dir_name.split('_')[1])
        print(f"Processing {dir_name}...")
        
        # Path to validation directory and existing assignment results
        validation_dir = os.path.join(iteration_dir, 'model_validation')
        assign_results_dir = os.path.join(validation_dir, 'assign_results')
        
        # Check if assignment results exist
        if not os.path.exists(assign_results_dir):
            print(f"No assignment results found for {dir_name}, skipping...")
            continue
            
        # Get the best cutoff from modeling results
        metrics_file = os.path.join(iteration_dir, "modeling_results/model_performance/model_performance_metrics.csv")
        if not os.path.exists(metrics_file):
            print(f"No metrics file found for {dir_name}, skipping...")
            continue
            
        # Get phage feature table path for this iteration
        phage_feature_table_path = os.path.join(iteration_dir, 'phage', 'features', 'feature_table.csv')
        if not os.path.exists(phage_feature_table_path):
            print(f"No phage feature table found for {dir_name}, skipping...")
            continue
            
        metrics_df = pd.read_csv(metrics_file)
        metrics_df = metrics_df.sort_values(['MCC', 'cut_off'], ascending=[False, False])
        best_cutoff = metrics_df['cut_off'].values[0]
        
        # Path to model directory with the best cutoff
        model_dir = os.path.join(iteration_dir, f'modeling_results', f'{best_cutoff}')
        
        # Create new prediction output directory with _new suffix
        new_predict_dir = os.path.join(validation_dir, 'predict_results_new')
        
        # Safety check - if directory exists and contains files, generate unique name
        if os.path.exists(new_predict_dir) and os.listdir(new_predict_dir):
            base_dir = new_predict_dir
            counter = 1
            while os.path.exists(new_predict_dir) and os.listdir(new_predict_dir):
                new_predict_dir = f"{base_dir}_{counter}"
                counter += 1
            print(f"Previous new predictions found, using directory: {new_predict_dir}")
            
        os.makedirs(new_predict_dir, exist_ok=True)
        
        # Run prediction workflow
        try:
            run_prediction_workflow(
                input_dir=assign_results_dir,
                phage_feature_table_path=phage_feature_table_path,
                model_dir=model_dir,
                output_dir=new_predict_dir,
                threads=threads
            )
            
            # Load and append new predictions
            new_predictions = pd.read_csv(os.path.join(new_predict_dir, 'strain_median_predictions.csv'))
            new_predictions['iteration'] = iteration_num
            final_predictions = pd.concat([final_predictions, new_predictions], ignore_index=True)
            
            print(f"Successfully processed {dir_name}")
            
        except Exception as e:
            print(f"Error processing {dir_name}: {str(e)}")
            continue
    
    # Save final concatenated predictions
    if not final_predictions.empty:
        final_predictions.to_csv(os.path.join(base_output_dir, 'final_predictions_new.csv'), index=False)
        print("Completed all iterations. New final predictions saved.")
    else:
        print("No predictions were generated.")

def main():
    parser = argparse.ArgumentParser(description="Re-run predictions using existing feature assignments.")
    parser.add_argument('--base_output_dir', type=str, required=True, 
                      help="Base directory containing iteration results")
    parser.add_argument('--threads', type=int, default=4,
                      help="Number of threads to use")

    args = parser.parse_args()

    rerun_predictions(
        base_output_dir=args.base_output_dir,
        threads=args.threads
    )

if __name__ == "__main__":
    main()