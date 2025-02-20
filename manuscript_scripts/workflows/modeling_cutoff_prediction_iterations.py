import os
import logging
from argparse import ArgumentParser
from phage_modeling.workflows.prediction_workflow import run_prediction_workflow

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def validate_paths(iteration_dir):
    """Validate required paths exist."""
    required_paths = [
        os.path.join(iteration_dir, 'phage', 'features', 'feature_table.csv'),
        os.path.join(iteration_dir, 'modeling_results'),
        os.path.join(iteration_dir, 'model_validation')
    ]
    for path in required_paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required path not found: {path}")

def cutoff_predict_workflow(
    input_dir,
    threads=4
):
    """
    Runs prediction workflow for each iteration and cutoff.
    
    Args:
        input_dir (str): Directory containing iteration folders
        threads (int): Number of threads for processing
    """
    logging.info("Starting prediction workflow...")
    
    # Get and validate iterations
    iterations = sorted([i for i in os.listdir(input_dir) if 'iteration' in i])
    if not iterations:
        raise ValueError(f"No iteration directories found in {input_dir}")

    for iteration in iterations:
        logging.info(f"Processing iteration: {iteration}")
        iteration_dir = os.path.join(input_dir, iteration)
        
        try:
            validate_paths(iteration_dir)
            
            output_base_dir = os.path.join(iteration_dir, 'model_validation', 'model_testing')
            os.makedirs(output_base_dir, exist_ok=True)

            phage_feature_table_path = os.path.join(iteration_dir, 'phage', 'features', 'feature_table.csv')
            modeling_results_dir = os.path.join(iteration_dir, 'modeling_results')
            assigned_features_dir = os.path.join(iteration_dir, 'model_validation', 'assign_results')
            
            # Get and validate cutoffs
            modeling_cutoffs = sorted([i for i in os.listdir(modeling_results_dir) if 'cutoff' in i])
            if not modeling_cutoffs:
                logging.warning(f"No cutoff directories found in {modeling_results_dir}")
                continue

            for cutoff in modeling_cutoffs:
                logging.info(f"Processing cutoff: {cutoff}")
                cutoff_dir = os.path.join(modeling_results_dir, cutoff)
                predict_output_dir = os.path.join(output_base_dir, cutoff)

                try:
                    run_prediction_workflow(
                        input_dir=assigned_features_dir,
                        phage_feature_table_path=phage_feature_table_path,
                        model_dir=cutoff_dir,
                        output_dir=predict_output_dir,
                        threads=threads
                    )
                    logging.info(f"Successfully completed prediction for {cutoff}")
                except Exception as e:
                    logging.error(f"Error processing cutoff {cutoff}: {str(e)}")
                    continue

        except Exception as e:
            logging.error(f"Error processing iteration {iteration}: {str(e)}")
            continue

    logging.info("Prediction workflow completed.")

def main():
    parser = ArgumentParser(description="Run prediction workflow for multiple iterations and cutoffs.")
    parser.add_argument('--input_dir', type=str, required=True, 
                       help="Directory containing iteration folders")
    parser.add_argument('--threads', type=int, default=4, 
                       help="Number of threads for processing")

    args = parser.parse_args()
    
    setup_logging()
    
    try:
        cutoff_predict_workflow(
            input_dir=args.input_dir,
            threads=args.threads
        )
    except Exception as e:
        logging.error(f"Workflow failed: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()