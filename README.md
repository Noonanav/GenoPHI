# Phage Modeling

Phage Modeling is a Python package designed to facilitate the machine-learning (ML)-based prediction of phage-host interactions based on gene content of phage and host genomes from amino acid (AA) sequences of protein coding genes and an interaction matrix of experimentally validate phage-host interactions. `phage_modeling` include the following functionalities in its workflows:

1. __Clustering with MMSeqs2__: Cluster AA sequences of phages and hosts into gene families.
2. __Gene Family Presence-Absence Matrices__: Generate presence-absence matrices of gene families in phage and host genomes
3. __Feature Table Construction__: Build feature tables based on co-occurrence patterns of gene families. Gene families with identical occurrence patterns across all genomes are collapsed into a single feature.
4. __Feature Selection__: Identification of protein family features that are predictive of phage-host interactions in a given dataset. Available methods include: recursive feature elimination (RFE), ANOVA using `SelectKBest`, chi-squared test, lasso linear regression, and SHAP feature importances.
5. __CatBoost ML Pipelines__: Execute CatBoost-based pipelines for interaction prediction, including hyperparameter tuning, performance evaluation, and performance vizualization.
6. __Prediction of New Interactions__: Assign gene family features to new genomes and predict their interactions with existing datasets.

## Installation

### Virtual Environment (_optional_)

We recommend installing this package in a virtual environment. To generate a `conda` environment for this package, run the following:

```bash
conda create -n phage_modeling python=3
conda activate phage_modeling
```

To install the `phage_modeling` package, clone the repository and install by running the following:

```bash
git clone https://github.com/Noonanav/phage_modeling.git
cd phage_modeling
pip install -e .
```

## External Dependencies
This package requires `MMseqs2` for clustering and sequence assignment. You can install it via conda or mamba: 
__With conda__:

```bash
conda install -c bioconda mmseqs2
```

__With Mamba__:

```bash
mamba install -c bioconda mmseqs2
```

For other installation methods, see the [MMSeqs2 Wiki](https://github.com/soedinglab/MMseqs2/wiki#installation)

## Workflow Overview

### Protein Family Feature Construction

The MMseqs2-based feature table generation process involves clustering protein sequences and creating a presence-absence matrix that represents the presence or absence of clusters (features) across genomes. First, sequences are processed into an MMseqs2 database and clustered based on sequence identity and coverage thresholds. For new sequences, they are assigned to the nearest existing clusters using MMseqs2 search. A presence-absence matrix is then generated to indicate the occurrence of clusters in each genome. This matrix is further used for feature selection and assignment, where unique clusters are identified for downstream analysis. The resulting feature tables can be merged with phenotype data for machine learning applications.

### Feature Selection

The feature selection workflow involves identifying the most predictive features from the dataset using various selection methods. The process allows users to select from Recursive Feature Elimination (RFE), SelectKBest (ANOVA F-test), Chi-Squared, Lasso, or SHAP-based methods. Data is split into training and testing sets, and feature selection is performed iteratively to identify top predictors of the target variable. Selected features are then used to train models, and their importance is ranked. After multiple runs, occurrence counts of selected features are compiled, and feature tables are generated for downstream analysis, providing a refined set of features for model training and evaluation.

### Modeling Workflow

The select feature modeling process involves training CatBoost models on filtered feature tables and assessing their performance. Multiple iterations are run for each feature table, where the data is split into training and testing sets. CatBoost models are trained using grid search to identify optimal hyperparameters. SHAP values are calculated for each model to determine feature importances. Performance metrics such as MCC, AUC, accuracy, precision, recall, and F1 are computed. Results are aggregated across runs and cutoffs, and SHAP summary plots visualize feature importances. Finally, performance metrics and SHAP values are saved for further analysis.

## Performance Metrics

- **AUC (Area Under the ROC Curve)**: The AUC represents the probability that a classifier will rank a randomly chosen positive instance higher than a randomly chosen negative one. It is calculated from the ROC curve as the area under the curve.

    $\text{AUC} = \int_{0}^{1} TPR(FPR) \, dFPR$

    where TPR is the True Positive Rate and FPR is the False Positive Rate.

- **Accuracy**: The proportion of true results (both true positives and true negatives) among the total number of cases.

    $\text{Accuracy} = \frac{TP + TN}{TP + TN + FP + FN}$

    where TP, TN, FP, and FN represent True Positives, True Negatives, False Positives, and False Negatives, respectively.

- **Precision**: The proportion of true positives among all instances that were predicted as positive.

    $\text{Precision} = \frac{TP}{TP + FP}$

- **Recall (Sensitivity or True Positive Rate)**: The proportion of true positives among all actual positive instances.

    $\text{Recall} = \frac{TP}{TP + FN}$

- **F1 Score**: The harmonic mean of precision and recall. It balances the two metrics, especially useful when the class distribution is imbalanced.

    $F1 = 2 \times \frac{\text{Precision} \times \text{Recall}}{\text{Precision} + \text{Recall}}$

- **MCC (Matthews Correlation Coefficient)**: A measure of the quality of binary classifications. It takes into account true and false positives and negatives, and is generally regarded as a balanced metric.

    $\text{MCC} = \frac{TP \times TN - FP \times FN}{\sqrt{(TP + FP)(TP + FN)(TN + FP)(TN + FN)}}$

#### Performance Plots

The following performance plots can be found in each `run_*` directory in the `<output_directory>/modeling_results/cutoff_*` directories for each modeling run:

 - Confusion matrices
 - Precision-Recall Curves
 - ROC AUC Curves
 - SHAP feature importance bar plots
 - SHAP value jitter plots

The following performance summary plots can also be found in `<output_directory>/modeling_results/model_performance/` directory, with an overview of model performance across feature selection cutoffs:

- SHAP value jitter plots showing SHAP values for all modeling runs:

![SHAP Summary](images/cutoff_5_shap_beeswarm.png)

- ROC AUC, precision-eecall, hit-rate, hit-ratio curves comparing feature selection cutoffs:

![ROC AUC](images/roc_curve.png)

![PR](images/pr_curve.png)

![Hit Rate](images/hit_rate_curve.png)

![Hit Ratio](images/hit_ratio_curve.png)

## Usage

### Full Workflow

#### Workflow Overview:

1. **Clustering of Genomes**: 
   The first step involves clustering the gene sequences of both host (strain) and phage genomes using MMseqs2. The clustering results in gene families, and each genome's genes are assigned to clusters. This step generates the cluster assignments for both the host and phage genomes.

2. **Feature Table Generation**: 
   Feature tables are created for the strain and phage genomes based on the gene families they are assigned to. These tables represent the presence or absence of specific gene families across all genomes.

3. **Merging Feature Tables**: 
   The strain and phage feature tables are merged using the provided interaction matrix, which contains validated interactions between strains and phages.

4. **Feature Selection**: 
   Feature selection methods, including Recursive Feature Elimination (RFE), ANOVA with `SelectKBest`, chi-squared test, lasso regression, and SHAP values, are applied to the merged feature table to identify the most predictive gene families (features) for host-phage interactions. These methods are run across multiple iterations, ensuring robust identification of relevant features. Results from each feature selection method are saved for easy comparison and further analysis.

5. **Modeling**: 
   Using the selected features, machine learning models (CatBoost) are trained to predict host-phage interactions. The model is evaluated using performance metrics like ROC curves and Precision-Recall curves, and the trained models are saved for future predictions.

6. **Interaction Prediction**: 
   Once the model is trained, it is used to predict interactions between new genomes and phages. The predictions are outputted with confidence scores for each potential interaction.

To run the entire workflow from feature table generation to modeling, you can use the `run-full-workflow` CLI command or call it directly in Python

__CLI__:

Help message:

```bash
usage: run-full-workflow [-h] -ih INPUT_STRAIN -ip INPUT_PHAGE -im INTERACTION_MATRIX [--suffix SUFFIX] [--strain_list STRAIN_LIST] [--phage_list PHAGE_LIST] [--strain_column STRAIN_COLUMN] [--phage_column PHAGE_COLUMN]
                         [--source_strain SOURCE_STRAIN] [--source_phage SOURCE_PHAGE] [--sample_column SAMPLE_COLUMN] [--phenotype_column PHENOTYPE_COLUMN] -o OUTPUT [--tmp TMP] [--min_seq_id MIN_SEQ_ID] [--coverage COVERAGE]
                         [--sensitivity SENSITIVITY] [--compare] [--filter_type FILTER_TYPE] [--method {rfe,select_k_best,chi_squared,lasso,shap}] [--num_features NUM_FEATURES] [--num_runs_fs NUM_RUNS_FS]
                         [--num_runs_modeling NUM_RUNS_MODELING] [--threads THREADS]

Complete workflow: Feature table generation, feature selection, and modeling.

options:
  -h, --help            show this help message and exit

Input data:
  -ih INPUT_STRAIN, --input_strain INPUT_STRAIN
                        Path to the input directory or file for strain clustering.
  -ip INPUT_PHAGE, --input_phage INPUT_PHAGE
                        Path to the input directory or file for phage clustering.
  -im INTERACTION_MATRIX, --interaction_matrix INTERACTION_MATRIX
                        Path to the interaction matrix.

Optional input arguments:
  --suffix SUFFIX       Suffix for input FASTA files (default: faa).
  --strain_list STRAIN_LIST
                        Path to a strain list file for filtering (default: none).
  --phage_list PHAGE_LIST
                        Path to a phage list file for filtering (default: none).
  --strain_column STRAIN_COLUMN
                        Column in the strain list containing strain names (default: strain).
  --phage_column PHAGE_COLUMN
                        Column in the phage list containing phage names (default: phage).
  --source_strain SOURCE_STRAIN
                        Prefix for naming selected features for strain in the assignment step (default: strain).
  --source_phage SOURCE_PHAGE
                        Prefix for naming selected features for phage in the assignment step (default: phage).
  --sample_column SAMPLE_COLUMN
                        Column name for the sample identifier (optional).
  --phenotype_column PHENOTYPE_COLUMN
                        Column name for the phenotype (optional).

Output arguments:
  -o OUTPUT, --output OUTPUT
                        Output directory to save results.
  --tmp TMP             Temporary directory for intermediate files (default: tmp).

Clustering:
  --min_seq_id MIN_SEQ_ID
                        Minimum sequence identity for clustering (default: 0.6).
  --coverage COVERAGE   Minimum coverage for clustering (default: 0.8).
  --sensitivity SENSITIVITY
                        Sensitivity for clustering (default: 7.5).
  --compare             Compare original clusters with assigned clusters.

Feature selection and modeling:
  --filter_type FILTER_TYPE
                        Filter type for the input data ('none', 'strain', 'phage', 'dataset'; default: none).
  --method {rfe,select_k_best,chi_squared,lasso,shap}
                        Feature selection method ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap'; default: rfe).
  --num_features NUM_FEATURES
                        Number of features to select (default: 100).
  --num_runs_fs NUM_RUNS_FS
                        Number of feature selection iterations to run (default: 10).
  --num_runs_modeling NUM_RUNS_MODELING
                        Number of runs per feature table for modeling (default: 10).

General:
  --threads THREADS     Number of threads to use (default: 4).
```

Example usage:

```bash
run-full-workflow --input_host /path/to/host/fasta --input_phage /path/to/phage/fasta --interaction_matrix /path/to/interaction_matrix.csv --output /path/to/output_dir --threads 4 --num_features 100 --num_runs_fs 10 --num_runs_modeling 20
```

__Python__:

```python
from phage_modeling.workflows.full_workflow import run_full_workflow

run_full_workflow(
    input_path_strain="/path/to/host/fasta",
    input_path_phage="/path/to/phage/fasta",
    interaction_matrix="/path/to/interaction_matrix.csv",
    output_dir="/path/to/output_dir",
    threads=4,
    num_features=100,
    num_runs_fs=10,
    num_runs_modeling=20
)
```

#### Inputs:
- `input_host`: Directory or file containing FASTA sequences of host gene AA sequences.
- `input_phage`: Directory or file containing FASTA sequences of phage gene AA sequences.
- `interaction_matrix`: CSV file with validated phage-host interactions (default headers: ['strain', 'phage', 'interaction']).
- `strain_list (optional)`: CSV file with list of strains in the input directory to process. Values in the list must match strain filenames in `<strain>.faa` or `<strain>.<suffix>`. If not provided, all genomes in input directory will be used.
- `phage_list (optional)`: CSV file with list of phages in the input directory to process. Values in the list must match phage filenames in `<phage>.faa` or `<phage>.<suffix>`. If not provided, all genomes in input directory will be used.
- `output_dir`: Directory to save the generated outputs.
- `threads`: Number of threads to use for processing.
- `tmp_dir`: Directory for temporary files during the workflow.
- `num_features`: Number of features to select during feature selection.
- `num_runs_fs`: Number of runs for feature selection to identify optimal features.
- `num_runs_modeling`: Number of runs for the modeling step, each with different random training/testing splits.
- `method`: Feature selection method (options: 'rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap')

#### Outputs:
- `presence_absence_matrix.csv`: Feature matrix indicating gene family presence/absence.
- `strain_feature_table.csv`: Feature table with collapsed features representing strain protein families.
- `phage_feature_table.csv`: Feature table with collapsed features representing phage protein families.
- `best_model.pkl`: Trained model after feature selection and modeling (`num_runs_modeling` models will be generated from random training/testing splits).
- `all_predictions.csv`: Predictions of host-phage interactions, with confidence scores, from all models.
- `mean_predictions.csv`: Aggregated predictions of phage-host interactions across all models, showing mean prediction scores.
- `roc_curve.png`: Receiver Operating Characteristic (ROC) curve of model performance.
- `precision_recall_curve.png`: Precision-Recall curve of model performance.

__Output Directory Structure__:

```bash
<output_directory>/
├── feature_selection/
│   ├── features_occurrence.csv          # Occurrence of features across genomes
│   ├── filtered_feature_tables/         # Contains filtered feature tables after feature selection
│   ├── run_X/                           # Individual feature selection runs (X indicates run number)
│   └── catboost_info/                   # Directory with CatBoost model logs and details
├── merged/
│   └── full_feature_table.csv           # Combined feature table for phage and strain genomes
├── modeling_results/
│   ├── cutoff_X/                        # Modeling results at different feature selection cutoffs (X indicates cutoff value)
│   ├── model_performance/               # Contains performance evaluation files of the models (figures and performance metrics)
│   ├── select_features_model_performance.csv # Overall performance metrics of all models generated with selected features
│   └── select_features_model_predictions.csv # Predictions made by the final model on the test set
├── phage/
│   ├── assigned_clusters.tsv            # Cluster assignments for phage genomes
│   ├── best_hits.tsv                    # Best cluster hits for phage sequences
│   ├── clusters.*                       # MMseqs2 clustering files for phage genomes
│   ├── clusters.tsv                     # Phage cluster results in TSV format
│   ├── features/                        # Directory containing feature-related files for phages
│   └── presence_absence_matrix.csv      # Presence-absence matrix for phage genomes
├── strain/
│   ├── assigned_clusters.tsv            # Cluster assignments for strain genomes
│   ├── best_hits.tsv                    # Best cluster hits for strain sequences
│   ├── clusters.*                       # MMseqs2 clustering files for strain genomes
│   ├── clusters.tsv                     # Strain cluster results in TSV format
│   ├── features/                        # Directory containing feature-related files for strains
│   └── presence_absence_matrix.csv      # Presence-absence matrix for strain genomes
├── tmp/
│   ├── phage/                           # Temporary files used during processing of phage genomes
│   └── strain/                          # Temporary files used during processing of strain genomes
```

### Assign and Predict Workflow

The `assign_predict` workflow involves assigning new sequences to pre-existing clusters, generating feature tables for new genomes, and making predictions of phage-host interactions using pre-trained models.

#### Workflow Overview:

1. __Assign Sequences to Clusters__: New genome sequences are assigned to existing MMseqs2 clusters based on similarity to previously clustered sequences. This step generates the cluster assignments for the new genomes.

2. __Feature Table Generation__: A feature table is created for each new genome based on the clusters they were assigned to. This table is used as input for interaction prediction models.

3. __Interaction Prediction__: Pre-trained models are used to predict interactions between new genomes (e.g., strains) and phages. The predictions are saved with confidence scores for each interaction.

__CLI__:

Help message:

```bash
usage: run-assign-and-predict-workflow [-h] --input_dir INPUT_DIR [--genome_list GENOME_LIST] [--genome_type {strain,phage}] [--genome_column GENOME_COLUMN] --mmseqs_db MMSEQS_DB --clusters_tsv CLUSTERS_TSV --feature_map FEATURE_MAP --tmp_dir
                                       TMP_DIR [--suffix SUFFIX] --model_dir MODEL_DIR --phage_feature_table PHAGE_FEATURE_TABLE --output_dir OUTPUT_DIR [--sensitivity SENSITIVITY] [--coverage COVERAGE] [--min_seq_id MIN_SEQ_ID]
                                       [--threads THREADS]

Assign new genes to existing clusters and predict interactions.

options:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR
                        Directory containing new genome FASTA files.
  --genome_list GENOME_LIST
                        CSV file with genome names.
  --genome_type {strain,phage}
                        Type of genome to process.
  --genome_column GENOME_COLUMN
                        Column name for genome identifiers in genome_list.
  --mmseqs_db MMSEQS_DB
                        Path to the existing MMseqs2 database.
  --clusters_tsv CLUSTERS_TSV
                        Path to the clusters TSV file.
  --feature_map FEATURE_MAP
                        Path to the feature map (selected_features.csv).
  --tmp_dir TMP_DIR     Temporary directory for intermediate files.
  --suffix SUFFIX       Suffix for FASTA files.
  --model_dir MODEL_DIR
                        Directory with trained models.
  --phage_feature_table PHAGE_FEATURE_TABLE
                        Path to the phage feature table.
  --output_dir OUTPUT_DIR
                        Directory to save results.
  --sensitivity SENSITIVITY
                        Sensitivity for MMseqs2 search.
  --coverage COVERAGE   Minimum coverage for assignment.
  --min_seq_id MIN_SEQ_ID
                        Minimum sequence identity for assignment.
  --threads THREADS     Number of threads for MMseqs2.
```

Example usage:

```bash
run-assign-and-predict-workflow --input_dir /path/to/genome/fasta --genome_list /path/to/genome_list.csv --genome_type strain --mmseqs_db /path/to/mmseqs_db --clusters_tsv /path/to/clusters.tsv --feature_map /path/to/selected_features.csv --output_dir /path/to/output_dir --tmp_dir /path/to/tmp_dir --model_dir /path/to/model_dir --phage_feature_table /path/to/phage_feature_table.csv --threads 4
```

__Python__:

```python
from phage_modeling.workflows.assign_predict_workflow import run_assign_and_predict_workflow

run_assign_and_predict_workflow(
    input_dir="/path/to/genome/fasta",
    genome_list="/path/to/genome_list.csv",
    genome_type="strain",  # or "phage"
    genome_column=None,  # If your genome list file has a specific column for genomes
    mmseqs_db="/path/to/mmseqs_db",
    clusters_tsv="/path/to/clusters.tsv",
    feature_map="/path/to/selected_features.csv",
    tmp_dir="/path/to/tmp_dir",
    model_dir="/path/to/model_dir",
    phage_feature_table="/path/to/phage_feature_table.csv",
    output_dir="/path/to/output_dir",
    threads=4
)
```

#### Inputs:
- `input_dir`: Directory containing FASTA files for new genomes.
- `genome_list`: CSV file with the list of genomes to process (e.g., phages or strains). 
- `genome_type`: Type of genome to process (`strain` or `phage`).
- `genome_column`: Optional column name from `genome_list` to use for genome identification.
- `mmseqs_db`: Path to the existing MMseqs2 database used for assigning clusters.
- `clusters_tsv`: TSV file of cluster assignments.
- `feature_map`: CSV file with selected features for the genome.
- `tmp_dir`: Directory for temporary files during the workflow.
- `model_dir`: Directory containing pre-trained models for interaction prediction.
- `phage_feature_table`: CSV file with feature table of phage genomes.
- `output_dir`: Directory to save the generated feature tables and predictions.
- `threads`: Number of threads to use for processing.

#### Outputs:
- `*_feature_table.csv`: Feature tables for each genome, indicating the presence/absence of gene families based on clustering.
- `all_predictions.csv`: Predictions of interactions between phages and hosts with confidence scores.
- `mean_predictions.csv`: Aggregated predictions showing the mean prediction score for each phage-host pair.

### Clustering and Feature Table Generation

The clustering and feature table generation step involves clustering protein sequences from phage and host genomes using MMseqs2. It generates a presence-absence matrix representing the gene families identified from both phages and hosts and prepares the feature tables for use in machine learning models.

__CLI__:

Help message:

```bash
usage: run-clustering-workflow [-h] -ih INPUT_STRAIN -ip INPUT_PHAGE -im INTERACTION_MATRIX -o OUTPUT [--tmp TMP] [--min_seq_id MIN_SEQ_ID] [--coverage COVERAGE] [--sensitivity SENSITIVITY] [--suffix SUFFIX] [--threads THREADS]
                               [--strain_list STRAIN_LIST] [--strain_column STRAIN_COLUMN] [--phage_list PHAGE_LIST] [--phage_column PHAGE_COLUMN] [--compare] [--source_strain SOURCE_STRAIN] [--source_phage SOURCE_PHAGE]

Run full feature table generation and merging workflow.

options:
  -h, --help            show this help message and exit
  -ih INPUT_STRAIN, --input_strain INPUT_STRAIN
                        Input path for strain clustering (directory or file).
  -ip INPUT_PHAGE, --input_phage INPUT_PHAGE
                        Input path for phage clustering (directory or file).
  -im INTERACTION_MATRIX, --interaction_matrix INTERACTION_MATRIX
                        Path to the interaction matrix.
  -o OUTPUT, --output OUTPUT
                        Output directory to save results.
  --tmp TMP             Temporary directory for intermediate files.
  --min_seq_id MIN_SEQ_ID
                        Minimum sequence identity for clustering.
  --coverage COVERAGE   Minimum coverage for clustering.
  --sensitivity SENSITIVITY
                        Sensitivity for clustering.
  --suffix SUFFIX       Suffix for input FASTA files.
  --threads THREADS     Number of threads to use.
  --strain_list STRAIN_LIST
                        Path to a strain list file for filtering.
  --strain_column STRAIN_COLUMN
                        Column in the strain list containing strain names.
  --phage_list PHAGE_LIST
                        Path to a phage list file for filtering.
  --phage_column PHAGE_COLUMN
                        Column in the phage list containing phage names.
  --compare             Compare original clusters with assigned clusters.
  --source_strain SOURCE_STRAIN
                        Prefix for naming selected features for strain in the assignment step.
  --source_phage SOURCE_PHAGE
                        Prefix for naming selected features for phage in the assignment step.
```

Example usage:

```bash
run-clustering-workflow --input_host /path/to/host/fasta --input_phage /path/to/phage/fasta --interaction_matrix /path/to/interaction_matrix.csv --output /path/to/output_dir --threads 4
```

__Python__:

```python
from phage_modeling.workflows.feature_table_workflow import run_full_feature_workflow

run_full_feature_workflow(
    input_path_strain="/path/to/host/fasta",
    input_path_phage="/path/to/phage/fasta",
    interaction_matrix="/path/to/interaction_matrix.csv",
    output_dir="/path/to/output_dir",
    threads=4
)
```

### Feature Selection

Perform feature selection on your generated feature table using recursive feature elimination (RFE) with the following command:

__CLI__:

Help message:

```bash
usage: run-feature-selection-workflow [-h] -i INPUT -o OUTPUT [--threads THREADS] [--num_features NUM_FEATURES] [--filter_type FILTER_TYPE] [--num_runs NUM_RUNS] [--method {rfe,select_k_best,chi_squared,lasso,shap}]

Run feature selection workflow.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input path for the full feature table.
  -o OUTPUT, --output OUTPUT
                        Base output directory for the results.
  --threads THREADS     Number of threads to use.
  --num_features NUM_FEATURES
                        Number of features to select during feature selection.
  --filter_type FILTER_TYPE
                        Type of filtering to use ('none', 'strain', 'phage').
  --num_runs NUM_RUNS   Number of feature selection iterations to run.
  --method {rfe,select_k_best,chi_squared,lasso,shap}
                        Feature selection method to use ('rfe', 'select_k_best', 'chi_squared', 'lasso', 'shap').
```

Usage example:

```bash
run-feature-selection-workflow --input /path/to/merged_feature_table.csv --output /path/to/output_dir --threads 4 --num_features 100 --num_runs_fs 50
```

__Python__:

```python
from phage_modeling.workflows.feature_selection_workflow import run_feature_selection_workflow

run_feature_selection_workflow(
    input_path="/path/to/merged_feature_table.csv",
    base_output_dir="/path/to/output_dir",
    threads=4,
    num_features=100,
    num_runs_fs=10
)
```

### Modeling Workflow

Run the modeling workflow using selected features and evaluate model performance using metrics such as MCC, AUC, F1-Score, Precision and Recall. Generates performance plots to evaluate and compare feature selection thresholds.

__CLI__:

Help message:

```bash
usage: run-modeling-workflow [-h] -i INPUT_DIR -o OUTPUT_DIR [--threads THREADS] [--num_runs NUM_RUNS] [--set_filter SET_FILTER] [--sample_column SAMPLE_COLUMN] [--phenotype_column PHENOTYPE_COLUMN]

Run modeling workflow on selected feature tables.

options:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Directory containing selected feature tables.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to save results of the experiments.
  --threads THREADS     Number of threads to use.
  --num_runs NUM_RUNS   Number of runs per feature table.
  --set_filter SET_FILTER
                        Filter for dataset ('none', 'strain', 'phage', 'dataset').
  --sample_column SAMPLE_COLUMN
                        Column name for the sample identifier (optional).
  --phenotype_column PHENOTYPE_COLUMN
                        Column name for the phenotype (optional).
```

Example usage:

```bash
run-modeling-workflow --input_dir /path/to/filtered_feature_tables --output_dir /path/to/output_dir --threads 4 --num_runs 100
```

__Python__:

```python
from phage_modeling.workflows.modeling_workflow import run_modeling_workflow

run_modeling_workflow(
    input_dir="/path/to/filtered_feature_tables",
    base_output_dir="/path/to/output_dir",
    threads=4,
    num_runs=10
)
```

## License and Copyright

This software is available under the MIT License. See the [LICENSE](LICENSE) file for details.

This software is also subject to Lawrence Berkeley National Laboratory copyright. 
See the [COPYRIGHT](COPYRIGHT) file for details.

This software was developed under funding from the U.S. Department of Energy and 
the U.S. Government consequently retains certain rights.
