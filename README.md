# Phage Modeling

Phage Modeling is a Python package designed to facilitate clustering, feature assignment, and machine learning-based analysis of phage genomes. This package leverages MMseqs2 for clustering and includes functions for processing biological sequence data, creating presence-absence matrices, and running machine learning workflows.

## Key Features

- **Clustering with MMseqs2**: Seamlessly run MMseqs2 to cluster genomic sequences and assign sequences to clusters. Supports genome-wide clustering for both phage and host genomes.
- **Feature Table Generation**: Create presence-absence matrices from clustered genomic sequences, with support for custom interaction matrices.
- **Feature Selection**: Perform recursive feature elimination (RFE) and other feature selection methods to optimize selected features based on clustering data.
- **Modeling Workflow**: Execute machine learning pipelines using CatBoost, grid search for hyperparameter tuning, and model evaluation with multiple metrics (e.g., F1-score, MCC).
- **Parallel Processing**: Support for multithreaded execution, allowing for efficient processing of large datasets.
- **Visualization**: Generate performance plots such as ROC curves, Precision-Recall curves, and hit rate/hit ratio curves to assess model performance.

## Installation

To install the package, clone the repository and run the following:

```bash
git clone https://github.com/Noonanav/phage_modeling.git
cd phage_modeling
pip install -e .
```

## External Dependencies
This package requires `MMseqs2` for clustering and sequence assignment. Install it via conda:

```bash
mamba install -c bioconda mmseqs2
```

## Usage

### Full Workflow

To run the entire workflow from feature table generation to modeling, you can use the `run-full-workflow` CLI command or call it directly in Python:

__CLI__:

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

### Clustering and Feature Table Generation

__CLI__:

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

Run the modeling workflow using selected features and evaluate model performance using metrics such as MCC and F1-Score.

CLI:

```bash
run-modeling-workflow --input_dir /path/to/filtered_feature_tables --output_dir /path/to/output_dir --threads 4 --num_runs 100
```

Python:

```python
from phage_modeling.workflows.modeling_workflow import run_modeling_workflow

run_modeling_workflow(
    input_dir="/path/to/filtered_feature_tables",
    base_output_dir="/path/to/output_dir",
    threads=4,
    num_runs=10
)
```


## License

This project is licensed under the MIT License. See the LICENSE file for more details.
