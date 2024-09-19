# Phage Modeling

Phage Modeling is a Python package designed to facilitate clustering, feature assignment, and machine learning-based analysis of phage genomes. This package leverages MMseqs2 for clustering and includes functions for processing biological sequence data, creating presence-absence matrices, and running machine learning workflows.

## Features

- **Clustering with MMseqs2**: Run clustering on genomic sequences and assign sequences to clusters.
- **Feature Selection**: Optimize and assign features for genomes based on presence-absence data.
- **Machine Learning Integration**: Includes support for CatBoost classifiers and evaluation using various metrics such as F1-score and Matthews Correlation Coefficient.

## Installation

To install the package, clone the repository and run the following:

```bash
git clone https://github.com/yourusername/phage_modeling.git
cd phage_modeling
pip install -e .
```

## Dependencies
This package requires `MMseqs2` for clustering and sequence assignment. Install it via conda:

```bash
conda install -c bioconda mmseqs2
```
