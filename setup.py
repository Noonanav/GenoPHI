from setuptools import setup, find_packages

# Read the contents of the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='phage_modeling',
    version='0.1.0',
    description='A package for phage modeling, clustering, and feature assignment using MMseqs2',
    long_description=long_description,
    long_description_content_type="text/markdown",  # This ensures README.md is correctly processed as markdown
    author='Avery Noonan',
    packages=find_packages(),  # Automatically discovers all sub-packages
    install_requires=[
        'pandas', 
        'biopython', 
        'scikit-learn',  # For train_test_split and various metrics
        'catboost',  # For CatBoostClassifier
        'matplotlib',  # For plotting
        'seaborn',  # For advanced plotting
        'numpy',  # For numerical computations
        'tqdm',  # For progress bars
        'joblib',  # For saving and loading models
        'plotnine',  # For ggplot-style plotting
        'shap',  # For SHAP-based feature selection
    ],
    entry_points={
        'console_scripts': [
            'run-clustering-workflow=phage_modeling.workflows.feature_table_workflow:main',
            'run-feature-selection-workflow=phage_modeling.workflows.feature_selection_workflow:main',
            'run-modeling-workflow=phage_modeling.workflows.modeling_workflow:main',
            'run-full-workflow=phage_modeling.workflows.full_workflow:main',  # Full workflow
            'run-assign-features-workflow=phage_modeling.workflows.assign_features_workflow:main',  # Optional feature assignment workflow
            'run-prediction-workflow=phage_modeling.workflows.prediction_workflow:main',
            'run-assign-and-predict-workflow=phage_modeling.workflows.assign_predict_workflow:main',
            'run-modeling-from-feature-table=phage_modeling.workflows.select_and_model_workflow:main',
            'run-predictive-proteins-workflow=phage_modeling.workflows.feature_annotations_workflow:main',
            'run-kmer-table-workflow=phage_modeling.workflows.kmer_table_workflow:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.7',
)
