from setuptools import setup, find_packages

setup(
    name='phage_modeling',
    version='0.1.0',
    description='A package for phage modeling, clustering, and feature assignment using MMseqs2',
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
    ],
    entry_points={
        'console_scripts': [
            # CLI entry points for different workflows
            'run-clustering-workflow=phage_modeling.workflows.feature_table_workflow:main',
            'run-feature-selection-workflow=phage_modeling.workflows.feature_selection_workflow:main',
            'run-modeling-workflow=phage_modeling.workflows.modeling_workflow:main',
            'run-full-workflow=phage_modeling.workflows.full_workflow:main',  # Full workflow
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',  # Adjust based on maturity
        'Programming Language :: Python :: 3.8',  # Adjust based on Python versions you support
        'License :: OSI Approved :: MIT License',  # Change if you're using a different license
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.7',
)
