I've completed the PRS (Polygenic Risk Score) optimization module implementation. This script provides a comprehensive framework for optimizing PRS parameters based on GWAS results. Let me explain the key components:

PRS Optimization Module Overview
The gwas_step3_polygenic_risk_score.py script is designed to optimize parameters for constructing polygenic risk scores using Optuna. It builds on the prior QC and statistical model optimization steps in our GWAS pipeline.
Key Functions and Features:

Data Preparation

Loads and standardizes GWAS summary statistics
Handles various input formats from different GWAS tools
Prepares input files for PRS calculation


Parameter Optimization

Optimizes key PRS parameters:

P-value thresholds for SNP inclusion
LD clumping parameters (r² and distance)
Optional MAF filtering
Scoring methods


Uses AUC (for binary traits) or R² (for continuous traits) as optimization metrics


PRS Calculation

Performs LD clumping with PLINK
Creates scoring files for PRS calculation
Calculates PRS on target data
Evaluates performance against phenotype


Evaluation and Visualization

Generates comprehensive evaluation metrics
Creates visualizations of PRS performance:

ROC curves for binary traits
Regression plots for continuous traits
PRS distribution plots
Odds ratio plots by PRS quantile


Visualizes optimization process


Output Generation

Generates optimal PRS model
Creates application script for applying the model to new data
Produces detailed reports and visualizations



Workflow:

The script first loads GWAS summary statistics and validates inputs
It then runs an optimization loop using Optuna to find the best parameters
For each trial, it:

Performs LD clumping with the current parameters
Calculates PRS using the clumped SNPs
Evaluates performance against phenotype data


After optimization, it:

Generates the final PRS model with optimal parameters
Evaluates the model thoroughly
Creates visualizations
Generates an application script


```
Usage:
bashpython prs_optimization.py \
  --gwas-results my_gwas_results.txt \
  --target-data-prefix validation_data \
  --target-pheno phenotype.txt \
  --pheno-name DISEASE \
  --output-dir prs_results \
  --n-trials 100
```

This implementation provides a flexible and comprehensive tool for optimizing PRS parameters, which is a critical step in the GWAS pipeline. It complements the QC and statistical analysis optimization modules we've already implemented, creating a full end-to-end pipeline for GWAS analysis.

