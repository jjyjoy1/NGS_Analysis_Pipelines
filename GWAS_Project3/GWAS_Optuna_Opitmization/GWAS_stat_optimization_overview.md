```
# GWAS Statistical Optimizer

This tool optimizes parameters for genome-wide association studies (GWAS) using Optuna, a hyperparameter optimization framework.

## Optimization Targets
The optimizer finds the best combination of:
- Analysis methods (basic association, logistic regression, or mixed linear model)
- Covariate handling and principal component selection
- Multiple testing correction methods
- Minor allele frequency (MAF) thresholds
- Missing data thresholds
- Method-specific parameters

## Key Features
- **Parameter Optimization**: Uses Optuna to find the optimal parameter combination
- **Cross-Validation**: Supports k-fold cross-validation for more robust results
- **Multiple Methods**: Supports PLINK and GCTA methods
- **Validation**: Can validate against known associations
- **Visualizations**: Generates Manhattan plots, QQ plots, power analysis, and more
- **Reporting**: Creates detailed summary reports and ready-to-use GWAS scripts

## Metrics Used for Optimization
The optimizer balances multiple evaluation metrics:
- **Genomic Inflation Factor**: Penalizes λGC that deviates from 1.0
- **Significant Hits**: Rewards higher numbers of significant associations
- **Power Score**: Evaluates power across MAF ranges
- **Recovery Rate**: Measures how well known associations are recovered (if provided)

## How to Use the Tool

### Basic Usage
```bash
python gwas_stat_optimizer.py \
  --plink-prefix /path/to/data/prefix \
  --output-dir results \
  --pheno-file phenotypes.txt \
  --pheno-name PHENOTYPE
```

### Advanced Options
```bash
python gwas_stat_optimizer.py \
  --plink-prefix /path/to/data/prefix \
  --output-dir results \
  --pheno-file phenotypes.txt \
  --pheno-name PHENOTYPE \
  --covar-file covariates.txt \
  --max-pcs 10 \
  --n-trials 200 \
  --n-jobs 4 \
  --cv-folds 5 \
  --known-associations known_snps.csv
```

## Requirements
The code requires the following dependencies:
- **Python**: 3.6+
- **GWAS Software**: PLINK (or PLINK2) and GCTA installed and available in PATH
- **Python Packages**: 
  - numpy
  - pandas
  - scipy
  - matplotlib
  - seaborn
  - optuna
  - sklearn

## Output
After running the optimization, the tool generates:
- **Optimization Database**: SQLite database of all trials
- **Power Analysis**: CSV file and plots showing power across MAF ranges
- **Visualizations**: Manhattan plot, QQ plot, parameter importance, etc.
- **Optimized Script**: Ready-to-use shell script with the best parameters
- **Summary Report**: Text report of optimization results

## Implementation Notes
- The code now includes all required methods and the main function
- Command-line parser is complete with reasonable defaults
- Error handling and logging are properly implemented
- Visualization functions generate useful plots for analysis

This tool helps researchers find optimal GWAS parameters to maximize statistical power while controlling for false positives and genomic inflation.
```

Key improvements:
1. Consistent header hierarchy
2. Proper code block formatting for commands
3. Clear section organization
4. Improved readability with bullet points
5. Maintained all technical details while making them more accessible
6. Added emphasis to important components
7. Consistent formatting for software names and metrics

