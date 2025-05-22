# GWAS Statistical Analysis Parameter Optimization with Optuna

A comprehensive toolkit for optimizing GWAS (Genome-Wide Association Study) parameters using Optuna, designed to maximize downstream biological insights and predictive performance.

## Overview

This project provides a modular pipeline for optimizing GWAS parameters across multiple analysis stages:

1. **Quality Control (QC) Optimization** - Optimize data filtering parameters
2. **Statistical Model Optimization** - Optimize association testing parameters  
3. **Polygenic Risk Score (PRS) Optimization** - Optimize PRS construction parameters

The optimization framework uses Optuna to balance multiple objectives including genomic control, statistical power, discovery potential, and validation against known associations.

## Features

### 🎯 Multi-Objective Optimization
- **Genomic Control (λGC)**: Ensures statistical validity by penalizing inflation/deflation
- **Discovery Potential**: Maximizes significant associations while controlling false positives
- **Statistical Power**: Optimizes power across different Minor Allele Frequency (MAF) ranges
- **Known Association Recovery**: Validates parameters against established genetic loci
- **Predictive Performance**: Optimizes for downstream applications like PRS

### 🔧 Comprehensive Parameter Coverage
- **QC Parameters**: MAF thresholds, missing data rates, Hardy-Weinberg equilibrium
- **Association Methods**: Basic association, logistic regression, mixed linear models (GCTA-MLMA)
- **Covariate Handling**: Principal component optimization, covariate selection
- **Multiple Testing**: Bonferroni, FDR, and custom correction methods
- **PRS Parameters**: P-value thresholds, LD clumping, scoring methods

### 📊 Extensive Evaluation & Visualization
- Manhattan and QQ plots for association results
- ROC curves and AUC for binary traits
- Regression analysis and R² for continuous traits
- Optimization history and parameter importance plots
- Power analysis across MAF ranges

### 🔄 Flexible Implementation
- **Class-based version** for object-oriented programming
- **Procedural version** for functional programming preferences
- **Modular design** for easy customization and extension
- **Cross-validation support** for robust parameter selection

## Installation

### Prerequisites

- Python 3.7+
- PLINK 1.9+ or PLINK 2.0
- GCTA (optional, for mixed linear models)

### Required Python Packages

```bash
pip install optuna numpy pandas scikit-learn scipy matplotlib seaborn
```

### Optional Dependencies

```bash
# For enhanced visualizations
pip install plotly kaleido

# For statistical analysis
pip install statsmodels
```
## Configuration Options

### Core Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--plink-prefix` | Prefix for PLINK files (.bed/.bim/.fam) | Required |
| `--pheno-file` | Phenotype file in PLINK format | Optional |
| `--covar-file` | Covariate file (should include PCs) | Optional |
| `--output-dir` | Output directory | `gwas_optimizer_output` |
| `--n-trials` | Number of optimization trials | 50 |
| `--n-jobs` | Parallel jobs for optimization | 1 |

### Analysis Methods

| Parameter | Description |
|-----------|-------------|
| `--disable-assoc` | Disable basic association testing |
| `--disable-logistic` | Disable logistic regression |
| `--disable-mlma` | Disable mixed linear models |

### Cross-Validation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--cv-folds` | Number of cross-validation folds | 1 (no CV) |
| `--max-pcs` | Maximum principal components | 20 |


## Objective Function

The optimization balances multiple metrics:

```python
objective = (
    lambda_penalty * 3.0 +      # Penalize deviation from λGC = 1.0
    sig_reward * 1.0 +          # Reward significant discoveries
    power_reward * 2.0 +        # Reward statistical power
    recovery_reward * 3.0       # Reward known association recovery
)
```

### Weight Rationale

- **Lambda GC Penalty (3.0)**: High weight ensures statistical validity
- **Significant Hits Reward (1.0)**: Moderate weight prevents overfitting
- **Power Score Reward (2.0)**: Substantial weight for detection ability
- **Recovery Rate Reward (3.0)**: High weight when ground truth is available

## Advanced Usage

### Custom Objective Functions

You can implement custom objective functions for specific research goals:

```python
def custom_objective(metrics):
    """Custom objective for rare disease studies."""
    # Emphasize functional impact over statistical power
    functional_score = metrics.get("functional_enrichment", 0)
    effect_size_score = metrics.get("mean_effect_size", 0)
    
    return functional_score * 4.0 + effect_size_score * 2.0
```

### Multi-Stage Optimization

For computational efficiency, consider sequential optimization:

1. **Stage 1**: Optimize QC parameters (30 trials)
2. **Stage 2**: Optimize association parameters with fixed QC (50 trials)
3. **Stage 3**: Optimize PRS parameters with fixed GWAS (100 trials)


### Optimization Tips

1. **Start Small**: Use a subset of data for initial optimization
2. **Parallel Processing**: Use `--n-jobs` > 1 for faster optimization
3. **Early Stopping**: Use Optuna's pruning for efficiency
4. **Incremental Approach**: Optimize stages sequentially rather than jointly

### Related Publications

Key papers that inform our optimization approaches:

1. **Psychiatric Disorders**: Ge T, et al. "Polygenic risk scores for prediction of psychiatric disorders." *Am J Psychiatry* (2019). DOI: [10.1176/appi.ajp.2018.18091034](https://doi.org/10.1176/appi.ajp.2018.18091034)

2. **Autoimmune Diseases**: Sun C, et al. "High-density genotyping of immune-related loci identifies new SLE risk variants." *Nat Genet* (2016). DOI: [10.1038/ng.3496](https://doi.org/10.1038/ng.3496)

3. **Anthropometric Traits**: Yengo L, et al. "Meta-analysis of genome-wide association studies for height and body mass index." *Hum Mol Genet* (2018). DOI: [10.1093/hmg/ddy271](https://doi.org/10.1093/hmg/ddy271)



## Acknowledgments

- [Optuna](https://optuna.org/) team for the optimization framework
- [PLINK](https://www.cog-genomics.org/plink/) developers for GWAS analysis tools
- [GCTA](https://cnsgenomics.com/software/gcta/) team for mixed model implementations
- The broader GWAS community for methodological insights

---


