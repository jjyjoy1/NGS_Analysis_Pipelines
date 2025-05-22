# GWAS QC Parameter Optimization

This toolkit provides an automated way to optimize quality control (QC) parameters for Genome-Wide Association Studies (GWAS) using Optuna, a hyperparameter optimization framework.

## Overview

The QC stage of GWAS analysis involves filtering out low-quality samples and variants that could introduce bias or errors. The choice of thresholds and parameters can significantly impact the results, especially in case-control studies for complex diseases like cancer.

This toolkit optimizes key QC parameters to balance:
- Data quality
- Sample and variant retention
- Case-control balance
- Covariate distribution preservation
- Genomic inflation control

## Requirements

- Python 3.6+
- PLINK (1.9 or 2.0)
- GCTA
- GWAS data in PLINK binary format (.bed, .bim, .fam)
- Phenotype file
- Covariate files (optional)

## Quick Start

1. Place the scripts in your working directory:
   - `gwas_qc_optimizer.py` (main optimization script)
   - `run_gwas_optimizer.sh` (setup and execution script)

2. Make the setup script executable:
   ```bash
   chmod +x run_gwas_optimizer.sh
   ```

3. Run the setup script:
   ```bash
   ./run_gwas_optimizer.sh
   ```

4. Follow the prompts to configure the optimization

5. After optimization completes, find results in the `optuna_results/` directory:
   - `best_params.txt`: The optimal QC parameters
   - `optimized_qc_pipeline.sh`: Ready-to-use PLINK commands with optimized parameters
   - Various plots and metrics files

## Parameters Optimized

This toolkit optimizes the following QC parameters:

### Sample QC
- Sample call rate threshold (`--mind`)
- Heterozygosity outlier threshold (±SD from mean)
- IBD relatedness threshold (`--min`)
- LD pruning window size, step size, and r² threshold for relatedness pruning

### Variant QC
- Variant call rate threshold (`--geno`)
- Minor Allele Frequency (MAF) threshold (`--maf`)
- Hardy-Weinberg Equilibrium p-value threshold (`--hwe`)
- Differential missingness p-value threshold

## How It Works

1. **Subsampling**: Creates a smaller, representative subset of your data for faster optimization
2. **Iterative Trials**: Tests different parameter combinations using Optuna
3. **Evaluation Metrics**: Scores each trial using:
   - Sample and SNP retention rates
   - Case-control ratio preservation
   - Genomic inflation factor (λ)
   - Covariate distribution preservation
4. **Optimization**: Identifies the optimal parameter set that maximizes the overall score
5. **Result Generation**: Creates an optimized pipeline script with the best parameters

## Cancer-Specific Considerations

For cancer GWAS, this optimizer pays special attention to:
- Preserving case-control balance
- Maintaining clinical covariate distributions (tumor type, diagnosis age, etc.)
- Balancing demographic covariates (age, sex)
- Controlling genomic inflation for reliable association results

## Customization

You can modify the `gwas_qc_optimizer.py` script to:
- Adjust the scoring weights in the `objective` function
- Change the parameter search ranges
- Add additional cancer-specific metrics
- Integrate with downstream analysis tools

## Advanced Usage

### Running on High-Performance Computing

To run on an HPC cluster, modify the setup script to create a job submission script instead of directly running the optimizer.

### Expanding to Multi-Stage Optimization

After optimizing QC parameters, you can extend the approach to optimize:
- PCA components for population stratification
- GWAS model covariates
- Fine-mapping and credible set parameters


