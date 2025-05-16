# GWAS Analysis Pipeline Using GCTA

A comprehensive pipeline for performing Genome-Wide Association Studies (GWAS) with GCTA, including quality control, imputation, heritability estimation, association testing, and downstream analyses.

## Overview

This pipeline performs a complete GWAS analysis workflow using GCTA (Genome-wide Complex Trait Analysis) and PLINK, designed for case-control studies such as breast cancer genetics research. The pipeline implements best practices including:

- Rigorous quality control of samples and variants
- Imputation to increase marker density
- Genomic relatedness matrix construction
- Heritability estimation
- Mixed linear model association analysis with LOCO approach
- Conditional and joint analysis for independent signals
- Credible set analysis to identify potential causal variants
- Colocalization analysis with eQTL data
- Polygenic risk score development

## Requirements

- GCTA (v1.94.1 or later) 
- PLINK (v1.9 and v2.0)
- R (v4.0.0 or later) with the following packages:
  - data.table
  - ggplot2
  - coloc
- Python (v3.6 or later)
- Workflow manager (Snakemake or Nextflow)
- High-Performance Computing environment with SLURM

## Installation

```bash
# Clone this repository
git clone https://github.com/yourusername/gcta-gwas-pipeline.git
cd gcta-gwas-pipeline

# Download and install GCTA
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_1.94.1.zip
unzip gcta_1.94.1.zip
mv gcta_1.94.1/gcta64 /usr/local/bin/

# Install PLINK
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
unzip plink_linux_x86_64_20231211.zip
mv plink /usr/local/bin/

# Install PLINK2
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20231212.zip
unzip plink2_linux_x86_64_20231212.zip
mv plink2 /usr/local/bin/

# Install R packages
Rscript -e 'install.packages(c("data.table", "ggplot2", "coloc"), repos="https://cran.rstudio.com/")'
```

## Pipeline Workflow

### 1. Sample and Variant Quality Control

```bash
# Run sample QC
./scripts/01_sample_qc.sh --bfile your_raw_data --out qc_output

# Run variant QC
./scripts/02_variant_qc.sh --bfile qc_output_sample_qc --out qc_output_complete
```

This performs:
- Sample call rate filtering (--mind 0.05)
- Sex check and discrepancy removal
- Heterozygosity outlier removal
- Identity-by-descent analysis to remove related individuals
- Principal Component Analysis for population stratification
- Variant call rate filtering (--geno 0.05)
- Minor Allele Frequency filtering (--maf 0.01)
- Hardy-Weinberg Equilibrium testing (--hwe 1e-6)
- Differential missingness test between cases and controls

### 2. Imputation Preparation

```bash
# Prepare files for imputation
./scripts/03_imputation_prep.sh --bfile qc_output_complete --out imputation_ready
```

This splits data by chromosome and converts to VCF format for submission to an imputation server (e.g., Michigan Imputation Server or TOPMed).

### 3. Imputation Post-processing

```bash
# Process imputed data
./scripts/04_imputation_postproc.sh --vcfdir imputed_vcfs --out imputed_merged
```

This converts imputed VCFs back to PLINK format and merges chromosomes.

### 4. GCTA Analysis

```bash
# Create Genomic Relationship Matrix (GRM)
./scripts/05_make_grm.sh --bfile imputed_merged --out project_grm

# Estimate heritability
./scripts/06_heritability.sh --grm project_grm --pheno phenotypes.txt --out h2_results

# Run GWAS with mixed linear model
./scripts/07_gwas_mlma_loco.sh --bfile imputed_merged --grm project_grm --pheno phenotypes.txt --out gwas_results

# Identify independent signals
./scripts/08_cojo_analysis.sh --bfile imputed_merged --gwas-summary gwas_results.loco.mlma --out independent_signals

# Generate credible sets
./scripts/09_credible_sets.sh --bfile imputed_merged --gwas-summary gwas_results.loco.mlma --independent-signals independent_signals.jma.cojo --out credible_sets
```

### 5. Colocalization Analysis

```bash
# Run colocalization with eQTL data
./scripts/10_colocalization.sh --gwas-summary gwas_results.loco.mlma --eqtl-data eqtl_data.txt --out coloc_results
```

### 6. Polygenic Risk Score

```bash
# Generate PRS
./scripts/11_create_prs.sh --bfile test_dataset --gwas-summary gwas_results.loco.mlma --out prs_results
```

## Key Features

### GCTA-MLMA-LOCO Approach

The "leave one chromosome out" (LOCO) approach prevents the SNP being tested from contributing to both fixed and random effects, reducing statistical bias. This is particularly important for complex traits like breast cancer where effect sizes may be modest.

### SNP-Heritability Estimation

GCTA-GREML quantifies the proportion of phenotypic variance explained by common genetic variants, providing insight into genetic architecture.

### GCTA-COJO Analysis

Conditional and joint analysis identifies independent signals, distinguishing true causal variants from those in linkage disequilibrium.

### Integration of Covariates

Throughout the pipeline, important covariates (sex, age, principal components, study site) are incorporated to control for population stratification and other confounding factors.

## Example Usage

```bash
# Full pipeline execution
./run_pipeline.sh --bfile raw_data --pheno phenotypes.txt --covar covariates.txt --out project_name
```

## Advanced Configuration

The pipeline can be customized by editing `config.yaml` to specify:
- Quality control thresholds
- Covariates to include
- HPC resource requirements
- Imputation parameters
- Analysis-specific options

## Output Files

- `*_h2.hsq`: Heritability estimates
- `*.loco.mlma`: GWAS summary statistics
- `*.jma.cojo`: Independent significant variants
- `*_credible_set.txt`: Credible sets for each locus
- `*_coloc_results.txt`: Colocalization results
- `*_prs.profile`: Polygenic risk scores

## Troubleshooting

See `docs/troubleshooting.md` for common issues and solutions.

## Citations

If you use this pipeline in your research, please cite:

- Yang J, Lee SH, Goddard ME, Visscher PM. GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet. 2011;88(1):76-82.
- Purcell S, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet. 2007;81(3):559-575.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributors

- Your Name (jiyang.jiang@gmail.com)

## Acknowledgments

- Acknowledgments for any funding sources or collaborators
