#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 8: Run GWAS with GCTA-MLMA-LOCO
$GCTA --mlma-loco \
    --bfile ${OUT_PREFIX}_imputed \
    --grm ${OUT_PREFIX}_grm \
    --pheno phenotype.txt \
    --covar covariates.txt \
    --qcovar qcovariates.txt \
    --out ${OUT_PREFIX}_gwas



