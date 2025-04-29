#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 7: Estimate SNP heritability with GCTA-GREML
$GCTA --grm ${OUT_PREFIX}_grm \
    --pheno phenotype.txt \
    --covar covariates.txt \
    --qcovar qcovariates.txt \
    --reml \
    --out ${OUT_PREFIX}_h2


