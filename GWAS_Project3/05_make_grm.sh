#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 5: Prune SNPs for GRM
$PLINK --bfile ${OUT_PREFIX}_imputed \
    --indep-pairwise 50 5 0.2 \
    --out ${OUT_PREFIX}_prune


