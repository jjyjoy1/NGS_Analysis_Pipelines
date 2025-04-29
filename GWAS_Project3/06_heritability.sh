#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 6: Create GRM with GCTA
$GCTA --bfile ${OUT_PREFIX}_imputed \
    --extract ${OUT_PREFIX}_prune.prune.in \
    --make-grm \
    --out ${OUT_PREFIX}_grm





