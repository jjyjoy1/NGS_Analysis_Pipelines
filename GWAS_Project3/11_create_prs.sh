#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 11: Generate PRS
# Get effect sizes from GWAS
awk '{print $2,$8}' ${OUT_PREFIX}_gwas.loco.mlma > ${OUT_PREFIX}_effects.txt

# Create PRS
$PLINK --bfile ${OUT_PREFIX}_imputed \
    --score ${OUT_PREFIX}_effects.txt 1 2 header \
    --out ${OUT_PREFIX}_prs




