#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 9: Identify independent signals with GCTA-COJO
$GCTA --bfile ${OUT_PREFIX}_imputed \
    --cojo-file ${OUT_PREFIX}_gwas.loco.mlma \
    --cojo-slct \
    --cojo-p 5e-8 \
    --out ${OUT_PREFIX}_cojo




