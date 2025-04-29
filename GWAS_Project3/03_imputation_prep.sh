#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# Continue from QC_complete dataset
QC_DATA="${OUT_PREFIX}_QC_complete"

# STEP 3: IMPUTATION PREPARATION
# Split by chromosome and convert to VCF for imputation
for CHR in {1..22}
do
    $PLINK --bfile $QC_DATA \
        --chr $CHR \
        --recode vcf \
        --out ${OUT_PREFIX}_chr${CHR}
    
    # Typically you'd submit these VCFs to an imputation server
    # like Michigan Imputation Server or TOPMed
done





