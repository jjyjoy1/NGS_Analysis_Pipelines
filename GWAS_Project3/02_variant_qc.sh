#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 2: VARIANT QC WITH PLINK
echo "Starting Variant QC..."

# Check variant call rates and remove variants with low call rates
$PLINK --bfile ${OUT_PREFIX}_sample_qc \
    --geno 0.05 \
    --make-bed \
    --out ${OUT_PREFIX}_geno

# Filter by Minor Allele Frequency (MAF)
$PLINK --bfile ${OUT_PREFIX}_geno \
    --maf 0.01 \
    --make-bed \
    --out ${OUT_PREFIX}_maf

# Filter by Hardy-Weinberg Equilibrium (HWE) - apply only to controls
$PLINK --bfile ${OUT_PREFIX}_maf \
    --hwe 1e-6 \
    --keep-if-cluster --filter controls.txt \
    --make-bed \
    --out ${OUT_PREFIX}_hwe

# Check for SNPs with differential missingness between cases and controls
$PLINK --bfile ${OUT_PREFIX}_hwe \
    --test-missing \
    --out ${OUT_PREFIX}_diffmiss

# Remove SNPs with differential missingness
awk '{ if ($5 < 1e-4) print $2 }' ${OUT_PREFIX}_diffmiss.missing > ${OUT_PREFIX}_diffmiss_snps.txt

$PLINK --bfile ${OUT_PREFIX}_hwe \
    --exclude ${OUT_PREFIX}_diffmiss_snps.txt \
    --make-bed \
    --out ${OUT_PREFIX}_QC_complete

echo "Sample and Variant QC completed. Final cleaned dataset: ${OUT_PREFIX}_QC_complete"



