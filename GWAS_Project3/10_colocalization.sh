#!/bin/bash

# STEP 10: Fine-mapping (example for one significant region)
# Extract SNPs in region (e.g., 1Mb window around top SNP)
TOP_SNP="rs12345"
CHR=$(grep $TOP_SNP ${OUT_PREFIX}_gwas.loco.mlma | cut -f 1)
POS=$(grep $TOP_SNP ${OUT_PREFIX}_gwas.loco.mlma | cut -f 3)
START=$((POS - 500000))
END=$((POS + 500000))

$PLINK --bfile ${OUT_PREFIX}_imputed \
    --chr $CHR \
    --from-bp $START \
    --to-bp $END \
    --make-bed \
    --out ${OUT_PREFIX}_finemapping_region

# Fine-mapping with GCTA-COJO conditional analysis
$GCTA --bfile ${OUT_PREFIX}_finemapping_region \
    --cojo-file ${OUT_PREFIX}_gwas.loco.mlma \
    --cojo-cond $TOP_SNP \
    --out ${OUT_PREFIX}_finemapping




