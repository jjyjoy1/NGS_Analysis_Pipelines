#!/bin/bash

# STEP 5: Prune SNPs for GRM
$PLINK --bfile ${OUT_PREFIX}_imputed \
    --indep-pairwise 50 5 0.2 \
    --out ${OUT_PREFIX}_prune


