#!/bin/bash

# STEP 8: Run GWAS with GCTA-MLMA-LOCO
$GCTA --mlma-loco \
    --bfile ${OUT_PREFIX}_imputed \
    --grm ${OUT_PREFIX}_grm \
    --pheno phenotype.txt \
    --covar covariates.txt \
    --qcovar qcovariates.txt \
    --out ${OUT_PREFIX}_gwas



