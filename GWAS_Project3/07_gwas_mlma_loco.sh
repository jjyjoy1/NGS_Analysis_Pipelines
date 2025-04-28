#!/bin/bash

# STEP 7: Estimate SNP heritability with GCTA-GREML
$GCTA --grm ${OUT_PREFIX}_grm \
    --pheno phenotype.txt \
    --covar covariates.txt \
    --qcovar qcovariates.txt \
    --reml \
    --out ${OUT_PREFIX}_h2


