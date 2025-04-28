#!/bin/bash

# STEP 9: Identify independent signals with GCTA-COJO
$GCTA --bfile ${OUT_PREFIX}_imputed \
    --cojo-file ${OUT_PREFIX}_gwas.loco.mlma \
    --cojo-slct \
    --cojo-p 5e-8 \
    --out ${OUT_PREFIX}_cojo




