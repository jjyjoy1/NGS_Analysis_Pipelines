#!/bin/bash

# STEP 6: Create GRM with GCTA
$GCTA --bfile ${OUT_PREFIX}_imputed \
    --extract ${OUT_PREFIX}_prune.prune.in \
    --make-grm \
    --out ${OUT_PREFIX}_grm





