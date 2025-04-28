#!/bin/bash

# STEP 4: After imputation, convert back to PLINK format
# Assuming imputed files have been downloaded as VCFs
for CHR in {1..22}
do
    $PLINK --vcf imputed_chr${CHR}.vcf.gz \
        --make-bed \
        --out imputed_chr${CHR}
done

# Merge imputed files
echo "imputed_chr1" > merge_list.txt
for CHR in {2..22}
do
    echo "imputed_chr${CHR}" >> merge_list.txt
done

$PLINK --merge-list merge_list.txt \
    --make-bed \
    --out ${OUT_PREFIX}_imputed



