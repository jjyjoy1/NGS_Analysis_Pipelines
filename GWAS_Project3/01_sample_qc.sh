#!/bin/bash

# Set paths to your tools and data
PLINK="/path/to/plink"
GCTA="/path/to/gcta"
RAW_DATA="your_raw_genotype_data"
OUT_PREFIX="breast_cancer_gwas"

# STEP 1: SAMPLE QC WITH PLINK
echo "Starting Sample QC..."

# Check sample call rates and remove samples with low call rates
$PLINK --bfile $RAW_DATA \
    --mind 0.05 \
    --make-bed \
    --out ${OUT_PREFIX}_mind

# Sex check - identify mismatches
$PLINK --bfile ${OUT_PREFIX}_mind \
    --check-sex \
    --out ${OUT_PREFIX}_sexcheck

# Remove samples with sex discrepancies
awk '{ if ($5 == "PROBLEM") print $1,$2 }' ${OUT_PREFIX}_sexcheck.sexcheck > ${OUT_PREFIX}_sex_discrepancy.txt

$PLINK --bfile ${OUT_PREFIX}_mind \
    --remove ${OUT_PREFIX}_sex_discrepancy.txt \
    --make-bed \
    --out ${OUT_PREFIX}_sexchecked

# Check for excessive heterozygosity
$PLINK --bfile ${OUT_PREFIX}_sexchecked \
    --het \
    --out ${OUT_PREFIX}_het

# Identify samples with outlier heterozygosity (±3 SD from mean)
Rscript -e '
  het <- read.table("'${OUT_PREFIX}'_het.het", header=TRUE);
  het$F <- (het$O.HOM - het$E.HOM)/het$E.HOM;
  mean_f <- mean(het$F);
  sd_f <- sd(het$F);
  outliers <- subset(het, F < mean_f-3*sd_f | F > mean_f+3*sd_f);
  write.table(outliers[,c(1,2)], "'${OUT_PREFIX}'_het_outliers.txt", 
              row.names=FALSE, col.names=FALSE, quote=FALSE);'

# Remove heterozygosity outliers
$PLINK --bfile ${OUT_PREFIX}_sexchecked \
    --remove ${OUT_PREFIX}_het_outliers.txt \
    --make-bed \
    --out ${OUT_PREFIX}_het_filtered

# Check for relatedness using Identity-By-Descent (IBD)
$PLINK --bfile ${OUT_PREFIX}_het_filtered \
    --indep-pairwise 50 5 0.2 \
    --out ${OUT_PREFIX}_pruned

$PLINK --bfile ${OUT_PREFIX}_het_filtered \
    --extract ${OUT_PREFIX}_pruned.prune.in \
    --genome \
    --min 0.185 \
    --out ${OUT_PREFIX}_ibd

# Keep one individual from each related pair (IBD > 0.185, approximately second-degree relatives)
awk '{ if ($9 > 0.185) print $1,$2,$3,$4,$9,$10 }' ${OUT_PREFIX}_ibd.genome > ${OUT_PREFIX}_related_pairs.txt

# Custom script to select one from each related cluster - keeping the one with higher call rate
Rscript -e '
  related <- read.table("'${OUT_PREFIX}'_related_pairs.txt", header=FALSE);
  samples_to_remove <- c();
  for(i in 1:nrow(related)) {
    id1 <- paste(related[i,1], related[i,2], sep=":");
    id2 <- paste(related[i,3], related[i,4], sep=":");
    if(!(id1 %in% samples_to_remove) && !(id2 %in% samples_to_remove)) {
      samples_to_remove <- c(samples_to_remove, id2);
    }
  }
  ids <- do.call(rbind, strsplit(samples_to_remove, ":"));
  write.table(ids, "'${OUT_PREFIX}'_related_to_remove.txt", 
              row.names=FALSE, col.names=FALSE, quote=FALSE);'

$PLINK --bfile ${OUT_PREFIX}_het_filtered \
    --remove ${OUT_PREFIX}_related_to_remove.txt \
    --make-bed \
    --out ${OUT_PREFIX}_unrelated

# Run PCA to identify population outliers using GCTA
$GCTA --bfile ${OUT_PREFIX}_unrelated \
    --make-grm \
    --out ${OUT_PREFIX}_grm

$GCTA --grm ${OUT_PREFIX}_grm \
    --pca 20 \
    --out ${OUT_PREFIX}_pca

# Identify population outliers (manual inspection needed)
# After visual inspection of PCs, create a file with outliers to remove
# Assuming you've identified outliers and saved them in a file called pca_outliers.txt
$PLINK --bfile ${OUT_PREFIX}_unrelated \
    --remove pca_outliers.txt \
    --make-bed \
    --out ${OUT_PREFIX}_sample_qc


