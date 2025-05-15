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

# STEP 5: Prune SNPs for GRM
$PLINK --bfile ${OUT_PREFIX}_imputed \
    --indep-pairwise 50 5 0.2 \
    --out ${OUT_PREFIX}_prune

# STEP 6: Create GRM with GCTA
$GCTA --bfile ${OUT_PREFIX}_imputed \
    --extract ${OUT_PREFIX}_prune.prune.in \
    --make-grm \
    --out ${OUT_PREFIX}_grm

# STEP 7: Estimate SNP heritability with GCTA-GREML
$GCTA --grm ${OUT_PREFIX}_grm \
    --pheno phenotype.txt \
    --covar covariates.txt \
    --qcovar qcovariates.txt \
    --reml \
    --out ${OUT_PREFIX}_h2

# STEP 8: Run GWAS with GCTA-MLMA-LOCO
$GCTA --mlma-loco \
    --bfile ${OUT_PREFIX}_imputed \
    --grm ${OUT_PREFIX}_grm \
    --pheno phenotype.txt \
    --covar covariates.txt \
    --qcovar qcovariates.txt \
    --out ${OUT_PREFIX}_gwas

# STEP 9: Identify independent signals with GCTA-COJO
$GCTA --bfile ${OUT_PREFIX}_imputed \
    --cojo-file ${OUT_PREFIX}_gwas.loco.mlma \
    --cojo-slct \
    --cojo-p 5e-8 \
    --out ${OUT_PREFIX}_cojo

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

# STEP 11: Generate PRS
# Get effect sizes from GWAS
awk '{print $2,$8}' ${OUT_PREFIX}_gwas.loco.mlma > ${OUT_PREFIX}_effects.txt

# Create PRS
$PLINK --bfile ${OUT_PREFIX}_imputed \
    --score ${OUT_PREFIX}_effects.txt 1 2 header \
    --out ${OUT_PREFIX}_prs


# Continue from previous workflow
# Assuming we've already completed GWAS and identified significant loci

# STEP 12: CREDIBLE SET ANALYSIS
# This identifies sets of SNPs that are likely to contain the causal variant
echo "Starting credible set analysis..."

# First identify significant regions from GWAS results
awk '$13 < 5e-8 {print $1, $3, $3, $2}' ${OUT_PREFIX}_gwas.loco.mlma | \
    sort -k1,1 -k2,2n > ${OUT_PREFIX}_significant_snps.txt

# Merge nearby significant SNPs into regions (within 500kb)
Rscript -e '
  snps <- read.table("'${OUT_PREFIX}'_significant_snps.txt", header=FALSE);
  colnames(snps) <- c("chr", "pos", "pos_end", "rsid");
  snps <- snps[order(snps$chr, snps$pos),];
  regions <- list();
  current_region <- c(snps[1,]);
  for(i in 2:nrow(snps)) {
    if(snps[i,]$chr == current_region[nrow(current_region),]$chr && 
       snps[i,]$pos - current_region[nrow(current_region),]$pos < 500000) {
      current_region <- rbind(current_region, snps[i,]);
    } else {
      regions[[length(regions) + 1]] <- current_region;
      current_region <- c(snps[i,]);
    }
  }
  regions[[length(regions) + 1]] <- current_region;
  
  # Output regions
  region_file <- data.frame(
    chr=numeric(), start=numeric(), end=numeric(), index_snp=character()
  );
  
  for(i in 1:length(regions)) {
    region <- regions[[i]];
    # Find index SNP (most significant in region)
    index_snp <- region[which.min(region$pval),]$rsid;
    region_file[i,] <- c(
      region[1,]$chr, 
      min(region$pos) - 50000,
      max(region$pos) + 50000,
      index_snp
    );
  }
  
  write.table(region_file, "'${OUT_PREFIX}'_regions.txt", 
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t");'

# For each region, extract SNPs and calculate posterior probabilities
mkdir -p credible_sets
cat ${OUT_PREFIX}_regions.txt | while read chr start end index_snp
do
  region_name="chr${chr}_${start}_${end}"
  
  # Extract SNPs in the region
  $PLINK --bfile ${OUT_PREFIX}_imputed \
      --chr $chr \
      --from-bp $start \
      --to-bp $end \
      --make-bed \
      --out "credible_sets/${region_name}"
  
  # Extract GWAS statistics for the region
  awk -v chr=$chr -v start=$start -v end=$end \
      '$1 == chr && $3 >= start && $3 <= end' ${OUT_PREFIX}_gwas.loco.mlma \
      > "credible_sets/${region_name}.stats"
  
  # Calculate Bayes factors and posterior probabilities
  Rscript -e '
    stats <- read.table("credible_sets/'${region_name}'.stats", header=TRUE);
    
    # Calculate Approximate Bayes Factors
    # Using Wakefield approximation
    W <- 0.04;  # Prior variance (standard choice for GWAS)
    V <- stats$se^2;
    z <- stats$b / stats$se;
    r <- W / (W + V);
    BF <- sqrt(1 - r) * exp(z^2 * r / 2);
    
    # Calculate posterior probabilities
    PP <- BF / sum(BF);
    
    # Sort by posterior probability
    result <- data.frame(
      SNP = stats$SNP,
      chromosome = stats$Chr,
      position = stats$bp,
      allele1 = stats$A1,
      allele2 = stats$A2,
      beta = stats$b,
      se = stats$se,
      pvalue = stats$p,
      BF = BF,
      PP = PP
    );
    
    result <- result[order(-result$PP),];
    
    # Identify 95% credible set
    cumsum_PP <- cumsum(result$PP);
    result$in_credible_set <- cumsum_PP <= 0.95;
    
    # Output full results and credible set
    write.table(result, "credible_sets/'${region_name}'_full_results.txt", 
                row.names=FALSE, quote=FALSE, sep="\t");
    write.table(result[result$in_credible_set,], "credible_sets/'${region_name}'_credible_set.txt", 
                row.names=FALSE, quote=FALSE, sep="\t");
    
    # Output summary
    cat("Region: ", "'${region_name}'", "\n");
    cat("Index SNP: ", "'${index_snp}'", "\n");
    cat("Number of SNPs in 95% credible set: ", sum(result$in_credible_set), "\n");
    cat("Top SNP in credible set: ", result$SNP[1], " (PP = ", result$PP[1], ")\n");
  ' > "credible_sets/${region_name}_summary.txt"
done

# Create combined report of all credible sets
cat credible_sets/*_summary.txt > ${OUT_PREFIX}_credible_sets_summary.txt

# STEP 13: COLOCALIZATION ANALYSIS
# This tests whether two traits share the same causal variants
echo "Starting colocalization analysis..."

# For breast cancer GWAS, we'll compare with eQTL data to identify potential mechanisms
# Assuming you have eQTL summary statistics files for relevant tissues

# For example, GTEx breast tissue eQTL data
EQTL_FILE="path/to/breast_tissue_eQTL.txt"

# Create directory for colocalization results
mkdir -p colocalization

# For each significant region, perform colocalization with eQTL data
cat ${OUT_PREFIX}_regions.txt | while read chr start end index_snp
do
  region_name="chr${chr}_${start}_${end}"
  
  # Extract GWAS statistics for the region
  awk -v chr=$chr -v start=$start -v end=$end \
      '$1 == chr && $3 >= start && $3 <= end' ${OUT_PREFIX}_gwas.loco.mlma \
      > "colocalization/${region_name}_gwas.txt"
  
  # Extract eQTL statistics for the region
  awk -v chr=$chr -v start=$start -v end=$end \
      '$1 == chr && $2 >= start && $2 <= end' $EQTL_FILE \
      > "colocalization/${region_name}_eqtl_all.txt"
  
  # Identify all genes tested in this region
  cut -f4 "colocalization/${region_name}_eqtl_all.txt" | sort | uniq > "colocalization/${region_name}_genes.txt"
  
  # Perform colocalization for each gene using coloc R package
  cat "colocalization/${region_name}_genes.txt" | while read gene
  do
    # Extract eQTL data for this gene
    grep $gene "colocalization/${region_name}_eqtl_all.txt" > "colocalization/${region_name}_${gene}_eqtl.txt"
    
    # Run coloc analysis using R
    Rscript -e '
      library(coloc);
      
      # Read datasets
      gwas <- read.table("colocalization/'${region_name}'_gwas.txt", header=TRUE);
      eqtl <- read.table("colocalization/'${region_name}'_'${gene}'_eqtl.txt", header=TRUE);
      
      # Prepare GWAS dataset
      dataset1 <- list(
        snp = gwas$SNP,
        position = gwas$bp,
        beta = gwas$b,
        varbeta = gwas$se^2,
        type = "cc",
        s = 0.4,  # Approximate case proportion in study
        N = 20000  # Total sample size - adjust as needed
      );
      
      # Prepare eQTL dataset (assuming standard format - adjust column names as needed)
      dataset2 <- list(
        snp = eqtl$SNP,
        position = eqtl$position,
        beta = eqtl$beta,
        varbeta = eqtl$se^2,
        type = "quant",
        N = 500  # GTEx sample size - adjust as needed
      );
      
      # Run colocalization analysis
      result <- coloc.abf(dataset1, dataset2);
      
      # Output results
      sink("colocalization/'${region_name}'_'${gene}'_coloc_results.txt");
      cat("Trait 1: Breast Cancer GWAS\n");
      cat("Trait 2: eQTL for '${gene}'\n\n");
      cat("Posterior probabilities:\n");
      cat("H0 (no association with either trait): ", result$summary[1], "\n");
      cat("H1 (association with trait 1 only): ", result$summary[2], "\n");
      cat("H2 (association with trait 2 only): ", result$summary[3], "\n");
      cat("H3 (association with both traits, different causal variants): ", result$summary[4], "\n");
      cat("H4 (association with both traits, shared causal variant): ", result$summary[5], "\n\n");
      
      if (result$summary[5] > 0.75) {
        cat("STRONG EVIDENCE for colocalization (PP.H4 > 0.75)\n");
      } else if (result$summary[5] > 0.5) {
        cat("SUGGESTIVE EVIDENCE for colocalization (PP.H4 > 0.5)\n");
      } else {
        cat("WEAK or NO EVIDENCE for colocalization (PP.H4 <= 0.5)\n");
      }
      sink();
      
      # Create simplified output for summary
      coloc_summary <- data.frame(
        region = "'${region_name}'",
        gene = "'${gene}'",
        PP_H0 = result$summary[1],
        PP_H1 = result$summary[2],
        PP_H2 = result$summary[3],
        PP_H3 = result$summary[4],
        PP_H4 = result$summary[5]
      );
      
      write.table(coloc_summary, "colocalization/'${region_name}'_'${gene}'_summary.txt", 
                  row.names=FALSE, quote=FALSE, sep="\t");
    '
  done
done

# Combine all colocalization results into a single file
echo -e "Region\tGene\tPP_H0\tPP_H1\tPP_H2\tPP_H3\tPP_H4" > ${OUT_PREFIX}_colocalization_results.txt
cat colocalization/*_summary.txt >> ${OUT_PREFIX}_colocalization_results.txt

# Generate report of strong colocalizations (PP.H4 > 0.75)
awk '$7 > 0.75 {print $0}' ${OUT_PREFIX}_colocalization_results.txt > ${OUT_PREFIX}_strong_colocalizations.txt

echo "Completed credible set and colocalization analyses."
