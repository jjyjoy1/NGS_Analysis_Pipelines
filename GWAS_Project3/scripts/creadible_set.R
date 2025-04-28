#!/usr/bin/env Rscript

# Credible Set Analysis Script
# This script calculates 95% credible sets for each locus identified in the COJO analysis

# Load required libraries
library(data.table)
library(ggplot2)

# Parse Snakemake parameters
bfile_prefix <- snakemake@params[["bfile_prefix"]]
out_prefix <- snakemake@params[["out_prefix"]]
window_size <- as.numeric(snakemake@params[["window_size"]])

# Input files
gwas_file <- snakemake@input[["gwas"]]
indep_file <- snakemake@input[["indep"]]

# Output file
credible_sets_file <- snakemake@output[["credible_sets"]]

# Set up logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "w")
sink(log_con, type = "message")
sink(log_con, type = "output")

cat("Starting credible set analysis...\n")

# Read GWAS results
cat("Reading GWAS results from:", gwas_file, "\n")
gwas_data <- fread(gwas_file)
colnames(gwas_data) <- c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")

# Read independent signals from COJO
cat("Reading independent signals from:", indep_file, "\n")
indep_signals <- fread(indep_file, skip = 1)
colnames(indep_signals) <- c("Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se", "pJ", "LD_r")

# Create directory for credible set results
dir.create(file.path(dirname(out_prefix)), showWarnings = FALSE, recursive = TRUE)

# Function to calculate Bayes factors and posterior probabilities
calculate_posteriors <- function(stats) {
  # Calculate Approximate Bayes Factors using Wakefield approximation
  W <- 0.04  # Prior variance
  V <- stats$se^2
  z <- stats$b / stats$se
  r <- W / (W + V)
  BF <- sqrt(1 - r) * exp(z^2 * r / 2)
  
  # Calculate posterior probabilities
  PP <- BF / sum(BF)
  
  # Create results data frame
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
  )
  
  # Sort by posterior probability
  result <- result[order(-result$PP), ]
  
  # Identify 95% credible set
  result$cumsum_PP <- cumsum(result$PP)
  result$in_credible_set <- result$cumsum_PP <= 0.95
  
  return(result)
}

# Process each independent signal
cat("Number of independent signals:", nrow(indep_signals), "\n")
all_credible_sets <- data.frame()

for (i in 1:nrow(indep_signals)) {
  index_snp <- indep_signals$SNP[i]
  chr <- indep_signals$Chr[i]
  pos <- indep_signals$bp[i]
  
  cat(sprintf("Processing signal %d/%d: %s (chr%s:%s)\n", 
              i, nrow(indep_signals), index_snp, chr, pos))
  
  # Define region boundaries
  start_pos <- max(1, pos - window_size)
  end_pos <- pos + window_size
  
  # Extract GWAS statistics for the region
  region_stats <- gwas_data[gwas_data$Chr == chr & 
                           gwas_data$bp >= start_pos & 
                           gwas_data$bp <= end_pos, ]
  
  if (nrow(region_stats) == 0) {
    cat("No SNPs found in region around", index_snp, "\n")
    next
  }
  
  cat("Number of SNPs in region:", nrow(region_stats), "\n")
  
  # Calculate posteriors
  credible_set <- calculate_posteriors(region_stats)
  
  # Add locus identifier
  credible_set$locus <- paste0("locus_", i, "_", index_snp)
  credible_set$index_snp <- index_snp
  
  # Save credible set for this locus
  locus_file <- paste0(out_prefix, "_locus_", i, ".txt")
  fwrite(credible_set, locus_file, sep = "\t", quote = FALSE)
  
  # Extract 95% credible set
  cs_95 <- credible_set[credible_set$in_credible_set, ]
  cat("Number of SNPs in 95% credible set:", nrow(cs_95), "\n")
  
  # Add to combined results
  all_credible_sets <- rbind(all_credible_sets, cs_95)
  
  # Optional: Create plot for this locus
  plot_file <- paste0(out_prefix, "_locus_", i, "_plot.pdf")
  pdf(plot_file, width = 10, height = 6)
  p <- ggplot(credible_set, aes(x = position, y = -log10(pvalue), color = in_credible_set)) +
    geom_point(aes(size = PP)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
    ggtitle(paste0("Locus ", i, ": ", index_snp)) +
    xlab(paste0("Position on chromosome ", chr)) +
    ylab("-log10(p-value)") +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "blue")) +
    scale_size_continuous(range = c(1, 5)) +
    theme_bw() +
    theme(legend.position = "bottom")
  print(p)
  dev.off()
}

# Write combined credible sets to file
cat("Writing combined credible sets to:", credible_sets_file, "\n")
fwrite(all_credible_sets, credible_sets_file, sep = "\t", quote = FALSE)

# Create summary table

