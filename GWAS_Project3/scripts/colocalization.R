#!/usr/bin/env Rscript

# Colocalization Analysis Script
# This script performs colocalization analysis between GWAS and eQTL data

# Load required libraries
library(data.table)
library(coloc)

# Parse Snakemake parameters
out_prefix <- snakemake@params[["out_prefix"]]
window_size <- as.numeric(snakemake@params[["window_size"]])

# Input files
gwas_file <- snakemake@input[["gwas"]]
indep_file <- snakemake@input[["indep"]]
eqtl_file <- snakemake@input[["eqtl"]]

# Output file
coloc_results_file <- snakemake@output[["coloc"]]

# Set up logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "w")
sink(log_con, type = "message")
sink(log_con, type = "output")

cat("Starting colocalization analysis...\n")

# Read GWAS results
cat("Reading GWAS results from:", gwas_file, "\n")
gwas_data <- fread(gwas_file)
colnames(gwas_data) <- c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")

# Read independent signals from COJO
cat("Reading independent signals from:", indep_file, "\n")
indep_signals <- fread(indep_file, skip = 1)
colnames(indep_signals) <- c("Chr", "SNP", "bp", "refA", "freq", "b", "se", "p", "n", "freq_geno", "bJ", "bJ_se", "pJ", "LD_r")

# Read eQTL data
cat("Reading eQTL data from:", eqtl_file, "\n")
eqtl_data <- fread(eqtl_file)

# Assuming eQTL data has columns: chr, pos, SNP, gene, beta, se, p
if (!all(c("chr", "pos", "SNP", "gene", "beta", "se", "p") %in% colnames(eqtl_data))) {
  stop("eQTL data must have columns: chr, pos, SNP, gene, beta, se, p")
}

# Create directory for colocalization results
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

# Initialize results data frame
coloc_results <- data.frame(
  region = character(),
  index_snp = character(),
  gene = character(),
  PP_H0 = numeric(),
  PP_H1 = numeric(),
  PP_H2 = numeric(),
  PP_H3 = numeric(),
  PP_H4 = numeric(),
  stringsAsFactors = FALSE
)

# Process each independent signal
cat("Number of independent signals:", nrow(indep_signals), "\n")

for (i in 1:nrow(indep_signals)) {
  index_snp <- indep_signals$SNP[i]
  chr <- indep_signals$Chr[i]
  pos <- indep_signals$bp[i]
  
  cat(sprintf("Processing signal %d/%d: %s (chr%s:%s)\n", 
              i, nrow(indep_signals), index_snp, chr, pos))
  
  # Define region boundaries
  start_pos <- max(1, pos - window_size)
  end_pos <- pos + window_size
  region_name <- paste0("chr", chr, "_", start_pos, "_", end_pos)
  
  # Extract GWAS statistics for the region
  region_gwas <- gwas_data[gwas_data$Chr == chr & 
                           gwas_data$bp >= start_pos & 
                           gwas_data$bp <= end_pos, ]
  
  if (nrow(region_gwas) == 0) {
    cat("No GWAS SNPs found in region around", index_snp, "\n")
    next
  }
  
  # Extract eQTL data for the region
  region_eqtl <- eqtl_data[eqtl_data$chr == chr & 
                           eqtl_data$pos >= start_pos & 
                           eqtl_data$pos <= end_pos, ]
  
  if (nrow(region_eqtl) == 0) {
    cat("No eQTL SNPs found in region around", index_snp, "\n")
    next
  }
  
  # Get unique genes in the region
  genes <- unique(region_eqtl$gene)
  cat("Number of genes in region:", length(genes), "\n")
  
  # Perform colocalization for each gene
  for (gene in genes) {
    cat("Processing gene:", gene, "\n")
    
    # Extract eQTL data for this gene
    gene_eqtl <- region_eqtl[region_eqtl$gene == gene, ]
    
    # Find common SNPs between GWAS and eQTL data
    common_snps <- intersect(region_gwas$SNP, gene_eqtl$SNP)
    
    if (length(common_snps) < 5) {
      cat("Too few common SNPs for gene", gene, ": ", length(common_snps), "\n")
      next
    }
    
    cat("Number of common SNPs:", length(common_snps), "\n")
    
    # Filter to common SNPs
    gwas_subset <- region_gwas[region_gwas$SNP %in% common_snps, ]
    eqtl_subset <- gene_eqtl[gene_eqtl$SNP %in% common_snps, ]
    
    # Ensure SNPs are in the same order
    setkey(gwas_subset, SNP)
    setkey(eqtl_subset, SNP)
    
    # Prepare datasets for coloc
    dataset1 <- list(
      snp = gwas_subset$SNP,
      position = gwas_subset$bp,
      beta = gwas_subset$b,
      varbeta = gwas_subset$se^2,
      type = "cc",         # Case-control study
      s = 0.4,             # Approximate case proportion
      N = 10000            # Sample size (adjust as needed)
    )
    
    dataset2 <- list(
      snp = eqtl_subset$SNP,
      position = eqtl_subset$pos,
      beta = eqtl_subset$beta,
      varbeta = eqtl_subset$se^2,
      type = "quant",      # Quantitative trait
      N = 500              # Sample size (adjust as needed)
    )
    
    # Run colocalization analysis
    cat("Running colocalization analysis...\n")
    result <- tryCatch({
      coloc.abf(dataset1, dataset2)
    }, error = function(e) {
      cat("Error in colocalization for gene", gene, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(result)) {
      next
    }
    
    # Extract results
    pp <- result$summary
    
    # Add to results data frame
    coloc_results <- rbind(coloc_results, data.frame(
      region = region_name,
      index_snp = index_snp,
      gene = gene,
      PP_H0 = pp[1],  # No association with either trait
      PP_H1 = pp[2],  # Association with trait 1 only
      PP_H2 = pp[3],  # Association with trait 2 only
      PP_H3 = pp[4],  # Association with both traits, different causal variants
      PP_H4 = pp[5],  # Association with both traits, shared causal variant
      stringsAsFactors = FALSE
    ))
    
    # Create detailed output for this gene
    gene_file <- paste0(out_prefix, "_", region_name, "_", gene, ".txt")
    sink(gene_file)
    cat("Trait 1: GWAS\n")
    cat("Trait 2: eQTL for", gene, "\n\n")
    cat("Posterior probabilities:\n")
    cat("H0 (no association with either trait):", pp[1], "\n")
    cat("H1 (association with trait 1 only):", pp[2], "\n")
    cat("H2 (association with trait 2 only):", pp[3], "\n")
    cat("H3 (association with both traits, different causal variants):", pp[4], "\n")
    cat("H4 (association with both traits, shared causal variant):", pp[5], "\n\n")
    
    if (pp[5] > 0.75) {
      cat("STRONG EVIDENCE for colocalization (PP.H4 > 0.75)\n")
    } else if (pp[5] > 0.5) {
      cat("SUGGESTIVE EVIDENCE for colocalization (PP.H4 > 0.5)\n")
    } else {
      cat("WEAK or NO EVIDENCE for colocalization (PP.H4 <= 0.5)\n")
    }
    sink()
  }
}

# Write combined colocalization results to file
cat("Writing combined colocalization results to:", coloc_results_file, "\n")
fwrite(coloc_results, coloc_results_file, sep = "\t", quote = FALSE)

# Create summary of strong colocalizations
strong_coloc <- coloc_results[coloc_results$PP_H4 > 0.75, ]
if (nrow(strong_coloc) > 0) {
  strong_file <- paste0(out_prefix, "_strong_colocalizations.txt")
  cat("Number of strong colocalizations:", nrow(strong_coloc), "\n")
  fwrite(strong_coloc, strong_file, sep = "\t", quote = FALSE)
} else {
  cat("No strong colocalizations found.\n")
}

cat("Colocalization analysis complete.\n")
