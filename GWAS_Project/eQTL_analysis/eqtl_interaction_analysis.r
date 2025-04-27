# eQTL Interaction Analysis Pipeline
# This script tests for SNP effects on expression that depend on covariates

# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("VariantAnnotation", quietly = TRUE)) BiocManager::install("VariantAnnotation")
if (!requireNamespace("MatrixEQTL", quietly = TRUE)) BiocManager::install("MatrixEQTL")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Load required libraries
library(lme4)
library(data.table)
library(VariantAnnotation)
library(MatrixEQTL)
library(ggplot2)

# Set file paths
vcf_file <- "genotypes.vcf"
expression_file <- "expression_counts.csv"
metadata_file <- "sample_metadata.csv"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

# Step 1: Load and format data
message("Loading data...")

# Read sample metadata
metadata <- fread(metadata_file)
print("Sample metadata summary:")
print(summary(metadata))

# Read expression data
expression_data <- fread(expression_file)
# Convert to matrix format for MatrixEQTL
gene_ids <- expression_data$gene_id  # Assuming first column has gene IDs
expr_matrix <- as.matrix(expression_data[, -1, with=FALSE])  # Remove gene_id column
rownames(expr_matrix) <- gene_ids
sample_ids <- colnames(expr_matrix)

# Read VCF file with VariantAnnotation
vcf <- readVcf(vcf_file)
# Extract genotype data and convert to numeric matrix
genotypes <- geno(vcf)$GT
# Convert genotypes to numeric (0, 1, 2 for ref/ref, ref/alt, alt/alt)
geno_matrix <- matrix(NA, nrow=nrow(genotypes), ncol=ncol(genotypes))
rownames(geno_matrix) <- rownames(genotypes)
colnames(geno_matrix) <- colnames(genotypes)
for (i in 1:nrow(genotypes)) {
  for (j in 1:ncol(genotypes)) {
    if (genotypes[i,j] == "0/0") geno_matrix[i,j] <- 0
    else if (genotypes[i,j] == "0/1" || genotypes[i,j] == "1/0") geno_matrix[i,j] <- 1
    else if (genotypes[i,j] == "1/1") geno_matrix[i,j] <- 2
    # Handle missing values
    else geno_matrix[i,j] <- NA
  }
}

# Ensure that sample IDs match across datasets
common_samples <- Reduce(intersect, list(metadata$sample_id, colnames(expr_matrix), colnames(geno_matrix)))
message(paste("Using", length(common_samples), "samples found in all datasets"))

metadata <- metadata[match(common_samples, metadata$sample_id),]
expr_matrix <- expr_matrix[, common_samples]
geno_matrix <- geno_matrix[, common_samples]

# Step 2: Prepare covariates of interest
# For this example, we'll use cell_type and disease_status as interaction covariates
covariates_of_interest <- c("cell_type", "disease_status")
print(paste("Using covariates for interaction:", paste(covariates_of_interest, collapse=", ")))

# Step 3: Function to test SNP-covariate interactions using lme4
test_interaction_lme4 <- function(gene_expr, snp_geno, metadata, covariate) {
  # Create data frame for modeling
  model_data <- data.frame(
    expression = gene_expr,
    genotype = snp_geno,
    covariate = metadata[[covariate]]
  )
  
  # Remove missing values
  model_data <- na.omit(model_data)
  
  # Fit models with and without interaction
  model_main <- lm(expression ~ genotype + covariate, data=model_data)
  model_interaction <- lm(expression ~ genotype * covariate, data=model_data)
  
  # Compare models
  anova_result <- anova(model_main, model_interaction)
  pval <- anova_result$`Pr(>F)`[2]
  
  # Get interaction coefficient
  interaction_coef <- summary(model_interaction)$coefficients[grep(":", rownames(summary(model_interaction)$coefficients)), ]
  
  return(list(
    pvalue = pval,
    coefficients = interaction_coef
  ))
}

# Step 4: Run interaction analysis for each gene-SNP pair and covariate
# For demonstration, we'll limit to top genes and SNPs
top_genes <- head(rownames(expr_matrix), 100)
top_snps <- head(rownames(geno_matrix), 100)

results_list <- list()

for (covariate in covariates_of_interest) {
  message(paste("Testing interactions with covariate:", covariate))
  
  result_table <- data.frame()
  
  for (gene in top_genes) {
    for (snp in top_snps) {
      gene_expr <- expr_matrix[gene, ]
      snp_geno <- geno_matrix[snp, ]
      
      # Test interaction
      test_result <- test_interaction_lme4(gene_expr, snp_geno, metadata, covariate)
      
      # Add to results
      result_row <- data.frame(
        gene = gene,
        snp = snp,
        covariate = covariate,
        pvalue = test_result$pvalue,
        stringsAsFactors = FALSE
      )
      
      result_table <- rbind(result_table, result_row)
    }
  }
  
  # Adjust p-values for multiple testing
  result_table$padj <- p.adjust(result_table$pvalue, method = "BH")
  
  # Sort by adjusted p-value
  result_table <- result_table[order(result_table$padj), ]
  
  # Save results
  results_list[[covariate]] <- result_table
  output_file <- file.path(output_dir, paste0("interaction_results_", covariate, ".csv"))
  write.csv(result_table, output_file, row.names = FALSE)
  message(paste("Results saved to", output_file))
}

# Step 5: Visualize top interactions
plot_top_interactions <- function(results_list, expr_matrix, geno_matrix, metadata, n_top=5) {
  for (covariate in names(results_list)) {
    result_table <- results_list[[covariate]]
    top_results <- head(result_table, n_top)
    
    message(paste("Plotting top", n_top, "interactions for", covariate))
    
    for (i in 1:nrow(top_results)) {
      gene <- top_results$gene[i]
      snp <- top_results$snp[i]
      
      # Extract data
      plot_data <- data.frame(
        expression = expr_matrix[gene, ],
        genotype = as.factor(geno_matrix[snp, ]),
        covariate = metadata[[covariate]]
      )
      
      # Remove missing values
      plot_data <- na.omit(plot_data)
      
      # Create plot
      p <- ggplot(plot_data, aes(x=genotype, y=expression, color=covariate)) +
        geom_boxplot() +
        geom_jitter(width=0.2, alpha=0.5) +
        labs(
          title=paste("Interaction between", snp, "and", gene),
          subtitle=paste("p-value =", signif(top_results$padj[i], 3)),
          x="Genotype",
          y="Expression",
          color=covariate
        ) +
        theme_minimal()
      
      # Save plot
      output_file <- file.path(output_dir, paste0("interaction_plot_", gene, "_", snp, "_", covariate, ".pdf"))
      ggsave(output_file, p, width=8, height=6)
    }
  }
}

# Plot top interactions
plot_top_interactions(results_list, expr_matrix, geno_matrix, metadata)

message("eQTL interaction analysis completed!")
