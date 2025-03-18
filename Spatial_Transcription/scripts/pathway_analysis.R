#!/usr/bin/env Rscript

# Pathway Analysis Script for 10x Genomics Visium Spatial Transcriptomics Data
# This script performs gene set enrichment analysis on spatially variable genes
# and marker genes identified in the spatial clusters

# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)  # Change to appropriate organism if not human
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(RColorBrewer)
library(DOSE)

# Set up logging
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Get input and output paths from Snakemake
markers_file <- snakemake@input[["markers"]]
output_pathways <- snakemake@output[["pathways"]]
output_plots <- snakemake@output[["plots"]]

cat("Starting pathway analysis\n")

# Load marker genes
cat("Loading marker genes\n")
markers <- read.csv(markers_file, stringsAsFactors = FALSE)

# Check if organism is human (default) - adjust as needed for your organism
organism <- "human"
if (exists("snakemake@params[['organism']]")) {
  organism <- snakemake@params[["organism"]]
}

# Set the correct organism database
if (organism == "human") {
  org_db <- org.Hs.eg.db
} else if (organism == "mouse") {
  org_db <- org.Mm.eg.db
} else if (organism == "rat") {
  org_db <- org.Rn.eg.db
} else {
  stop("Unsupported organism. Please use human, mouse, or rat.")
}

# Function to convert gene symbols to Entrez IDs
symbol_to_entrez <- function(gene_symbols) {
  cat(paste0("Converting ", length(gene_symbols), " gene symbols to Entrez IDs\n"))
  entrez_ids <- bitr(gene_symbols, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org_db)
  return(entrez_ids)
}

# Function to run GO enrichment analysis
run_go_enrichment <- function(gene_list, gene_universe = NULL, ont = "BP", 
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  cat(paste0("Running GO enrichment analysis (", ont, ") for ", length(gene_list), " genes\n"))
  
  # Convert gene symbols to Entrez IDs
  genes_entrez <- symbol_to_entrez(gene_list)
  
  # Check if we have a gene universe
  if (!is.null(gene_universe)) {
    universe_entrez <- symbol_to_entrez(gene_universe)
    universe <- universe_entrez$ENTREZID
  } else {
    universe <- NULL
  }
  
  # Run enrichment analysis
  ego <- enrichGO(gene = genes_entrez$ENTREZID,
                  universe = universe,
                  OrgDb = org_db,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = pvalueCutoff,
                  qvalueCutoff = qvalueCutoff,
                  readable = TRUE)
  
  return(ego)
}

# Function to run KEGG pathway enrichment analysis
run_kegg_enrichment <- function(gene_list, gene_universe = NULL, 
                               pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  cat(paste0("Running KEGG pathway enrichment analysis for ", length(gene_list), " genes\n"))
  
  # Convert gene symbols to Entrez IDs
  genes_entrez <- symbol_to_entrez(gene_list)
  
  # Check if we have a gene universe
  if (!is.null(gene_universe)) {
    universe_entrez <- symbol_to_entrez(gene_universe)
    universe <- universe_entrez$ENTREZID
  } else {
    universe <- NULL
  }
  
  # Run KEGG enrichment
  kk <- enrichKEGG(gene = genes_entrez$ENTREZID,
                  universe = universe,
                  organism = ifelse(organism == "human", "hsa", 
                                   ifelse(organism == "mouse", "mmu", "rno")),
                  pAdjustMethod = "BH",
                  pvalueCutoff = pvalueCutoff,
                  qvalueCutoff = qvalueCutoff)
  
  return(kk)
}

# Function to run Reactome pathway analysis
run_reactome_enrichment <- function(gene_list, gene_universe = NULL, 
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  cat(paste0("Running Reactome pathway enrichment analysis for ", length(gene_list), " genes\n"))
  
  # Convert gene symbols to Entrez IDs
  genes_entrez <- symbol_to_entrez(gene_list)
  
  # Check if we have a gene universe
  if (!is.null(gene_universe)) {
    universe_entrez <- symbol_to_entrez(gene_universe)
    universe <- universe_entrez$ENTREZID
  } else {
    universe <- NULL
  }
  
  # Run Reactome enrichment
  x <- enrichPathway(gene = genes_entrez$ENTREZID,
                    universe = universe,
                    organism = ifelse(organism == "human", "human", 
                                     ifelse(organism == "mouse", "mouse", "rat")),
                    pAdjustMethod = "BH",
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    readable = TRUE)
  
  return(x)
}

# Process marker genes by cluster
cat("Processing marker genes by cluster\n")
clusters <- unique(markers$cluster)
all_results <- list()

# Create a PDF file for plots
pdf(output_plots, width = 12, height = 10)

# Process each cluster
for (cluster in clusters) {
  cat(paste0("Processing cluster ", cluster, "\n"))
  
  # Get marker genes for this cluster
  cluster_markers <- markers %>%
    filter(cluster == !!cluster & p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))
  
  # Skip if no significant markers
  if (nrow(cluster_markers) == 0) {
    cat(paste0("No significant markers for cluster ", cluster, "\n"))
    next
  }
  
  # Get top markers (top 100 or all if less)
  top_genes <- cluster_markers$gene[1:min(100, nrow(cluster_markers))]
  
  # Run GO Biological Process enrichment
  go_bp <- run_go_enrichment(top_genes, ont = "BP")
  
  # Run GO Molecular Function enrichment
  go_mf <- run_go_enrichment(top_genes, ont = "MF")
  
  # Run GO Cellular Component enrichment
  go_cc <- run_go_enrichment(top_genes, ont = "CC")
  
  # Run KEGG pathway enrichment
  kegg <- try(run_kegg_enrichment(top_genes), silent = TRUE)
  
  # Run Reactome pathway enrichment
  reactome <- try(run_reactome_enrichment(top_genes), silent = TRUE)
  
  # Store results
  cluster_results <- list(
    cluster = cluster,
    go_bp = go_bp,
    go_mf = go_mf,
    go_cc = go_cc,
    kegg = if (class(kegg) != "try-error") kegg else NULL,
    reactome = if (class(reactome) != "try-error") reactome else NULL
  )
  
  all_results[[paste0("cluster_", cluster)]] <- cluster_results
  
  # Generate plots
  
  # Plot title
  plot_title <- paste0("Pathway Enrichment for Cluster ", cluster)
  
  # 1. GO Biological Process dotplot
  if (nrow(as.data.frame(go_bp)) > 0) {
    p1 <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process") + 
      theme(axis.text.y = element_text(size = 8))
    print(p1)
  }
  
  # 2. GO Molecular Function dotplot
  if (nrow(as.data.frame(go_mf)) > 0) {
    p2 <- dotplot(go_mf, showCategory = 15, title = "GO Molecular Function") + 
      theme(axis.text.y = element_text(size = 8))
    print(p2)
  }
  
  # 3. KEGG pathway dotplot
  if (!is.null(cluster_results$kegg) && nrow(as.data.frame(cluster_results$kegg)) > 0) {
    p3 <- dotplot(cluster_results$kegg, showCategory = 15, title = "KEGG Pathways") + 
      theme(axis.text.y = element_text(size = 8))
    print(p3)
  }
  
  # 4. Reactome pathway dotplot
  if (!is.null(cluster_results$reactome) && nrow(as.data.frame(cluster_results$reactome)) > 0) {
    p4 <- dotplot(cluster_results$reactome, showCategory = 15, title = "Reactome Pathways") + 
      theme(axis.text.y = element_text(size = 8))
    print(p4)
  }
  
  # 5. GO Biological Process network
  if (nrow(as.data.frame(go_bp)) > 5) {
    p5 <- emapplot(pairwise_termsim(go_bp), showCategory = 20, 
                   cex_label_category = 0.6)
    print(p5)
  }
  
  # 6. GSEA plot for top pathways
  if (nrow(as.data.frame(go_bp)) > 0) {
    for (i in 1:min(5, nrow(as.data.frame(go_bp)))) {
      tryCatch({
        category <- go_bp$Description[i]
        p6 <- cnetplot(go_bp, categorySize = "pvalue", 
                      showCategory = category,
                      cex_gene = 0.5, cex_category = 0.8)
        print(p6)
      }, error = function(e) {
        cat(paste0("Error generating cnetplot for category ", i, ": ", e$message, "\n"))
      })
    }
  }
}

# Generate combined pathway results table
cat("Generating combined pathway results table\n")
all_pathways <- data.frame()

for (cluster_name in names(all_results)) {
  cluster_data <- all_results[[cluster_name]]
  cluster_id <- cluster_data$cluster
  
  # Process GO BP results
  if (!is.null(cluster_data$go_bp) && nrow(as.data.frame(cluster_data$go_bp)) > 0) {
    go_bp_df <- as.data.frame(cluster_data$go_bp)
    go_bp_df$cluster <- cluster_id
    go_bp_df$category <- "GO_BP"
    go_bp_df <- go_bp_df %>% 
      select(cluster, category, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    all_pathways <- rbind(all_pathways, go_bp_df)
  }
  
  # Process GO MF results
  if (!is.null(cluster_data$go_mf) && nrow(as.data.frame(cluster_data$go_mf)) > 0) {
    go_mf_df <- as.data.frame(cluster_data$go_mf)
    go_mf_df$cluster <- cluster_id
    go_mf_df$category <- "GO_MF"
    go_mf_df <- go_mf_df %>% 
      select(cluster, category, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    all_pathways <- rbind(all_pathways, go_mf_df)
  }
  
  # Process GO CC results
  if (!is.null(cluster_data$go_cc) && nrow(as.data.frame(cluster_data$go_cc)) > 0) {
    go_cc_df <- as.data.frame(cluster_data$go_cc)
    go_cc_df$cluster <- cluster_id
    go_cc_df$category <- "GO_CC"
    go_cc_df <- go_cc_df %>% 
      select(cluster, category, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    all_pathways <- rbind(all_pathways, go_cc_df)
  }
  
  # Process KEGG results
  if (!is.null(cluster_data$kegg) && nrow(as.data.frame(cluster_data$kegg)) > 0) {
    kegg_df <- as.data.frame(cluster_data$kegg)
    kegg_df$cluster <- cluster_id
    kegg_df$category <- "KEGG"
    kegg_df <- kegg_df %>% 
      select(cluster, category, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    all_pathways <- rbind(all_pathways, kegg_df)
  }
  
  # Process Reactome results
  if (!is.null(cluster_data$reactome) && nrow(as.data.frame(cluster_data$reactome)) > 0) {
    reactome_df <- as.data.frame(cluster_data$reactome)
    reactome_df$cluster <- cluster_id
    reactome_df$category <- "Reactome"
    reactome_df <- reactome_df %>% 
      select(cluster, category, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    all_pathways <- rbind(all_pathways, reactome_df)
  }
}

# Plot pathways heatmap across clusters
if (nrow(all_pathways) > 0) {
  # Get top pathways by adjusted p-value
  top_pathways <- all_pathways %>%
    group_by(category, Description) %>%
    summarize(min_padj = min(p.adjust), .groups = "drop") %>%
    arrange(min_padj) %>%
    head(30)
  
  # Filter to these top pathways
  plot_data <- all_pathways %>%
    filter(paste0(category, "_", Description) %in% 
             paste0(top_pathways$category, "_", top_pathways$Description))
  
  # Create matrix of -log10(p.adjust) values
  plot_matrix <- plot_data %>%
    mutate(neg_log_padj = -log10(p.adjust),
           pathway = paste0(category, ": ", Description)) %>%
    select(cluster, pathway, neg_log_padj) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = neg_log_padj, values_fill = 0)
  
  # Convert to matrix
  pathway_matrix <- as.matrix(plot_matrix[, -1])
  rownames(pathway_matrix) <- plot_matrix$pathway
  
  # Plot heatmap
  pheatmap::pheatmap(
    pathway_matrix,
    color = colorRampPalette(c("white", "red"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Top Pathways Across Clusters (-log10 adjusted p-value)",
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45
  )
}

# Close PDF
dev.off()

# Save all pathway results
cat("Saving pathway results\n")
write.csv(all_pathways, file = output_pathways, row.names = FALSE)

cat("Pathway analysis completed\n")

# Close log file
sink(type = "message")
sink(type = "output")
close(log)
