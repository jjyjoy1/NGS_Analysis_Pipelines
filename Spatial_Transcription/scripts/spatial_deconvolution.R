#!/usr/bin/env Rscript

# Spatial Deconvolution Script for 10x Genomics Visium data
# This script uses SPOTlight for spatial deconvolution

# Load libraries
library(Seurat)
library(SPOTlight)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(scater)
library(scran)
library(BayesSpace)
library(viridis)

# Set up logging
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Starting spatial deconvolution analysis\n")

# Get input and output paths from Snakemake
seurat_object_path <- snakemake@input[["seurat"]]
output_results <- snakemake@output[["results"]]
output_plots <- snakemake@output[["plots"]]

# Load the Seurat object
cat("Loading Seurat object\n")
visium_data <- readRDS(seurat_object_path)

# Reference single-cell dataset
# For demonstration, we'll use a mock reference. In a real analysis,
# you would load your own single-cell reference data or use a public dataset
cat("Creating/loading reference single-cell data\n")

# Check if a reference scRNA-seq dataset is available
ref_path <- file.path(dirname(dirname(seurat_object_path)), "reference", "reference_scrna.rds")

if (file.exists(ref_path)) {
  cat("Loading existing reference scRNA-seq dataset\n")
  reference <- readRDS(ref_path)
} else {
  cat("No reference found. Using mock reference data\n")
  # Create a mock reference for demonstration
  # In a real analysis, you would replace this with your own reference data
  
  # For demonstration only: simulating reference data
  # This would be replaced by your actual reference scRNA-seq data
  set.seed(123)
  n_cells <- 5000
  n_genes <- 2000
  n_cell_types <- 8
  
  # Get gene names from the Visium data
  visium_genes <- rownames(visium_data)
  genes_to_use <- visium_genes[1:min(n_genes, length(visium_genes))]
  
  # Create expression matrix
  expr_matrix <- matrix(rpois(n_cells * n_genes, lambda = 0.5), nrow = n_genes)
  rownames(expr_matrix) <- genes_to_use
  colnames(expr_matrix) <- paste0("cell_", 1:n_cells)
  
  # Assign cell types
  cell_types <- sample(paste0("CellType", 1:n_cell_types), n_cells, replace = TRUE)
  names(cell_types) <- colnames(expr_matrix)
  
  # Simulate expression patterns for each cell type
  for (ct in unique(cell_types)) {
    ct_cells <- which(cell_types == ct)
    # Select random marker genes for this cell type
    marker_genes <- sample(1:n_genes, 50)
    # Increase expression for these genes in this cell type
    expr_matrix[marker_genes, ct_cells] <- expr_matrix[marker_genes, ct_cells] + 
                                           rpois(length(marker_genes) * length(ct_cells), lambda = 5)
  }
  
  # Create Seurat object
  reference <- CreateSeuratObject(counts = expr_matrix)
  reference$cell_type <- cell_types
  
  # Process reference data
  reference <- NormalizeData(reference)
  reference <- FindVariableFeatures(reference)
  reference <- ScaleData(reference)
  reference <- RunPCA(reference)
  reference <- FindNeighbors(reference)
  reference <- FindClusters(reference, resolution = 0.8)
  reference <- RunUMAP(reference, dims = 1:30)
  
  # Use assigned cell types
  Idents(reference) <- reference$cell_type
  
  cat("Mock reference created\n")
}

# Find marker genes in the reference data
cat("Finding marker genes in the reference data\n")
markers <- FindAllMarkers(reference, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)

# Prepare data for SPOTlight
cat("Preparing data for SPOTlight deconvolution\n")

# Check if reference genes match visium genes
genes_to_use <- intersect(rownames(reference), rownames(visium_data))
cat(paste0("Number of common genes between reference and spatial data: ", length(genes_to_use), "\n"))

if (length(genes_to_use) < 100) {
  stop("Too few common genes between reference and spatial data. Check gene naming.")
}

# Subset data to common genes
reference <- reference[genes_to_use, ]
visium_data <- visium_data[genes_to_use, ]

# Prepare reference expression matrix
ref_counts <- GetAssayData(reference, slot = "counts")
ref_metadata <- data.frame(cell_type = reference$cell_type)
rownames(ref_metadata) <- colnames(ref_counts)

# Prepare spatial expression matrix
spatial_counts <- GetAssayData(visium_data, assay = "Spatial", slot = "counts")

# Run SPOTlight deconvolution
cat("Running SPOTlight deconvolution\n")
set.seed(123)

# Train the model with reference data
spotlight_model <- SPOTlight(
  x = ref_counts,
  y = ref_metadata$cell_type,
  min_prop = 0.01,
  n_top_genes = 200
)

# Predict cell type proportions
cat("Predicting cell type proportions\n")
deconvolution_results <- SPOTlight::predictCellTypes(
  spotlight_model,
  spatial_counts
)

# Format results
deconv_df <- as.data.frame(deconvolution_results$mat)
rownames(deconv_df) <- colnames(spatial_counts)

# Add spot coordinates
spot_coords <- visium_data@images[[1]]@coordinates
deconv_df$row <- spot_coords[rownames(deconv_df), "row"]
deconv_df$col <- spot_coords[rownames(deconv_df), "col"]
deconv_df$imagerow <- spot_coords[rownames(deconv_df), "imagerow"]
deconv_df$imagecol <- spot_coords[rownames(deconv_df), "imagecol"]

# Save deconvolution results
cat("Saving deconvolution results\n")
write.csv(deconv_df, file = output_results, row.names = TRUE)

# Add deconvolution results to the Seurat object
cat("Adding deconvolution results to Seurat object\n")
for (cell_type in colnames(deconv_df)[1:(ncol(deconv_df)-4)]) {
  visium_data[[cell_type]] <- deconv_df[colnames(visium_data), cell_type]
}

# Generate visualizations
cat("Generating visualizations\n")
pdf(output_plots, width = 12, height = 12)

# Plot cell type proportions on tissue
for (cell_type in colnames(deconv_df)[1:(ncol(deconv_df)-4)]) {
  p <- SpatialFeaturePlot(
    visium_data,
    features = cell_type,
    pt.size.factor = 1.5,
    alpha = c(0.1, 1)
  ) + 
  scale_fill_viridis(option = "plasma") +
  ggtitle(paste("Proportion of", cell_type)) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right")
  
  print(p)
}

# Stacked barplot of proportions by cluster
cell_props_by_cluster <- data.frame()
clusters <- levels(visium_data$seurat_clusters)

for (cluster in clusters) {
  cells_in_cluster <- WhichCells(visium_data, idents = cluster)
  cluster_props <- colMeans(deconv_df[cells_in_cluster, 1:(ncol(deconv_df)-4), drop = FALSE])
  
  cluster_df <- data.frame(
    cluster = cluster,
    cell_type = names(cluster_props),
    proportion = as.numeric(cluster_props)
  )
  
  cell_props_by_cluster <- rbind(cell_props_by_cluster, cluster_df)
}

# Plot stacked barplot
p <- ggplot(cell_props_by_cluster, aes(x = cluster, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Cell Type Composition by Spatial Cluster",
       x = "Spatial Cluster",
       y = "Cell Type Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)

# Heatmap of cell type proportions by cluster
cell_props_matrix <- reshape2::dcast(cell_props_by_cluster, 
                                   cluster ~ cell_type, 
                                   value.var = "proportion")
rownames(cell_props_matrix) <- cell_props_matrix$cluster
cell_props_matrix$cluster <- NULL

# Plot heatmap
pheatmap::pheatmap(
  cell_props_matrix,
  color = viridis(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Cell Type Proportions by Cluster",
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)

# Close PDF
dev.off()

cat("Spatial deconvolution analysis completed\n")

# Close log file
sink(type = "message")
sink(type = "output")
close(log)
