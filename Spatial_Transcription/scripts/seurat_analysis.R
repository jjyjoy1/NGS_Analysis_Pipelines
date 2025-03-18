#!/usr/bin/env Rscript

# Seurat Analysis Script for 10x Genomics Visium Spatial Transcriptomics
# This script processes Space Ranger output using Seurat for spatial analysis

# Load required libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(future)
library(SpatialExperiment)
library(scater)
library(BiocParallel)

# Enable parallel processing
plan("multiprocess", workers = snakemake@threads)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# Get input and output paths from Snakemake
sample_dirs <- snakemake@input[["matrices"]]
spatial_dirs <- snakemake@input[["spatial"]]
output_rds <- snakemake@output[["rds"]]
output_clusters <- snakemake@output[["clusters_pdf"]]
output_markers <- snakemake@output[["markers_csv"]]
output_svg <- snakemake@output[["svg_genes_csv"]]
log_file <- snakemake@log[[1]]

# Setup logging
log <- file(log_file, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Starting Seurat analysis for Visium data\n")

# Function to process a single sample
process_sample <- function(data_dir, spatial_dir, sample_name) {
  cat(paste0("Processing sample: ", sample_name, "\n"))
  
  # Load the data
  visium_data <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample_name,
    filter.matrix = TRUE,
    image = NULL
  )
  
  # Add the image
  image_path <- file.path(spatial_dir, "tissue_hires_image.png")
  scale_factors_path <- file.path(spatial_dir, "scalefactors_json.json")
  visium_data <- AddImage(
    object = visium_data,
    image = image_path,
    scale.factors = jsonlite::read_json(scale_factors_path)
  )
  
  # Quality control
  visium_data[["percent.mt"]] <- PercentageFeatureSet(visium_data, pattern = "^MT-")
  cat("QC metrics calculated\n")
  
  # Filter cells based on QC metrics
  visium_data <- subset(
    visium_data,
    subset = nFeature_Spatial > 200 & nFeature_Spatial < 6000 & percent.mt < 20
  )
  
  # Normalize data
  visium_data <- SCTransform(visium_data, assay = "Spatial", verbose = FALSE)
  cat("Data normalized with SCTransform\n")
  
  return(visium_data)
}

# Process all samples and create a list of Seurat objects
sample_names <- basename(dirname(dirname(sample_dirs)))
seurat_objects <- list()

for (i in seq_along(sample_dirs)) {
  sample_name <- sample_names[i]
  data_dir <- sample_dirs[i]
  spatial_dir <- spatial_dirs[i]
  
  seurat_objects[[sample_name]] <- process_sample(data_dir, spatial_dir, sample_name)
}

# If multiple samples, integrate them
if (length(seurat_objects) > 1) {
  cat("Integrating multiple samples\n")
  
  # Select features for integration
  features <- SelectIntegrationFeatures(
    object.list = seurat_objects,
    nfeatures = 3000
  )
  
  # Prepare for integration
  seurat_objects <- PrepSCTIntegration(
    object.list = seurat_objects,
    anchor.features = features
  )
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(
    object.list = seurat_objects,
    normalization.method = "SCT",
    anchor.features = features
  )
  
  # Integrate data
  visium_integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT"
  )
  
  # Use integrated data for downstream analysis
  visium_data <- visium_integrated
  DefaultAssay(visium_data) <- "integrated"
  
} else {
  # If only one sample, use it directly
  visium_data <- seurat_objects[[1]]
  DefaultAssay(visium_data) <- "SCT"
}

# Dimensionality reduction
cat("Running PCA\n")
visium_data <- RunPCA(visium_data, assay = "SCT", verbose = FALSE)

# Determine dimensionality
pc_dims <- min(30, ncol(visium_data))
cat(paste0("Using ", pc_dims, " PCs for downstream analysis\n"))

# UMAP for visualization
cat("Running UMAP\n")
visium_data <- RunUMAP(visium_data, dims = 1:pc_dims, verbose = FALSE)

# Clustering
cat("Finding clusters\n")
visium_data <- FindNeighbors(visium_data, dims = 1:pc_dims, verbose = FALSE)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (res in resolutions) {
  visium_data <- FindClusters(visium_data, resolution = res, verbose = FALSE)
}

# Find spatially variable features
cat("Finding spatially variable features\n")
visium_data <- FindSpatiallyVariableFeatures(
  visium_data,
  assay = "SCT",
  features = VariableFeatures(visium_data),
  selection.method = "markvariogram"
)

top_spatial_features <- SpatiallyVariableFeatures(
  visium_data,
  selection.method = "markvariogram",
  nfeatures = 100
)

# Find marker genes for each cluster
cat("Finding marker genes\n")
Idents(visium_data) <- "seurat_clusters"
all_markers <- FindAllMarkers(
  visium_data,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save outputs
cat("Saving outputs\n")

# Save markers
write.csv(all_markers, file = output_markers, row.names = FALSE)

# Save spatially variable genes
svg_results <- as.data.frame(SpatiallyVariableFeatures(
  visium_data,
  selection.method = "markvariogram",
  nfeatures = 1000
))
colnames(svg_results) <- "gene"
write.csv(svg_results, file = output_svg, row.names = FALSE)

# Create spatial visualizations
cat("Creating visualizations\n")
pdf(output_clusters, width = 12, height = 10)

# Plot spatial clusters for each original sample
for (sample_name in names(seurat_objects)) {
  if (length(seurat_objects) > 1) {
    # For integrated dataset, extract cells from this sample
    sample_cells <- WhichCells(visium_data, idents = NULL, cells = colnames(seurat_objects[[sample_name]]))
    p1 <- SpatialDimPlot(visium_data, cells = sample_cells, 
                        group.by = "seurat_clusters", 
                        images = sample_name, 
                        pt.size.factor = 1.6) + 
          ggtitle(paste0("Clusters - ", sample_name))
  } else {
    p1 <- SpatialDimPlot(visium_data, 
                        group.by = "seurat_clusters", 
                        pt.size.factor = 1.6) + 
          ggtitle("Spatial Clusters")
  }
  
  print(p1)
  
  # Plot top spatially variable genes
  for (i in 1:min(5, length(top_spatial_features))) {
    gene <- top_spatial_features[i]
    if (length(seurat_objects) > 1) {
      p2 <- SpatialFeaturePlot(visium_data, 
                              features = gene, 
                              images = sample_name, 
                              cells = sample_cells) + 
            ggtitle(paste0(gene, " - ", sample_name))
    } else {
      p2 <- SpatialFeaturePlot(visium_data, 
                              features = gene) + 
            ggtitle(gene)
    }
    print(p2)
  }
}

# Add UMAP plot
p3 <- DimPlot(visium_data, reduction = "umap", group.by = "seurat_clusters")
print(p3)

dev.off()

# Save the Seurat object
cat("Saving Seurat object\n")
saveRDS(visium_data, file = output_rds)

cat("Seurat analysis completed successfully\n")

# Close log file
sink(type = "message")
sink(type = "output")
close(log)
