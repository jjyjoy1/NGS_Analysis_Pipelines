#!/usr/bin/env python
# Scanpy Analysis Script for 10x Genomics Visium Spatial Transcriptomics

import os
import sys
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import squidpy as sq
from datetime import datetime

# Set up logging
def setup_logger(log_file):
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# Get Snakemake parameters
matrices_dirs = snakemake.input.matrices
spatial_dirs = snakemake.input.spatial
output_h5ad = snakemake.output.h5ad
output_clusters_pdf = snakemake.output.clusters_pdf
output_markers_csv = snakemake.output.markers_csv
output_umap_pdf = snakemake.output.umap_pdf
log_file = snakemake.log[0]
n_threads = snakemake.threads

# Setup logger
logger = setup_logger(log_file)
logger.info(f"Starting Scanpy analysis for Visium data with {n_threads} threads")

# Set the number of threads for parallelization
sc.settings.n_jobs = n_threads

# Set scanpy settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 8))

# Function to process a single sample
def process_sample(matrix_dir, spatial_dir, sample_name):
    logger.info(f"Processing sample: {sample_name}")
    
    # Read the count matrix
    adata = sc.read_10x_h5(os.path.join(matrix_dir, "filtered_feature_bc_matrix.h5"))
    
    # Read spatial coordinates
    spatial_path = os.path.join(spatial_dir, "tissue_positions_list.csv")
    positions = pd.read_csv(spatial_path, header=None)
    positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    positions.index = positions['barcode']
    
    # Filter barcodes that are in the tissue
    positions_in_tissue = positions[positions['in_tissue'] == 1]
    adata = adata[adata.obs.index.isin(positions_in_tissue.index)]
    
    # Add spatial coordinates to AnnData
    adata.obsm['spatial'] = positions_in_tissue.loc[adata.obs.index, ['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
    
    # Add image to AnnData
    img_path = os.path.join(spatial_dir, "tissue_hires_image.png")
    adata.uns['spatial'] = {sample_name: {}}
    adata.uns['spatial'][sample_name]['images'] = {'hires': plt.imread(img_path)}
    
    # Read scale factors
    import json
    with open(os.path.join(spatial_dir, "scalefactors_json.json"), 'r') as f:
        scale_factors = json.load(f)
    
    adata.uns['spatial'][sample_name]['scalefactors'] = scale_factors
    
    # Add sample name to observations
    adata.obs['sample'] = sample_name
    
    return adata

# Process all samples
sample_names = [os.path.basename(os.path.dirname(os.path.dirname(d))) for d in matrices_dirs]
adatas = []

for i, sample_name in enumerate(sample_names):
    matrix_dir = matrices_dirs[i]
    spatial_dir = spatial_dirs[i]
    
    adata = process_sample(matrix_dir, spatial_dir, sample_name)
    adatas.append(adata)

# Concatenate all samples if more than one
if len(adatas) > 1:
    logger.info(f"Concatenating {len(adatas)} samples")
    adata = ad.concat(adatas, join='outer', label='sample', index_unique='-')
else:
    adata = adatas[0]

# Basic preprocessing
logger.info("Performing quality control")
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Filter cells
logger.info("Filtering cells")
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Filter genes
logger.info("Filtering genes")
sc.pp.filter_genes(adata, min_cells=3)

# Normalization
logger.info("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Feature selection
logger.info("Finding highly variable genes")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', batch_key='sample' if len(adatas) > 1 else None)

# Scale data
logger.info("Scaling data")
sc.pp.scale(adata, max_value=10)

# PCA
logger.info("Running PCA")
sc.tl.pca(adata, svd_solver='arpack')

# UMAP
logger.info("Running UMAP")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

# Clustering at different resolutions
logger.info("Finding clusters")
resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_res{res}')

# Default clustering for downstream analysis
adata.obs['clusters'] = adata.obs['leiden_res0.6']

# Find marker genes
logger.info("Finding marker genes")
sc.tl.rank_genes_groups(adata, 'clusters', method='wilcoxon')
marker_genes = sc.get.rank_genes_groups_df(adata, group=None)

# Save marker genes
logger.info("Saving marker genes")
marker_genes.to_csv(output_markers_csv, index=False)

# Save the AnnData object
logger.info("Saving AnnData object")
adata.write(output_h5ad)

# Generate visualizations
logger.info("Creating visualizations")

# UMAP plot
plt.figure(figsize=(10, 8))
sc.pl.umap(adata, color=['clusters'], show=False)
plt.savefig(output_umap_pdf, bbox_inches='tight')
plt.close()

# Spatial plots
with plt.rc_context({'figure.figsize': (12, 12)}):
    # Create PDF for spatial plots
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(output_clusters_pdf) as pdf:
        # Plot clusters for each sample
        for sample in adata.obs['sample'].unique():
            sample_adata = adata[adata.obs['sample'] == sample]
            
            # Cluster plot
            plt.figure(figsize=(10, 10))
            sc.pl.spatial(sample_adata, img_key="hires", color='clusters', 
                         title=f"Clusters - {sample}", show=False)
            pdf.savefig()
            plt.close()
            
            # Top marker genes for spatial visualization
            for group in adata.obs['clusters'].unique():
                markers = sc.get.rank_genes_groups_df(adata, group=group)
                if len(markers) > 0:
                    top_gene = markers.iloc[0]['names']
                    plt.figure(figsize=(10, 10))
                    sc.pl.spatial(sample_adata, img_key="hires", color=top_gene, 
                                 title=f"Top marker for cluster {group}: {top_gene} - {sample}", 
                                 show=False)
                    pdf.savefig()
                    plt.close()
        
        # Add UMAP visualization
        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color=['clusters'], show=False)
        pdf.savefig()
        plt.close()

# Run spatial neighborhood analysis with Squidpy
logger.info("Running spatial analysis with Squidpy")
try:
    # Calculate spatial graph
    sq.gr.spatial_neighbors(adata, coord_type='generic')
    
    # Calculate spatial autocorrelation (Moran's I)
    sq.gr.spatial_autocorr(
        adata,
        mode='moran',
        genes=adata.var_names[adata.var.highly_variable],
        n_perms=100,
        n_jobs=n_threads
    )
    
    # Identify spatially variable genes
    spatial_genes = sq.gr.spatial_neighbors(adata, coord_type='generic')
    
    # Add to the AnnData object
    adata.uns['spatial_variable_genes'] = spatial_genes
    
    # Update the saved object
    adata.write(output_h5ad)
except Exception as e:
    logger.error(f"Error in Squidpy analysis: {e}")

logger.info("Scanpy analysis completed successfully")
