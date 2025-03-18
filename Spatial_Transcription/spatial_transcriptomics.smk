# Snakefile for 10x Genomics Visium Spatial Transcriptomics Analysis
# This pipeline processes data from raw FASTQ files to final visualization and analysis

import os
from os.path import join
import glob
import pandas as pd
from snakemake.utils import validate, min_version

# Set minimum Snakemake version
min_version("7.0.0")

# Load configuration
configfile: "config.yaml"

# Define sample information
SAMPLES = config["samples"]
REFERENCE = config["reference_genome"]
TRANSCRIPTOME = config["reference_transcriptome"]
RESULTS_DIR = config["results_dir"]
SPACERANGER_PATH = config["spaceranger_path"]
TISSUE_POSITIONS = config["tissue_positions_file"]
SLIDE_FILE = config["slide_file"]
IMAGE_FILE = config["image_file"]

# Create directories
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(join(RESULTS_DIR, "qc"), exist_ok=True)
os.makedirs(join(RESULTS_DIR, "spaceranger"), exist_ok=True)
os.makedirs(join(RESULTS_DIR, "seurat"), exist_ok=True)
os.makedirs(join(RESULTS_DIR, "scanpy"), exist_ok=True)
os.makedirs(join(RESULTS_DIR, "final_results"), exist_ok=True)

# Target rule to define the desired outputs
rule all:
    input:
        # QC reports
        expand(join(RESULTS_DIR, "qc", "{sample}_fastqc.html"), sample=SAMPLES),
        # Space Ranger outputs
        expand(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "web_summary.html"), sample=SAMPLES),
        expand(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "filtered_feature_bc_matrix"), sample=SAMPLES),
        # Seurat analysis
        join(RESULTS_DIR, "seurat", "seurat_analysis.rds"),
        # Scanpy analysis
        join(RESULTS_DIR, "scanpy", "scanpy_analysis.h5ad"),
        # Final visualizations and results
        join(RESULTS_DIR, "final_results", "spatial_clusters.pdf"),
        join(RESULTS_DIR, "final_results", "spatially_variable_genes.csv"),
        join(RESULTS_DIR, "final_results", "deconvolution_results.csv"),
        join(RESULTS_DIR, "final_results", "pathway_analysis.csv"),
        join(RESULTS_DIR, "final_results", "final_report.html")

# QC of raw FASTQ files
rule fastqc:
    input:
        lambda wildcards: config["samples"][wildcards.sample]["fastq_path"]
    output:
        html=join(RESULTS_DIR, "qc", "{sample}_fastqc.html"),
        zip=join(RESULTS_DIR, "qc", "{sample}_fastqc.zip")
    log:
        join(RESULTS_DIR, "logs", "fastqc_{sample}.log")
    threads: 4
    shell:
        """
        fastqc --threads {threads} --outdir $(dirname {output.html}) {input} &> {log}
        """

# MultiQC report for all samples
rule multiqc:
    input:
        expand(join(RESULTS_DIR, "qc", "{sample}_fastqc.zip"), sample=SAMPLES)
    output:
        join(RESULTS_DIR, "qc", "multiqc_report.html")
    log:
        join(RESULTS_DIR, "logs", "multiqc.log")
    shell:
        """
        multiqc {input} -o $(dirname {output}) &> {log}
        """

# Run Space Ranger for each sample
rule spaceranger_count:
    input:
        fastq=lambda wildcards: config["samples"][wildcards.sample]["fastq_path"],
        image=lambda wildcards: config["samples"][wildcards.sample]["image_path"],
        tissue=lambda wildcards: config["samples"][wildcards.sample]["tissue_positions_path"]
    output:
        summary=join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "web_summary.html"),
        matrix_dir=directory(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "filtered_feature_bc_matrix")),
        spatial_dir=directory(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "spatial"))
    params:
        sample_id="{sample}",
        transcriptome=TRANSCRIPTOME,
        slide=SLIDE_FILE,
        area=lambda wildcards: config["samples"][wildcards.sample]["tissue_area"]
    log:
        join(RESULTS_DIR, "logs", "spaceranger_{sample}.log")
    threads: 16
    resources:
        mem_mb=64000,
        time="24:00:00"
    shell:
        """
        {SPACERANGER_PATH} count \
            --id={params.sample_id} \
            --transcriptome={params.transcriptome} \
            --fastqs={input.fastq} \
            --image={input.image} \
            --slide={params.slide} \
            --area={params.area} \
            --localcores={threads} \
            --localmem=64 \
            --tissue-positions-file={input.tissue} \
            --output-dir={RESULTS_DIR}/spaceranger &> {log}
        """

# Seurat analysis (in R)
rule seurat_analysis:
    input:
        matrices=expand(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "filtered_feature_bc_matrix"), sample=SAMPLES),
        spatial=expand(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "spatial"), sample=SAMPLES)
    output:
        rds=join(RESULTS_DIR, "seurat", "seurat_analysis.rds"),
        clusters_pdf=join(RESULTS_DIR, "seurat", "spatial_clusters.pdf"),
        markers_csv=join(RESULTS_DIR, "seurat", "spatial_markers.csv"),
        svg_genes_csv=join(RESULTS_DIR, "seurat", "spatially_variable_genes.csv")
    log:
        join(RESULTS_DIR, "logs", "seurat_analysis.log")
    threads: 8
    resources:
        mem_mb=32000
    script:
        "scripts/seurat_analysis.R"

# Scanpy analysis (in Python)
rule scanpy_analysis:
    input:
        matrices=expand(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "filtered_feature_bc_matrix"), sample=SAMPLES),
        spatial=expand(join(RESULTS_DIR, "spaceranger", "{sample}", "outs", "spatial"), sample=SAMPLES)
    output:
        h5ad=join(RESULTS_DIR, "scanpy", "scanpy_analysis.h5ad"),
        clusters_pdf=join(RESULTS_DIR, "scanpy", "spatial_clusters.pdf"),
        markers_csv=join(RESULTS_DIR, "scanpy", "spatial_markers.csv"),
        umap_pdf=join(RESULTS_DIR, "scanpy", "umap_embedding.pdf")
    log:
        join(RESULTS_DIR, "logs", "scanpy_analysis.log")
    threads: 8
    resources:
        mem_mb=32000
    script:
        "scripts/scanpy_analysis.py"

# Spatial deconvolution analysis
rule cell_type_deconvolution:
    input:
        seurat=join(RESULTS_DIR, "seurat", "seurat_analysis.rds")
    output:
        results=join(RESULTS_DIR, "final_results", "deconvolution_results.csv"),
        plots=join(RESULTS_DIR, "final_results", "deconvolution_plots.pdf")
    log:
        join(RESULTS_DIR, "logs", "deconvolution.log")
    script:
        "scripts/spatial_deconvolution.R"

# Pathway analysis
rule pathway_analysis:
    input:
        markers=join(RESULTS_DIR, "seurat", "spatial_markers.csv")
    output:
        pathways=join(RESULTS_DIR, "final_results", "pathway_analysis.csv"),
        plots=join(RESULTS_DIR, "final_results", "pathway_plots.pdf")
    log:
        join(RESULTS_DIR, "logs", "pathway_analysis.log")
    script:
        "scripts/pathway_analysis.R"

# Generate final report
rule final_report:
    input:
        seurat_rds=join(RESULTS_DIR, "seurat", "seurat_analysis.rds"),
        scanpy_h5ad=join(RESULTS_DIR, "scanpy", "scanpy_analysis.h5ad"),
        deconvolution=join(RESULTS_DIR, "final_results", "deconvolution_results.csv"),
        pathways=join(RESULTS_DIR, "final_results", "pathway_analysis.csv"),
        spatial_markers=join(RESULTS_DIR, "seurat", "spatial_markers.csv")
    output:
        report=join(RESULTS_DIR, "final_results", "final_report.html"),
        clusters=join(RESULTS_DIR, "final_results", "spatial_clusters.pdf"),
        svg=join(RESULTS_DIR, "final_results", "spatially_variable_genes.csv")
    log:
        join(RESULTS_DIR, "logs", "final_report.log")
    script:
        "scripts/generate_report.Rmd"
