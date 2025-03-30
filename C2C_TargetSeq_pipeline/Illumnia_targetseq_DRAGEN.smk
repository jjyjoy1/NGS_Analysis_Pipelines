# Snakemake pipeline for Illumina DRAGEN bioinformatics analysis

import os
import glob
from os.path import join

# Configuration
configfile: "config.yaml"

# Define sample information from config
SAMPLES = config["samples"]
REFERENCE = config["reference"]["genome"]
OUTPUT_DIR = config["output_dir"]
DRAGEN_PATH = config["dragen_path"]

# Define output directories
FASTQ_DIR = join(OUTPUT_DIR, "fastq")
ALIGN_DIR = join(OUTPUT_DIR, "aligned")
VARIANT_DIR = join(OUTPUT_DIR, "variants")
QC_DIR = join(OUTPUT_DIR, "qc")

# Final outputs to generate
rule all:
    input:
        # Aligned BAM files
        expand(join(ALIGN_DIR, "{sample}.bam"), sample=SAMPLES),
        # Variant calls
        expand(join(VARIANT_DIR, "{sample}.vcf.gz"), sample=SAMPLES),
        # QC reports
        expand(join(QC_DIR, "{sample}_report.html"), sample=SAMPLES),
        # MultiQC report
        join(QC_DIR, "multiqc_report.html")

# C2C data preprocessing
rule preprocess_c2c:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        r1 = join(FASTQ_DIR, "{sample}_preprocessed_R1.fastq.gz"),
        r2 = join(FASTQ_DIR, "{sample}_preprocessed_R2.fastq.gz")
    log:
        join(OUTPUT_DIR, "logs", "preprocess", "{sample}.log")
    threads: 8
    shell:
        """
        # DRAGEN command for C2C preprocessing
        {DRAGEN_PATH}/dragen \
            --c2c-preprocess \
            --input-r1 {input.r1} \
            --input-r2 {input.r2} \
            --output-r1 {output.r1} \
            --output-r2 {output.r2} \
            --threads {threads} \
            2> {log}
        """

# DRAGEN alignment
rule dragen_align:
    input:
        r1 = join(FASTQ_DIR, "{sample}_preprocessed_R1.fastq.gz"),
        r2 = join(FASTQ_DIR, "{sample}_preprocessed_R2.fastq.gz")
    output:
        bam = join(ALIGN_DIR, "{sample}.bam"),
        bai = join(ALIGN_DIR, "{sample}.bam.bai")
    params:
        sample_name = "{sample}",
        output_prefix = join(ALIGN_DIR, "{sample}")
    log:
        join(OUTPUT_DIR, "logs", "align", "{sample}.log")
    threads: 16
    shell:
        """
        # DRAGEN alignment command
        {DRAGEN_PATH}/dragen \
            --ref-dir {REFERENCE} \
            --fastq-file1 {input.r1} \
            --fastq-file2 {input.r2} \
            --output-directory $(dirname {params.output_prefix}) \
            --output-file-prefix $(basename {params.output_prefix}) \
            --RGID {params.sample_name} \
            --RGSM {params.sample_name} \
            --enable-map-align true \
            --enable-sort true \
            --enable-map-align-output true \
            --enable-bam-indexing true \
            --c2c-processing true \
            --threads {threads} \
            2> {log}
        """

# DRAGEN variant calling
rule dragen_variant_call:
    input:
        bam = join(ALIGN_DIR, "{sample}.bam"),
        bai = join(ALIGN_DIR, "{sample}.bam.bai")
    output:
        vcf = join(VARIANT_DIR, "{sample}.vcf.gz"),
        vcf_idx = join(VARIANT_DIR, "{sample}.vcf.gz.tbi")
    params:
        sample_name = "{sample}",
        output_prefix = join(VARIANT_DIR, "{sample}")
    log:
        join(OUTPUT_DIR, "logs", "variant", "{sample}.log")
    threads: 16
    shell:
        """
        # DRAGEN variant calling command
        {DRAGEN_PATH}/dragen \
            --ref-dir {REFERENCE} \
            --bam-input {input.bam} \
            --output-directory $(dirname {params.output_prefix}) \
            --output-file-prefix $(basename {params.output_prefix}) \
            --RGID {params.sample_name} \
            --RGSM {params.sample_name} \
            --enable-variant-caller true \
            --enable-sv true \
            --sv-detection-mode normal \
            --vc-emit-ref-confidence GVCF \
            --c2c-variant-caller true \
            --threads {threads} \
            2> {log}
        """

# QC for individual samples
rule dragen_qc:
    input:
        bam = join(ALIGN_DIR, "{sample}.bam"),
        vcf = join(VARIANT_DIR, "{sample}.vcf.gz")
    output:
        report = join(QC_DIR, "{sample}_report.html")
    params:
        output_prefix = join(QC_DIR, "{sample}")
    log:
        join(OUTPUT_DIR, "logs", "qc", "{sample}.log")
    shell:
        """
        # DRAGEN QC command
        {DRAGEN_PATH}/dragen_qc \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --output-prefix {params.output_prefix} \
            2> {log}
        """

# MultiQC for aggregated project reports
rule multiqc:
    input:
        reports = expand(join(QC_DIR, "{sample}_report.html"), sample=SAMPLES)
    output:
        report = join(QC_DIR, "multiqc_report.html")
    log:
        join(OUTPUT_DIR, "logs", "multiqc.log")
    shell:
        """
        # Run MultiQC to aggregate QC reports
        multiqc {QC_DIR} -o {QC_DIR} -f -n multiqc_report.html 2> {log}
        """


