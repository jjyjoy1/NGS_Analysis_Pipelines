# Snakemake workflow for Tumor Mutational Burden (TMB) analysis
# This pipeline handles both tumor-only and tumor-normal paired samples
# From FASTQ to annotation

import os
import glob
from snakemake.utils import min_version

# Set minimum Snakemake version
min_version("6.0")

# Configuration file
configfile: "config.yaml"

# Define sample information
SAMPLES = config["samples"]
TUMOR_SAMPLES = config["tumor_samples"]
NORMAL_SAMPLES = config.get("normal_samples", {})
PAIRED_ANALYSIS = config.get("paired_analysis", False)

# Define reference files
REFERENCE = config["reference"]["genome"]
KNOWN_VARIANTS = config["reference"]["known_variants"]
DBSNP = config["reference"]["dbsnp"]
CAPTURE_BED = config["reference"]["capture_bed"]
REFSEQ = config["reference"]["refseq"]

# Define output directories
RESULT_DIR = config["output_dir"]
FASTQC_DIR = os.path.join(RESULT_DIR, "fastqc")
ALIGN_DIR = os.path.join(RESULT_DIR, "aligned")
CALL_DIR = os.path.join(RESULT_DIR, "variants")
FILTER_DIR = os.path.join(RESULT_DIR, "filtered")
ANNOTATION_DIR = os.path.join(RESULT_DIR, "annotated")
TMB_DIR = os.path.join(RESULT_DIR, "tmb")

# Target rule
rule all:
    input:
        # FastQC reports
        expand(os.path.join(FASTQC_DIR, "{sample}_R{read}_fastqc.html"), 
               sample=SAMPLES.keys(), read=[1, 2]),
        
        # Aligned BAMs
        expand(os.path.join(ALIGN_DIR, "{sample}.recal.bam"), 
               sample=SAMPLES.keys()),
        
        # Variants
        expand(os.path.join(CALL_DIR, "{tumor}.vcf.gz"), 
               tumor=TUMOR_SAMPLES.keys()) if not PAIRED_ANALYSIS else [],
        
        # Somatic variants (for paired analysis)
        expand(os.path.join(CALL_DIR, "{tumor}_vs_{normal}.vcf.gz"), 
               zip, tumor=config.get("pairs", {}).keys(), 
               normal=config.get("pairs", {}).values()) if PAIRED_ANALYSIS else [],
        
        # Funcotator annotations
        expand(os.path.join(ANNOTATION_DIR, "{tumor}.funcotator.vcf.gz"), 
               tumor=TUMOR_SAMPLES.keys()) if not PAIRED_ANALYSIS else [],
        expand(os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.funcotator.vcf.gz"), 
               zip, tumor=config.get("pairs", {}).keys(), 
               normal=config.get("pairs", {}).values()) if PAIRED_ANALYSIS else [],
        
        # Oncotator annotations
        expand(os.path.join(ANNOTATION_DIR, "{tumor}.oncotator.maf.gz"), 
               tumor=TUMOR_SAMPLES.keys()) if not PAIRED_ANALYSIS else [],
        expand(os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.oncotator.maf.gz"), 
               zip, tumor=config.get("pairs", {}).keys(), 
               normal=config.get("pairs", {}).values()) if PAIRED_ANALYSIS else [],
        
        # Combined annotated variants
        expand(os.path.join(ANNOTATION_DIR, "{tumor}.annotated.vcf.gz"), 
               tumor=TUMOR_SAMPLES.keys()) if not PAIRED_ANALYSIS else [],
        expand(os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.annotated.vcf.gz"), 
               zip, tumor=config.get("pairs", {}).keys(), 
               normal=config.get("pairs", {}).values()) if PAIRED_ANALYSIS else [],
        
        # TMB reports
        expand(os.path.join(TMB_DIR, "{tumor}.tmb.txt"), 
               tumor=TUMOR_SAMPLES.keys()) if not PAIRED_ANALYSIS else [],
        expand(os.path.join(TMB_DIR, "{tumor}_vs_{normal}.tmb.txt"), 
               zip, tumor=config.get("pairs", {}).keys(), 
               normal=config.get("pairs", {}).values()) if PAIRED_ANALYSIS else []

# Quality control
rule fastqc:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["r1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["r2"]
    output:
        html_r1 = os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.html"),
        html_r2 = os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.html"),
        zip_r1 = os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.zip"),
        zip_r2 = os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.zip")
    log:
        os.path.join(FASTQC_DIR, "logs", "{sample}.fastqc.log")
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc --outdir {FASTQC_DIR} \
            --threads {threads} \
            {input.r1} {input.r2} &> {log}
        """

# Adapter trimming
rule trim_adapters:
    input:
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["r1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["r2"]
    output:
        r1 = temp(os.path.join(RESULT_DIR, "trimmed", "{sample}_R1_trimmed.fastq.gz")),
        r2 = temp(os.path.join(RESULT_DIR, "trimmed", "{sample}_R2_trimmed.fastq.gz")),
        html = os.path.join(RESULT_DIR, "trimmed", "{sample}.trimming_report.html"),
        json = os.path.join(RESULT_DIR, "trimmed", "{sample}.trimming_report.json")
    log:
        os.path.join(RESULT_DIR, "trimmed", "logs", "{sample}.trimmomatic.log")
    threads: 8
    conda:
        "envs/trim.yaml"
    params:
        adapters = config["adapters"],
        extra = config.get("trimmomatic_params", "")
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} /dev/null \
            {output.r2} /dev/null \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            {params.extra} \
            &> {log}
        """

# Alignment
rule bwa_mem:
    input:
        r1 = rules.trim_adapters.output.r1,
        r2 = rules.trim_adapters.output.r2,
        ref = REFERENCE
    output:
        temp(os.path.join(ALIGN_DIR, "{sample}.raw.bam"))
    log:
        os.path.join(ALIGN_DIR, "logs", "{sample}.bwa.log")
    threads: 16
    conda:
        "envs/align.yaml"
    params:
        read_group = lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"
    shell:
        """
        bwa mem -M -t {threads} \
            -R '{params.read_group}' \
            {input.ref} {input.r1} {input.r2} 2> {log} | \
            samtools sort -@ {threads} -o {output} -
        """

# Mark duplicates
rule mark_duplicates:
    input:
        rules.bwa_mem.output
    output:
        bam = temp(os.path.join(ALIGN_DIR, "{sample}.dedup.bam")),
        metrics = os.path.join(ALIGN_DIR, "{sample}.dedup.metrics.txt")
    log:
        os.path.join(ALIGN_DIR, "logs", "{sample}.markdup.log")
    conda:
        "envs/align.yaml"
    shell:
        """
        picard MarkDuplicates \
            INPUT={input} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            VALIDATION_STRINGENCY=LENIENT \
            &> {log}
        """

# Base quality recalibration
rule base_recalibration:
    input:
        bam = rules.mark_duplicates.output.bam,
        ref = REFERENCE,
        known = KNOWN_VARIANTS,
        dbsnp = DBSNP
    output:
        recal = temp(os.path.join(ALIGN_DIR, "{sample}.recal.table")),
        bam = os.path.join(ALIGN_DIR, "{sample}.recal.bam")
    log:
        os.path.join(ALIGN_DIR, "logs", "{sample}.bqsr.log")
    conda:
        "envs/align.yaml"
    shell:
        """
        # Generate recalibration table
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.known} \
            --known-sites {input.dbsnp} \
            -O {output.recal} &> {log}
        
        # Apply recalibration
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {output.recal} \
            -O {output.bam} &>> {log}
        """

# Index final BAM
rule index_bam:
    input:
        rules.base_recalibration.output.bam
    output:
        os.path.join(ALIGN_DIR, "{sample}.recal.bam.bai")
    log:
        os.path.join(ALIGN_DIR, "logs", "{sample}.index.log")
    conda:
        "envs/align.yaml"
    shell:
        """
        samtools index {input} &> {log}
        """

# Tumor-only variant calling with Mutect2
rule call_variants_tumor_only:
    input:
        bam = os.path.join(ALIGN_DIR, "{tumor}.recal.bam"),
        bai = os.path.join(ALIGN_DIR, "{tumor}.recal.bam.bai"),
        ref = REFERENCE,
        intervals = CAPTURE_BED
    output:
        vcf = os.path.join(CALL_DIR, "{tumor}.vcf.gz"),
        idx = os.path.join(CALL_DIR, "{tumor}.vcf.gz.tbi"),
        stats = os.path.join(CALL_DIR, "{tumor}.vcf.gz.stats")
    log:
        os.path.join(CALL_DIR, "logs", "{tumor}.mutect2.log")
    conda:
        "envs/variants.yaml"
    params:
        extra = "--max-mnp-distance 0"  # Disable MNP detection for TMB calculation
    shell:
        """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.intervals} \
            --germline-resource {DBSNP} \
            {params.extra} \
            -O {output.vcf} &> {log}
        """

# Somatic variant calling with Mutect2 (paired tumor-normal)
rule call_variants_paired:
    input:
        tumor_bam = os.path.join(ALIGN_DIR, "{tumor}.recal.bam"),
        tumor_bai = os.path.join(ALIGN_DIR, "{tumor}.recal.bam.bai"),
        normal_bam = os.path.join(ALIGN_DIR, "{normal}.recal.bam"),
        normal_bai = os.path.join(ALIGN_DIR, "{normal}.recal.bam.bai"),
        ref = REFERENCE,
        intervals = CAPTURE_BED
    output:
        vcf = os.path.join(CALL_DIR, "{tumor}_vs_{normal}.vcf.gz"),
        idx = os.path.join(CALL_DIR, "{tumor}_vs_{normal}.vcf.gz.tbi"),
        stats = os.path.join(CALL_DIR, "{tumor}_vs_{normal}.vcf.gz.stats")
    log:
        os.path.join(CALL_DIR, "logs", "{tumor}_vs_{normal}.mutect2.log")
    conda:
        "envs/variants.yaml"
    params:
        extra = "--max-mnp-distance 0"  # Disable MNP detection for TMB calculation
    shell:
        """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -normal {wildcards.normal} \
            -L {input.intervals} \
            --germline-resource {DBSNP} \
            {params.extra} \
            -O {output.vcf} &> {log}
        """

# Filter variants - tumor only
rule filter_variants_tumor_only:
    input:
        vcf = rules.call_variants_tumor_only.output.vcf,
        stats = rules.call_variants_tumor_only.output.stats,
        ref = REFERENCE
    output:
        vcf = os.path.join(FILTER_DIR, "{tumor}.filtered.vcf.gz"),
        idx = os.path.join(FILTER_DIR, "{tumor}.filtered.vcf.gz.tbi")
    log:
        os.path.join(FILTER_DIR, "logs", "{tumor}.filter.log")
    conda:
        "envs/variants.yaml"
    params:
        extra = "--stats {input.stats}"
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            {params.extra} \
            -O {output.vcf} &> {log}
        """

# Filter variants - paired
rule filter_variants_paired:
    input:
        vcf = rules.call_variants_paired.output.vcf,
        stats = rules.call_variants_paired.output.stats,
        ref = REFERENCE
    output:
        vcf = os.path.join(FILTER_DIR, "{tumor}_vs_{normal}.filtered.vcf.gz"),
        idx = os.path.join(FILTER_DIR, "{tumor}_vs_{normal}.filtered.vcf.gz.tbi")
    log:
        os.path.join(FILTER_DIR, "logs", "{tumor}_vs_{normal}.filter.log")
    conda:
        "envs/variants.yaml"
    params:
        extra = "--stats {input.stats}"
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            {params.extra} \
            -O {output.vcf} &> {log}
        """

# Annotation with Funcotator - tumor only
rule funcotator_tumor_only:
    input:
        vcf = rules.filter_variants_tumor_only.output.vcf,
        ref = REFERENCE,
        refseq = REFSEQ
    output:
        vcf = os.path.join(ANNOTATION_DIR, "{tumor}.funcotator.vcf.gz"),
        idx = os.path.join(ANNOTATION_DIR, "{tumor}.funcotator.vcf.gz.tbi")
    log:
        os.path.join(ANNOTATION_DIR, "logs", "{tumor}.funcotator.log")
    conda:
        "envs/annotation.yaml"
    params:
        data_sources = config["funcotator_data_sources"],
        output_format = "VCF"
    shell:
        """
        gatk Funcotator \
            --variant {input.vcf} \
            --reference {input.ref} \
            --ref-version hg38 \
            --data-sources-path {params.data_sources} \
            --output {output.vcf} \
            --output-file-format {params.output_format} &> {log}
        """

# Annotation with Funcotator - paired
rule funcotator_paired:
    input:
        vcf = rules.filter_variants_paired.output.vcf,
        ref = REFERENCE,
        refseq = REFSEQ
    output:
        vcf = os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.funcotator.vcf.gz"),
        idx = os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.funcotator.vcf.gz.tbi")
    log:
        os.path.join(ANNOTATION_DIR, "logs", "{tumor}_vs_{normal}.funcotator.log")
    conda:
        "envs/annotation.yaml"
    params:
        data_sources = config["funcotator_data_sources"],
        output_format = "VCF"
    shell:
        """
        gatk Funcotator \
            --variant {input.vcf} \
            --reference {input.ref} \
            --ref-version hg38 \
            --data-sources-path {params.data_sources} \
            --output {output.vcf} \
            --output-file-format {params.output_format} &> {log}
        """

# Annotation with Oncotator - tumor only
rule oncotator_tumor_only:
    input:
        vcf = rules.filter_variants_tumor_only.output.vcf,
        ref = REFERENCE
    output:
        maf = os.path.join(ANNOTATION_DIR, "{tumor}.oncotator.maf.gz")
    log:
        os.path.join(ANNOTATION_DIR, "logs", "{tumor}.oncotator.log")
    conda:
        "envs/oncotator.yaml"
    params:
        db_dir = config["oncotator_db_dir"],
        tx_mode = "TRANSCRIPT"
    shell:
        """
        # First, extract and convert VCF to intermediate format
        zcat {input.vcf} | grep -v "^#" | awk '{{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t.\t.\t.\tvariant\ttumor"}}' > {wildcards.tumor}.tsv
        
        # Run Oncotator
        oncotator -v --db-dir {params.db_dir} \
            --input_format MAFLITE \
            --output_format TCGAMAF \
            --tx-mode {params.tx_mode} \
            {wildcards.tumor}.tsv {output.maf} {input.ref} \
            &> {log}
            
        # Compress output MAF file
        gzip -f {output.maf}
        
        # Clean up intermediate file
        rm {wildcards.tumor}.tsv
        """

# Annotation with Oncotator - paired
rule oncotator_paired:
    input:
        vcf = rules.filter_variants_paired.output.vcf,
        ref = REFERENCE
    output:
        maf = os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.oncotator.maf.gz")
    log:
        os.path.join(ANNOTATION_DIR, "logs", "{tumor}_vs_{normal}.oncotator.log")
    conda:
        "envs/oncotator.yaml"
    params:
        db_dir = config["oncotator_db_dir"],
        tx_mode = "TRANSCRIPT"
    shell:
        """
        # First, extract and convert VCF to intermediate format
        zcat {input.vcf} | grep -v "^#" | awk '{{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t.\t.\t.\tvariant\t{wildcards.tumor}"}}' > {wildcards.tumor}_vs_{wildcards.normal}.tsv
        
        # Run Oncotator
        oncotator -v --db-dir {params.db_dir} \
            --input_format MAFLITE \
            --output_format TCGAMAF \
            --tx-mode {params.tx_mode} \
            {wildcards.tumor}_vs_{wildcards.normal}.tsv {output.maf} {input.ref} \
            &> {log}
            
        # Compress output MAF file
        gzip -f {output.maf}
        
        # Clean up intermediate file
        rm {wildcards.tumor}_vs_{wildcards.normal}.tsv
        """

# Combine annotations - tumor only
rule combine_annotations_tumor_only:
    input:
        funcotator = rules.funcotator_tumor_only.output.vcf,
        oncotator = rules.oncotator_tumor_only.output.maf
    output:
        vcf = os.path.join(ANNOTATION_DIR, "{tumor}.annotated.vcf.gz"),
        idx = os.path.join(ANNOTATION_DIR, "{tumor}.annotated.vcf.gz.tbi")
    log:
        os.path.join(ANNOTATION_DIR, "logs", "{tumor}.combine_annotations.log")
    conda:
        "envs/annotation.yaml"
    script:
        "scripts/combine_annotations.py"

# Combine annotations - paired
rule combine_annotations_paired:
    input:
        funcotator = rules.funcotator_paired.output.vcf,
        oncotator = rules.oncotator_paired.output.maf
    output:
        vcf = os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.annotated.vcf.gz"),
        idx = os.path.join(ANNOTATION_DIR, "{tumor}_vs_{normal}.annotated.vcf.gz.tbi")
    log:
        os.path.join(ANNOTATION_DIR, "logs", "{tumor}_vs_{normal}.combine_annotations.log")
    conda:
        "envs/annotation.yaml"
    script:
        "scripts/combine_annotations.py"

# Calculate TMB - tumor only
rule calculate_tmb_tumor_only:
    input:
        vcf = rules.annotate_variants_tumor_only.output.vcf,
        bed = CAPTURE_BED
    output:
        report = os.path.join(TMB_DIR, "{tumor}.tmb.txt")
    log:
        os.path.join(TMB_DIR, "logs", "{tumor}.tmb.log")
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/calculate_tmb.py"

# Calculate TMB - paired
rule calculate_tmb_paired:
    input:
        vcf = rules.annotate_variants_paired.output.vcf,
        bed = CAPTURE_BED
    output:
        report = os.path.join(TMB_DIR, "{tumor}_vs_{normal}.tmb.txt")
    log:
        os.path.join(TMB_DIR, "logs", "{tumor}_vs_{normal}.tmb.log")
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/calculate_tmb.py"

