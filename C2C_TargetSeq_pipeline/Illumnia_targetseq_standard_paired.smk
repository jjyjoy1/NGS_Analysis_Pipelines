# Snakemake pipeline for Illumina C2C sequencing using standard bioinformatics tools
# File: Snakefile

import os
from os.path import join

# Configuration
configfile: "config_standard.yaml"

# Define sample information from config
SAMPLES = config["samples"]
REFERENCE = config["reference"]["genome"]
THREADS = config["threads"]
OUTPUT_DIR = config["output_dir"]

# Define output directories
FASTQ_DIR = join(OUTPUT_DIR, "fastq")
ALIGN_DIR = join(OUTPUT_DIR, "aligned")
VARIANT_DIR = join(OUTPUT_DIR, "variants")
QC_DIR = join(OUTPUT_DIR, "qc")
LOG_DIR = join(OUTPUT_DIR, "logs")

# Create directories if they don't exist
for dir_path in [FASTQ_DIR, ALIGN_DIR, VARIANT_DIR, QC_DIR, LOG_DIR]:
    os.makedirs(dir_path, exist_ok=True)

# Get tumor samples
def get_tumor_samples():
    return [s for s, v in config["samples"].items() if v.get("type") == "tumor"]

# Get tumor-normal pairs
def get_pairs():
    pairs = []
    for sample, values in config["samples"].items():
        if values.get("type") == "tumor" and "matched_normal" in values:
            pairs.append((sample, values["matched_normal"]))
    return pairs

# Final outputs to generate
rule all:
    input:
        # Aligned BAM files
        expand(join(ALIGN_DIR, "{sample}.sorted.bam"), sample=SAMPLES),
        expand(join(ALIGN_DIR, "{sample}.sorted.bam.bai"), sample=SAMPLES),
        # Variant calls
        expand(join(VARIANT_DIR, "{sample}.snps.vcf.gz"), sample=SAMPLES),
        # Paired somatic variant calls
        expand(join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.vcf.gz"), 
               zip, tumor=[t for t, n in get_pairs()], normal=[n for t, n in get_pairs()]),
        # Paired structural variant calls
        expand(join(VARIANT_DIR, "{tumor}_vs_{normal}.sv.vcf.gz"), 
               zip, tumor=[t for t, n in get_pairs()], normal=[n for t, n in get_pairs()]),
        # QC reports
        expand(join(QC_DIR, "{sample}_fastqc.zip"), sample=SAMPLES),
        expand(join(QC_DIR, "{tumor}_vs_{normal}.somatic_report.html"), 
               zip, tumor=[t for t, n in get_pairs()], normal=[n for t, n in get_pairs()]),
        join(QC_DIR, "multiqc_report.html")

# Quality control for raw reads
rule fastqc:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        join(QC_DIR, "{sample}_fastqc.zip")
    log:
        join(LOG_DIR, "fastqc", "{sample}.log")
    threads: 4
    shell:
        """
        # FastQC for quality assessment
        mkdir -p {QC_DIR}/{wildcards.sample}
        fastqc -o {QC_DIR}/{wildcards.sample} -t {threads} {input.r1} {input.r2} 2> {log}
        mv {QC_DIR}/{wildcards.sample}/*_fastqc.zip {output}
        """

# Preprocessing for C2C data using standard tools
rule preprocess_reads:
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        r1 = join(FASTQ_DIR, "{sample}_trimmed_R1.fastq.gz"),
        r2 = join(FASTQ_DIR, "{sample}_trimmed_R2.fastq.gz")
    params:
        quality = config["params"]["preprocessing"]["quality_threshold"],
        min_length = config["params"]["preprocessing"]["min_length"],
        adapter_r1 = config["params"]["preprocessing"]["adapter_r1"],
        adapter_r2 = config["params"]["preprocessing"]["adapter_r2"]
    log:
        join(LOG_DIR, "trimming", "{sample}.log")
    threads: 8
    shell:
        """
        # Trimmomatic for adapter removal and quality trimming
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {FASTQ_DIR}/{wildcards.sample}_unpaired_R1.fastq.gz \
            {output.r2} {FASTQ_DIR}/{wildcards.sample}_unpaired_R2.fastq.gz \
            ILLUMINACLIP:{params.adapter_r1}:{params.adapter_r2}:2:30:10 \
            SLIDINGWINDOW:4:{params.quality} \
            MINLEN:{params.min_length} \
            2> {log}
        """

# Generate consensus reads for C2C data
rule generate_consensus:
    input:
        r1 = join(FASTQ_DIR, "{sample}_trimmed_R1.fastq.gz"),
        r2 = join(FASTQ_DIR, "{sample}_trimmed_R2.fastq.gz")
    output:
        consensus = join(FASTQ_DIR, "{sample}_consensus.fastq.gz")
    params:
        min_quality = config["params"]["consensus"]["min_quality"],
        min_coverage = config["params"]["consensus"]["min_coverage"]
    log:
        join(LOG_DIR, "consensus", "{sample}.log")
    threads: 16
    shell:
        """
        # Using ccs (Circular Consensus Sequencing) from pbbioconda
        # This is an approximation for C2C processing
        # In a real pipeline, you might need custom scripts for C2C specifics
        ccs --min-rq {params.min_quality} \
            --min-passes {params.min_coverage} \
            --num-threads {threads} \
            --report-file {FASTQ_DIR}/{wildcards.sample}.ccs_report.txt \
            --hifi-kinetics \
            {input.r1} {input.r2} \
            {output.consensus} \
            2> {log}
        """

# BWA-MEM alignment
rule bwa_align:
    input:
        reads = join(FASTQ_DIR, "{sample}_consensus.fastq.gz"),
        reference = REFERENCE
    output:
        bam = temp(join(ALIGN_DIR, "{sample}.bam"))
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    log:
        join(LOG_DIR, "align", "{sample}.log")
    threads: 16
    shell:
        """
        # BWA-MEM alignment with read group information
        bwa mem -t {threads} -R '{params.rg}' {input.reference} {input.reads} | \
        samtools view -Sb - > {output.bam} 2> {log}
        """

# Sort and index BAM files
rule sort_bam:
    input:
        bam = join(ALIGN_DIR, "{sample}.bam")
    output:
        sorted_bam = join(ALIGN_DIR, "{sample}.sorted.bam"),
        index = join(ALIGN_DIR, "{sample}.sorted.bam.bai")
    log:
        join(LOG_DIR, "sort", "{sample}.log")
    threads: 8
    shell:
        """
        # Sort and index BAM file
        samtools sort -@ {threads} -o {output.sorted_bam} {input.bam} 2> {log}
        samtools index {output.sorted_bam} 2>> {log}
        """

# Mark duplicates
rule mark_duplicates:
    input:
        bam = join(ALIGN_DIR, "{sample}.sorted.bam")
    output:
        bam = join(ALIGN_DIR, "{sample}.dedup.bam"),
        metrics = join(QC_DIR, "{sample}.dedup_metrics.txt")
    log:
        join(LOG_DIR, "dedup", "{sample}.log")
    shell:
        """
        # Mark duplicates with Picard
        picard MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            ASSUME_SORTED=true \
            CREATE_INDEX=true \
            2> {log}
        """

# Base quality score recalibration
rule bqsr:
    input:
        bam = join(ALIGN_DIR, "{sample}.dedup.bam"),
        reference = REFERENCE,
        known_sites = config["reference"]["known_sites"]
    output:
        recal_table = join(ALIGN_DIR, "{sample}.recal_data.table"),
        bam = join(ALIGN_DIR, "{sample}.recal.bam"),
        bai = join(ALIGN_DIR, "{sample}.recal.bam.bai")  # Ensure index is created
    log:
        join(LOG_DIR, "bqsr", "{sample}.log")
    shell:
        """
        # Generate recalibration table
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.reference} \
            --known-sites {input.known_sites} \
            -O {output.recal_table} \
            2> {log}
            
        # Apply recalibration
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.reference} \
            --bqsr-recal-file {output.recal_table} \
            -O {output.bam} \
            --create-output-bam-index true \
            2>> {log}
        """

# Call SNPs and small indels with GATK
rule call_variants_gatk:
    input:
        bam = join(ALIGN_DIR, "{sample}.recal.bam"),
        reference = REFERENCE
    output:
        vcf = join(VARIANT_DIR, "{sample}.snps.vcf.gz"),
        index = join(VARIANT_DIR, "{sample}.snps.vcf.gz.tbi")
    log:
        join(LOG_DIR, "variant_calling", "{sample}.gatk.log")
    params:
        ploidy = config["params"]["variant_calling"]["ploidy"]
    shell:
        """
        # Call variants with HaplotypeCaller
        gatk HaplotypeCaller \
            -I {input.bam} \
            -R {input.reference} \
            -O {output.vcf} \
            --sample-ploidy {params.ploidy} \
            2> {log}
        """

# Call structural variants with Manta for single samples
rule call_sv_manta_single:
    input:
        bam = join(ALIGN_DIR, "{sample}.recal.bam"),
        reference = REFERENCE
    output:
        vcf = join(VARIANT_DIR, "{sample}.structural.vcf.gz"),
        index = join(VARIANT_DIR, "{sample}.structural.vcf.gz.tbi")
    log:
        join(LOG_DIR, "sv_calling", "{sample}.manta.log")
    params:
        outdir = join(VARIANT_DIR, "{sample}_manta")
    threads: 16
    shell:
        """
        # Configure Manta
        configManta.py \
            --bam {input.bam} \
            --referenceFasta {input.reference} \
            --runDir {params.outdir} \
            2> {log}
            
        # Run Manta
        {params.outdir}/runWorkflow.py -j {threads} 2>> {log}
        
        # Copy results
        cp {params.outdir}/results/variants/diploidSV.vcf.gz {output.vcf}
        cp {params.outdir}/results/variants/diploidSV.vcf.gz.tbi {output.index}
        """

# Call structural variants with Manta for tumor-normal pairs
rule call_sv_manta_paired:
    input:
        tumor_bam = join(ALIGN_DIR, "{tumor}.recal.bam"),
        normal_bam = join(ALIGN_DIR, "{normal}.recal.bam"),
        tumor_bai = join(ALIGN_DIR, "{tumor}.recal.bam.bai"),
        normal_bai = join(ALIGN_DIR, "{normal}.recal.bam.bai"),
        reference = REFERENCE
    output:
        vcf = join(VARIANT_DIR, "{tumor}_vs_{normal}.sv.vcf.gz"),
        index = join(VARIANT_DIR, "{tumor}_vs_{normal}.sv.vcf.gz.tbi")
    log:
        join(LOG_DIR, "sv_calling", "{tumor}_vs_{normal}.manta.log")
    params:
        outdir = join(VARIANT_DIR, "{tumor}_vs_{normal}_manta")
    threads: 16
    shell:
        """
        # Configure Manta for somatic analysis
        configManta.py \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam} \
            --referenceFasta {input.reference} \
            --runDir {params.outdir} \
            2> {log}
            
        # Run Manta
        {params.outdir}/runWorkflow.py -j {threads} 2>> {log}
        
        # Copy results - for tumor-normal, we use the somaticSV file
        cp {params.outdir}/results/variants/somaticSV.vcf.gz {output.vcf}
        cp {params.outdir}/results/variants/somaticSV.vcf.gz.tbi {output.index}
        """

# Call somatic variants with VarDict for tumor-normal pairs
rule call_paired_somatic_vardict:
    input:
        tumor_bam = join(ALIGN_DIR, "{tumor}.recal.bam"),
        normal_bam = join(ALIGN_DIR, "{normal}.recal.bam"),
        tumor_bai = join(ALIGN_DIR, "{tumor}.recal.bam.bai"),
        normal_bai = join(ALIGN_DIR, "{normal}.recal.bam.bai"),
        reference = REFERENCE,
        bed = config["reference"]["regions_of_interest"]
    output:
        vcf_raw = temp(join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.raw.vcf")),
        vcf = join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.vcf.gz"),
        index = join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.vcf.gz.tbi")
    params:
        tumor_name = "{tumor}",
        normal_name = "{normal}",
        min_af = config["params"]["somatic_calling"]["min_allele_freq"],
        min_qval = config["params"]["somatic_calling"]["min_qvalue"],
        max_pval = config["params"]["somatic_calling"]["max_pvalue"],
        vardict_path = config["tools"]["vardict_path"]
    log:
        join(LOG_DIR, "somatic_calling", "{tumor}_vs_{normal}.vardict.log")
    threads: 16
    shell:
        """
        # Run VarDict for tumor-normal paired somatic variant calling
        {params.vardict_path}/VarDict \
            -G {input.reference} \
            -f {params.min_af} \
            -N {params.tumor_name} \
            -b "{input.tumor_bam}|{input.normal_bam}" \
            -z -c 1 -S 2 -E 3 -g 4 \
            -th {threads} \
            {input.bed} | \
        {params.vardict_path}/testsomatic.R | \
        {params.vardict_path}/var2vcf_paired.pl \
            -N "{params.tumor_name}|{params.normal_name}" \
            -f {params.min_af} \
            -Q {params.min_qval} \
            -P {params.max_pval} > {output.vcf_raw} 2> {log}
            
        # Sort, compress, and index the VCF file
        cat {output.vcf_raw} | vcf-sort | \
        bgzip -c > {output.vcf} && \
        tabix -p vcf {output.vcf} 2>> {log}
        """

# Filter and annotate paired somatic variants
rule filter_annotate_paired_somatic:
    input:
        vcf = join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.vcf.gz"),
        reference = REFERENCE
    output:
        filtered_vcf = join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.filtered.vcf.gz"),
        filtered_index = join(VARIANT_DIR, "{tumor}_vs_{normal}.somatic.filtered.vcf.gz.tbi"),
        html = join(QC_DIR, "{tumor}_vs_{normal}.somatic_report.html")
    params:
        min_af = config["params"]["somatic_calling"]["min_allele_freq"],
        min_depth = config["params"]["somatic_calling"]["min_coverage"],
        somatic_status = "'SOMATIC'" # Filter for somatic-only variants
    log:
        join(LOG_DIR, "somatic_filtering", "{tumor}_vs_{normal}.log")
    shell:
        """
        # Filter somatic variants - for paired analysis, we use STATUS=SOMATIC
        bcftools filter \
            -i "STATUS={params.somatic_status} && AF >= {params.min_af} && DP >= {params.min_depth}" \
            -O z -o {output.filtered_vcf} \
            {input.vcf} 2> {log}
            
        # Index filtered VCF
        tabix -p vcf {output.filtered_vcf} 2>> {log}
        
        # Generate HTML report for somatic variants
        bcftools stats {output.filtered_vcf} > {VARIANT_DIR}/{wildcards.tumor}_vs_{wildcards.normal}.somatic.stats
        plot-vcfstats -p {QC_DIR}/{wildcards.tumor}_vs_{wildcards.normal}_somatic_plots {VARIANT_DIR}/{wildcards.tumor}_vs_{wildcards.normal}.somatic.stats
        
        # Create HTML report
        echo "<html><head><title>Somatic Variant Report - {wildcards.tumor} vs {wildcards.normal}</title></head>" > {output.html}
        echo "<body><h1>Somatic Variant Report for {wildcards.tumor} vs {wildcards.normal}</h1>" >> {output.html}
        echo "<p>Tumor sample: {wildcards.tumor}</p>" >> {output.html}
        echo "<p>Normal sample: {wildcards.normal}</p>" >> {output.html}
        echo "<p>Total somatic variants: $(bcftools view -H {output.filtered_vcf} | wc -l)</p>" >> {output.html}
        echo "<p>Filter criteria: Somatic status, Allele frequency >= {params.min_af}, Read depth >= {params.min_depth}</p>" >> {output.html}
        
        # Add variant type counts
        echo "<h2>Variant Type Counts</h2>" >> {output.html}
        echo "<pre>" >> {output.html}
        echo "SNVs: $(bcftools view -H {output.filtered_vcf} | grep -v "INDEL" | wc -l)" >> {output.html}
        echo "Indels: $(bcftools view -H {output.filtered_vcf} | grep "INDEL" | wc -l)" >> {output.html}
        echo "</pre>" >> {output.html}
        
        echo "<h2>Variant Statistics</h2>" >> {output.html}
        echo "<img src='{wildcards.tumor}_vs_{wildcards.normal}_somatic_plots/indels.0.png'>" >> {output.html}
        echo "<img src='{wildcards.tumor}_vs_{wildcards.normal}_somatic_plots/SNPs.0.png'>" >> {output.html}
        echo "</body></html>" >> {output.html}
        """

# Run MultiQC to aggregate QC metrics
rule multiqc:
    input:
        fastqc = expand(join(QC_DIR, "{sample}_fastqc.zip"), sample=SAMPLES),
        dedup = expand(join(QC_DIR, "{sample}.dedup_metrics.txt"), sample=SAMPLES),
        somatic = expand(join(QC_DIR, "{tumor}_vs_{normal}.somatic_report.html"), 
                         zip, tumor=[t for t, n in get_pairs()], normal=[n for t, n in get_pairs()])
    output:
        report = join(QC_DIR, "multiqc_report.html")
    log:
        join(LOG_DIR, "multiqc.log")
    shell:
        """
        # Run MultiQC to aggregate QC reports
        multiqc {OUTPUT_DIR} -o {QC_DIR} -f -n multiqc_report.html 2> {log}
        """


