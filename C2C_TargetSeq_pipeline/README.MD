# C2C Sequencing Analysis Pipeline

A comprehensive Snakemake pipeline for processing Illumina Circle-to-Circle Amplification (C2C) sequencing data with support for paired tumor-normal somatic variant detection.

## Overview

This pipeline provides a complete workflow for processing high-accuracy C2C sequencing data, with special focus on somatic variant detection in cancer samples. It includes quality control, read processing, alignment, and variant calling using both standard open-source tools and optionally Illumina's DRAGEN platform.

### Key Features

- Complete workflow from FASTQ to variant calls
- Support for paired tumor-normal somatic variant detection
- Optimized for C2C sequencing data characteristics
- Options for both DRAGEN and standard bioinformatics tools
- Comprehensive somatic variant detection:
  - SNVs and small indels with VarDict
  - Structural variants with Manta
- Detailed QC reporting and visualization

## Requirements

### Standard Pipeline

- **Software Environment**:
  - Python 3.8+ with Snakemake 7.18.2+
  - Conda or Mamba for environment management

- **Core Dependencies**:
  - BWA (0.7.17+)
  - Samtools (1.16+)
  - GATK (4.3.0+)
  - VarDict (1.8.3+)
  - Manta (1.6.0+)
  - Trimmomatic (0.39+)
  - FastQC (0.11.9+)
  - MultiQC (1.13+)
  - bcftools (1.15+)
  - Python BioConda packages for consensus generation

### DRAGEN Pipeline

- **Hardware**:
  - DRAGEN server or access to DRAGEN in the cloud
  
- **Software**:
  - Illumina DRAGEN Bio-IT Platform
  - Snakemake 7.18.2+

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/c2c-sequencing-pipeline.git
cd c2c-sequencing-pipeline
```

2. Create and activate the environment:
```bash
# Using conda
conda env create -f environment.yml
conda activate c2c_pipeline

# Alternatively, using mamba (faster)
mamba env create -f environment.yml
mamba activate c2c_pipeline
```

3. Install VarDict (if not using conda):
```bash
git clone https://github.com/AstraZeneca-NGS/VarDictJava.git
cd VarDictJava
./gradlew clean installDist
# Add the path to config.yaml
```

## Usage

### 1. Configure Your Analysis

Edit the appropriate configuration file:

- For standard tools: `config_standard.yaml`
- For DRAGEN: `config.yaml`

Update with your sample information and paths:
```yaml
# Example for paired tumor-normal samples
samples:
  patient1_tumor:
    type: "tumor"
    r1: "/path/to/patient1_tumor_R1.fastq.gz"
    r2: "/path/to/patient1_tumor_R2.fastq.gz"
    matched_normal: "patient1_normal"
    
  patient1_normal:
    type: "normal"
    r1: "/path/to/patient1_normal_R1.fastq.gz"
    r2: "/path/to/patient1_normal_R2.fastq.gz"
```

### 2. Run the Pipeline

For the standard pipeline:
```bash
snakemake -s Snakefile --configfile config_standard.yaml --cores 32 --use-conda
```

For the DRAGEN pipeline:
```bash
snakemake -s Snakefile_dragen --configfile config.yaml --cores 32
```

### 3. Pipeline Output

The pipeline creates the following directory structure:
```
output_dir/
├── fastq/         # Processed FASTQ files
├── aligned/       # BAM files and indices
├── variants/      # VCF files with variants
│   ├── *.snps.vcf.gz              # Germline variants
│   ├── *.structural.vcf.gz        # Structural variants
│   ├── *_vs_*.somatic.vcf.gz      # Somatic variants (paired analysis)
│   └── *_vs_*.sv.vcf.gz           # Somatic structural variants
├── qc/            # Quality control reports
│   ├── *_fastqc.zip               # FastQC reports
│   ├── *_vs_*.somatic_report.html # Somatic variant reports
│   └── multiqc_report.html        # Aggregated QC report
└── logs/          # Log files from each step
```

## C2C Sequencing Background

Illumina's Circle-to-Circle Amplification (C2C) sequencing technology offers:

- Higher accuracy through consensus sequencing of circular templates
- Reduced systematic errors with multiple passes over the same DNA
- Improved detection of low-frequency variants
- Better performance in repetitive and complex regions

This pipeline is optimized to leverage these advantages for somatic variant detection.

## Pipeline Details

### Standard Pipeline Workflow

1. **Initial QC**: FastQC assessment of raw reads
2. **Read Preprocessing**: Adapter removal and quality trimming
3. **Consensus Generation**: Processing C2C-specific read characteristics
4. **Alignment**: Mapping to reference genome with BWA-MEM
5. **BAM Processing**: Sorting, duplicate marking, and base quality recalibration
6. **Variant Calling**:
   - Germline variants with GATK
   - Somatic variants with VarDict (paired tumor-normal mode)
   - Structural variants with Manta (both single sample and paired modes)
7. **Variant Filtering**: Post-processing and filtering based on quality metrics
8. **Reporting**: Comprehensive QC and variant reports

### DRAGEN Pipeline Workflow

1. **C2C Preprocessing**: DRAGEN-specific processing of C2C reads
2. **Alignment**: Using DRAGEN's optimized mapper
3. **Variant Calling**: DRAGEN variant callers with C2C-specific parameters
4. **QC**: DRAGEN QC metrics and MultiQC aggregation

## Advanced Configuration

### Customizing Resources

Adjust threads and memory in the config files:
```yaml
# Example resource configuration
threads: 32
params:
  memory_per_thread: "4G"
```

### Tuning Variant Calling Parameters

For somatic variant detection sensitivity:
```yaml
params:
  somatic_calling:
    min_allele_freq: 0.01  # Detect variants at 1% frequency
    min_coverage: 20       # Minimum read depth
```

## Troubleshooting

### Common Issues

1. **Missing dependencies**:
   - Ensure all tools are properly installed or use `--use-conda`
   
2. **Resource limitations**:
   - Adjust thread count with `--cores` parameter
   - For memory issues, reduce threads or use a cluster

3. **VarDict errors**:
   - Check that paths in `config_standard.yaml` point to the correct VarDict installation

### Getting Help

Create an issue on GitHub with:
- Error messages and logs
- Configuration file (with sensitive paths removed)
- Command line used to run the pipeline

## Citation

If you use this pipeline in your research, please cite:

```
Author, A. (2024). C2C Sequencing Analysis Pipeline. GitHub Repository: https://github.com/yourusername/c2c-sequencing-pipeline
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- Illumina for C2C sequencing technology
- The developers of all open-source tools used in this pipeline
- The Snakemake community for workflow management

