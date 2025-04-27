# TMB Pipeline Overview

I've created a comprehensive Snakemake pipeline for Tumor Mutational Burden (TMB) analysis that can handle both tumor-only and tumor-normal paired sample analysis.

## Pipeline Structure

### Main Snakemake Workflow (`tmb-pipeline.smk`)
- Handles all steps from FASTQ to TMB calculation
- Supports both tumor-only and tumor-normal paired sample analysis

### TMB Calculation Script (`calculate_tmb.py`)
- Processes annotated VCF files to calculate TMB
- Filters for relevant somatic mutations
- Generates detailed TMB reports with mutation types and impacts
- Calculates mutations per megabase of targeted regions

### Configuration System (`config.yaml`)
- Flexible sample definition for both tumor-only and paired analysis
- Reference genome and annotation resources
- Customizable parameters for each pipeline step
- Ability to specify which samples have matched normals

### Conda Environments
- Separate environments for different pipeline stages
- Ensures reproducibility and proper dependency management
- Includes all necessary tools for each analysis step

## Key Features

### Quality Control and Preprocessing
- FastQC for raw read quality assessment
- Adapter trimming with Trimmomatic
- Read alignment with BWA-MEM
- Duplicate marking with Picard
- Base quality recalibration with GATK

### Variant Calling
- Tumor-only mode using Mutect2 with germline resource filtering
- Paired tumor-normal mode for better somatic variant detection
- Variant filtering with FilterMutectCalls

### Annotation and TMB Calculation
- Variant annotation using Funcotator for functional impact
- Filtering for non-synonymous mutations
- TMB calculation based on the number of somatic mutations per megabase
- Detailed reporting on mutation types and impacts

## How to Use the Pipeline

### Setup
1. Clone the repository containing the pipeline files
2. Create conda environments using the provided YAML files
3. Modify the `config.yaml` file to specify your samples and reference files

### Running the Pipeline
```bash
# Dry run to check workflow
snakemake -n

# Run the pipeline with 16 cores
snakemake --use-conda --cores 16 -s tumor_mutational_burden.smk

# Run specific targets
snakemake --use-conda --cores 16 -s tumor_mutational_burden.smk results/tmb/tumor1.tmb.txt


