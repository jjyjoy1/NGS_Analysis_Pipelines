There are many other important public bioinformatics datasets beyond TCGA. Here's a comprehensive overview of major public datasets in bioinformatics:
Genomics Datasets. I try to created a comprehensive python and R scripts to download expected data. 

1. GTEx (Genotype-Tissue Expression)
Multi-tissue gene expression and genetic variation data from non-diseased individuals
54 tissue types from ~1,000 individuals
Accessible via: https://gtexportal.org/

2. CCLE (Cancer Cell Line Encyclopedia)
Genomic data from ~1,000 cancer cell lines
Gene expression, mutations, drug sensitivity
Accessible via: https://portals.broadinstitute.org/ccle

3. ENCODE (Encyclopedia of DNA Elements)
Comprehensive parts list of functional elements in the human genome
ChIP-seq, RNA-seq, DNase-seq data
Accessible via: https://www.encodeproject.org/

4. GEO (Gene Expression Omnibus)
Repository of high-throughput gene expression data
Microarray and sequencing data from many studies
Accessible via: https://www.ncbi.nlm.nih.gov/geo/

5. UK Biobank
Prospective cohort study with genetic data from 500,000 individuals
Genotyping, imaging, health records
Requires application: https://www.ukbiobank.ac.uk/

Proteomics Datasets
1. PRIDE (PRoteomics IDEntifications Database)
Repository of mass spectrometry proteomics data
Accessible via: https://www.ebi.ac.uk/pride/

2. CPTAC (Clinical Proteomic Tumor Analysis Consortium)
Comprehensive proteomic characterization of tumors
Complements TCGA with proteomic data
Accessible via: https://proteomics.cancer.gov/data-portal

Single-Cell Datasets
1. Human Cell Atlas
Comprehensive reference maps of all human cells
Multiple single-cell technologies
Accessible via: https://www.humancellatlas.org/

2. CellxGene
Single-cell gene expression datasets
Interactive data explorer
Accessible via: https://cellxgene.cziscience.com/

3. Single Cell Portal (Broad Institute)
Collection of single-cell datasets
Accessible via: https://singlecell.broadinstitute.org/

Microbiome Datasets
1. Human Microbiome Project (HMP)
Characterization of the human microbiome
16S rRNA and metagenomic sequencing data
Accessible via: https://hmpdacc.org/

2. MGnify (formerly EBI Metagenomics)
Archive of metagenomic and metatranscriptomic data
Accessible via: https://www.ebi.ac.uk/metagenomics/

Drug/Compound Datasets
1. GDSC (Genomics of Drug Sensitivity in Cancer)
Drug response data for cancer cell lines
Accessible via: https://www.cancerrxgene.org/

2. LINCS (Library of Integrated Network-Based Cellular Signatures)
Gene expression and cellular response signatures to perturbations
Accessible via: https://lincsproject.org/

3. DrugBank
Detailed drug data with comprehensive drug target information
Accessible via: https://go.drugbank.com/

Tools for Accessing Multiple Datasets
1. Bioconductor
ExperimentHub package provides access to many datasets
Install with: BiocManager::install("ExperimentHub")

2. recount3
R package for accessing uniformly processed RNA-seq datasets
Install with: BiocManager::install("recount3")

3. cBioPortal
Web resource for exploring cancer genomics datasets
Accessible via: https://www.cbioportal.org/

