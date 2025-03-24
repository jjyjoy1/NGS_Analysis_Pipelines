The HMP data is hosted at the HMP Data Analysis and Coordination Center (DACC):

Main portal: https://hmpdacc.org/
Data browser: https://portal.hmpdacc.org/

For most data types, we can:

Navigate to the data browser
Select data types and samples of interest
Generate a manifest file
Download using the manifest


Understanding HMP Manifest Files
HMP manifest files are typically tab-separated values (TSV) files with these columns:

id: Unique identifier for the file
md5: MD5 checksum for file integrity verification
size: File size in bytes
urls: URL(s) where the file can be downloaded
tags: Metadata tags for the file
Additional metadata columns (varies by data type)

Types of HMP Data Available

16S rRNA Sequencing Data

Raw reads (FASTQ)
OTU tables
Taxonomic profiles


Metagenomic WGS Data

Raw reads (FASTQ)
Assembled contigs
Gene catalogs


Metabolomics Data

LC-MS data
NMR data


Viromics Data

Viral sequences
Phage annotations
