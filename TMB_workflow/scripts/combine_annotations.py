#!/usr/bin/env python3
"""
Script to combine annotations from Funcotator (VCF) and Oncotator (MAF)
into a single annotated VCF file.
"""

import os
import gzip
import pandas as pd
import pysam
import logging
import sys

# Get inputs from Snakemake
funcotator_vcf = snakemake.input.funcotator
oncotator_maf = snakemake.input.oncotator
output_vcf = snakemake.output.vcf
output_idx = snakemake.output.idx
log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger()

def load_oncotator_maf(maf_file):
    """
    Load Oncotator MAF file into a pandas DataFrame.
    Create a chromosome:position:ref:alt key for merging.
    """
    try:
        # Handle compressed MAF files
        if maf_file.endswith('.gz'):
            with gzip.open(maf_file, 'rt') as f:
                # Skip comment lines starting with #
                comment_lines = 0
                for line in f:
                    if line.startswith('#'):
                        comment_lines += 1
                    else:
                        break
            
            # Read the MAF file
            maf_df = pd.read_csv(maf_file, sep='\t', comment='#', skiprows=comment_lines)
        else:
            # Uncompressed MAF file
            maf_df = pd.read_csv(maf_file, sep='\t', comment='#')
        
        # Create a unique key for each variant
        maf_df['variant_key'] = maf_df.apply(
            lambda x: f"{x['Chromosome']}:{x['Start_Position']}:{x['Reference_Allele']}:{x['Tumor_Seq_Allele2']}",
            axis=1
        )
        
        # Select relevant columns for annotation
        relevant_columns = [
            'variant_key',
            'Hugo_Symbol',
            'Entrez_Gene_Id',
            'Variant_Classification',
            'Variant_Type',
            'Transcript_ID',
            'Exon_Number',
            'COSMIC_overlapping',
            'COSMIC_n_overlapping',
            'COSMIC_fusion_genes',
            'dbSNP_RS',
            'dbSNP_Val_Status',
            'Protein_Change',
            'Tumor_Sample_Barcode'
        ]
        
        # Filter to only columns that exist in the MAF file
        available_columns = [col for col in relevant_columns if col in maf_df.columns]
        
        return maf_df[available_columns].set_index('variant_key')
        
    except Exception as e:
        logger.error(f"Error loading MAF file: {e}")
        return pd.DataFrame()

def combine_annotations(vcf_file, maf_df, output_file):
    """
    Add Oncotator annotations to the Funcotator-annotated VCF.
    """
    try:
        # Open input and output VCF files
        vcf = pysam.VariantFile(vcf_file)
        
        # Add Oncotator header fields to the VCF
        oncotator_headers = {
            'ONCOTATOR_GENE': {'Number': '1', 'Type': 'String', 'Description': 'Gene symbol from Oncotator'},
            'ONCOTATOR_VARIANT_CLASS': {'Number': '1', 'Type': 'String', 'Description': 'Variant classification from Oncotator'},
            'ONCOTATOR_PROTEIN': {'Number': '1', 'Type': 'String', 'Description': 'Protein change from Oncotator'},
            'ONCOTATOR_COSMIC': {'Number': '1', 'Type': 'String', 'Description': 'COSMIC IDs from Oncotator'},
            'ONCOTATOR_DBSNP': {'Number': '1', 'Type': 'String', 'Description': 'dbSNP ID from Oncotator'}
        }
        
        for header_id, header_info in oncotator_headers.items():
            vcf.header.add_meta('INFO', items=[
                ('ID', header_id),
                ('Number', header_info['Number']),
                ('Type', header_info['Type']),
                ('Description', header_info['Description'])
            ])
        
        # Open output VCF file
        with pysam.VariantFile(output_file, 'w', header=vcf.header) as out_vcf:
            # Process each variant
            for variant in vcf.fetch():
                # Create a variant key
                variant_key = f"{variant.chrom}:{variant.pos}:{variant.ref}:{variant.alts[0]}"
                
                # Add Oncotator annotations if the variant is in the MAF file
                if variant_key in maf_df.index:
                    oncotator_data = maf_df.loc[variant_key]
                    
                    # Add Oncotator gene
                    if 'Hugo_Symbol' in oncotator_data:
                        variant.info['ONCOTATOR_GENE'] = oncotator_data['Hugo_Symbol']
                    
                    # Add variant classification
                    if 'Variant_Classification' in oncotator_data:
                        variant.info['ONCOTATOR_VARIANT_CLASS'] = oncotator_data['Variant_Classification']
                    
                    # Add protein change
                    if 'Protein_Change' in oncotator_data:
                        variant.info['ONCOTATOR_PROTEIN'] = oncotator_data['Protein_Change']
                    
                    # Add COSMIC info
                    if 'COSMIC_overlapping' in oncotator_data and not pd.isna(oncotator_data['COSMIC_overlapping']):
                        variant.info['ONCOTATOR_COSMIC'] = oncotator_data['COSMIC_overlapping']
                    
                    # Add dbSNP info
                    if 'dbSNP_RS' in oncotator_data and not pd.isna(oncotator_data['dbSNP_RS']):
                        variant.info['ONCOTATOR_DBSNP'] = oncotator_data['dbSNP_RS']
                
                # Write the variant to the output
                out_vcf.write(variant)
        
        # Create tabix index for the output VCF
        pysam.tabix_index(output_file, preset='vcf', force=True)
        
        logger.info(f"Successfully combined annotations from {vcf_file} and {maf_df.shape[0]} variants from Oncotator")
        
    except Exception as e:
        logger.error(f"Error combining annotations: {e}")
        sys.exit(1)

def main():
    try:
        # Load Oncotator MAF file
        logger.info(f"Loading Oncotator MAF file: {oncotator_maf}")
        maf_df = load_oncotator_maf(oncotator_maf)
        logger.info(f"Loaded {maf_df.shape[0]} variants from Oncotator")
        
        # Combine annotations
        logger.info(f"Combining annotations from Funcotator VCF: {funcotator_vcf}")
        combine_annotations(funcotator_vcf, maf_df, output_vcf)
        
    except Exception as e:
        logger.error(f"Error in annotation combining: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
else:
    # Being run from Snakemake
    main()

