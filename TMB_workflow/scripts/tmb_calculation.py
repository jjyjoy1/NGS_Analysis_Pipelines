#!/usr/bin/env python3
"""
Script to calculate Tumor Mutational Burden (TMB) from a VCF file.
TMB is defined as the number of somatic mutations per megabase of
the target region covered.
"""

import argparse
import pysam
import pandas as pd
import numpy as np
import os
import sys
from pybedtools import BedTool

# Get inputs from Snakemake
vcf_file = snakemake.input.vcf
bed_file = snakemake.input.bed
output_file = snakemake.output.report
log_file = snakemake.log[0]

# Configure logging
import logging
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger()

def get_target_region_size(bed_file):
    """Calculate the total size of the target regions in megabases."""
    bed = BedTool(bed_file)
    total_size = sum(interval.length for interval in bed)
    return total_size / 1e6  # Convert to megabases

def filter_somatic_variants(vcf_file):
    """
    Filter somatic variants from the VCF file.
    Keep only the following types:
    - SNVs (single nucleotide variants)
    - Small indels (< 50bp)
    - Non-synonymous mutations
    
    Exclude:
    - Germline variants
    - Synonymous mutations
    - Variants in introns/intergenic regions
    - Low quality variants
    """
    try:
        vcf = pysam.VariantFile(vcf_file)
        somatic_variants = []
        
        for variant in vcf.fetch():
            # Skip low-quality variants
            if 'FILTER' in variant.info and variant.info['FILTER'] != 'PASS':
                continue
                
            # Skip germline variants (for tumor-normal mode)
            if 'germline' in variant.info:
                continue
                
            # Check variant type
            if len(variant.ref) == 1 and len(variant.alts[0]) == 1:
                # SNV
                variant_type = 'SNV'
            elif len(variant.ref) > len(variant.alts[0]):
                # Deletion
                if len(variant.ref) - len(variant.alts[0]) >= 50:
                    continue  # Skip large deletions
                variant_type = 'DEL'
            elif len(variant.ref) < len(variant.alts[0]):
                # Insertion
                if len(variant.alts[0]) - len(variant.ref) >= 50:
                    continue  # Skip large insertions
                variant_type = 'INS'
            else:
                # Other type, skip
                continue
                
            # Check for functional impact (if available)
            impact = 'UNKNOWN'
            if 'FUNCOTATION' in variant.info:
                func_fields = variant.info['FUNCOTATION'][0].split('|')
                variant_classification = func_fields[5]  # Adjust index based on Funcotator output
                
                # Skip synonymous variants and intronic/intergenic variants
                if variant_classification in ['SILENT', 'INTRON', 'IGR', '3_UTR', '5_UTR']:
                    continue
                    
                impact = variant_classification
            
            # Add to somatic variants list
            somatic_variants.append({
                'CHROM': variant.chrom,
                'POS': variant.pos,
                'REF': variant.ref,
                'ALT': variant.alts[0],
                'TYPE': variant_type,
                'IMPACT': impact
            })
            
        return somatic_variants
    except Exception as e:
        logger.error(f"Error processing VCF file: {e}")
        sys.exit(1)

def calculate_tmb(somatic_variants, target_size_mb):
    """Calculate TMB as mutations per megabase."""
    if target_size_mb == 0:
        logger.error("Target region size is 0. Cannot calculate TMB.")
        return 0
        
    return len(somatic_variants) / target_size_mb

def generate_report(somatic_variants, tmb, target_size_mb, output_file):
    """Generate a TMB report."""
    df = pd.DataFrame(somatic_variants)
    
    # Variant type summary
    type_counts = df['TYPE'].value_counts().to_dict()
    
    # Impact summary (if available)
    impact_counts = {}
    if 'IMPACT' in df.columns:
        impact_counts = df['IMPACT'].value_counts().to_dict()
    
    # Write report
    with open(output_file, 'w') as f:
        f.write("# Tumor Mutational Burden (TMB) Report\n\n")
        f.write(f"Total target region size: {target_size_mb:.2f} Mb\n")
        f.write(f"Total somatic mutations: {len(somatic_variants)}\n")
        f.write(f"TMB: {tmb:.2f} mutations/Mb\n\n")
        
        f.write("## Mutation Type Summary\n\n")
        for mutation_type, count in type_counts.items():
            f.write(f"{mutation_type}: {count}\n")
        
        f.write("\n## Functional Impact Summary\n\n")
        for impact, count in impact_counts.items():
            f.write(f"{impact}: {count}\n")
        
        # Save variant details to a separate file
        detail_file = os.path.splitext(output_file)[0] + "_details.tsv"
        df.to_csv(detail_file, sep='\t', index=False)
        f.write(f"\nDetailed variant information saved to: {os.path.basename(detail_file)}\n")

def main():
    try:
        # Calculate target region size
        logger.info(f"Calculating target region size from: {bed_file}")
        target_size_mb = get_target_region_size(bed_file)
        logger.info(f"Target region size: {target_size_mb:.2f} Mb")
        
        # Filter somatic variants
        logger.info(f"Filtering somatic variants from: {vcf_file}")
        somatic_variants = filter_somatic_variants(vcf_file)
        logger.info(f"Found {len(somatic_variants)} somatic mutations")
        
        # Calculate TMB
        tmb = calculate_tmb(somatic_variants, target_size_mb)
        logger.info(f"TMB: {tmb:.2f} mutations/Mb")
        
        # Generate report
        logger.info(f"Generating report: {output_file}")
        generate_report(somatic_variants, tmb, target_size_mb, output_file)
        
    except Exception as e:
        logger.error(f"Error in TMB calculation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
else:
    # Being run from Snakemake
    main()


