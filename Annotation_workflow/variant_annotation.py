#!/usr/bin/env python3
"""
Variant Annotation, Curation, and Classification Pipeline

This script processes variant data (VCF format) and performs:
1. Annotation using public databases (ClinVar, COSMIC, gnomAD)
2. Annotation using proprietary databases (OncoKB)
3. Curation of annotations to resolve conflicts
4. Classification according to variant significance guidelines
"""

import os
import pandas as pd
import numpy as np
import argparse
import logging
from pysam import VariantFile
import vcfpy
import json
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class VariantProcessor:
    """Main class to handle variant processing, annotation and classification"""
    
    def __init__(self, 
                 vcf_file, 
                 clinvar_db, 
                 cosmic_db, 
                 gnomad_db, 
                 oncokb_db,
                 output_dir):
        """
        Initialize with paths to input VCF and reference databases
        
        Args:
            vcf_file (str): Path to input VCF file with variants
            clinvar_db (str): Path to ClinVar database
            cosmic_db (str): Path to COSMIC database
            gnomad_db (str): Path to gnomAD database
            oncokb_db (str): Path to OncoKB database
            output_dir (str): Directory to save results
        """
        self.vcf_file = vcf_file
        self.clinvar_db = clinvar_db
        self.cosmic_db = cosmic_db
        self.gnomad_db = gnomad_db
        self.oncokb_db = oncokb_db
        self.output_dir = output_dir
        
        # Initialize variant dataframe
        self.variants_df = None
        
        # Initialize classification rules
        self.classification_rules = {
            'pathogenic': {
                'min_criteria': 2,
                'criteria': [
                    'clinvar_pathogenic', 
                    'cosmic_high_count', 
                    'low_gnomad_freq',
                    'oncokb_oncogenic'
                ]
            },
            'likely_pathogenic': {
                'min_criteria': 1,
                'criteria': [
                    'clinvar_likely_pathogenic',
                    'cosmic_medium_count',
                    'below_gnomad_threshold',
                    'oncokb_likely_oncogenic'
                ]
            },
            'vus': {
                'min_criteria': 1, 
                'criteria': [
                    'clinvar_uncertain',
                    'cosmic_low_count',
                    'rare_in_gnomad',
                    'oncokb_unknown'
                ]
            },
            'likely_benign': {
                'min_criteria': 1,
                'criteria': [
                    'clinvar_likely_benign',
                    'not_in_cosmic',
                    'common_in_gnomad',
                    'oncokb_likely_neutral'
                ]
            },
            'benign': {
                'min_criteria': 2,
                'criteria': [
                    'clinvar_benign',
                    'high_gnomad_freq',
                    'oncokb_neutral'
                ]
            }
        }
        
    def load_vcf(self):
        """Load VCF file and convert to DataFrame"""
        logger.info(f"Loading VCF file: {self.vcf_file}")
        
        variants = []
        vcf_reader = VariantFile(self.vcf_file)
        
        for record in vcf_reader:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            
            for alt in record.alts:
                variant_data = {
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'id': record.id if record.id else f"{chrom}-{pos}-{ref}-{alt}"
                }
                
                # Extract any existing annotations from INFO field
                for key, value in record.info.items():
                    if isinstance(value, tuple) and len(value) == 1:
                        variant_data[key] = value[0]
                    else:
                        variant_data[key] = value
                
                variants.append(variant_data)
        
        self.variants_df = pd.DataFrame(variants)
        logger.info(f"Loaded {len(self.variants_df)} variants")
        
    def annotate_with_clinvar(self):
        """Annotate variants with ClinVar data"""
        logger.info("Annotating with ClinVar database")
        
        # Load ClinVar database (assuming tab-delimited format)
        clinvar_df = pd.read_csv(self.clinvar_db, sep='\t', low_memory=False)
        
        # Create lookup keys for joining
        self.variants_df['clinvar_key'] = self.variants_df.apply(
            lambda row: f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}", 
            axis=1
        )
        clinvar_df['clinvar_key'] = clinvar_df.apply(
            lambda row: f"{row['CHROM']}-{row['POS']}-{row['REF']}-{row['ALT']}", 
            axis=1
        )
        
        # Merge with variant dataframe
        self.variants_df = self.variants_df.merge(
            clinvar_df[['clinvar_key', 'CLNSIG', 'CLNREVSTAT', 'CLNDN']],
            how='left',
            on='clinvar_key'
        )
        
        # Add classification flags
        self.variants_df['clinvar_pathogenic'] = self.variants_df['CLNSIG'].str.contains('Pathogenic', na=False) & ~self.variants_df['CLNSIG'].str.contains('Likely', na=False)
        self.variants_df['clinvar_likely_pathogenic'] = self.variants_df['CLNSIG'].str.contains('Likely_pathogenic', na=False)
        self.variants_df['clinvar_uncertain'] = self.variants_df['CLNSIG'].str.contains('Uncertain', na=False)
        self.variants_df['clinvar_likely_benign'] = self.variants_df['CLNSIG'].str.contains('Likely_benign', na=False)
        self.variants_df['clinvar_benign'] = self.variants_df['CLNSIG'].str.contains('Benign', na=False) & ~self.variants_df['CLNSIG'].str.contains('Likely', na=False)
        
        # Drop temporary key
        self.variants_df.drop('clinvar_key', axis=1, inplace=True)
        
        logger.info(f"ClinVar annotations added: {self.variants_df['CLNSIG'].notna().sum()} variants annotated")
        
    def annotate_with_cosmic(self):
        """Annotate variants with COSMIC data"""
        logger.info("Annotating with COSMIC database")
        
        # Load COSMIC database (assuming tab-delimited format)
        cosmic_df = pd.read_csv(self.cosmic_db, sep='\t', low_memory=False)
        
        # Create lookup keys for joining
        self.variants_df['cosmic_key'] = self.variants_df.apply(
            lambda row: f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}", 
            axis=1
        )
        cosmic_df['cosmic_key'] = cosmic_df.apply(
            lambda row: f"{row['CHROM']}-{row['POS']}-{row['REF']}-{row['ALT']}", 
            axis=1
        )
        
        # Merge with variant dataframe - select relevant columns
        self.variants_df = self.variants_df.merge(
            cosmic_df[['cosmic_key', 'ID', 'CNT', 'GENE', 'FATHMM']],
            how='left',
            on='cosmic_key'
        )
        
        # Rename columns to avoid confusion
        self.variants_df.rename(columns={
            'ID': 'COSMIC_ID',
            'CNT': 'COSMIC_COUNT',
            'GENE': 'COSMIC_GENE',
            'FATHMM': 'COSMIC_FATHMM'
        }, inplace=True)
        
        # Add classification flags based on COSMIC count
        self.variants_df['not_in_cosmic'] = self.variants_df['COSMIC_COUNT'].isna()
        self.variants_df['cosmic_low_count'] = (self.variants_df['COSMIC_COUNT'] > 0) & (self.variants_df['COSMIC_COUNT'] < 5)
        self.variants_df['cosmic_medium_count'] = (self.variants_df['COSMIC_COUNT'] >= 5) & (self.variants_df['COSMIC_COUNT'] < 10)
        self.variants_df['cosmic_high_count'] = self.variants_df['COSMIC_COUNT'] >= 10
        
        # Drop temporary key
        self.variants_df.drop('cosmic_key', axis=1, inplace=True)
        
        logger.info(f"COSMIC annotations added: {self.variants_df['COSMIC_ID'].notna().sum()} variants annotated")
        
    def annotate_with_gnomad(self):
        """Annotate variants with gnomAD data"""
        logger.info("Annotating with gnomAD database")
        
        # Load gnomAD database (assuming tab-delimited format)
        gnomad_df = pd.read_csv(self.gnomad_db, sep='\t', low_memory=False)
        
        # Create lookup keys for joining
        self.variants_df['gnomad_key'] = self.variants_df.apply(
            lambda row: f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}", 
            axis=1
        )
        gnomad_df['gnomad_key'] = gnomad_df.apply(
            lambda row: f"{row['CHROM']}-{row['POS']}-{row['REF']}-{row['ALT']}", 
            axis=1
        )
        
        # Merge with variant dataframe - select relevant columns
        self.variants_df = self.variants_df.merge(
            gnomad_df[['gnomad_key', 'AF', 'AF_popmax', 'AC', 'AN']],
            how='left',
            on='gnomad_key'
        )
        
        # Rename columns to avoid confusion
        self.variants_df.rename(columns={
            'AF': 'GNOMAD_AF',
            'AF_popmax': 'GNOMAD_AF_POPMAX',
            'AC': 'GNOMAD_AC',
            'AN': 'GNOMAD_AN'
        }, inplace=True)
        
        # Add classification flags based on allele frequency
        self.variants_df['high_gnomad_freq'] = self.variants_df['GNOMAD_AF'] > 0.05
        self.variants_df['common_in_gnomad'] = (self.variants_df['GNOMAD_AF'] > 0.01) & (self.variants_df['GNOMAD_AF'] <= 0.05)
        self.variants_df['below_gnomad_threshold'] = (self.variants_df['GNOMAD_AF'] > 0.001) & (self.variants_df['GNOMAD_AF'] <= 0.01)
        self.variants_df['rare_in_gnomad'] = (self.variants_df['GNOMAD_AF'] > 0) & (self.variants_df['GNOMAD_AF'] <= 0.001)
        self.variants_df['low_gnomad_freq'] = (self.variants_df['GNOMAD_AF'] < 0.0001) | self.variants_df['GNOMAD_AF'].isna()
        
        # Drop temporary key
        self.variants_df.drop('gnomad_key', axis=1, inplace=True)
        
        logger.info(f"gnomAD annotations added: {self.variants_df['GNOMAD_AF'].notna().sum()} variants annotated")
        
    def annotate_with_oncokb(self):
        """Annotate variants with OncoKB data"""
        logger.info("Annotating with OncoKB database")
        
        # Load OncoKB database (assuming JSON format)
        with open(self.oncokb_db, 'r') as f:
            oncokb_data = json.load(f)
        
        # Convert to dataframe for easier merging
        oncokb_records = []
        for variant_id, variant_data in oncokb_data.items():
            # Parse variant ID to get components
            components = variant_id.split('-')
            if len(components) >= 4:
                record = {
                    'chrom': components[0],
                    'pos': int(components[1]),
                    'ref': components[2],
                    'alt': components[3],
                    'ONCOKB_ONCOGENIC': variant_data.get('oncogenic', ''),
                    'ONCOKB_MUTATION_EFFECT': variant_data.get('mutationEffect', ''),
                    'ONCOKB_DRUGS': ','.join(variant_data.get('drugs', [])),
                    'ONCOKB_GENE': variant_data.get('gene', '')
                }
                oncokb_records.append(record)
        
        oncokb_df = pd.DataFrame(oncokb_records)
        
        # Create lookup keys for joining
        oncokb_df['oncokb_key'] = oncokb_df.apply(
            lambda row: f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}", 
            axis=1
        )
        self.variants_df['oncokb_key'] = self.variants_df.apply(
            lambda row: f"{row['chrom']}-{row['pos']}-{row['ref']}-{row['alt']}", 
            axis=1
        )
        
        # Merge with variant dataframe
        self.variants_df = self.variants_df.merge(
            oncokb_df[['oncokb_key', 'ONCOKB_ONCOGENIC', 'ONCOKB_MUTATION_EFFECT', 'ONCOKB_DRUGS', 'ONCOKB_GENE']],
            how='left',
            on='oncokb_key'
        )
        
        # Add classification flags based on OncoKB oncogenicity
        self.variants_df['oncokb_oncogenic'] = self.variants_df['ONCOKB_ONCOGENIC'] == 'Oncogenic'
        self.variants_df['oncokb_likely_oncogenic'] = self.variants_df['ONCOKB_ONCOGENIC'] == 'Likely Oncogenic'
        self.variants_df['oncokb_unknown'] = (self.variants_df['ONCOKB_ONCOGENIC'] == 'Unknown') | self.variants_df['ONCOKB_ONCOGENIC'].isna()
        self.variants_df['oncokb_likely_neutral'] = self.variants_df['ONCOKB_ONCOGENIC'] == 'Likely Neutral'
        self.variants_df['oncokb_neutral'] = self.variants_df['ONCOKB_ONCOGENIC'] == 'Neutral'
        
        # Drop temporary key
        self.variants_df.drop('oncokb_key', axis=1, inplace=True)
        
        logger.info(f"OncoKB annotations added: {self.variants_df['ONCOKB_ONCOGENIC'].notna().sum()} variants annotated")
        
    def curate_annotations(self):
        """Curate annotations to resolve conflicts between sources"""
        logger.info("Curating annotations to resolve conflicts")
        
        # Example: Prefer ClinVar over OncoKB when both are available and conflicting
        conflict_mask = (
            (self.variants_df['CLNSIG'].notna()) & 
            (self.variants_df['ONCOKB_ONCOGENIC'].notna()) &
            (
                (self.variants_df['clinvar_pathogenic'] & self.variants_df['oncokb_neutral']) |
                (self.variants_df['clinvar_benign'] & self.variants_df['oncokb_oncogenic'])
            )
        )
        
        conflict_count = conflict_mask.sum()
        logger.info(f"Found {conflict_count} variants with conflicting annotations")
        
        if conflict_count > 0:
            # Create a conflict resolution column
            self.variants_df['CONFLICT_RESOLUTION'] = None
            
            # Resolve in favor of ClinVar when it has high review status
            high_confidence_clinvar = (
                conflict_mask & 
                self.variants_df['CLNREVSTAT'].str.contains('practice|expert', na=False)
            )
            
            self.variants_df.loc[high_confidence_clinvar, 'CONFLICT_RESOLUTION'] = 'Resolved in favor of ClinVar (high confidence)'
            
            # For remaining conflicts, mark for manual review
            remaining_conflicts = conflict_mask & ~high_confidence_clinvar
            self.variants_df.loc[remaining_conflicts, 'CONFLICT_RESOLUTION'] = 'Requires manual review'
            
            logger.info(f"Auto-resolved {high_confidence_clinvar.sum()} conflicts; {remaining_conflicts.sum()} require manual review")
        
    def classify_variants(self):
        """Classify variants based on annotation data"""
        logger.info("Classifying variants based on combined annotations")
        
        # Initialize classification column
        self.variants_df['CLASSIFICATION'] = 'Unknown'
        
        # Apply classification rules
        for class_name, rule in self.classification_rules.items():
            criteria = rule['criteria']
            min_criteria = rule['min_criteria']
            
            # Count how many criteria are met for each variant
            criteria_met = self.variants_df[criteria].sum(axis=1)
            
            # Assign classification if minimum criteria are met
            self.variants_df.loc[criteria_met >= min_criteria, 'CLASSIFICATION'] = class_name
        
        # Count variants in each classification
        classification_counts = self.variants_df['CLASSIFICATION'].value_counts()
        logger.info("Variant classification complete:")
        for class_name, count in classification_counts.items():
            logger.info(f"  {class_name}: {count} variants")
    
    def save_results(self):
        """Save processed variants to output files"""
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Save full annotated dataset
        output_file = os.path.join(self.output_dir, 'annotated_variants.tsv')
        self.variants_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Full annotated dataset saved to {output_file}")
        
        # Save variants by classification
        for classification in self.variants_df['CLASSIFICATION'].unique():
            class_file = os.path.join(self.output_dir, f"{classification}_variants.tsv")
            class_df = self.variants_df[self.variants_df['CLASSIFICATION'] == classification]
            class_df.to_csv(class_file, sep='\t', index=False)
            logger.info(f"{classification} variants ({len(class_df)}) saved to {class_file}")
        
        # Save variants requiring manual review
        if 'CONFLICT_RESOLUTION' in self.variants_df.columns:
            manual_review = self.variants_df[self.variants_df['CONFLICT_RESOLUTION'] == 'Requires manual review']
            if len(manual_review) > 0:
                review_file = os.path.join(self.output_dir, 'variants_for_manual_review.tsv')
                manual_review.to_csv(review_file, sep='\t', index=False)
                logger.info(f"{len(manual_review)} variants requiring manual review saved to {review_file}")
    
    def run_pipeline(self):
        """Run the complete variant annotation and classification pipeline"""
        logger.info("Starting variant annotation and classification pipeline")
        
        self.load_vcf()
        self.annotate_with_clinvar()
        self.annotate_with_cosmic()
        self.annotate_with_gnomad()
        self.annotate_with_oncokb()
        self.curate_annotations()
        self.classify_variants()
        self.save_results()
        
        logger.info("Pipeline completed successfully")
        return self.variants_df

def main():
    """Main function to parse arguments and run the pipeline"""
    parser = argparse.ArgumentParser(description='Variant Annotation, Curation, and Classification Pipeline')
    
    parser.add_argument('--vcf', required=True, help='Input VCF file with variants')
    parser.add_argument('--clinvar', required=True, help='Path to ClinVar database')
    parser.add_argument('--cosmic', required=True, help='Path to COSMIC database')
    parser.add_argument('--gnomad', required=True, help='Path to gnomAD database')
    parser.add_argument('--oncokb', required=True, help='Path to OncoKB database')
    parser.add_argument('--output', required=True, help='Output directory for results')
    
    args = parser.parse_args()
    
    # Create and run the pipeline
    processor = VariantProcessor(
        vcf_file=args.vcf,
        clinvar_db=args.clinvar,
        cosmic_db=args.cosmic,
        gnomad_db=args.gnomad,
        oncokb_db=args.oncokb,
        output_dir=args.output
    )
    
    processor.run_pipeline()

if __name__ == '__main__':
    main()


