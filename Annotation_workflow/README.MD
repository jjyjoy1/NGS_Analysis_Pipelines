# Variant Annotation Script

This script performs the following key functions:

## Features

### Loading and Parsing VCF
- Reads variant data in standard VCF format

### Annotation with Multiple Databases
Annotates variants using:
- **ClinVar** for clinical significance
- **COSMIC** for cancer variant data
- **gnomAD** for population frequency
- **OncoKB** for oncogenic potential

### Curation
- Resolves conflicts between different annotation sources
- Gives preference to high-confidence ClinVar annotations
- Flags other conflicts for manual review

### Classification
Classifies variants using a rule-based system into categories:
1. Pathogenic
2. Likely pathogenic
3. Variant of uncertain significance (VUS)
4. Likely benign
5. Benign

### Output
Generates:
- A complete annotated dataset
- Separate files for each classification category
- A list of variants requiring manual review


