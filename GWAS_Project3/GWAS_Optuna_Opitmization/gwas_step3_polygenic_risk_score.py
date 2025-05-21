#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRS Optimization Module for GWAS Pipeline

This script optimizes PRS parameters based on the results from the GWAS statistical analysis.
It takes optimized GWAS results as input and finds the best parameters for constructing PRS models.
"""

import os
import sys
import argparse
import logging
import subprocess
import numpy as np
import pandas as pd
import optuna
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, r2_score
from datetime import datetime
import json
import shutil

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("prs_optimizer.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("PRS-Optimizer")

# Global variables to store state
global_state = {
    'args': None,
    'study': None,
    'current_trial_number': 0,
    'gwas_results': None,
    'base_prefix': None,
    'output_dir': None,
    'plink_path': None,
    'target_data_prefix': None,
    'db_path': None,
    'temp_dir': None
}

def setup_paths(args):
    """
    Set up paths for input and output files.
    
    Args:
        args: Command line arguments
    """
    # Create output directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    global_state['output_dir'] = os.path.join(args.output_dir, f"prs_opt_{timestamp}")
    os.makedirs(global_state['output_dir'], exist_ok=True)
    
    # Create temp directory for intermediate files
    global_state['temp_dir'] = os.path.join(global_state['output_dir'], "temp")
    os.makedirs(global_state['temp_dir'], exist_ok=True)
    
    # Set PLINK path
    global_state['plink_path'] = args.plink_path if args.plink_path else "plink"
    
    # Set target data prefix
    global_state['target_data_prefix'] = args.target_data_prefix
    
    # Set GWAS results path
    global_state['gwas_results'] = args.gwas_results
    
    # Extract base prefix for naming
    global_state['base_prefix'] = os.path.basename(args.target_data_prefix)
    
    # Set Optuna database path
    global_state['db_path'] = os.path.join(global_state['output_dir'], "optuna_trials.db")
    
    logger.info(f"Output directory: {global_state['output_dir']}")
    logger.info(f"Temporary files will be stored in: {global_state['temp_dir']}")

def validate_inputs(args):
    """
    Validate input files and dependencies.
    
    Args:
        args: Command line arguments
    """
    # Check GWAS results file
    if not os.path.exists(args.gwas_results):
        raise FileNotFoundError(f"GWAS results file not found: {args.gwas_results}")
    
    # Check target data files
    required_target_files = [
        f"{args.target_data_prefix}.bed",
        f"{args.target_data_prefix}.bim",
        f"{args.target_data_prefix}.fam"
    ]
    
    for file_path in required_target_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Required target data file not found: {file_path}")
    
    # Check phenotype file if provided
    if args.target_pheno and not os.path.exists(args.target_pheno):
        raise FileNotFoundError(f"Target phenotype file not found: {args.target_pheno}")
    
    # Check PLINK installation
    try:
        subprocess.run([global_state['plink_path'], "--help"], 
                      stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE)
    except:
        logger.warning("PLINK executable not found. Please check path.")
    
    logger.info("All input validations passed")

def setup_optuna_study(args):
    """
    Set up Optuna study for parameter optimization.
    
    Args:
        args: Command line arguments
    """
    storage_url = f"sqlite:///{global_state['db_path']}"
    
    # Check if we need to load existing study
    if args.resume_study and os.path.exists(global_state['db_path']):
        global_state['study'] = optuna.load_study(
            study_name=args.study_name,
            storage=storage_url
        )
        logger.info(f"Resuming existing study '{args.study_name}' with {len(global_state['study'].trials)} previous trials")
    else:
        # Create new study
        global_state['study'] = optuna.create_study(
            study_name=args.study_name,
            storage=storage_url,
            direction="maximize",  # We'll define a custom objective that should be maximized
            load_if_exists=False
        )
        logger.info(f"Created new optimization study: {args.study_name}")

def load_gwas_summary_stats(file_path, args):
    """
    Load and preprocess GWAS summary statistics.
    
    Args:
        file_path: Path to GWAS summary statistics file
        args: Command line arguments
        
    Returns:
        DataFrame with processed summary statistics
    """
    logger.info(f"Loading GWAS summary statistics from {file_path}")
    
    # Try to detect file format (space, tab, or comma-separated)
    with open(file_path, 'r') as f:
        header = f.readline().strip()
    
    if '\t' in header:
        sep = '\t'
    elif ',' in header:
        sep = ','
    else:
        sep = None  # Let pandas infer (usually whitespace)
    
    # Load the file
    gwas_data = pd.read_csv(file_path, sep=sep)
    
    # Standardize column names to lowercase for easier matching
    gwas_data.columns = [col.lower() for col in gwas_data.columns]
    
    # Try to identify key columns
    snp_col = None
    p_col = None
    effect_col = None
    other_col = None
    effect_size_col = None
    se_col = None
    
    # SNP ID column
    for col in ['snp', 'rsid', 'rs', 'marker', 'markername', 'snpid']:
        if col in gwas_data.columns:
            snp_col = col
            break
    
    # P-value column
    for col in ['p', 'pval', 'pvalue', 'p_value', 'p-value', 'p.value', 'pval_nominal']:
        if col in gwas_data.columns:
            p_col = col
            break
    
    # Effect allele
    for col in ['a1', 'allele1', 'effect_allele', 'effectallele', 'alt', 'effect']:
        if col in gwas_data.columns:
            effect_col = col
            break
    
    # Other allele
    for col in ['a2', 'allele2', 'other_allele', 'otherallele', 'ref', 'non_effect']:
        if col in gwas_data.columns:
            other_col = col
            break
    
    # Effect size
    for col in ['beta', 'b', 'effect', 'effects', 'effect_size', 'or', 'odds_ratio']:
        if col in gwas_data.columns:
            effect_size_col = col
            break
    
    # Standard error
    for col in ['se', 'std_err', 'stderr', 'standard_error']:
        if col in gwas_data.columns:
            se_col = col
            break
    
    # Check if we have the necessary columns
    missing_cols = []
    if snp_col is None:
        missing_cols.append("SNP")
    if p_col is None:
        missing_cols.append("P-value")
    
    if missing_cols:
        raise ValueError(f"Missing required columns in GWAS summary stats: {', '.join(missing_cols)}")
    
    # Report found columns
    logger.info(f"Detected columns in GWAS summary stats:")
    logger.info(f"  SNP column: {snp_col}")
    logger.info(f"  P-value column: {p_col}")
    logger.info(f"  Effect allele column: {effect_col}")
    logger.info(f"  Other allele column: {other_col}")
    logger.info(f"  Effect size column: {effect_size_col}")
    logger.info(f"  Standard error column: {se_col}")
    
    # Create standardized dataframe with required columns
    standardized_data = pd.DataFrame()
    standardized_data['SNP'] = gwas_data[snp_col]
    standardized_data['P'] = gwas_data[p_col]
    
    # Add effect allele if available
    if effect_col is not None:
        standardized_data['A1'] = gwas_data[effect_col]
    
    # Add other allele if available
    if other_col is not None:
        standardized_data['A2'] = gwas_data[other_col]
    
    # Add effect size if available
    if effect_size_col is not None:
        # Check if this is odds ratio and convert to beta if needed
        if effect_size_col.lower() in ['or', 'odds_ratio']:
            standardized_data['BETA'] = np.log(gwas_data[effect_size_col])
            logger.info("Converting odds ratio to beta by taking natural log")
        else:
            standardized_data['BETA'] = gwas_data[effect_size_col]
    
    # Add standard error if available
    if se_col is not None:
        standardized_data['SE'] = gwas_data[se_col]
    
    # Filter out missing values
    standardized_data = standardized_data.dropna(subset=['SNP', 'P'])
    
    # Convert p-values to numeric, handling any string values
    standardized_data['P'] = pd.to_numeric(standardized_data['P'], errors='coerce')
    standardized_data = standardized_data.dropna(subset=['P'])
    
    # Filter out invalid p-values
    standardized_data = standardized_data[(standardized_data['P'] >= 0) & (standardized_data['P'] <= 1)]
    
    logger.info(f"Processed GWAS summary stats: {len(standardized_data)} valid variants")
    
    return standardized_data

def prepare_prs_input_files(gwas_data, args):
    """
    Prepare input files for PRS calculation.
    
    Args:
        gwas_data: DataFrame with GWAS summary statistics
        args: Command line arguments
        
    Returns:
        Path to the base file for PRS calculation
    """
    # Create a base file for PRS calculation
    base_file = os.path.join(global_state['temp_dir'], "gwas_base.txt")
    
    # Check which columns we have available
    has_effect = 'A1' in gwas_data.columns and 'A2' in gwas_data.columns
    has_beta = 'BETA' in gwas_data.columns
    
    # Minimum required columns for PRS
    required_cols = ['SNP', 'P']
    
    # Add allele columns if available
    if has_effect:
        required_cols.extend(['A1', 'A2'])
    
    # Add effect size if available
    if has_beta:
        required_cols.append('BETA')
    
    # Save to file
    gwas_data[required_cols].to_csv(base_file, sep='\t', index=False)
    
    logger.info(f"Prepared PRS base file: {base_file}")
    
    return base_file

def run_prs_clumping(base_file, params, args):
    """
    Run LD clumping for PRS calculation.
    
    Args:
        base_file: Path to the base file with GWAS summary statistics
        params: Dictionary with PRS parameters
        args: Command line arguments
        
    Returns:
        Path to the clumped file
    """
    # Create output prefix for clumping
    clump_prefix = os.path.join(global_state['temp_dir'], 
                              f"clump_r2_{params['clump_r2']}_kb_{params['clump_kb']}")
    
    # Construct PLINK command for clumping
    cmd = [
        global_state['plink_path'],
        "--bfile", global_state['target_data_prefix'],
        "--clump", base_file,
        "--clump-p1", str(params['p_threshold']),
        "--clump-p2", str(params['p_threshold']),
        "--clump-r2", str(params['clump_r2']),
        "--clump-kb", str(params['clump_kb']),
        "--out", clump_prefix
    ]
    
    # Add reference population if specified
    if args.ld_reference_prefix:
        cmd = [
            global_state['plink_path'],
            "--bfile", args.ld_reference_prefix,
            "--clump", base_file,
            "--clump-p1", str(params['p_threshold']),
            "--clump-p2", str(params['p_threshold']),
            "--clump-r2", str(params['clump_r2']),
            "--clump-kb", str(params['clump_kb']),
            "--out", clump_prefix
        ]
    
    # Run clumping
    logger.info(f"Running LD clumping with r2={params['clump_r2']}, kb={params['clump_kb']}, p={params['p_threshold']}")
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if process.returncode != 0:
        logger.warning(f"Clumping failed: {process.stderr.decode()}")
        return None
    
    # Check if clumped file was created
    clumped_file = f"{clump_prefix}.clumped"
    if not os.path.exists(clumped_file):
        logger.warning(f"Clumped file not found: {clumped_file}")
        return None
    
    return clumped_file

def create_score_file(clumped_file, gwas_data, params, args):
    """
    Create a score file for PRS calculation.
    
    Args:
        clumped_file: Path to the clumped SNP file
        gwas_data: DataFrame with GWAS summary statistics
        params: Dictionary with PRS parameters
        args: Command line arguments
        
    Returns:
        Path to the score file
    """
    # Load clumped SNPs
    if clumped_file is not None and os.path.exists(clumped_file):
        clumped_snps = pd.read_csv(clumped_file, sep='\s+')
        
        # Get the SNP column name (can be 'SNP' or 'ID' depending on PLINK version)
        snp_col = 'SNP' if 'SNP' in clumped_snps.columns else 'ID'
        
        if snp_col in clumped_snps.columns:
            selected_snps = set(clumped_snps[snp_col].tolist())
            logger.info(f"Loaded {len(selected_snps)} clumped SNPs")
        else:
            logger.warning(f"SNP column not found in clumped file. Available columns: {clumped_snps.columns}")
            selected_snps = set()
    else:
        # If clumping failed or was not requested, select SNPs based on p-value threshold only
        selected_snps = set(gwas_data.loc[gwas_data['P'] <= params['p_threshold'], 'SNP'].tolist())
        logger.info(f"Selected {len(selected_snps)} SNPs based on p-value threshold {params['p_threshold']}")
    
    # Create a new dataframe with selected SNPs
    score_data = gwas_data[gwas_data['SNP'].isin(selected_snps)].copy()
    
    # Apply additional filters if specified
    if params.get('maf_filter', 0) > 0:
        # If we have MAF information in the GWAS data
        if 'MAF' in gwas_data.columns:
            score_data = score_data[score_data['MAF'] >= params['maf_filter']]
            logger.info(f"Applied MAF filter: {len(score_data)} SNPs remaining")
    
    # Check if we have effect sizes
    has_beta = 'BETA' in score_data.columns
    
    # Create score file with optimal parameters
    clumped_file = f"{clump_prefix}.clumped"
    score_file = os.path.join(final_dir, "optimal_score.txt")
    
    if os.path.exists(clumped_file):
        # Load clumped SNPs
        clumped_snps = pd.read_csv(clumped_file, sep='\s+')
        
        # Get the SNP column name
        snp_col = 'SNP' if 'SNP' in clumped_snps.columns else 'ID'
        
        if snp_col in clumped_snps.columns:
            selected_snps = set(clumped_snps[snp_col].tolist())
            logger.info(f"Selected {len(selected_snps)} SNPs for optimal PRS")
            
            # Create score file
            score_data = gwas_data[gwas_data['SNP'].isin(selected_snps)].copy()
            
            if 'BETA' in score_data.columns:
                # Use effect sizes if available
                score_data[['SNP', 'A1', 'BETA']].to_csv(score_file, sep='\t', index=False, header=False)
            else:
                # Use binary scoring if no effect sizes
                binary_score = pd.DataFrame({
                    'SNP': score_data['SNP'],
                    'A1': score_data['A1'] if 'A1' in score_data.columns else 'A',
                    'SCORE': 1
                })
                binary_score.to_csv(score_file, sep='\t', index=False, header=False)
        else:
            logger.warning(f"SNP column not found in clumped file. Available columns: {clumped_snps.columns}")
            return None
    else:
        # If clumping failed, select SNPs based on p-value threshold only
        selected_snps = gwas_data.loc[gwas_data['P'] <= best_params['p_threshold'], 'SNP'].tolist()
        logger.info(f"Selected {len(selected_snps)} SNPs based on p-value threshold {best_params['p_threshold']}")
        
        # Create score file
        score_data = gwas_data[gwas_data['SNP'].isin(selected_snps)].copy()
        
        if 'BETA' in score_data.columns:
            # Use effect sizes if available
            score_data[['SNP', 'A1', 'BETA']].to_csv(score_file, sep='\t', index=False, header=False)
        else:
            # Use binary scoring if no effect sizes
            binary_score = pd.DataFrame({
                'SNP': score_data['SNP'],
                'A1': score_data['A1'] if 'A1' in score_data.columns else 'A',
                'SCORE': 1
            })
            binary_score.to_csv(score_file, sep='\t', index=False, header=False)
    
    # Calculate final PRS
    prs_prefix = os.path.join(final_dir, "optimal_prs")
    
    # Determine score parameters based on file format
    with open(score_file, 'r') as f:
        first_line = f.readline().strip().split()
        has_beta = len(first_line) >= 3  # Check if we have effect size
    
    # Construct PLINK command for PRS calculation
    cmd = [
        global_state['plink_path'],
        "--bfile", global_state['target_data_prefix'],
        "--score", score_file
    ]
    
    # Add score parameters based on format
    if has_beta:
        cmd.extend(["1", "2", "3"])  # SNP, A1, BETA columns
    
    # Add output prefix
    cmd.extend(["--out", prs_prefix])
    
    # Run PRS calculation
    logger.info(f"Calculating final PRS with optimal parameters")
    subprocess.run(cmd)
    
    # Check if PRS file was created
    prs_file = f"{prs_prefix}.profile"
    if not os.path.exists(prs_file):
        logger.warning(f"Final PRS file not found: {prs_file}")
        return None
    
    # Copy the final PRS file to output directory
    final_prs_file = os.path.join(global_state['output_dir'], "optimal_prs.profile")
    shutil.copy(prs_file, final_prs_file)
    
    logger.info(f"Generated optimal PRS: {final_prs_file}")
    
    return final_prs_file

def evaluate_final_prs(prs_file, args):
    """
    Evaluate the final PRS and generate visualizations.
    
    Args:
        prs_file: Path to the final PRS file
        args: Command line arguments
    """
    # Create output directory for evaluation results
    eval_dir = os.path.join(global_state['output_dir'], "evaluation")
    os.makedirs(eval_dir, exist_ok=True)
    
    # Load PRS results
    prs_results = pd.read_csv(prs_file, sep='\s+')
    
    # Basic distribution plot of PRS
    plt.figure(figsize=(10, 6))
    sns.histplot(prs_results['SCORE'], kde=True)
    plt.title('Distribution of Polygenic Risk Scores')
    plt.xlabel('PRS')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(eval_dir, "prs_distribution.png"), dpi=300)
    plt.close()
    
    # If phenotype data is available, create additional plots
    if args.target_pheno:
        # Load phenotype data
        pheno_data = pd.read_csv(args.target_pheno, sep='\s+')
        
        # Get phenotype column name
        if args.pheno_name and args.pheno_name in pheno_data.columns:
            pheno_col = args.pheno_name
        else:
            # Assume phenotype is in the first non-ID column
            id_cols = ['FID', 'IID']
            pheno_cols = [col for col in pheno_data.columns if col not in id_cols]
            
            if not pheno_cols:
                logger.warning("No phenotype columns found in phenotype file")
                return
            
            pheno_col = pheno_cols[0]
        
        # Set index for joining
        pheno_data = pheno_data.set_index(['FID', 'IID'])
        prs_results = prs_results.set_index(['FID', 'IID'])
        
        # Join phenotype with PRS
        combined_data = prs_results.join(pheno_data[[pheno_col]])
        
        # Check if we have missing values
        if combined_data[pheno_col].isna().any():
            logger.warning(f"Found {combined_data[pheno_col].isna().sum()} samples with missing phenotype")
            combined_data = combined_data.dropna(subset=[pheno_col])
        
        # Determine if phenotype is binary or continuous
        unique_values = combined_data[pheno_col].unique()
        is_binary = len(unique_values) <= 2
        
        if is_binary:
            # For binary traits
            
            # Convert phenotype to binary (0/1)
            combined_data['PHENO'] = combined_data[pheno_col].astype(int)
            
            # Box plot of PRS by phenotype status
            plt.figure(figsize=(10, 6))
            sns.boxplot(x='PHENO', y='SCORE', data=combined_data)
            plt.title('PRS Distribution by Phenotype Status')
            plt.xlabel('Phenotype (Case/Control)')
            plt.ylabel('PRS')
            plt.savefig(os.path.join(eval_dir, "prs_by_phenotype.png"), dpi=300)
            plt.close()
            
            # Calculate ROC curve and AUC
            y_true = combined_data['PHENO']
            y_score = combined_data['SCORE']
            
            fpr, tpr, thresholds = roc_curve(y_true, y_score)
            roc_auc = auc(fpr, tpr)
            
            # Plot ROC curve
            plt.figure(figsize=(10, 6))
            plt.plot(fpr, tpr, label=f'AUC = {roc_auc:.3f}')
            plt.plot([0, 1], [0, 1], 'k--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('Receiver Operating Characteristic')
            plt.legend(loc='lower right')
            plt.savefig(os.path.join(eval_dir, "roc_curve.png"), dpi=300)
            plt.close()
            
            # Calculate and plot odds ratios by PRS quantile
            combined_data['PRS_quantile'] = pd.qcut(combined_data['SCORE'], 10, labels=False)
            
            # Use the lowest quantile as reference
            reference_quantile = 0
            
            # Calculate odds ratios
            odds_ratios = []
            ci_lower = []
            ci_upper = []
            quantiles = []
            
            from scipy.stats import fisher_exact
            from statsmodels.stats.proportion import proportion_confint
            
            for q in range(10):
                if q == reference_quantile:
                    # Skip reference quantile
                    odds_ratios.append(1.0)
                    ci_lower.append(1.0)
                    ci_upper.append(1.0)
                    quantiles.append(q)
                    continue
                
                # Create 2x2 contingency table
                ref_cases = sum((combined_data['PRS_quantile'] == reference_quantile) & (combined_data['PHENO'] == 1))
                ref_controls = sum((combined_data['PRS_quantile'] == reference_quantile) & (combined_data['PHENO'] == 0))
                q_cases = sum((combined_data['PRS_quantile'] == q) & (combined_data['PHENO'] == 1))
                q_controls = sum((combined_data['PRS_quantile'] == q) & (combined_data['PHENO'] == 0))
                
                # Calculate odds ratio
                table = [[q_cases, q_controls], [ref_cases, ref_controls]]
                odds_ratio, p_value = fisher_exact(table)
                
                # Calculate confidence intervals
                ci = proportion_confint(q_cases, q_cases + q_controls, method='wilson')
                
                odds_ratios.append(odds_ratio)
                ci_lower.append(ci[0])
                ci_upper.append(ci[1])
                quantiles.append(q)
            
            # Plot odds ratios
            plt.figure(figsize=(12, 6))
            plt.errorbar(quantiles, odds_ratios, yerr=[
                [or_val - lower for or_val, lower in zip(odds_ratios, ci_lower)],
                [upper - or_val for or_val, upper in zip(odds_ratios, ci_upper)]
            ], fmt='o')
            plt.axhline(y=1.0, color='r', linestyle='-')
            plt.xlabel('PRS Decile')
            plt.ylabel('Odds Ratio (relative to lowest decile)')
            plt.title('Odds Ratios by PRS Decile')
            plt.xticks(range(10))
            plt.savefig(os.path.join(eval_dir, "odds_ratios.png"), dpi=300)
            plt.close()
            
            # Save evaluation metrics
            metrics = {
                'auc': roc_auc,
                'odds_ratios': {
                    'quantiles': quantiles,
                    'values': odds_ratios,
                    'ci_lower': ci_lower,
                    'ci_upper': ci_upper
                }
            }
            
        else:
            # For continuous traits
            
            # Scatter plot of PRS vs phenotype
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x='SCORE', y=pheno_col, data=combined_data, alpha=0.5)
            plt.title('PRS vs Phenotype')
            plt.xlabel('PRS')
            plt.ylabel(pheno_col)
            plt.savefig(os.path.join(eval_dir, "prs_vs_phenotype.png"), dpi=300)
            plt.close()
            
            # Calculate R²
            from sklearn.linear_model import LinearRegression
            
            X = combined_data[['SCORE']].values
            y = combined_data[pheno_col].values
            
            model = LinearRegression()
            model.fit(X, y)
            
            r2 = r2_score(y, model.predict(X))
            
            # Plot regression line
            plt.figure(figsize=(10, 6))
            sns.regplot(x='SCORE', y=pheno_col, data=combined_data)
            plt.title(f'PRS vs Phenotype (R² = {r2:.3f})')
            plt.xlabel('PRS')
            plt.ylabel(pheno_col)
            plt.savefig(os.path.join(eval_dir, "prs_regression.png"), dpi=300)
            plt.close()
            
            # Save evaluation metrics
            metrics = {
                'r2': r2,
                'coefficient': model.coef_[0],
                'intercept': model.intercept_
            }
        
        # Save metrics to file
        with open(os.path.join(eval_dir, "evaluation_metrics.json"), 'w') as f:
            json.dump(metrics, f, indent=2)
        
        logger.info(f"Evaluation completed: {eval_dir}")

def create_prs_script(best_params, args):
    """
    Create a script to apply the optimal PRS to new data.
    
    Args:
        best_params: Best parameters from optimization
        args: Command line arguments
        
    Returns:
        Path to the script
    """
    script_dir = os.path.join(global_state['output_dir'], "scripts")
    os.makedirs(script_dir, exist_ok=True)
    
    script_path = os.path.join(script_dir, "apply_optimal_prs.sh")
    
    with open(script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Script to apply optimal PRS to new data\n")
        f.write(f"# Generated by prs_optimization.py on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("# Input parameters\n")
        f.write("TARGET_DATA=$1\n")
        f.write("OUTPUT_PREFIX=$2\n\n")
        
        f.write("if [ -z \"$TARGET_DATA\" ] || [ -z \"$OUTPUT_PREFIX\" ]; then\n")
        f.write("    echo \"Usage: $0 <target_data_prefix> <output_prefix>\"\n")
        f.write("    echo \"Example: $0 new_data optimal_prs\"\n")
        f.write("    exit 1\n")
        f.write("fi\n\n")
        
        f.write("# Check if required files exist\n")
        f.write("if [ ! -f \"${TARGET_DATA}.bed\" ] || [ ! -f \"${TARGET_DATA}.bim\" ] || [ ! -f \"${TARGET_DATA}.fam\" ]; then\n")
        f.write("    echo \"Error: Required target data files not found.\"\n")
        f.write("    exit 1\n")
        f.write("fi\n\n")
        
        # Copy optimal score file to script directory
        score_file = os.path.join(global_state['output_dir'], "final_model", "optimal_score.txt")
        script_score_file = os.path.join(script_dir, "optimal_score.txt")
        
        if os.path.exists(score_file):
            shutil.copy(score_file, script_score_file)
            
            # Determine if we have effect sizes
            with open(score_file, 'r') as sf:
                first_line = sf.readline().strip().split()
                has_beta = len(first_line) >= 3
            
            f.write("# Calculate PRS using optimal parameters\n")
            f.write(f"PLINK_PATH=\"{global_state['plink_path']}\"\n")
            f.write("SCORE_FILE=\"optimal_score.txt\"\n\n")
            
            f.write("# Run PRS calculation\n")
            f.write("$PLINK_PATH --bfile $TARGET_DATA --score $SCORE_FILE")
            
            if has_beta:
                f.write(" 1 2 3")  # SNP, A1, BETA columns
            
            f.write(" --out $OUTPUT_PREFIX\n\n")
            
            f.write("echo \"PRS calculation complete: ${OUTPUT_PREFIX}.profile\"\n")
        else:
            f.write("# Error: Optimal score file not found\n")
            f.write("echo \"Error: Optimal score file not found. Run optimization first.\"\n")
            f.write("exit 1\n")
    
    # Make script executable
    os.chmod(script_path, 0o755)
    
    logger.info(f"Created PRS application script: {script_path}")
    
    return script_path

def run_optimization(args):
    """
    Run the optimization process.
    
    Args:
        args: Command line arguments
        
    Returns:
        Tuple of (best_params, best_value)
    """
    logger.info(f"Starting optimization with {args.n_trials} trials")
    
    # Load GWAS summary statistics
    gwas_data = load_gwas_summary_stats(args.gwas_results, args)
    global_state['gwas_data'] = gwas_data
    
    # Prepare input files for PRS
    base_file = prepare_prs_input_files(gwas_data, args)
    global_state['base_file'] = base_file
    
    # Run optimization
    global_state['study'].optimize(
        objective, 
        n_trials=args.n_trials,
        n_jobs=args.n_jobs
    )
    
    # Get best parameters
    best_params = global_state['study'].best_params
    best_value = global_state['study'].best_value
    
    logger.info(f"Best parameters: {best_params}")
    if global_state['study'].best_trial.user_attrs.get('auc') is not None:
        logger.info(f"Best AUC: {global_state['study'].best_trial.user_attrs['auc']:.4f}")
    elif global_state['study'].best_trial.user_attrs.get('r2') is not None:
        logger.info(f"Best R²: {global_state['study'].best_trial.user_attrs['r2']:.4f}")
    
    # Generate optimal PRS
    final_prs = generate_optimal_prs(best_params, args)
    
    # Evaluate final PRS
    if final_prs:
        evaluate_final_prs(final_prs, args)
    
    # Create PRS application script
    create_prs_script(best_params, args)
    
    # Create visualization of optimization results
    vis_dir = os.path.join(global_state['output_dir'], "visualization")
    os.makedirs(vis_dir, exist_ok=True)
    
    try:
        # Plot optimization history
        fig = optuna.visualization.plot_optimization_history(global_state['study'])
        fig.write_image(os.path.join(vis_dir, "optimization_history.png"))
        
        # Plot parameter importance
        fig = optuna.visualization.plot_param_importances(global_state['study'])
        fig.write_image(os.path.join(vis_dir, "param_importance.png"))
        
        # Plot parameter contour plots
        fig = optuna.visualization.plot_contour(global_state['study'])
        fig.write_image(os.path.join(vis_dir, "param_contour.png"))
        
        # Plot individual parameter effects
        fig = optuna.visualization.plot_slice(global_state['study'])
        fig.write_image(os.path.join(vis_dir, "param_slice.png"))
        
        logger.info(f"Created optimization visualizations: {vis_dir}")
    except Exception as e:
        logger.error(f"Failed to create visualizations: {str(e)}")
    
    # Save optimization results to JSON
    results = {
        'best_params': best_params,
        'best_value': best_value,
        'n_trials': len(global_state['study'].trials),
        'optimization_metric': 'auc' if global_state['study'].best_trial.user_attrs.get('auc') is not None else 'r2',
        'timestamp': datetime.now().isoformat()
    }
    
    results_file = os.path.join(global_state['output_dir'], "optimization_results.json")
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Saved optimization results to: {results_file}")
    
    return best_params, best_value

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="PRS Parameter Optimization with Optuna",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    parser.add_argument("--gwas-results", required=True, 
                       help="Path to GWAS summary statistics file")
    parser.add_argument("--target-data-prefix", required=True, 
                       help="Prefix for target PLINK files (.bed, .bim, .fam)")
    parser.add_argument("--target-pheno", 
                       help="Phenotype file for target data in PLINK format")
    parser.add_argument("--pheno-name", 
                       help="Column name in phenotype file to use")
    parser.add_argument("--ld-reference-prefix", 
                       help="Prefix for LD reference PLINK files (if different from target)")
    
    # Output options
    parser.add_argument("--output-dir", default="prs_optimizer_output",
                       help="Directory for output files")
    parser.add_argument("--study-name", default="prs_optimization",
                       help="Name for the Optuna study")
    parser.add_argument("--resume-study", action="store_true",
                       help="Resume an existing study with the same name")
    
    # Tool paths
    parser.add_argument("--plink-path", default="plink",
                       help="Path to PLINK executable")
    
    # Optimization settings
    parser.add_argument("--n-trials", type=int, default=50,
                       help="Number of optimization trials to run")
    parser.add_argument("--n-jobs", type=int, default=1,
                       help="Number of parallel jobs for optimization")
    parser.add_argument("--optimize-maf", action="store_true",
                       help="Optimize MAF thresholds")
    parser.add_argument("--optimize-scoring", action="store_true",
                       help="Optimize scoring methods")
    
    return parser.parse_args()

def main():
    """Main function to run the PRS parameter optimization."""
    try:
        # Parse command line arguments
        args = parse_arguments()
        global_state['args'] = args
        
        # Display configuration summary
        print("\n" + "=" * 60)
        print("PRS Parameter Optimization")
        print("=" * 60)
        print(f"GWAS Results: {args.gwas_results}")
        print(f"Target Data: {args.target_data_prefix}")
        print(f"Output Directory: {args.output_dir}")
        print(f"Number of Trials: {args.n_trials}")
        
        if args.target_pheno:
            print(f"Target Phenotype: {args.target_pheno}")
            if args.pheno_name:
                print(f"Phenotype Column: {args.pheno_name}")
        
        if args.ld_reference_prefix:
            print(f"LD Reference: {args.ld_reference_prefix}")
        
        print("\nOptimization Parameters:")
        print(f"- Optimize MAF thresholds: {'Yes' if args.optimize_maf else 'No'}")
        print(f"- Optimize scoring methods: {'Yes' if args.optimize_scoring else 'No'}")
        
        print("=" * 60 + "\n")
        
        # Initialize paths and directories
        setup_paths(args)
        
        # Validate input files
        validate_inputs(args)
        
        # Setup Optuna study
        setup_optuna_study(args)
        
        # Run optimization
        best_params, best_value = run_optimization(args)
        
        # Print final results
        print("\n" + "=" * 60)
        print("Optimization Complete!")
        print("=" * 60)
        print(f"Best parameters: {best_params}")
        print(f"Best objective value: {best_value:.4f}")
        print(f"Results saved to: {global_state['output_dir']}")
        
        print("\nOutput Files:")
        print(f"- Optimal PRS: {os.path.join(global_state['output_dir'], 'optimal_prs.profile')}")
        print(f"- Evaluation: {os.path.join(global_state['output_dir'], 'evaluation')}")
        print(f"- Visualization: {os.path.join(global_state['output_dir'], 'visualization')}")
        print(f"- Application Script: {os.path.join(global_state['output_dir'], 'scripts', 'apply_optimal_prs.sh')}")
        
        print("=" * 60)
        
        return 0  # Success
    
    except FileNotFoundError as e:
        print(f"\nError: {e}", file=sys.stderr)
        print("Please check that all input files exist and are accessible.", file=sys.stderr)
        return 1  # File not found error
    
    except ValueError as e:
        print(f"\nError: {e}", file=sys.stderr)
        print("Please check your input parameters and try again.", file=sys.stderr)
        return 2  # Invalid parameter error
    
    except KeyboardInterrupt:
        print("\nOptimization interrupted by user.", file=sys.stderr)
        return 130  # Interrupted by user
    
    except Exception as e:
        print(f"\nUnexpected error: {e}", file=sys.stderr)
        print("For detailed error information, check the log file.", file=sys.stderr)
        logging.error(f"Unexpected error: {e}", exc_info=True)
        return 3  # Other error


if __name__ == "__main__":
    sys.exit(main()) file
    score_file = os.path.join(global_state['temp_dir'], 
                             f"score_p{params['p_threshold']}_snps{len(score_data)}.txt")
    
    if has_beta:
        # Use effect sizes if available
        score_data[['SNP', 'A1', 'BETA']].to_csv(score_file, sep='\t', index=False, header=False)
    else:
        # Use binary scoring if no effect sizes
        # Create a dataframe with SNP, A1, and a score of 1
        binary_score = pd.DataFrame({
            'SNP': score_data['SNP'],
            'A1': score_data['A1'] if 'A1' in score_data.columns else 'A',  # Default to 'A' if missing
            'SCORE': 1
        })
        binary_score.to_csv(score_file, sep='\t', index=False, header=False)
    
    logger.info(f"Created score file with {len(score_data)} SNPs: {score_file}")
    
    return score_file

def calculate_prs(score_file, params, args):
    """
    Calculate PRS using PLINK.
    
    Args:
        score_file: Path to the score file
        params: Dictionary with PRS parameters
        args: Command line arguments
        
    Returns:
        Path to the PRS results file
    """
    # Create output prefix for PRS
    prs_prefix = os.path.join(global_state['temp_dir'], 
                            f"prs_p{params['p_threshold']}_r2{params['clump_r2']}_kb{params['clump_kb']}")
    
    # Determine score parameters based on file format
    with open(score_file, 'r') as f:
        first_line = f.readline().strip().split()
        has_beta = len(first_line) >= 3  # Check if we have effect size in the score file
    
    # Construct PLINK command for PRS calculation
    cmd = [
        global_state['plink_path'],
        "--bfile", global_state['target_data_prefix'],
        "--score", score_file
    ]
    
    # Add score parameters based on format
    if has_beta:
        cmd.extend(["1", "2", "3"])  # SNP, A1, BETA columns
    
    # Add output prefix
    cmd.extend(["--out", prs_prefix])
    
    # Run PRS calculation
    logger.info(f"Calculating PRS with {os.path.basename(score_file)}")
    process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if process.returncode != 0:
        logger.warning(f"PRS calculation failed: {process.stderr.decode()}")
        return None
    
    # Check if PRS file was created
    prs_file = f"{prs_prefix}.profile"
    if not os.path.exists(prs_file):
        logger.warning(f"PRS file not found: {prs_file}")
        return None
    
    return prs_file

def evaluate_prs_performance(prs_file, args):
    """
    Evaluate PRS performance against phenotype.
    
    Args:
        prs_file: Path to the PRS results file
        args: Command line arguments
        
    Returns:
        Dictionary with performance metrics
    """
    # Load PRS results
    prs_results = pd.read_csv(prs_file, sep='\s+')
    
    # Initialize metrics
    metrics = {
        'auc': None,
        'r2': None,
        'n_samples': len(prs_results),
        'n_snps_used': prs_results['CNT'].iloc[0] if 'CNT' in prs_results.columns else None
    }
    
    # Load phenotype data if provided
    if args.target_pheno:
        pheno_data = pd.read_csv(args.target_pheno, sep='\s+')
        
        # Get phenotype column name
        if args.pheno_name and args.pheno_name in pheno_data.columns:
            pheno_col = args.pheno_name
        else:
            # Assume phenotype is in the first non-ID column
            id_cols = ['FID', 'IID']
            pheno_cols = [col for col in pheno_data.columns if col not in id_cols]
            
            if not pheno_cols:
                logger.warning("No phenotype columns found in phenotype file")
                return metrics
            
            pheno_col = pheno_cols[0]
            logger.info(f"Using {pheno_col} as phenotype column")
        
        # Set index for joining
        pheno_data = pheno_data.set_index(['FID', 'IID'])
        prs_results = prs_results.set_index(['FID', 'IID'])
        
        # Join phenotype with PRS
        combined_data = prs_results.join(pheno_data[[pheno_col]])
        
        # Check if we have missing values
        if combined_data[pheno_col].isna().any():
            logger.warning(f"Found {combined_data[pheno_col].isna().sum()} samples with missing phenotype")
            combined_data = combined_data.dropna(subset=[pheno_col])
        
        # Determine if phenotype is binary or continuous
        unique_values = combined_data[pheno_col].unique()
        is_binary = len(unique_values) <= 2
        
        if is_binary:
            # For binary traits, calculate AUC
            try:
                # Ensure binary coding (0/1)
                y_true = combined_data[pheno_col].astype(int)
                
                # Use SCORE column for PRS
                y_score = combined_data['SCORE']
                
                # Calculate ROC curve and AUC
                fpr, tpr, _ = roc_curve(y_true, y_score)
                roc_auc = auc(fpr, tpr)
                
                metrics['auc'] = roc_auc
                logger.info(f"Calculated AUC: {roc_auc:.4f}")
                
                # Save ROC curve data for visualization
                roc_data = pd.DataFrame({
                    'fpr': fpr,
                    'tpr': tpr
                })
                roc_file = os.path.join(os.path.dirname(prs_file), f"{os.path.basename(prs_file)}_roc.csv")
                roc_data.to_csv(roc_file, index=False)
                
            except Exception as e:
                logger.warning(f"Failed to calculate AUC: {str(e)}")
        else:
            # For continuous traits, calculate R²
            try:
                from sklearn.linear_model import LinearRegression
                
                # Prepare data
                X = combined_data[['SCORE']].values
                y = combined_data[pheno_col].values
                
                # Fit linear regression
                model = LinearRegression()
                model.fit(X, y)
                
                # Calculate R²
                r2 = r2_score(y, model.predict(X))
                
                metrics['r2'] = r2
                logger.info(f"Calculated R²: {r2:.4f}")
                
            except Exception as e:
                logger.warning(f"Failed to calculate R²: {str(e)}")
    
    return metrics

def define_parameter_space(trial, args):
    """
    Define the parameter space for Optuna to search.
    
    Args:
        trial: Optuna trial object
        args: Command line arguments
        
    Returns:
        Dictionary of parameters for this trial
    """
    params = {}
    
    # 1. P-value threshold for SNP inclusion
    params["p_threshold"] = trial.suggest_float("p_threshold", 1e-8, 0.5, log=True)
    
    # 2. LD clumping parameters
    params["clump_r2"] = trial.suggest_float("clump_r2", 0.1, 0.9)
    params["clump_kb"] = trial.suggest_int("clump_kb", 100, 1000)
    
    # 3. Optional MAF filter
    if args.optimize_maf:
        params["maf_filter"] = trial.suggest_float("maf_filter", 0.001, 0.1, log=True)
    else:
        params["maf_filter"] = 0
    
    # 4. Scoring method (if multiple methods are available)
    if args.optimize_scoring:
        params["scoring_method"] = trial.suggest_categorical(
            "scoring_method", 
            ["standard", "center", "weighted"]
        )
    else:
        params["scoring_method"] = "standard"
    
    return params

def objective(trial):
    """
    Optuna objective function to maximize.
    
    Args:
        trial: Optuna trial object
    
    Returns:
        Objective value to maximize (AUC or R²)
    """
    args = global_state['args']
    
    # Keep track of current trial number for file naming
    global_state['current_trial_number'] = trial.number
    
    # Get parameters for this trial
    params = define_parameter_space(trial, args)
    
    # Log parameters
    logger.info(f"Trial {trial.number}: {params}")
    
    # Load previously prepared GWAS data
    gwas_data = global_state.get('gwas_data')
    if gwas_data is None:
        gwas_data = load_gwas_summary_stats(args.gwas_results, args)
        global_state['gwas_data'] = gwas_data
    
    # Prepare base file if not already done
    base_file = global_state.get('base_file')
    if base_file is None:
        base_file = prepare_prs_input_files(gwas_data, args)
        global_state['base_file'] = base_file
    
    # Run clumping
    clumped_file = run_prs_clumping(base_file, params, args)
    
    # Create score file
    score_file = create_score_file(clumped_file, gwas_data, params, args)
    
    # Calculate PRS
    prs_file = calculate_prs(score_file, params, args)
    
    if prs_file is None:
        logger.warning("PRS calculation failed")
        return -100.0  # Return a very low value for failed trials
    
    # Evaluate PRS performance
    metrics = evaluate_prs_performance(prs_file, args)
    
    # Store metrics in trial user attributes
    for key, value in metrics.items():
        if value is not None:
            trial.set_user_attr(key, value)
    
    # Determine objective value based on phenotype type
    if metrics['auc'] is not None:
        # Binary trait, maximize AUC
        objective_value = metrics['auc']
    elif metrics['r2'] is not None:
        # Continuous trait, maximize R²
        objective_value = metrics['r2']
    else:
        # If no metrics available, return a low value
        objective_value = -100.0
    
    return objective_value

def generate_optimal_prs(best_params, args):
    """
    Generate final PRS using optimal parameters.
    
    Args:
        best_params: Best parameters from optimization
        args: Command line arguments
        
    Returns:
        Path to the final PRS file
    """
    # Create output directory for final results
    final_dir = os.path.join(global_state['output_dir'], "final_model")
    os.makedirs(final_dir, exist_ok=True)
    
    # Load GWAS data
    gwas_data = global_state.get('gwas_data')
    if gwas_data is None:
        gwas_data = load_gwas_summary_stats(args.gwas_results, args)
    
    # Prepare base file
    base_file = os.path.join(final_dir, "gwas_base.txt")
    if 'BETA' in gwas_data.columns:
        gwas_data[['SNP', 'A1', 'A2', 'P', 'BETA']].to_csv(base_file, sep='\t', index=False)
    else:
        gwas_data[['SNP', 'A1', 'A2', 'P']].to_csv(base_file, sep='\t', index=False)
    
    # Run clumping with optimal parameters
    clump_prefix = os.path.join(final_dir, f"optimal_clump")
    cmd = [
        global_state['plink_path'],
        "--bfile", global_state['target_data_prefix'],
        "--clump", base_file,
        "--clump-p1", str(best_params['p_threshold']),
        "--clump-p2", str(best_params['p_threshold']),
        "--clump-r2", str(best_params['clump_r2']),
        "--clump-kb", str(best_params['clump_kb']),
        "--out", clump_prefix
    ]
    
    # Add reference population if specified
    if args.ld_reference_prefix:
        cmd = [
            global_state['plink_path'],
            "--bfile", args.ld_reference_prefix,
            "--clump", base_file,
            "--clump-p1", str(best_params['p_threshold']),
            "--clump-p2", str(best_params['p_threshold']),
            "--clump-r2", str(best_params['clump_r2']),
            "--clump-kb", str(best_params['clump_kb']),
            "--out", clump_prefix
        ]
    
    # Run clumping
    logger.info(f"Running final LD clumping with optimal parameters")
    subprocess.run(cmd)
    
    # Create score


