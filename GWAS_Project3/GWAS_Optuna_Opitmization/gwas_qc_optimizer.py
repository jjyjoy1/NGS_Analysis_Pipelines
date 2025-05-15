#!/usr/bin/env python3
# gwas_qc_optimizer.py
# Optimize GWAS QC parameters using Optuna for a cancer case-control study

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import optuna
import scipy.stats as stats
import logging
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedShuffleSplit
import sqlite3
import joblib

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"gwas_optim_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('gwas_optimizer')

# Configuration
class Config:
    # Path configurations
    PLINK_PATH = "/path/to/plink"  # Update with your plink path
    GCTA_PATH = "/path/to/gcta"    # Update with your gcta path
    RAW_DATA = "your_raw_genotype_data"  # Update with your raw data path
    OUT_PREFIX = "breast_cancer_gwas"
    TEMP_DIR = "optuna_runs"
    
    # Covariate file paths
    PHENOTYPE_FILE = "phenotype.txt"  # Format: FID IID PHENO (0=control, 1=case)
    COVAR_FILE = "covariates.txt"     # Categorical covariates (e.g., sex, cancer_subtype)
    QCOVAR_FILE = "qcovariates.txt"   # Quantitative covariates (e.g., age, PC1-10)
    
    # Subsampling parameters
    SUBSAMPLE_SIZE = 5000  # Number of samples to use for optimization
    SUBSAMPLE_SEED = 42    # For reproducibility
    
    # Optimization parameters
    N_TRIALS = 50          # Number of parameter combinations to try
    TIMEOUT = 72 * 3600    # Maximum runtime in seconds (72 hours)
    N_JOBS = 1             # Number of parallel jobs (careful with disk I/O)
    
    # Database settings
    DB_PATH = "sqlite:///gwas_optuna.db"  # SQLite database file
    STUDY_NAME = "gwas_qc_optimization"   # Study name in the database
    
    # Paths to store results
    RESULTS_DIR = "optuna_results"

# Utility functions for sample and SNP tracking
def count_samples(bfile_prefix):
    """Count number of samples in PLINK binary file"""
    cmd = f"{Config.PLINK_PATH} --bfile {bfile_prefix} --freq --out temp_count"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open("temp_count.log", 'r') as f:
        for line in f:
            if "samples" in line:
                return int(line.strip().split()[0])
    return 0

def count_snps(bfile_prefix):
    """Count number of SNPs in PLINK binary file"""
    cmd = f"{Config.PLINK_PATH} --bfile {bfile_prefix} --freq --out temp_count"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open("temp_count.log", 'r') as f:
        for line in f:
            if "variants" in line or "SNPs" in line:
                return int(line.strip().split()[0])
    return 0

def get_case_control_counts(bfile_prefix, phenotype_file):
    """Get counts of cases and controls"""
    # Extract IDs from PLINK file
    cmd = f"{Config.PLINK_PATH} --bfile {bfile_prefix} --write-samples --out temp_ids"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Read IDs and phenotype file
    try:
        ids = pd.read_csv("temp_ids.samples", sep='\t', header=0)
        pheno = pd.read_csv(phenotype_file, sep='\s+', header=None, names=['FID', 'IID', 'PHENO'])
        
        # Merge to get phenotypes of included samples
        merged = pd.merge(ids, pheno, left_on=['ID_1', 'ID_2'], right_on=['FID', 'IID'])
        
        # Count cases and controls
        cases = sum(merged['PHENO'] == 2)  # In PLINK, 2=case, 1=control, 0=missing
        controls = sum(merged['PHENO'] == 1)
        
        return cases, controls
    except Exception as e:
        logger.error(f"Error getting case-control counts: {e}")
        return 0, 0

def calculate_lambda_gc(gwas_results_file):
    """Calculate genomic inflation factor lambda"""
    try:
        # Read GWAS results
        results = pd.read_csv(gwas_results_file, sep='\s+')
        
        # Extract p-values
        pvals = results['p']
        
        # Calculate observed and expected chi-square values
        obs_chi2 = -2 * np.log(pvals)
        df = 1  # Degrees of freedom
        expected_chi2_median = stats.chi2.ppf(0.5, df)
        observed_chi2_median = np.median(obs_chi2)
        
        # Calculate lambda GC
        lambda_gc = observed_chi2_median / expected_chi2_median
        
        return lambda_gc
    except Exception as e:
        logger.error(f"Error calculating lambda GC: {e}")
        return None

def compare_covariate_distribution(original_samples, qc_samples, covariate_file, categorical_cols=None, quantitative_cols=None):
    """
    Compare covariate distributions before and after QC.
    Returns a score indicating the degree of change (lower is better).
    """
    try:
        # Read covariate file
        covars = pd.read_csv(covariate_file, sep='\s+')
        
        # Filter to original and QC samples
        orig_covars = covars[covars['IID'].isin(original_samples)]
        qc_covars = covars[covars['IID'].isin(qc_samples)]
        
        total_score = 0
        
        # Compare categorical variables using chi-square
        if categorical_cols:
            for col in categorical_cols:
                if col in covars.columns:
                    # Create contingency tables
                    orig_counts = orig_covars[col].value_counts(normalize=True)
                    qc_counts = qc_covars[col].value_counts(normalize=True)
                    
                    # Align indexes
                    all_cats = sorted(set(orig_counts.index) | set(qc_counts.index))
                    orig_props = [orig_counts.get(cat, 0) for cat in all_cats]
                    qc_props = [qc_counts.get(cat, 0) for cat in all_cats]
                    
                    # Chi-square test for goodness of fit
                    expected_counts = np.array(orig_props) * len(qc_covars)
                    observed_counts = np.array(qc_props) * len(qc_covars)
                    
                    # Replace zeros to avoid division errors
                    expected_counts = np.where(expected_counts == 0, 0.001, expected_counts)
                    
                    chi2 = sum((observed_counts - expected_counts)**2 / expected_counts)
                    total_score += chi2
        
        # Compare quantitative variables using KS test
        if quantitative_cols:
            for col in quantitative_cols:
                if col in covars.columns:
                    # Perform Kolmogorov-Smirnov test
                    orig_vals = orig_covars[col].dropna()
                    qc_vals = qc_covars[col].dropna()
                    
                    if len(orig_vals) > 0 and len(qc_vals) > 0:
                        ks_stat, _ = stats.ks_2samp(orig_vals, qc_vals)
                        total_score += ks_stat
        
        return total_score
    except Exception as e:
        logger.error(f"Error comparing covariate distributions: {e}")
        return 100  # Return high penalty if error occurs

def create_subsample(input_prefix, output_prefix, sample_size, seed):
    """Create a stratified subsample of the data for faster optimization"""
    try:
        # Get phenotype data
        pheno = pd.read_csv(Config.PHENOTYPE_FILE, sep='\s+', header=None, names=['FID', 'IID', 'PHENO'])
        
        # Create stratified split
        splitter = StratifiedShuffleSplit(n_splits=1, test_size=sample_size, random_state=seed)
        indices = list(range(len(pheno)))
        
        # Get indices for the subsample
        for _, subsample_indices in splitter.split(indices, pheno['PHENO']):
            subsample = pheno.iloc[subsample_indices]
        
        # Write subsample IDs to file
        subsample[['FID', 'IID']].to_csv('subsample.txt', sep=' ', index=False, header=False)
        
        # Extract subsample with PLINK
        cmd = f"{Config.PLINK_PATH} --bfile {input_prefix} --keep subsample.txt --make-bed --out {output_prefix}"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        logger.info(f"Created subsample with {len(subsample)} samples")
        return True
    except Exception as e:
        logger.error(f"Error creating subsample: {e}")
        return False

# Main function to run GWAS QC with given parameters
def run_gwas_qc(trial_id, params):
    """
    Run GWAS QC with specified parameters and return metrics
    
    Args:
        trial_id: Optuna trial ID
        params: Dictionary of QC parameters
    
    Returns:
        Dictionary of metrics
    """
    # Create trial directory
    trial_dir = os.path.join(Config.TEMP_DIR, f"trial_{trial_id}")
    os.makedirs(trial_dir, exist_ok=True)
    os.chdir(trial_dir)
    
    # Set paths for this trial
    trial_prefix = f"trial_{trial_id}"
    subsample_prefix = f"../subsample"
    
    # Create symbolic link to the subsample data
    for ext in ['.bed', '.bim', '.fam']:
        if os.path.exists(f"{subsample_prefix}{ext}"):
            if not os.path.exists(f"{trial_prefix}_raw{ext}"):
                os.symlink(f"{subsample_prefix}{ext}", f"{trial_prefix}_raw{ext}")
    
    # Extract original sample IDs for later comparison
    cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_raw --write-samples --out original_samples"
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    original_samples = pd.read_csv("original_samples.samples", sep='\t')['ID_2'].tolist()
    
    # Get initial counts
    initial_samples = count_samples(f"{trial_prefix}_raw")
    initial_snps = count_snps(f"{trial_prefix}_raw")
    initial_cases, initial_controls = get_case_control_counts(f"{trial_prefix}_raw", Config.PHENOTYPE_FILE)
    
    # Store metrics
    metrics = {
        'initial_samples': initial_samples,
        'initial_snps': initial_snps,
        'initial_cases': initial_cases,
        'initial_controls': initial_controls,
        'initial_case_control_ratio': initial_cases / initial_controls if initial_controls > 0 else 0
    }
    
    try:
        # STEP 1: Sample call rate
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_raw --mind {params['mind_threshold']} --make-bed --out {trial_prefix}_mind"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # STEP 2: Sex check (simplified for optimization)
        # In a real run, you'd handle sex mismatches properly
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_mind --make-bed --out {trial_prefix}_sexchecked"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # STEP 3: Heterozygosity check
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_sexchecked --het --out {trial_prefix}_het"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Identify het outliers with R script
        het_script = f'''
        het <- read.table("{trial_prefix}_het.het", header=TRUE)
        het$F <- (het$O.HOM - het$E.HOM)/het$E.HOM
        mean_f <- mean(het$F)
        sd_f <- sd(het$F)
        outliers <- subset(het, F < mean_f-{params['het_sd']}*sd_f | F > mean_f+{params['het_sd']}*sd_f)
        write.table(outliers[,c(1,2)], "{trial_prefix}_het_outliers.txt", 
                  row.names=FALSE, col.names=FALSE, quote=FALSE)
        '''
        
        with open("het_script.R", "w") as f:
            f.write(het_script)
        
        subprocess.run("Rscript het_script.R", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_sexchecked --remove {trial_prefix}_het_outliers.txt --make-bed --out {trial_prefix}_het_filtered"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # STEP 4: Relatedness check
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_het_filtered --indep-pairwise {params['ld_window']} {params['ld_step']} {params['ld_r2']} --out {trial_prefix}_pruned"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_het_filtered --extract {trial_prefix}_pruned.prune.in --genome --min {params['ibd_threshold']} --out {trial_prefix}_ibd"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Custom R script to select one from each related pair
        ibd_script = f'''
        if(file.exists("{trial_prefix}_ibd.genome")) {{
          related <- read.table("{trial_prefix}_ibd.genome", header=TRUE)
          if(nrow(related) > 0) {{
            samples_to_remove <- c()
            for(i in 1:nrow(related)) {{
              id1 <- paste(related[i,"FID1"], related[i,"IID1"], sep=":")
              id2 <- paste(related[i,"FID2"], related[i,"IID2"], sep=":")
              if(!(id1 %in% samples_to_remove) && !(id2 %in% samples_to_remove)) {{
                samples_to_remove <- c(samples_to_remove, id2)
              }}
            }}
            ids <- do.call(rbind, strsplit(samples_to_remove, ":"))
            write.table(ids, "{trial_prefix}_related_to_remove.txt", 
                      row.names=FALSE, col.names=FALSE, quote=FALSE)
          }} else {{
            # Create empty file if no related samples
            file.create("{trial_prefix}_related_to_remove.txt")
          }}
        }} else {{
          # Create empty file if no genome file
          file.create("{trial_prefix}_related_to_remove.txt")
        }}
        '''
        
        with open("ibd_script.R", "w") as f:
            f.write(ibd_script)
        
        subprocess.run("Rscript ibd_script.R", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Check if related_to_remove file exists and has content
        if os.path.exists(f"{trial_prefix}_related_to_remove.txt") and os.path.getsize(f"{trial_prefix}_related_to_remove.txt") > 0:
            cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_het_filtered --remove {trial_prefix}_related_to_remove.txt --make-bed --out {trial_prefix}_unrelated"
        else:
            cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_het_filtered --make-bed --out {trial_prefix}_unrelated"
        
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # For optimization purposes, skip PCA and population structure steps
        # In a real analysis, you would include these steps
        
        # STEP 5: Variant QC - call rate
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_unrelated --geno {params['geno_threshold']} --make-bed --out {trial_prefix}_geno"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # STEP 6: MAF filtering
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_geno --maf {params['maf_threshold']} --make-bed --out {trial_prefix}_maf"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # STEP 7: HWE filtering (in controls only)
        # Create controls file
        with open("controls.txt", "w") as f:
            f.write("CONTROL\n")
        
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_maf --hwe {params['hwe_threshold']} --make-bed --out {trial_prefix}_hwe"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # STEP 8: Differential missingness
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_hwe --test-missing --out {trial_prefix}_diffmiss"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Extract SNPs with differential missingness
        diff_miss_script = f'''
        if(file.exists("{trial_prefix}_diffmiss.missing")) {{
          diffmiss <- read.table("{trial_prefix}_diffmiss.missing", header=TRUE)
          sig_snps <- subset(diffmiss, P < {params['diffmiss_threshold']})
          write.table(sig_snps$SNP, "{trial_prefix}_diffmiss_snps.txt", 
                    row.names=FALSE, col.names=FALSE, quote=FALSE)
        }} else {{
          # Create empty file if missing file doesn't exist
          file.create("{trial_prefix}_diffmiss_snps.txt")
        }}
        '''
        
        with open("diff_miss_script.R", "w") as f:
            f.write(diff_miss_script)
        
        subprocess.run("Rscript diff_miss_script.R", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Remove SNPs with differential missingness
        if os.path.exists(f"{trial_prefix}_diffmiss_snps.txt") and os.path.getsize(f"{trial_prefix}_diffmiss_snps.txt") > 0:
            cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_hwe --exclude {trial_prefix}_diffmiss_snps.txt --make-bed --out {trial_prefix}_QC_complete"
        else:
            cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_hwe --make-bed --out {trial_prefix}_QC_complete"
        
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Run a quick GWAS to evaluate parameters
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_QC_complete --assoc --adjust --out {trial_prefix}_gwas"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Get final counts
        final_samples = count_samples(f"{trial_prefix}_QC_complete")
        final_snps = count_snps(f"{trial_prefix}_QC_complete")
        final_cases, final_controls = get_case_control_counts(f"{trial_prefix}_QC_complete", Config.PHENOTYPE_FILE)
        
        # Extract final sample IDs for comparison
        cmd = f"{Config.PLINK_PATH} --bfile {trial_prefix}_QC_complete --write-samples --out final_samples"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        final_sample_ids = pd.read_csv("final_samples.samples", sep='\t')['ID_2'].tolist()
        
        # Calculate metrics
        lambda_gc = calculate_lambda_gc(f"{trial_prefix}_gwas.assoc")
        
        # Calculate covariate distribution changes
        if os.path.exists(Config.COVAR_FILE):
            cat_cols = ['sex', 'cancer_subtype']  # Adjust based on your actual categorical covariates
            covar_change = compare_covariate_distribution(original_samples, final_sample_ids, Config.COVAR_FILE, categorical_cols=cat_cols)
        else:
            covar_change = 0
        
        if os.path.exists(Config.QCOVAR_FILE):
            quant_cols = ['age', 'diagnosis_age']  # Adjust based on your actual quantitative covariates
            qcovar_change = compare_covariate_distribution(original_samples, final_sample_ids, Config.QCOVAR_FILE, quantitative_cols=quant_cols)
        else:
            qcovar_change = 0
        
        # Update metrics
        metrics.update({
            'final_samples': final_samples,
            'final_snps': final_snps,
            'final_cases': final_cases,
            'final_controls': final_controls,
            'final_case_control_ratio': final_cases / final_controls if final_controls > 0 else 0,
            'sample_retention': final_samples / initial_samples if initial_samples > 0 else 0,
            'snp_retention': final_snps / initial_snps if initial_snps > 0 else 0,
            'case_retention': final_cases / initial_cases if initial_cases > 0 else 0,
            'control_retention': final_controls / initial_controls if initial_controls > 0 else 0,
            'case_control_ratio_change': abs((final_cases / final_controls) - (initial_cases / initial_controls)) / (initial_cases / initial_controls) if initial_controls > 0 and final_controls > 0 else float('inf'),
            'lambda_gc': lambda_gc,
            'covariate_change': covar_change,
            'qcovariate_change': qcovar_change
        })
        
        # Save trial metrics to file
        pd.DataFrame([metrics]).to_csv(f"{trial_prefix}_metrics.csv", index=False)
        
        # Create a plot of the QC process
        plot_qc_flowchart(trial_id, metrics)
        
        # Return to original directory
        os.chdir('../..')
        
        return metrics
    
    except Exception as e:
        logger.error(f"Error in trial {trial_id}: {e}")
        os.chdir('../..')
        return {'error': str(e)}

def plot_qc_flowchart(trial_id, metrics):
    """Create a flowchart of the QC process with sample/SNP counts"""
    plt.figure(figsize=(10, 8))
    
    # Create text labels for each step
    steps = [
        f"Raw data\n{metrics['initial_samples']} samples\n{metrics['initial_snps']} SNPs",
        f"After sample call rate filter\n(mind = {params['mind_threshold']})",
        f"After heterozygosity filter\n(±{params['het_sd']} SD)",
        f"After relatedness filter\n(IBD > {params['ibd_threshold']})",
        f"After variant call rate filter\n(geno = {params['geno_threshold']})",
        f"After MAF filter\n(MAF > {params['maf_threshold']})",
        f"After HWE filter\n(p > {params['hwe_threshold']})",
        f"Final QC'd data\n{metrics['final_samples']} samples\n{metrics['final_snps']} SNPs"
    ]
    
    # Create flow chart
    for i, step in enumerate(steps):
        plt.text(0.5, 1 - i*0.14, step, ha='center', va='center', bbox=dict(facecolor='lightblue', alpha=0.5))
        if i < len(steps) - 1:
            plt.arrow(0.5, 1 - i*0.14 - 0.05, 0, -0.04, head_width=0.02, head_length=0.01, fc='black', ec='black')
    
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    plt.title(f"QC Flowchart - Trial {trial_id}")
    plt.tight_layout()
    plt.savefig(f"trial_{trial_id}_flowchart.png", dpi=300)
    plt.close()

def objective(trial):
    """Optuna objective function for GWAS QC parameter optimization"""
    # Generate parameters for this trial
    params = {
        'mind_threshold': trial.suggest_float('mind_threshold', 0.01, 0.10),
        'het_sd': trial.suggest_float('het_sd', 2.0, 4.0),
        'ibd_threshold': trial.suggest_float('ibd_threshold', 0.1, 0.25),
        'ld_window': trial.suggest_int('ld_window', 25, 200),
        'ld_step': trial.suggest_int('ld_step', 2, 20),
        'ld_r2': trial.suggest_float('ld_r2', 0.1, 0.5),
        'geno_threshold': trial.suggest_float('geno_threshold', 0.01, 0.10),
        'maf_threshold': trial.suggest_float('maf_threshold', 0.001, 0.05),
        'hwe_threshold': trial.suggest_float('hwe_threshold', 1e-10, 1e-4, log=True),
        'diffmiss_threshold': trial.suggest_float('diffmiss_threshold', 1e-6, 1e-3, log=True)
    }
    
    # Log parameters
    logger.info(f"Trial {trial.number}: {params}")
    
    # Run GWAS QC with these parameters
    metrics = run_gwas_qc(trial.number, params)
    
    if 'error' in metrics:
        # If error occurred, return a very bad score
        return float('-inf')
    
    # Calculate optimization score - adjust weights based on your priorities
    score = (
        0.25 * metrics['sample_retention'] +                     # Higher retention is better (25%)
        0.25 * metrics['snp_retention'] +                        # Higher SNP retention is better (25%)
        -0.20 * abs(metrics['lambda_gc'] - 1.0 if metrics['lambda_gc'] else 10) +  # Lambda closer to 1 is better (20%)
        -0.15 * metrics['case_control_ratio_change'] +           # Minimize case/control ratio change (15%)
        -0.10 * (metrics['covariate_change'] + metrics['qcovariate_change']) +  # Minimize covariate distribution changes (10%)
        0.05 * min(metrics['case_retention'], metrics['control_retention'])     # Balance case/control retention (5%)
    )
    
    # Log score
    logger.info(f"Trial {trial.number} score: {score}")
    
    return score

def save_results(study):
    """Save and visualize optimization results"""
    os.makedirs(Config.RESULTS_DIR, exist_ok=True)
    
    # Save study statistics
    df = study.trials_dataframe()
    df.to_csv(os.path.join(Config.RESULTS_DIR, "all_trials.csv"), index=False)
    
    # Get best parameters
    best_params = study.best_params
    best_value = study.best_value
    
    # Save best parameters
    with open(os.path.join(Config.RESULTS_DIR, "best_params.txt"), "w") as f:
        f.write(f"Best Score: {best_value}\n\n")
        f.write("Best Parameters:\n")
        for param, value in best_params.items():
            f.write(f"{param}: {value}\n")
    
    # Create parameter importance plot
    try:
        param_importances = optuna.importance.get_param_importances(study)
        plt.figure(figsize=(10, 6))
        importance_df = pd.DataFrame(
            {'Parameter': list(param_importances.keys()), 
             'Importance': list(param_importances.values())}
        ).sort_values('Importance', ascending=False)
        
        sns.barplot(x='Importance', y='Parameter', data=importance_df)
        plt.title('Parameter Importance')
        plt.tight_layout()
        plt.savefig(os.path.join(Config.RESULTS_DIR, "parameter_importance.png"), dpi=300)
        plt.close()
    except:
        logger.warning("Could not generate parameter importance plot")
    
    # Create optimization history plot
    plt.figure(figsize=(10, 6))
    optuna.visualization.matplotlib.plot_optimization_history(study)
    plt.tight_layout()
    plt.savefig(os.path.join(Config.RESULTS_DIR, "optimization_history.png"), dpi=300)
    plt.close()
    
    # Create parallel coordinate plot
    plt.figure(figsize=(12, 8))
    optuna.visualization.matplotlib.plot_parallel_coordinate(study)
    plt.tight_layout()
    plt.savefig(os.path.join(Config.RESULTS_DIR, "parallel_coordinate.png"), dpi=300)
    plt.close()

    # Save best trial metrics
    best_trial_metrics_file = os.path.join(Config.TEMP_DIR, f"trial_{study.best_trial.number}", 
                                         f"trial_{study.best_trial.number}_metrics.csv")
    if os.path.exists(best_trial_metrics_file):
        best_metrics = pd.read_csv(best_trial_metrics_file)
        best_metrics.to_csv(os.path.join(Config.RESULTS_DIR, "best_trial_metrics.csv"), index=False)
    
    logger.info(f"Results saved to {Config.RESULTS_DIR}")
    logger.info(f"Best parameters: {best_params}")
    logger.info(f"Best score: {best_value}")

def create_optimized_script(best_params):
    """Create a bash script with the optimized parameters"""
    script_content = f"""#!/bin/bash
# Optimized GWAS QC parameters
# Generated by Optuna optimization on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

# Set paths to your tools and data
PLINK="{Config.PLINK_PATH}"
GCTA="{Config.GCTA_PATH}"
RAW_DATA="{Config.RAW_DATA}"
OUT_PREFIX="{Config.OUT_PREFIX}"

echo "Starting Sample QC with optimized parameters..."

# Step 1: Sample call rate filter (optimized)
$PLINK --bfile $RAW_DATA \\
    --mind {best_params['mind_threshold']:.6f} \\
    --make-bed \\
    --out ${{OUT_PREFIX}}_mind

# Step 2: Sex check - identify mismatches
$PLINK --bfile ${{OUT_PREFIX}}_mind \\
    --check-sex \\
    --out ${{OUT_PREFIX}}_sexcheck

# Remove samples with sex discrepancies
awk '{{ if ($5 == "PROBLEM") print $1,$2 }}' ${{OUT_PREFIX}}_sexcheck.sexcheck > ${{OUT_PREFIX}}_sex_discrepancy.txt

$PLINK --bfile ${{OUT_PREFIX}}_mind \\
    --remove ${{OUT_PREFIX}}_sex_discrepancy.txt \\
    --make-bed \\
    --out ${{OUT_PREFIX}}_sexchecked

# Step 3: Check for excessive heterozygosity (optimized threshold)
$PLINK --bfile ${{OUT_PREFIX}}_sexchecked \\
    --het \\
    --out ${{OUT_PREFIX}}_het

# Identify samples with outlier heterozygosity (optimized SD threshold)
Rscript -e '
  het <- read.table("'${{OUT_PREFIX}}'_het.het", header=TRUE);
  het$F <- (het$O.HOM - het$E.HOM)/het$E.HOM;
  mean_f <- mean(het$F);
  sd_f <- sd(het$F);
  outliers <- subset(het, F < mean_f-{best_params['het_sd']:.6f}*sd_f | F > mean_f+{best_params['het_sd']:.6f}*sd_f);
  write.table(outliers[,c(1,2)], "'${{OUT_PREFIX}}'_het_outliers.txt", 
              row.names=FALSE, col.names=FALSE, quote=FALSE);'

# Remove heterozygosity outliers
$PLINK --bfile ${{OUT_PREFIX}}_sexchecked \\
    --remove ${{OUT_PREFIX}}_het_outliers.txt \\
    --make-bed \\
    --out ${{OUT_PREFIX}}_het_filtered

# Step 4: Check for relatedness using Identity-By-Descent (optimized parameters)
$PLINK --bfile ${{OUT_PREFIX}}_het_filtered \\
    --indep-pairwise {best_params['ld_window']} {best_params['ld_step']} {best_params['ld_r2']:.6f} \\
    --out ${{OUT_PREFIX}}_pruned

$PLINK --bfile ${{OUT_PREFIX}}_het_filtered \\
    --extract ${{OUT_PREFIX}}_pruned.prune.in \\
    --genome \\
    --min {best_params['ibd_threshold']:.6f} \\
    --out ${{OUT_PREFIX}}_ibd

# Keep one individual from each related pair
awk '{{ if ($9 > {best_params['ibd_threshold']:.6f}) print $1,$2,$3,$4,$9,$10 }}' ${{OUT_PREFIX}}_ibd.genome > ${{OUT_PREFIX}}_related_pairs.txt

# Custom script to select one from each related cluster
Rscript -e '
  related <- read.table("'${{OUT_PREFIX}}'_related_pairs.txt", header=FALSE);
  samples_to_remove <- c();
  for(i in 1:nrow(related)) {{
    id1 <- paste(related[i,1], related[i,2], sep=":");
    id2 <- paste(related[i,3], related[i,4], sep=":");
    if(!(id1 %in% samples_to_remove) && !(id2 %in% samples_to_remove)) {{
      samples_to_remove <- c(samples_to_remove, id2);
    }}
  }}
  ids <- do.call(rbind, strsplit(samples_to_remove, ":"));
  write.table(ids, "'${{OUT_PREFIX}}'_related_to_remove.txt", 
              row.names=FALSE, col.names=FALSE, quote=FALSE);'

$PLINK --bfile ${{OUT_PREFIX}}_het_filtered \\
    --remove ${{OUT_PREFIX}}_related_to_remove.txt \\
    --make-bed \\
    --out ${{OUT_PREFIX}}_unrelated

# Run PCA to identify population outliers
$GCTA --bfile ${{OUT_PREFIX}}_unrelated \\
    --make-grm \\
    --out ${{OUT_PREFIX}}_grm

$GCTA --grm ${{OUT_PREFIX}}_grm \\
    --pca 20 \\
    --out ${{OUT_PREFIX}}_pca

echo "Completed sample QC. Review PCA results and create pca_outliers.txt before continuing."
echo "After creating pca_outliers.txt, run the next command:"
echo "$PLINK --bfile ${{OUT_PREFIX}}_unrelated --remove pca_outliers.txt --make-bed --out ${{OUT_PREFIX}}_sample_qc"

echo ""
echo "Then proceed with Variant QC:"

echo "# STEP 5-8: VARIANT QC WITH OPTIMIZED PARAMETERS"
echo "$PLINK --bfile ${{OUT_PREFIX}}_sample_qc \\"
echo "    --geno {best_params['geno_threshold']:.6f} \\"
echo "    --make-bed \\"
echo "    --out ${{OUT_PREFIX}}_geno"
echo ""
echo "$PLINK --bfile ${{OUT_PREFIX}}_geno \\"
echo "    --maf {best_params['maf_threshold']:.6f} \\"
echo "    --make-bed \\"
echo "    --out ${{OUT_PREFIX}}_maf"
echo ""
echo "$PLINK --bfile ${{OUT_PREFIX}}_maf \\"
echo "    --hwe {best_params['hwe_threshold']:.6e} \\"
echo "    --filter controls.txt \\"
echo "    --make-bed \\"
echo "    --out ${{OUT_PREFIX}}_hwe"
echo ""
echo "$PLINK --bfile ${{OUT_PREFIX}}_hwe \\"
echo "    --test-missing \\"
echo "    --out ${{OUT_PREFIX}}_diffmiss"
echo ""
echo "awk '{{ if ($5 < {best_params['diffmiss_threshold']:.6e}) print $2 }}' ${{OUT_PREFIX}}_diffmiss.missing > ${{OUT_PREFIX}}_diffmiss_snps.txt"
echo ""
echo "$PLINK --bfile ${{OUT_PREFIX}}_hwe \\"
echo "    --exclude ${{OUT_PREFIX}}_diffmiss_snps.txt \\"
echo "    --make-bed \\"
echo "    --out ${{OUT_PREFIX}}_QC_complete"
echo ""
echo "echo \"Sample and Variant QC completed with optimized parameters.\""
echo "echo \"Final cleaned dataset: ${{OUT_PREFIX}}_QC_complete\""
"""
    
    # Save script
    script_path = os.path.join(Config.RESULTS_DIR, "optimized_qc_pipeline.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(script_path, 0o755)
    
    logger.info(f"Optimized script saved to {script_path}")
    return script_path

def main():
    """Main function to run the optimization"""
    # Create directories
    os.makedirs(Config.TEMP_DIR, exist_ok=True)
    os.makedirs(Config.RESULTS_DIR, exist_ok=True)
    
    logger.info("Starting GWAS QC parameter optimization")
    logger.info(f"Using PLINK at {Config.PLINK_PATH}")
    logger.info(f"Using raw data: {Config.RAW_DATA}")
    
    # Create subsample for faster optimization
    logger.info(f"Creating subsample with {Config.SUBSAMPLE_SIZE} samples")
    success = create_subsample(Config.RAW_DATA, f"{Config.TEMP_DIR}/subsample", Config.SUBSAMPLE_SIZE, Config.SUBSAMPLE_SEED)
    
    if not success:
        logger.error("Failed to create subsample. Exiting.")
        return
    
    # Create or load Optuna study
    logger.info(f"Creating/loading Optuna study from database {Config.DB_PATH}")
    storage = optuna.storages.RDBStorage(
        url=Config.DB_PATH,
        engine_kwargs={"connect_args": {"timeout": 10}}
    )
    
    study = optuna.create_study(
        study_name=Config.STUDY_NAME,
        direction="maximize",  # Higher score is better
        sampler=optuna.samplers.TPESampler(seed=Config.SUBSAMPLE_SEED),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=5),
        storage=storage,
        load_if_exists=True  # Allow resuming previous study
    )
    
    # Check if study already has trials
    existing_trials = len(study.trials)
    if existing_trials > 0:
        logger.info(f"Found existing study with {existing_trials} completed trials")
        logger.info(f"Best score so far: {study.best_value}")
        
        # Ask if user wants to continue or start fresh
        if input("Continue existing study? (y/n): ").lower() != 'y':
            logger.info("Creating new study")
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            new_study_name = f"{Config.STUDY_NAME}_{timestamp}"
            study = optuna.create_study(
                study_name=new_study_name,
                direction="maximize",
                sampler=optuna.samplers.TPESampler(seed=Config.SUBSAMPLE_SEED),
                pruner=optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=5),
                storage=storage
            )
            logger.info(f"Created new study '{new_study_name}'")
    
    # Number of trials to run (remaining if continuing)
    n_remaining = max(0, Config.N_TRIALS - len(study.trials))
    if n_remaining == 0 and existing_trials > 0:
        logger.info("Target number of trials already completed. You can increase N_TRIALS to run more.")
        n_remaining = int(input("Enter additional trials to run (0 to skip optimization): "))
    
    # Optimize parameters
    if n_remaining > 0:
        logger.info(f"Starting optimization with {n_remaining} trials")
        
        # Use joblib for parallel execution if configured with n_jobs > 1
        if Config.N_JOBS > 1:
            logger.info(f"Using parallel execution with {Config.N_JOBS} jobs")
            study.optimize(
                objective, 
                n_trials=n_remaining, 
                timeout=Config.TIMEOUT,
                n_jobs=Config.N_JOBS,
                gc_after_trial=True  # Garbage collect to prevent memory issues
            )
        else:
            study.optimize(objective, n_trials=n_remaining, timeout=Config.TIMEOUT)
    
    # Save and visualize results
    logger.info("Optimization complete. Saving results.")
    save_results(study)
    
    # Create optimized script
    script_path = create_optimized_script(study.best_params)
    
    # Export study for easier analysis
    export_study_data(study)
    
    logger.info("Optimization pipeline completed successfully.")
    logger.info(f"Run the optimized pipeline with: {script_path}")
    
    # Summary
    print("\n" + "="*50)
    print("GWAS QC PARAMETER OPTIMIZATION COMPLETE")
    print("="*50)
    print(f"Best parameters saved to: {os.path.join(Config.RESULTS_DIR, 'best_params.txt')}")
    print(f"Optimized script saved to: {script_path}")
    print(f"Visualization plots saved to: {Config.RESULTS_DIR}")
    print(f"Database stored at: {Config.DB_PATH.replace('sqlite:///', '')}")
    print("\nTo analyze the database, use:")
    print("  python gwas_qc_optimizer.py analyze")
    print("="*50 + "\n")

def export_study_data(study):
    """Export study data to CSV and joblib for easier analysis"""
    # Export trials data to CSV
    trials_df = study.trials_dataframe()
    trials_df.to_csv(os.path.join(Config.RESULTS_DIR, "all_trials.csv"), index=False)
    
    # Export study object using joblib
    joblib.dump(study, os.path.join(Config.RESULTS_DIR, "optuna_study.joblib"))
    
    # Export trial metrics data if available
    metrics_data = []
    for trial in study.trials:
        if trial.state == optuna.trial.TrialState.COMPLETE:
            trial_metrics_file = os.path.join(Config.TEMP_DIR, f"trial_{trial.number}", 
                                            f"trial_{trial.number}_metrics.csv")
            if os.path.exists(trial_metrics_file):
                trial_metrics = pd.read_csv(trial_metrics_file)
                trial_metrics['trial_number'] = trial.number
                trial_metrics['value'] = trial.value
                metrics_data.append(trial_metrics)
    
    if metrics_data:
        all_metrics = pd.concat(metrics_data, ignore_index=True)
        all_metrics.to_csv(os.path.join(Config.RESULTS_DIR, "all_trial_metrics.csv"), index=False)

def analyze_database():
    """Analyze the Optuna database and generate insights"""
    # Check if database exists
    db_file = Config.DB_PATH.replace('sqlite:///', '')
    if not os.path.exists(db_file):
        print(f"Database file {db_file} not found.")
        return
    
    # Create connection to SQLite database
    conn = sqlite3.connect(db_file)
    
    # Get list of studies
    studies_df = pd.read_sql("SELECT study_name, COUNT(*) as trials FROM trials GROUP BY study_name", conn)
    
    print("\n" + "="*50)
    print("OPTUNA DATABASE ANALYSIS")
    print("="*50)
    
    print("\nStudies in database:")
    print(studies_df)
    
    # Let user select a study
    if len(studies_df) > 1:
        print("\nSelect a study to analyze:")
        for i, study_name in enumerate(studies_df['study_name']):
            print(f"{i+1}. {study_name} ({studies_df.iloc[i]['trials']} trials)")
        
        selection = int(input("\nEnter number: ")) - 1
        selected_study = studies_df.iloc[selection]['study_name']
    else:
        selected_study = studies_df.iloc[0]['study_name']
    
    print(f"\nAnalyzing study: {selected_study}")
    
    # Load study
    study = optuna.load_study(study_name=selected_study, storage=Config.DB_PATH)
    
    # Get basic statistics
    print("\nBasic Statistics:")
    print(f"Total trials: {len(study.trials)}")
    print(f"Completed trials: {sum(1 for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE)}")
    print(f"Failed trials: {sum(1 for t in study.trials if t.state == optuna.trial.TrialState.FAIL)}")
    print(f"Best value: {study.best_value}")
    
    # Show best parameters
    print("\nBest Parameters:")
    for param, value in study.best_params.items():
        print(f"  {param}: {value}")
    
    # Generate and save additional plots
    os.makedirs("db_analysis", exist_ok=True)
    
    print("\nGenerating analysis plots in 'db_analysis' directory...")
    
    # Optimization history
    plt.figure(figsize=(10, 6))
    optuna.visualization.matplotlib.plot_optimization_history(study)
    plt.tight_layout()
    plt.savefig("db_analysis/optimization_history.png", dpi=300)
    
    # Parameter importance
    try:
        plt.figure(figsize=(10, 6))
        optuna.visualization.matplotlib.plot_param_importances(study)
        plt.tight_layout()
        plt.savefig("db_analysis/param_importances.png", dpi=300)
    except:
        print("Could not generate parameter importance plot")
    
    # Parallel coordinate plot
    try:
        plt.figure(figsize=(12, 8))
        optuna.visualization.matplotlib.plot_parallel_coordinate(study)
        plt.tight_layout()
        plt.savefig("db_analysis/parallel_coordinate.png", dpi=300)
    except:
        print("Could not generate parallel coordinate plot")
    
    # Slice plot
    try:
        plt.figure(figsize=(12, 8))
        optuna.visualization.matplotlib.plot_slice(study)
        plt.tight_layout()
        plt.savefig("db_analysis/slice_plot.png", dpi=300)
    except:
        print("Could not generate slice plot")
    
    # Contour plot for top two important parameters
    try:
        param_importance = optuna.importance.get_param_importances(study)
        top_params = list(param_importance.keys())[:2]
        
        if len(top_params) >= 2:
            plt.figure(figsize=(10, 8))
            optuna.visualization.matplotlib.plot_contour(study, params=top_params)
            plt.tight_layout()
            plt.savefig("db_analysis/contour_plot.png", dpi=300)
    except:
        print("Could not generate contour plot")
    
    # Export trials data for manual analysis
    trials_df = study.trials_dataframe()
    trials_df.to_csv("db_analysis/trials_data.csv", index=False)
    
    print("\nAnalysis complete. Results saved to 'db_analysis' directory.")
    print("\nFor detailed analysis, you can use the trials data in 'db_analysis/trials_data.csv'")
    
    # Close connection
    conn.close()

if __name__ == "__main__":
    # Check for command line arguments
    if len(sys.argv) > 1 and sys.argv[1] == "analyze":
        analyze_database()
    else:
        main()