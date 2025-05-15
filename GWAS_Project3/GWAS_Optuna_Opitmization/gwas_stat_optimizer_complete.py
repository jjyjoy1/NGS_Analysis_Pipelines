#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GWAS Statistical Analysis Parameter Optimization with Optuna

This script optimizes GWAS statistical analysis parameters using Optuna,
balancing genomic control, statistical power, and discovery potential
while supporting cross-validation and known association validation.
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
from sklearn.model_selection import StratifiedKFold
from scipy import stats
from datetime import datetime
from optuna.visualization import plot_optimization_history, plot_param_importances

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("gwas_optimizer.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("GWAS-Optimizer")

class GWASOptimizer:
    """Class for optimizing GWAS statistical parameters with Optuna."""
    
    def __init__(self, args):
        """Initialize the optimizer with command line arguments."""
        self.args = args
        self.setup_paths()
        self.validate_inputs()
        
        # Create study directory if it doesn't exist
        os.makedirs(self.study_dir, exist_ok=True)
        
        # Initialize Optuna study
        self.setup_optuna_study()
        
        # Load known associations if provided
        self.known_associations = None
        if self.args.known_associations:
            self.known_associations = pd.read_csv(self.args.known_associations)
            logger.info(f"Loaded {len(self.known_associations)} known associations")
    
    def setup_paths(self):
        """Set up paths for input and output files."""
        # Base file names without extensions
        self.base_prefix = os.path.basename(self.args.plink_prefix)
        
        # Study directory for output
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.study_dir = os.path.join(self.args.output_dir, f"gwas_opt_{timestamp}")
        
        # Paths for PLINK and GCTA
        self.plink_path = self.args.plink_path if self.args.plink_path else "plink"
        self.gcta_path = self.args.gcta_path if self.args.gcta_path else "gcta64"
        
        # Full paths to input files
        self.plink_prefix = self.args.plink_prefix
        
        # Output database for Optuna
        self.db_path = os.path.join(self.study_dir, "optuna_trials.db")
    
    def validate_inputs(self):
        """Validate input files and dependencies."""
        required_plink_files = [f"{self.plink_prefix}.bed", 
                               f"{self.plink_prefix}.bim", 
                               f"{self.plink_prefix}.fam"]
        
        for file_path in required_plink_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Required PLINK file not found: {file_path}")
        
        # Validate PLINK and GCTA installations
        try:
            subprocess.run([self.plink_path, "--help"], 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE)
        except:
            logger.warning("PLINK executable not found. Please check path.")
        
        # Only check for GCTA if we're going to use it
        if not self.args.disable_mlma:
            try:
                subprocess.run([self.gcta_path, "--help"], 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE)
            except:
                logger.warning("GCTA executable not found. Some optimization options may be unavailable.")
        
        # Validate phenotype file if provided
        if self.args.pheno_file and not os.path.exists(self.args.pheno_file):
            raise FileNotFoundError(f"Phenotype file not found: {self.args.pheno_file}")
        
        # Validate covariate file if provided
        if self.args.covar_file and not os.path.exists(self.args.covar_file):
            raise FileNotFoundError(f"Covariate file not found: {self.args.covar_file}")
        
        # Validate known associations file if provided
        if self.args.known_associations and not os.path.exists(self.args.known_associations):
            raise FileNotFoundError(f"Known associations file not found: {self.args.known_associations}")
        
        logger.info("All input validations passed")
    
    def setup_optuna_study(self):
        """Set up Optuna study for parameter optimization."""
        storage_url = f"sqlite:///{self.db_path}"
        
        # Check if we need to load existing study
        if self.args.resume_study and os.path.exists(self.db_path):
            self.study = optuna.load_study(
                study_name=self.args.study_name,
                storage=storage_url
            )
            logger.info(f"Resuming existing study '{self.args.study_name}' with {len(self.study.trials)} previous trials")
        else:
            # Create new study
            self.study = optuna.create_study(
                study_name=self.args.study_name,
                storage=storage_url,
                direction="maximize",  # We'll define a custom objective that should be maximized
                load_if_exists=False
            )
            logger.info(f"Created new optimization study: {self.args.study_name}")
    
    def define_parameter_space(self, trial):
        """Define the parameter space for Optuna to search."""
        params = {}
        
        # 1. Analysis method selection
        available_methods = []
        if not self.args.disable_assoc:
            available_methods.append("assoc")
        if not self.args.disable_logistic:
            available_methods.append("logistic")
        if not self.args.disable_mlma:
            available_methods.append("mlma")
        
        # Make sure we have at least one method
        if not available_methods:
            raise ValueError("All analysis methods are disabled. Enable at least one.")
        
        params["analysis_method"] = trial.suggest_categorical("analysis_method", available_methods)
        
        # 2. Covariate handling
        if self.args.covar_file:
            params["use_covariates"] = trial.suggest_categorical("use_covariates", [True, False])
            
            # Only suggest PC count if we're using covariates
            if params["use_covariates"]:
                max_pcs = min(20, self.args.max_pcs)
                params["n_pcs"] = trial.suggest_int("n_pcs", 0, max_pcs)
            else:
                params["n_pcs"] = 0
        else:
            params["use_covariates"] = False
            params["n_pcs"] = 0
        
        # 3. Multiple testing correction
        params["multiple_testing"] = trial.suggest_categorical(
            "multiple_testing", 
            ["bonferroni", "fdr", "none"]
        )
        
        # 4. Method-specific parameters
        if params["analysis_method"] == "mlma":
            params["use_loco"] = trial.suggest_categorical("use_loco", [True, False])
        elif params["analysis_method"] == "logistic":
            params["genotype_certainty"] = trial.suggest_categorical(
                "genotype_certainty", 
                ["dosage", "hard_call"]
            )
        
        # 5. MAF threshold
        params["maf_threshold"] = trial.suggest_float("maf_threshold", 0.001, 0.05)
        
        # 6. Missing rate threshold
        params["geno_threshold"] = trial.suggest_float("geno_threshold", 0.01, 0.1)
        
        return params
    
    def run_gwas(self, params, fold=None):
        """
        Run GWAS analysis with the given parameters.
        
        Args:
            params: Dictionary of GWAS parameters
            fold: Optional fold index for cross-validation
        
        Returns:
            Path to the results file
        """
        # Generate output prefix for this run
        if fold is not None:
            output_prefix = os.path.join(self.study_dir, 
                                        f"trial_{self.current_trial_number}_fold_{fold}")
        else:
            output_prefix = os.path.join(self.study_dir, 
                                        f"trial_{self.current_trial_number}")
        
        # Create command for the selected analysis method
        if params["analysis_method"] == "assoc":
            cmd = self.create_plink_assoc_command(params, output_prefix)
        elif params["analysis_method"] == "logistic":
            cmd = self.create_plink_logistic_command(params, output_prefix)
        elif params["analysis_method"] == "mlma":
            cmd = self.create_gcta_mlma_command(params, output_prefix)
        else:
            raise ValueError(f"Unknown analysis method: {params['analysis_method']}")
        
        # Run the command
        logger.info(f"Running GWAS with command: {' '.join(cmd)}")
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        if process.returncode != 0:
            logger.error(f"GWAS run failed: {process.stderr.decode()}")
            return None
        
        # Determine output file based on method
        if params["analysis_method"] == "assoc":
            result_file = f"{output_prefix}.assoc"
        elif params["analysis_method"] == "logistic":
            result_file = f"{output_prefix}.assoc.logistic"
        elif params["analysis_method"] == "mlma":
            result_file = f"{output_prefix}.mlma"
        
        if not os.path.exists(result_file):
            logger.error(f"Expected result file not found: {result_file}")
            return None
        
        return result_file
    
    def create_plink_assoc_command(self, params, output_prefix):
        """Create PLINK command for basic association analysis."""
        cmd = [
            self.plink_path,
            "--bfile", self.plink_prefix,
            "--assoc",
            "--out", output_prefix,
            "--maf", str(params["maf_threshold"]),
            "--geno", str(params["geno_threshold"])
        ]
        
        # Add phenotype file if provided
        if self.args.pheno_file:
            cmd.extend(["--pheno", self.args.pheno_file])
            if self.args.pheno_name:
                cmd.extend(["--pheno-name", self.args.pheno_name])
        
        # Add covariates if enabled
        if params["use_covariates"] and self.args.covar_file:
            cmd.extend(["--covar", self.args.covar_file])
            
            # Add PCs if specified
            if params["n_pcs"] > 0:
                pc_names = ",".join([f"PC{i}" for i in range(1, params["n_pcs"]+1)])
                cmd.extend(["--covar-name", pc_names])
        
        return cmd
    
    def create_plink_logistic_command(self, params, output_prefix):
        """Create PLINK command for logistic regression analysis."""
        cmd = [
            self.plink_path,
            "--bfile", self.plink_prefix,
            "--logistic",
            "--out", output_prefix,
            "--maf", str(params["maf_threshold"]),
            "--geno", str(params["geno_threshold"])
        ]
        
        # Use dosage or hard calls
        if params["genotype_certainty"] == "dosage":
            cmd.append("--dosage")
        
        # Add phenotype file if provided
        if self.args.pheno_file:
            cmd.extend(["--pheno", self.args.pheno_file])
            if self.args.pheno_name:
                cmd.extend(["--pheno-name", self.args.pheno_name])
        
        # Add covariates if enabled
        if params["use_covariates"] and self.args.covar_file:
            cmd.extend(["--covar", self.args.covar_file])
            
            # Add PCs if specified
            if params["n_pcs"] > 0:
                pc_names = ",".join([f"PC{i}" for i in range(1, params["n_pcs"]+1)])
                cmd.extend(["--covar-name", pc_names])
        
        return cmd
    
    def create_gcta_mlma_command(self, params, output_prefix):
        """Create GCTA command for mixed linear model analysis."""
        # First, generate GRM
        grm_prefix = f"{output_prefix}_grm"
        grm_cmd = [
            self.gcta_path,
            "--bfile", self.plink_prefix,
            "--make-grm",
            "--out", grm_prefix,
            "--maf", str(params["maf_threshold"])
        ]
        
        subprocess.run(grm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Then run MLMA
        cmd = [
            self.gcta_path,
            "--mlma",
            "--bfile", self.plink_prefix,
            "--grm", grm_prefix,
            "--out", output_prefix,
            "--maf", str(params["maf_threshold"])
        ]
        
        # Add LOCO if enabled
        if params["use_loco"]:
            cmd.append("--loco")
        
        # Add phenotype file if provided
        if self.args.pheno_file:
            cmd.extend(["--pheno", self.args.pheno_file])
            if self.args.pheno_name:
                cmd.extend(["--pheno-name", self.args.pheno_name])
        
        # Add covariates if enabled
        if params["use_covariates"] and self.args.covar_file:
            cmd.extend(["--qcovar", self.args.covar_file])
            
            # Add PCs if specified
            if params["n_pcs"] > 0:
                # GCTA expects 1-based indexing for columns
                pc_indices = ",".join([str(i) for i in range(1, params["n_pcs"]+1)])
                cmd.extend(["--qcovar-col", pc_indices])
        
        return cmd
    
    def calculate_metrics(self, result_file, params):
        """
        Calculate evaluation metrics from GWAS results.
        
        Args:
            result_file: Path to GWAS results file
            params: Parameters used for the analysis
        
        Returns:
            Dictionary of metrics
        """
        metrics = {}
        
        # Load results file
        if params["analysis_method"] == "assoc":
            results = pd.read_table(result_file, delim_whitespace=True)
            p_col = "P"
        elif params["analysis_method"] == "logistic":
            results = pd.read_table(result_file, delim_whitespace=True)
            # Logistic regression outputs multiple test types, use ADD (additive)
            results = results[results["TEST"] == "ADD"]
            p_col = "P"
        elif params["analysis_method"] == "mlma":
            results = pd.read_table(result_file, delim_whitespace=True)
            p_col = "p"
        
        # Filter out missing p-values
        results = results.dropna(subset=[p_col])
        
        # 1. Calculate genomic inflation factor (lambda GC)
        observed_p = results[p_col].values
        expected_p = np.linspace(0, 1, len(observed_p) + 2)[1:-1]
        
        # Convert to chi-square
        observed_chi2 = stats.chi2.ppf(1 - observed_p, 1)
        expected_chi2 = stats.chi2.ppf(1 - expected_p, 1)
        
        # Calculate lambda (median observed / median expected)
        lambda_gc = np.median(observed_chi2) / stats.chi2.ppf(0.5, 1)
        metrics["lambda_gc"] = lambda_gc
        
        # 2. Count significant hits (applying chosen multiple testing correction)
        n_tests = len(results)
        if params["multiple_testing"] == "bonferroni":
            alpha = 0.05 / n_tests
        elif params["multiple_testing"] == "fdr":
            # Calculate FDR-adjusted p-values
            from statsmodels.stats.multitest import multipletests
            _, adj_p, _, _ = multipletests(results[p_col], method='fdr_bh')
            results["ADJ_P"] = adj_p
            alpha = 0.05  # Use adjusted p-values
            p_col = "ADJ_P"
        else:  # none
            alpha = 5e-8  # Standard GWAS threshold
        
        # Count significant hits
        n_sig = sum(results[p_col] < alpha)
        metrics["n_significant"] = n_sig
        
        # 3. Calculate power score across MAF ranges
        maf_ranges = [(0.01, 0.05), (0.05, 0.1), (0.1, 0.3), (0.3, 0.5)]
        power_scores = []
        
        for maf_min, maf_max in maf_ranges:
            # Filter results by MAF range
            if "MAF" in results.columns:
                maf_results = results[(results["MAF"] >= maf_min) & (results["MAF"] < maf_max)]
            elif "AF" in results.columns:
                maf_results = results[(results["AF"] >= maf_min) & (results["AF"] < maf_max)]
            else:
                # If MAF information is not available, skip this metric
                logger.warning("MAF information not found in results, skipping power score calculation")
                break
            
            # Calculate power for this MAF range
            if len(maf_results) > 0:
                maf_power = sum(maf_results[p_col] < alpha) / len(maf_results)
                power_scores.append(maf_power)
        
        # Calculate overall power score (mean across MAF ranges)
        if power_scores:
            metrics["power_score"] = np.mean(power_scores)
        else:
            # Assign a neutral power score if we couldn't calculate it
            metrics["power_score"] = 0.0
        
        # 4. Check for recovery of known associations (if provided)
        if self.known_associations is not None:
            n_recovered = 0
            total_known = len(self.known_associations)
            
            for _, known in self.known_associations.iterrows():
                # Find the closest SNP to the known one
                if "CHR" in results.columns and "BP" in results.columns:
                    # Filter to same chromosome
                    chr_results = results[results["CHR"] == known["CHR"]]
                    
                    if len(chr_results) > 0:
                        # Find SNPs within window
                        window = self.args.known_window  # Window size in base pairs
                        nearby = chr_results[
                            (chr_results["BP"] >= known["BP"] - window) & 
                            (chr_results["BP"] <= known["BP"] + window)
                        ]
                        
                        # Check if any pass significance threshold
                        if len(nearby) > 0 and any(nearby[p_col] < alpha):
                            n_recovered += 1
            
            if total_known > 0:
                metrics["recovery_rate"] = n_recovered / total_known
            else:
                metrics["recovery_rate"] = 0.0
        else:
            metrics["recovery_rate"] = None
        
        return metrics
    
    def calculate_objective(self, metrics):
        """
        Calculate the objective value to maximize.
        
        This function combines multiple metrics into a single objective
        value that Optuna will try to maximize.
        
        Args:
            metrics: Dictionary of evaluation metrics
        
        Returns:
            Objective value (higher is better)
        """
        # 1. Penalize lambda_gc that deviates from 1.0
        lambda_penalty = -abs(metrics["lambda_gc"] - 1.0) * 3.0
        
        # 2. Reward number of significant hits (normalized)
        if metrics["n_significant"] > 0:
            sig_reward = min(np.log10(metrics["n_significant"]) / 2.0, 1.0)
        else:
            sig_reward = 0.0
        
        # 3. Reward power score
        power_reward = metrics["power_score"] * 2.0
        
        # 4. Reward recovery of known associations (if available)
        if metrics["recovery_rate"] is not None:
            recovery_reward = metrics["recovery_rate"] * 3.0
        else:
            recovery_reward = 0.0
        
        # Combine components
        objective = lambda_penalty + sig_reward + power_reward + recovery_reward
        
        return objective
    
    def objective(self, trial):
        """
        Optuna objective function to maximize.
        
        Args:
            trial: Optuna trial object
        
        Returns:
            Objective value to maximize
        """
        # Keep track of current trial number for file naming
        self.current_trial_number = trial.number
        
        # Get parameters for this trial
        params = self.define_parameter_space(trial)
        
        # Log parameters
        logger.info(f"Trial {trial.number}: {params}")
        
        # Check if we're using cross-validation
        if self.args.cv_folds > 1:
            # Implement k-fold cross-validation
            metrics_list = []
            
            # Create stratified folds
            if self.args.pheno_file and self.args.pheno_name:
                # Load phenotype data for stratification
                pheno_data = pd.read_csv(self.args.pheno_file, sep="\s+")
                pheno_data.set_index(pheno_data.columns[0], inplace=True)
                
                # Match with FAM file to get the right order
                fam = pd.read_csv(f"{self.plink_prefix}.fam", sep="\s+", header=None)
                fam.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
                
                # Align phenotypes with FAM file
                y = []
                for iid in fam["IID"]:
                    if iid in pheno_data.index:
                        y.append(pheno_data.loc[iid, self.args.pheno_name])
                    else:
                        y.append(np.nan)
                
                y = np.array(y)
                
                # Remove missing phenotypes
                valid_idx = ~np.isnan(y)
                y = y[valid_idx]
                
                # Create stratified folds
                skf = StratifiedKFold(n_splits=self.args.cv_folds, shuffle=True, random_state=42)
                fold_indices = list(skf.split(np.zeros(len(y)), y))
            else:
                # If no phenotype file, just create random folds
                n_samples = sum(1 for _ in open(f"{self.plink_prefix}.fam"))
                indices = np.arange(n_samples)
                fold_size = n_samples // self.args.cv_folds
                fold_indices = [(np.setdiff1d(indices, indices[i*fold_size:(i+1)*fold_size]), 
                                indices[i*fold_size:(i+1)*fold_size]) 
                                for i in range(self.args.cv_folds)]
            
            # Run GWAS on each fold
            for fold, (_, test_idx) in enumerate(fold_indices):
                # Create fold-specific files
                fold_prefix = os.path.join(self.study_dir, f"fold_{fold}")
                
                # Create a list of samples to include for this fold
                sample_file = f"{fold_prefix}.sample_list"
                with open(sample_file, "w") as f:
                    for idx in test_idx:
                        f.write(f"{idx}\n")
                
                # Run GWAS on this fold
                result_file = self.run_gwas(params, fold)
                
                if result_file:
                    # Calculate metrics for this fold
                    fold_metrics = self.calculate_metrics(result_file, params)
                    metrics_list.append(fold_metrics)
            
            # Average metrics across folds
            if metrics_list:
                avg_metrics = {}
                for key in metrics_list[0].keys():
                    avg_metrics[key] = np.mean([m[key] for m in metrics_list if m[key] is not None])
                
                # Calculate objective from averaged metrics
                objective_value = self.calculate_objective(avg_metrics)
            else:
                # If all folds failed, return a very low value
                objective_value = -100.0
        else:
            # No cross-validation, just run once on all data
            result_file = self.run_gwas(params)
            
            if result_file:
                # Calculate metrics
                metrics = self.calculate_metrics(result_file, params)
                
                # Calculate objective
                objective_value = self.calculate_objective(metrics)
                
                # Store metrics in trial user attributes for later analysis
                trial.set_user_attr("lambda_gc", metrics["lambda_gc"])
                trial.set_user_attr("n_significant", metrics["n_significant"])
                trial.set_user_attr("power_score", metrics["power_score"])
                if metrics["recovery_rate"] is not None:
                    trial.set_user_attr("recovery_rate", metrics["recovery_rate"])
            else:
                # If GWAS failed, return a very low value
                objective_value = -100.0
        
        return objective_value
    
    def calculate_power(self, best_params):
        """
        Calculate statistical power across MAF ranges for the best parameters.
        
        Args:
            best_params: Best parameters from optimization
        
        Returns:
            DataFrame with power estimates by MAF
        """
        # Run GWAS with best parameters
        self.current_trial_number = "power_analysis"
        result_file = self.run_gwas(best_params)
        
        if not result_file:
            logger.error("Failed to run GWAS for power analysis")
            return None
        
        # Load results
        if best_params["analysis_method"] == "assoc":
            results = pd.read_table(result_file, delim_whitespace=True)
            p_col = "P"
        elif best_params["analysis_method"] == "logistic":
            results = pd.read_table(result_file, delim_whitespace=True)
            results = results[results["TEST"] == "ADD"]
            p_col = "P"
        elif best_params["analysis_method"] == "mlma":
            results = pd.read_table(result_file, delim_whitespace=True)
            p_col = "p"
        
        # Filter out missing p-values
        results = results.dropna(subset=[p_col])
        
        # Determine significance threshold
        n_tests = len(results)
        if best_params["multiple_testing"] == "bonferroni":
            alpha = 0.05 / n_tests
        elif best_params["multiple_testing"] == "fdr":
            from statsmodels.stats.multitest import multipletests
            _, adj_p, _, _ = multipletests(results[p_col], method='fdr_bh')
            results["ADJ_P"] = adj_p
            alpha = 0.05
            p_col = "ADJ_P"
        else:  # none
            alpha = 5e-8
        
        # Define MAF bins
        maf_bins = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        
        # Calculate power for each MAF bin
        power_df = pd.DataFrame(columns=["MAF_min", "MAF_max", "n_variants", "n_significant", "power"])
        
        for i in range(len(maf_bins)-1):
            maf_min = maf_bins[i]
            maf_max = maf_bins[i+1]
            
            # Filter results by MAF range
            if "MAF" in results.columns:
                maf_results = results[(results["MAF"] >= maf_min) & (results["MAF"] < maf_max)]
            elif "AF" in results.columns:
                maf_results = results[(results["AF"] >= maf_min) & (results["AF"] < maf_max)]
            else:
                logger.warning("MAF information not found in results, skipping power analysis")
                return None
            
            # Calculate power for this MAF range
            n_variants = len(maf_results)
            n_significant = sum(maf_results[p_col] < alpha)
            
            if n_variants > 0:
                power = n_significant / n_variants
            else:
                power = 0.0
            
            # Add to DataFrame
            power_df = pd.concat([
                power_df,
                pd.DataFrame({
                    "MAF_min": [maf_min],
                    "MAF_max": [maf_max],
                    "n_variants": [n_variants],
                    "n_significant": [n_significant],
                    "power": [power]
                })
            ])
        
        return power_df
    
    def generate_optimized_script(self, best_params):
        """
        Generate a script to run GWAS with the optimized parameters.
        
        Args:
            best_params: Best parameters from optimization
        
        Returns:
            Path to the generated script
        """
        script_path = os.path.join(self.study_dir, "run_optimized_gwas.sh")
        
        with open(script_path, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("# GWAS script with optimized parameters\n")
            f.write(f"# Generated by gwas_stat_optimizer.py on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Create command based on method
            if best_params["analysis_method"] == "assoc":
                cmd = self.create_plink_assoc_command(best_params, "${OUT_PREFIX}")
            elif best_params["analysis_method"] == "logistic":
                cmd = self.create_plink_logistic_command(best_params, "${OUT_PREFIX}")
            elif best_params["analysis_method"] == "mlma":
                cmd = self.create_gcta_mlma_command(best_params, "${OUT_PREFIX}")
            
            # Write command as variables
            f.write("# Input variables (modify as needed)\n")
            f.write(f"PLINK_PREFIX=\"{self.plink_prefix}\"\n")
            f.write(f"OUT_PREFIX=\"output/gwas_results\"\n\n")
            
            if self.args.pheno_file:
                f.write(f"PHENO_FILE=\"{self.args.pheno_file}\"\n")
            if self.args.covar_file:
                f.write(f"COVAR_FILE=\"{self.args.covar_file}\"\n")
            f.write("\n")
            
            # Write command
            f.write("# GWAS command with optimized parameters\n")
            cmd_str = " ".join(cmd).replace(self.plink_prefix, "${PLINK_PREFIX}")
            if self.args.pheno_file:
                cmd_str = cmd_str.replace(self.args.pheno_file, "${PHENO_FILE}")
            if self.args.covar_file:
                cmd_str = cmd_str.replace(self.args.covar_file, "${COVAR_FILE}")
            f.write(f"{cmd_str}\n\n")
            
            # Make executable
            os.chmod(script_path, 0o755)
            
            return script_path
            
    def run_optimization(self):
        """Run the optimization process."""
        logger.info(f"Starting optimization with {self.args.n_trials} trials")
        
        # Run optimization
        self.study.optimize(
            self.objective, 
            n_trials=self.args.n_trials,
            n_jobs=self.args.n_jobs
        )
        
        # Get best parameters
        best_params = self.study.best_params
        best_value = self.study.best_value
        
        logger.info(f"Best parameters: {best_params}")
        logger.info(f"Best objective value: {best_value}")
        
        # Run power analysis with best parameters
        logger.info("Running power analysis with best parameters")
        power_df = self.calculate_power(best_params)
        if power_df is not None:
            power_path = os.path.join(self.study_dir, "power_analysis.csv")
            power_df.to_csv(power_path, index=False)
            logger.info(f"Power analysis saved to {power_path}")
            
            # Generate power plot
            vis_dir = os.path.join(self.study_dir, "visualization")
            os.makedirs(vis_dir, exist_ok=True)
            self.generate_power_plot(power_df, vis_dir)
        
        # Generate optimized script
        script_path = self.generate_optimized_script(best_params)
        logger.info(f"Optimized GWAS script saved to {script_path}")
        
        # Generate visualization
        self.generate_visualization()
        
        # Create summary report
        report_path = self.create_summary_report()
        logger.info(f"Summary report saved to {report_path}")
        
        return best_params, best_value
    
    def generate_visualization(self):
        """Generate visualizations of the optimization process."""
        vis_dir = os.path.join(self.study_dir, "visualization")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Plot optimization history
        try:
            fig = plot_optimization_history(self.study)
            fig.write_image(os.path.join(vis_dir, "optimization_history.png"))
            
            # Plot parameter importance
            fig = plot_param_importances(self.study)
            fig.write_image(os.path.join(vis_dir, "param_importance.png"))
            
            # Create custom plots
            self.plot_best_trial_metrics(vis_dir)
            self.plot_lambda_gc_distribution(vis_dir)
            
            logger.info(f"Visualization saved to {vis_dir}")
        except Exception as e:
            logger.error(f"Error generating visualization: {e}")
    
    def plot_best_trial_metrics(self, vis_dir):
        """
        Plot metrics from the best trial.
        
        Args:
            vis_dir: Directory to save visualization
        """
        # Get best parameters
        best_params = self.study.best_params
        
        # Run GWAS with best parameters
        self.current_trial_number = "best"
        result_file = self.run_gwas(best_params)
        
        if not result_file:
            logger.error("Failed to run GWAS for best trial metrics plot")
            return
        
        # Load results
        if best_params["analysis_method"] == "assoc":
            results = pd.read_table(result_file, delim_whitespace=True)
            p_col = "P"
        elif best_params["analysis_method"] == "logistic":
            results = pd.read_table(result_file, delim_whitespace=True)
            results = results[results["TEST"] == "ADD"]
            p_col = "P"
        elif best_params["analysis_method"] == "mlma":
            results = pd.read_table(result_file, delim_whitespace=True)
            p_col = "p"
        
        # Filter out missing p-values
        results = results.dropna(subset=[p_col])
        
        # Create Manhattan plot
        try:
            plt.figure(figsize=(12, 6))
            
            # Get chromosome and position
            if "CHR" in results.columns and "BP" in results.columns:
                # Calculate x-axis position
                results = results.sort_values(["CHR", "BP"])
                results["pos"] = 0
                
                # Calculate cumulative position
                chrom_ends = {}
                pos = 0
                for chrom in sorted(results["CHR"].unique()):
                    chrom_data = results[results["CHR"] == chrom]
                    results.loc[results["CHR"] == chrom, "pos"] = pos + chrom_data["BP"]
                    pos = results.loc[results["CHR"] == chrom, "pos"].max() + 1e6  # 1Mb spacing
                    chrom_ends[chrom] = pos
                
                # Plot points by chromosome
                chroms = sorted(results["CHR"].unique())
                colors = ["#1f77b4", "#ff7f0e"]  # Alternating colors
                
                for i, chrom in enumerate(chroms):
                    chrom_data = results[results["CHR"] == chrom]
                    plt.scatter(
                        chrom_data["pos"], 
                        -np.log10(chrom_data[p_col]),
                        s=3,
                        color=colors[i % len(colors)],
                        alpha=0.8
                    )
                
                # Add chromosome labels
                chrom_mids = {}
                for chrom in chroms:
                    if chrom > 1:
                        prev_end = chrom_ends[chrom - 1]
                    else:
                        prev_end = 0
                    chrom_mids[chrom] = (prev_end + chrom_ends[chrom]) / 2
                
                plt.xticks(list(chrom_mids.values()), list(chrom_mids.keys()))
                plt.xlim(0, max(chrom_ends.values()))
            else:
                # If no chromosome info, just plot by index
                plt.scatter(
                    range(len(results)),
                    -np.log10(results[p_col]),
                    s=3,
                    alpha=0.8
                )
            
            # Add significance threshold
            n_tests = len(results)
            if best_params["multiple_testing"] == "bonferroni":
                alpha = 0.05 / n_tests
            elif best_params["multiple_testing"] == "fdr":
                alpha = 0.05  # Using adjusted p-values
            else:  # none
                alpha = 5e-8
            
            plt.axhline(-np.log10(alpha), color="red", linestyle="--")
            
            plt.xlabel("Chromosome")
            plt.ylabel("-log10(p-value)")
            plt.title("Manhattan Plot of Best Trial Results")
            
            plt.tight_layout()
            plt.savefig(os.path.join(vis_dir, "manhattan_plot.png"), dpi=300)
            plt.close()
            
            # Create QQ plot
            plt.figure(figsize=(8, 8))
            
            observed_p = results[p_col].values
            observed_p.sort()
            expected_p = np.linspace(0, 1, len(observed_p) + 2)[1:-1]
            
            plt.scatter(
                -np.log10(expected_p),
                -np.log10(observed_p),
                s=3,
                alpha=0.8
            )
            
            # Add diagonal line
            max_val = max(-np.log10(expected_p).max(), -np.log10(observed_p[observed_p > 0]).max())
            plt.plot([0, max_val], [0, max_val], "r--")
            
            plt.xlabel("Expected -log10(p)")
            plt.ylabel("Observed -log10(p)")
            plt.title(f"QQ Plot (λGC = {np.median(stats.chi2.ppf(1 - observed_p, 1)) / stats.chi2.ppf(0.5, 1):.3f})")
            
            plt.tight_layout()
            plt.savefig(os.path.join(vis_dir, "qq_plot.png"), dpi=300)
            plt.close()
            
        except Exception as e:
            logger.error(f"Error creating Manhattan/QQ plots: {e}")
    
    def plot_lambda_gc_distribution(self, vis_dir):
        """
        Plot distribution of lambda GC values across trials.
        
        Args:
            vis_dir: Directory to save visualization
        """
        try:
            # Extract lambda GC values from trials
            lambda_values = []
            for trial in self.study.trials:
                if trial.state == optuna.trial.TrialState.COMPLETE:
                    # Find lambda in trial user attributes
                    for key, value in trial.user_attrs.items():
                        if key == "lambda_gc":
                            lambda_values.append(value)
                            break
            
            if not lambda_values:
                logger.warning("No lambda GC values found in trials")
                return
            
            # Create histogram
            plt.figure(figsize=(10, 6))
            sns.histplot(lambda_values, bins=20, kde=True)
            
            plt.axvline(1.0, color="red", linestyle="--")
            plt.xlabel("Genomic Inflation Factor (λGC)")
            plt.ylabel("Number of Trials")
            plt.title("Distribution of Genomic Inflation Factors Across Trials")
            
            plt.tight_layout()
            plt.savefig(os.path.join(vis_dir, "lambda_distribution.png"), dpi=300)
            plt.close()
            
        except Exception as e:
            logger.error(f"Error creating lambda GC distribution plot: {e}")
    
    def generate_power_plot(self, power_df, vis_dir):
        """
        Generate power analysis plot.
        
        Args:
            power_df: DataFrame with power analysis results
            vis_dir: Directory to save visualization
        """
        try:
            if power_df is None or len(power_df) == 0:
                logger.warning("No power data available for plotting")
                return
            
            plt.figure(figsize=(10, 6))
            
            # Plot power by MAF
            plt.bar(
                range(len(power_df)),
                power_df["power"],
                alpha=0.7
            )
            
            # Add bar labels
            for i, p in enumerate(power_df["power"]):
                plt.text(
                    i,
                    p + 0.02,
                    f"{p:.3f}",
                    ha="center"
                )
            
            # Customize x-axis labels
            x_labels = [f"{row['MAF_min']:.3f}-{row['MAF_max']:.3f}" for _, row in power_df.iterrows()]
            plt.xticks(range(len(power_df)), x_labels, rotation=45)
            
            plt.xlabel("Minor Allele Frequency Range")
            plt.ylabel("Statistical Power")
            plt.title("Statistical Power by Minor Allele Frequency")
            plt.ylim(0, 1.1)
            
            plt.tight_layout()
            plt.savefig(os.path.join(vis_dir, "power_analysis.png"), dpi=300)
            plt.close()
            
            # Plot number of variants by MAF
            plt.figure(figsize=(10, 6))
            
            plt.bar(
                range(len(power_df)),
                power_df["n_variants"],
                alpha=0.7
            )
            
            # Add bar labels
            for i, n in enumerate(power_df["n_variants"]):
                plt.text(
                    i,
                    n + max(power_df["n_variants"]) * 0.02,
                    str(n),
                    ha="center"
                )
            
            plt.xticks(range(len(power_df)), x_labels, rotation=45)
            
            plt.xlabel("Minor Allele Frequency Range")
            plt.ylabel("Number of Variants")
            plt.title("Variant Count by Minor Allele Frequency")
            
            plt.tight_layout()
            plt.savefig(os.path.join(vis_dir, "variant_count.png"), dpi=300)
            plt.close()
            
        except Exception as e:
            logger.error(f"Error creating power analysis plot: {e}")
            
    def create_summary_report(self):
        """
        Create a summary report of the optimization results.
        
        Returns:
            Path to the summary report
        """
        report_path = os.path.join(self.study_dir, "optimization_summary.txt")
        
        with open(report_path, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("GWAS PARAMETER OPTIMIZATION SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            # Date and time
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Study Name: {self.args.study_name}\n")
            f.write(f"PLINK Files: {self.plink_prefix}\n\n")
            
            # Best parameters
            f.write("OPTIMAL PARAMETERS\n")
            f.write("-" * 50 + "\n")
            
            best_params = self.study.best_params
            best_value = self.study.best_value
            
            for param, value in best_params.items():
                f.write(f"{param}: {value}\n")
            
            f.write(f"\nObjective Value: {best_value:.4f}\n\n")
            
            # Best trial metrics
            best_trial = self.study.best_trial
            if best_trial.user_attrs:
                f.write("PERFORMANCE METRICS\n")
                f.write("-" * 50 + "\n")
                
                for key, value in best_trial.user_attrs.items():
                    f.write(f"{key}: {value:.4f}\n")
                
                f.write("\n")
            
            # Command to run GWAS with best parameters
            f.write("OPTIMIZED GWAS COMMAND\n")
            f.write("-" * 50 + "\n")
            
            if best_params["analysis_method"] == "assoc":
                cmd = self.create_plink_assoc_command(best_params, "output/gwas_results")
            elif best_params["analysis_method"] == "logistic":
                cmd = self.create_plink_logistic_command(best_params, "output/gwas_results")
            elif best_params["analysis_method"] == "mlma":
                cmd = self.create_gcta_mlma_command(best_params, "output/gwas_results")
            
            cmd_str = " ".join(cmd)
            f.write(f"{cmd_str}\n\n")
            
            # Trial statistics
            n_completed = len([t for t in self.study.trials if t.state == optuna.trial.TrialState.COMPLETE])
            n_pruned = len([t for t in self.study.trials if t.state == optuna.trial.TrialState.PRUNED])
            n_failed = len([t for t in self.study.trials if t.state == optuna.trial.TrialState.FAIL])
            
            f.write("OPTIMIZATION STATISTICS\n")
            f.write("-" * 50 + "\n")
            f.write(f"Total Trials: {len(self.study.trials)}\n")
            f.write(f"Completed Trials: {n_completed}\n")
            f.write(f"Pruned Trials: {n_pruned}\n")
            f.write(f"Failed Trials: {n_failed}\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("End of Report\n")
        
def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="GWAS Statistical Analysis Parameter Optimization with Optuna",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "--plink-prefix",
        required=True,
        help="Path prefix to PLINK files (.bed, .bim, .fam)"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to save optimization results"
    )
    
    # Optional arguments for GWAS parameters
    parser.add_argument(
        "--pheno-file",
        help="Path to phenotype file"
    )
    parser.add_argument(
        "--pheno-name",
        help="Name of phenotype column in phenotype file"
    )
    parser.add_argument(
        "--covar-file",
        help="Path to covariate file"
    )
    parser.add_argument(
        "--max-pcs",
        type=int,
        default=10,
        help="Maximum number of principal components to consider"
    )
    
    # Optuna parameters
    parser.add_argument(
        "--study-name",
        default="gwas_opt",
        help="Name of the Optuna study"
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=100,
        help="Number of optimization trials to run"
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel jobs for optimization"
    )
    parser.add_argument(
        "--resume-study",
        action="store_true",
        help="Resume an existing study with the same name"
    )
    
    # Cross-validation parameters
    parser.add_argument(
        "--cv-folds",
        type=int,
        default=1,
        help="Number of cross-validation folds (1 means no CV)"
    )
    
    # Method selection
    parser.add_argument(
        "--disable-assoc",
        action="store_true",
        help="Disable PLINK's basic association test"
    )
    parser.add_argument(
        "--disable-logistic",
        action="store_true",
        help="Disable PLINK's logistic regression"
    )
    parser.add_argument(
        "--disable-mlma",
        action="store_true",
        help="Disable GCTA's mixed linear model association"
    )
    
    # Path to executables
    parser.add_argument(
        "--plink-path",
        help="Path to PLINK executable"
    )
    parser.add_argument(
        "--gcta-path",
        help="Path to GCTA executable"
    )
    
    # Known associations for validation
    parser.add_argument(
        "--known-associations",
        help="Path to CSV file with known SNP associations (requires CHR, BP columns)"
    )
    parser.add_argument(
        "--known-window",
        type=int,
        default=100000,
        help="Window size in base pairs around known associations"
    )
    
    return parser.parse_args()

def main():
    """Main function to run the GWAS parameter optimization."""
    # Parse command line arguments
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Initialize optimizer
        optimizer = GWASOptimizer(args)
        
        # Run optimization
        best_params, best_value = optimizer.run_optimization()
        
        logger.info(f"Optimization complete. Best objective: {best_value:.4f}")
        logger.info(f"Best parameters: {best_params}")
        
        # Exit with success
        return 0
        
    except Exception as e:
        logger.error(f"Error during optimization: {e}")
        logger.exception(e)
        
        # Exit with error
        return 1

if __name__ == "__main__":
    sys.exit(main())