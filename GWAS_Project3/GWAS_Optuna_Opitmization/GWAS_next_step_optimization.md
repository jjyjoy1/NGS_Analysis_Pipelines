
# Additional GWAS Optimization Metrics

## 1. Heritability Estimation
**Description**: Estimate the proportion of phenotypic variance explained by genetic variants.

**Implementation**:
```python
def calculate_heritability(self, result_file, params):
    """Calculate heritability from GWAS results."""
    # Load results
    results = pd.read_table(result_file, delim_whitespace=True)
    
    # Identify significant SNPs
    if params["analysis_method"] == "assoc":
        p_col = "P"
    elif params["analysis_method"] == "logistic":
        results = results[results["TEST"] == "ADD"]
        p_col = "P"
    elif params["analysis_method"] == "mlma":
        p_col = "p"
    
    # Determine significance threshold
    n_tests = len(results)
    if params["multiple_testing"] == "bonferroni":
        alpha = 0.05 / n_tests
    elif params["multiple_testing"] == "fdr":
        alpha = 0.05  # Using adjusted p-values
    else:
        alpha = 5e-8  # Standard GWAS threshold
    
    # Filter significant SNPs
    sig_snps = results[results[p_col] < alpha]
    
    # Calculate heritability using LDSC or GCTA-GREML
    if self.args.heritability_method == "ldsc":
        # Implement LD Score Regression
        h2 = self._calculate_ldsc(sig_snps)
    elif self.args.heritability_method == "greml":
        # Implement GCTA-GREML
        h2 = self._calculate_greml(sig_snps)
    else:
        # Simple variance explained calculation
        h2 = self._calculate_simple_h2(sig_snps)
    
    return h2
```

**Why it's valuable**: Heritability estimates help assess how well the model captures genetic architecture and can identify parameter combinations that maximize explained variance.

---

## 2. Replication Consistency
**Description**: Split the dataset into discovery and replication sets and measure consistency of results.

**Implementation**:
```python
def calculate_replication_consistency(self, result_file, params):
    """Calculate consistency of results in discovery and replication sets."""
    # Split dataset
    discovery_results = self._run_on_subset(params, subset="discovery")
    replication_results = self._run_on_subset(params, subset="replication")
    
    # Load and process results
    disc_df = pd.read_table(discovery_results, delim_whitespace=True)
    repl_df = pd.read_table(replication_results, delim_whitespace=True)
    
    # Find significant SNPs in discovery
    p_col = "P" if params["analysis_method"] != "mlma" else "p"
    alpha = 5e-8  # Or use appropriate threshold
    disc_sig = disc_df[disc_df[p_col] < alpha]
    
    # Check if significant SNPs replicate
    replication_rate = 0
    direction_consistency = 0
    
    if len(disc_sig) > 0:
        # Match SNPs between datasets
        common_snps = pd.merge(disc_sig, repl_df, on="SNP")
        
        # Calculate replication at nominal significance
        replicated = sum(common_snps[f"{p_col}_y"] < 0.05)
        replication_rate = replicated / len(disc_sig) if len(disc_sig) > 0 else 0
        
        # Calculate effect direction consistency
        if "BETA" in common_snps.columns:
            same_direction = sum(np.sign(common_snps["BETA_x"]) == np.sign(common_snps["BETA_y"]))
        elif "OR" in common_snps.columns:
            same_direction = sum((common_snps["OR_x"] > 1) == (common_snps["OR_y"] > 1))
        else:
            same_direction = 0
            
        direction_consistency = same_direction / len(common_snps) if len(common_snps) > 0 else 0
    
    return {"replication_rate": replication_rate, "direction_consistency": direction_consistency}
```

**Why it's valuable**: Replication consistency indicates whether the parameter settings are robust and likely to generalize to new datasets, which is crucial for real-world applicability.

---

## 3. Enrichment of Functional Annotations
**Description**: Measure enrichment of significant SNPs in functional genomic regions.

**Implementation**:
```python
def calculate_functional_enrichment(self, result_file, params):
    """Calculate enrichment of significant SNPs in functional annotations."""
    # Load results
    results = pd.read_table(result_file, delim_whitespace=True)
    
    # Identify significant SNPs
    p_col = "P" if params["analysis_method"] != "mlma" else "p"
    alpha = 5e-8  # Or appropriate threshold
    sig_snps = results[results[p_col] < alpha]
    
    # Load functional annotations (from a BED file, for example)
    annotations = pd.read_csv(self.args.annotation_file, sep="\t", 
                             names=["chr", "start", "end", "feature"])
    
    # Map SNPs to genomic positions
    if all(x in results.columns for x in ["CHR", "BP"]):
        # Count SNPs in functional regions
        sig_in_func = 0
        all_in_func = 0
        
        for _, snp in sig_snps.iterrows():
            chr_annot = annotations[annotations["chr"] == f"chr{snp['CHR']}"]
            if any((chr_annot["start"] <= snp["BP"]) & (chr_annot["end"] >= snp["BP"])):
                sig_in_func += 1
        
        for _, snp in results.iterrows():
            chr_annot = annotations[annotations["chr"] == f"chr{snp['CHR']}"]
            if any((chr_annot["start"] <= snp["BP"]) & (chr_annot["end"] >= snp["BP"])):
                all_in_func += 1
        
        # Calculate enrichment
        sig_ratio = sig_in_func / len(sig_snps) if len(sig_snps) > 0 else 0
        all_ratio = all_in_func / len(results) if len(results) > 0 else 0
        
        enrichment = sig_ratio / all_ratio if all_ratio > 0 else 0
        return enrichment
    else:
        return 0  # No position information available
```

**Why it's valuable**: Functional enrichment helps identify parameter settings that prioritize biologically relevant variants, potentially increasing the interpretability of results.

---

## 4. LD Structure Preservation
**Description**: Evaluate how well the chosen parameters preserve the linkage disequilibrium (LD) structure.

**Implementation**:
```python
def calculate_ld_preservation(self, result_file, params):
    """Calculate how well the LD structure is preserved in significant results."""
    # Load GWAS results
    results = pd.read_table(result_file, delim_whitespace=True)
    
    # Identify significant SNPs
    p_col = "P" if params["analysis_method"] != "mlma" else "p"
    alpha = 5e-8  # Or appropriate threshold
    sig_snps = results[results[p_col] < alpha]
    
    # Check if we have significant hits
    if len(sig_snps) < 10:  # Need enough SNPs for meaningful LD analysis
        return 0
    
    # Calculate LD among significant SNPs using PLINK
    sig_snps_file = os.path.join(self.study_dir, f"trial_{self.current_trial_number}_sig_snps.txt")
    with open(sig_snps_file, "w") as f:
        for snp in sig_snps["SNP"]:
            f.write(f"{snp}\n")
    
    # Run PLINK to calculate LD
    ld_out = os.path.join(self.study_dir, f"trial_{self.current_trial_number}_ld")
    cmd = [
        self.plink_path,
        "--bfile", self.plink_prefix,
        "--extract", sig_snps_file,
        "--r2",
        "--ld-snp-list", sig_snps_file,
        "--out", ld_out
    ]
    
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Load LD results
    ld_file = f"{ld_out}.ld"
    if os.path.exists(ld_file):
        ld_data = pd.read_table(ld_file, delim_whitespace=True)
        
        # Calculate proportion of SNP pairs in expected LD patterns
        # For example, SNPs in close proximity should have higher LD
        valid_patterns = 0
        total_pairs = 0
        
        for _, row in ld_data.iterrows():
            bp_distance = abs(row["BP_A"] - row["BP_B"])
            r2 = row["R2"]
            
            # Define expected patterns (simplified example)
            if bp_distance < 10000 and r2 > 0.5:  # Close SNPs should have high LD
                valid_patterns += 1
            elif bp_distance > 100000 and r2 < 0.2:  # Distant SNPs should have low LD
                valid_patterns += 1
                
            total_pairs += 1
        
        ld_score = valid_patterns / total_pairs if total_pairs > 0 else 0
        return ld_score
    else:
        return 0  # LD calculation failed
```

**Why it's valuable**: This metric helps ensure that the chosen parameters don't disrupt the natural LD structure in the genome, which could lead to artificial signals or missed associations.

---

```


```
## 5. Polygenic Risk Score Performance
**Description**: Test how well the identified SNPs predict the phenotype through polygenic risk scoring.

**Implementation**:
```python
def calculate_prs_performance(self, result_file, params):
    """Calculate predictive performance of PRS using the GWAS results."""
    # Split dataset into training and testing
    training_results = self._run_on_subset(params, subset="training")
    
    # Extract significant SNPs and their effect sizes
    train_df = pd.read_table(training_results, delim_whitespace=True)
    
    # Determine p-value column and effect size column
    p_col = "P" if params["analysis_method"] != "mlma" else "p"
    
    if "BETA" in train_df.columns:
        effect_col = "BETA"
    elif "OR" in train_df.columns:
        effect_col = "OR"
        # Convert OR to BETA for PRS calculation
        train_df["BETA"] = np.log(train_df["OR"])
        effect_col = "BETA"
    else:
        return 0  # No effect size available
    
    # Select SNPs for PRS (using different p-value thresholds)
    prs_thresholds = [5e-8, 1e-5, 1e-3, 0.01, 0.05, 0.1, 0.5, 1.0]
    best_r2 = 0
    
    for threshold in prs_thresholds:
        # Extract SNPs and effect sizes
        prs_snps = train_df[train_df[p_col] <= threshold]
        
        if len(prs_snps) < 10:  # Need enough SNPs for meaningful PRS
            continue
        
        # Create SNP file for PLINK
        snp_file = os.path.join(self.study_dir, f"trial_{self.current_trial_number}_prs_snps.txt")
        prs_snps[["SNP", "BETA"]].to_csv(snp_file, sep="\t", index=False, header=False)
        
        # Calculate PRS using PLINK
        prs_out = os.path.join(self.study_dir, f"trial_{self.current_trial_number}_prs")
        cmd = [
            self.plink_path,
            "--bfile", self.plink_prefix,
            "--score", snp_file, "1", "2",  # Use SNP name and effect size
            "--out", prs_out
        ]
        
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Load PRS and phenotype data
        prs_file = f"{prs_out}.profile"
        if os.path.exists(prs_file):
            prs_data = pd.read_table(prs_file, delim_whitespace=True)
            
            # Load phenotype data
            pheno_data = pd.read_csv(self.args.pheno_file, sep="\s+")
            
            # Merge PRS with phenotype
            merged_data = pd.merge(prs_data, pheno_data, left_on="IID", right_on=pheno_data.columns[0])
            
            # Calculate correlation between PRS and phenotype
            corr = np.corrcoef(merged_data["SCORE"], merged_data[self.args.pheno_name])[0, 1]
            r2 = corr ** 2
            
            if r2 > best_r2:
                best_r2 = r2
    
    return best_r2
```

**Why it's valuable**: PRS performance directly measures the practical utility of the GWAS results for phenotype prediction, which is often a key goal of genetic studies.

---

## 6. Population Stratification Robustness
**Description**: Evaluate how well the parameters handle population stratification.

**Implementation**:
```python
def calculate_population_robustness(self, result_file, params):
    """Calculate robustness to population stratification."""
    # Run GWAS with and without population stratification control
    # This could mean running with different PC counts
    
    # First run: no PCs
    params_no_pcs = params.copy()
    params_no_pcs["n_pcs"] = 0
    result_no_pcs = self._run_with_params(params_no_pcs)
    
    # Second run: with optimal PCs
    params_with_pcs = params.copy()
    params_with_pcs["n_pcs"] = self.args.max_pcs
    result_with_pcs = self._run_with_params(params_with_pcs)
    
    # Load results
    no_pcs_df = pd.read_table(result_no_pcs, delim_whitespace=True)
    with_pcs_df = pd.read_table(result_with_pcs, delim_whitespace=True)
    
    # Determine p-value column
    p_col = "P" if params["analysis_method"] != "mlma" else "p"
    
    # Calculate lambda GC for both
    lambda_no_pcs = self._calculate_lambda_gc(no_pcs_df[p_col].values)
    lambda_with_pcs = self._calculate_lambda_gc(with_pcs_df[p_col].values)
    
    # Measure robustness as the difference in lambda
    robustness = 1.0 - abs(lambda_no_pcs - lambda_with_pcs)
    
    # Also consider overlap in top SNPs
    top_n = 100  # Number of top SNPs to compare
    
    # Get top SNPs from both analyses
    top_no_pcs = no_pcs_df.nsmallest(top_n, p_col)["SNP"].tolist()
    top_with_pcs = with_pcs_df.nsmallest(top_n, p_col)["SNP"].tolist()
    
    # Calculate overlap
    overlap = len(set(top_no_pcs).intersection(set(top_with_pcs))) / top_n
    
    # Combine metrics (higher is better)
    combined_score = (robustness + overlap) / 2
    
    return combined_score
```

**Why it's valuable**: This metric helps identify parameter settings that produce results robust to population structure, reducing the risk of spurious associations.

---

## 7. Genetic Correlation with Related Traits
**Description**: Measure genetic correlation between the GWAS results and those from related phenotypes.

**Implementation**:
```python
def calculate_genetic_correlation(self, result_file, params):
    """Calculate genetic correlation with related traits."""
    # This would typically use tools like LDSC
    # For simplicity, we'll use a direct comparison with known GWAS results
    
    # Load current GWAS results
    results = pd.read_table(result_file, delim_whitespace=True)
    
    # Load reference GWAS results for related traits
    ref_gwas_files = self.args.reference_gwas.split(",")
    
    # Calculate average correlation
    correlations = []
    
    for ref_file in ref_gwas_files:
        if os.path.exists(ref_file):
            ref_results = pd.read_table(ref_file, delim_whitespace=True)
            
            # Merge results by SNP
            merged = pd.merge(results, ref_results, on="SNP", suffixes=("_curr", "_ref"))
            
            # Get effect columns
            if "BETA_curr" in merged.columns and "BETA_ref" in merged.columns:
                corr = np.corrcoef(merged["BETA_curr"], merged["BETA_ref"])[0, 1]
                correlations.append(abs(corr))  # Use absolute correlation
            elif "OR_curr" in merged.columns and "OR_ref" in merged.columns:
                # Convert OR to log scale
                log_or_curr = np.log(merged["OR_curr"])
                log_or_ref = np.log(merged["OR_ref"])
                corr = np.corrcoef(log_or_curr, log_or_ref)[0, 1]
                correlations.append(abs(corr))
    
    # Return average correlation
    if correlations:
        return sum(correlations) / len(correlations)
    else:
        return 0
```

**Why it's valuable**: Genetic correlation helps validate findings against established results and can reveal shared genetic architecture between traits.

---

## Summary of Optimization Metrics
| Metric | Weight | Purpose |
|--------|--------|---------|
| Heritability Estimation | 20% | Assess genetic architecture capture |
| Replication Consistency | 20% | Validate robustness across subsets |
| Functional Enrichment | 15% | Prioritize biologically relevant variants |
| LD Structure Preservation | 15% | Maintain natural genomic patterns |
| PRS Performance | 15% | Measure predictive utility |
| Population Robustness | 10% | Control for stratification effects |
| Genetic Correlation | 5% | Validate against known associations |

These metrics collectively ensure the optimization process balances statistical power, biological relevance, and practical utility of the GWAS results.
```
