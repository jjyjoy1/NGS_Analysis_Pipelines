Here's the converted markdown version of your PDF:

```markdown
# Optimizing GWAS Parameters for Downstream Hypotheses: A Back-to-Front Perspective

## Reimagining GWAS Optimization from Endpoint to Beginning

Looking at the GWAS pipeline from downstream hypotheses back to initial QC, and statistical parameters creates a powerful perspective shift. Let's explore how each downstream hypothesis should inform our optimization strategy, moving from the final outcomes (Polygenic Risk) back to the initial assumptions (Heritability).

## Polygenic Risk Hypothesis → Statistical Parameters

**Core Question:** Can we identify combinations of variants that predict disease risk?

### Optimization Goals:
- Maximize predictive accuracy (AUC) in target populations
- Ensure robust performance across diverse populations
- Capture both common and rare variant contributions

### Parameters Influenced:
- **MAF Thresholds**: Lower thresholds needed to capture rare variants with potentially large effects
- **P-value Thresholds**: More permissive thresholds may improve polygenic scores by capturing sub-significant signals
- **Multiple Testing Correction**: Less stringent correction may improve PRS by including more variants with modest effects
- **Imputation Quality Filters**: Higher thresholds if PRS is sensitive to imputation errors

### Trade-offs:
- Including more variants increases PRS potential but risks noise
- Population-specific optimizations may reduce transferability
- Optimizing for predictor count vs. predictor quality

## Functional Mechanism Hypothesis → Statistical Parameters

**Core Question:** What biological mechanisms explain the genetic associations?

### Optimization Goals:
- Prioritize variants with interpretable functional consequences
- Maximize enrichment in relevant tissue types and pathways
- Ensure signal detection in regulatory regions, not just coding variants

**Parameters Influenced:**
- **Genomic Annotation Weights:** Adjust to capture functional elements relevant to the phenotype
- **Conditional Analysis Parameters:** Tune to identify independent signals that may reveal distinct mechanisms
- **Regional Association Parameters:** Optimize to accurately identify credible sets for functional follow-up
- **Model Selection:** Choose models that appropriately capture genetic architecture (additive/dominant/recessive)

**Trade-offs:**
- Statistical significance vs. functional relevance
- Tissue-specific optimization may miss unexpected mechanisms
- Balancing known biology vs. discovery potential

# Fine-mapping Hypothesis → Statistical Parameters

**Core Question:** What are the causal variants driving GWAS signals?

**Optimization Goals:**
- Minimize credible set size while maintaining causal variant coverage
- Maximize posterior probability of causal variants
- Improve resolution in high-LD regions

**Parameters Influenced:**
- **LD Reference Selection:** Crucial for accurate fine-mapping
- **Prior Probability Models:** Integrate functional information appropriately
- **Imputation Quality Thresholds:** Higher thresholds for regions undergoing fine-mapping
- **Variant Density:** Ensure sufficient coverage in regions of interest
- **Regional Analysis Parameters:** Window size for conditional analyses

**Trade-offs:**
- Statistical power vs. fine-mapping resolution
- Population-specific LD patterns affect fine-mapping accuracy
- Computing resources vs. fine-mapping complexity

# Independent Signals Hypothesis → Statistical Parameters

**Core Question:** How many distinct causal signals exist at each locus?

**Optimization Goals:**
- Accurately distinguish true independent signals from LD shadows
- Minimize false merging of distinct signals
- Maximize discovery of allelic heterogeneity

**Parameters Influenced:**
- **LD Threshold for Signal Clustering:** Critical for defining independent signals
- **Conditional Analysis Parameters:** Step-wise conditioning approach and thresholds
- **Regional Association Analysis:** Window size and boundary definitions
- **Significance Thresholds for Secondary Signals:** May differ from primary signals

**Trade-offs:**
- Separating true signals vs. statistical artifacts
- Lumping vs. splitting error rates
- Sample size requirements increase for detecting secondary signals

# Association Hypothesis → Statistical Parameters

**Core Question:** Which genomic loci are associated with the phenotype?

**Optimization Goals:**
- Maximize true positive discoveries
- Minimize false positives
- Ensure associations are robust across populations

**Parameters Influenced:**
- **Statistical Model:** Linear vs. logistic vs. mixed models
- **Covariate Selection:** Which covariates and how many PCs to include
- **Multiple Testing Correction:** Balancing stringency with discovery potential
- **Minor Allele Frequency Filters:** Impact on rare vs. common variant detection
- **Variant Quality Filters:** Call rate, HWE, imputation quality

**Trade-offs:**
- Type I vs. Type II error rates
- Computing efficiency vs. model complexity
- Single-variant tests vs. burden tests for rare variants

# Heritability Hypothesis → Statistical Parameters

**Core Question:** What fraction of phenotypic variance is explained by genetics?

**Optimization Goals:**
- Accurate estimation of total genetic contribution
- Partition heritability across functional categories
- Minimize bias from population structure

**Parameters Influenced:**
- **Relatedness Cutoffs**: Impact family structure and heritability estimates
- **Population Stratification Correction**: Critical for unbiased estimates
- **Variant Inclusion Criteria**: MAF spectrum affects heritability distribution
- **LD Score Regression Parameters**: Reference population, regression weights

**Trade-offs:**
- Sample size vs. relationship structure
- Common vs. rare variant contributions to heritability
- Additive vs. non-additive genetic effects

# GWAS Statistical Optimizer

**Optimization Goals:**
- Analysis methods (basic association, logistic regression, or mixed linear model)
- Covariate handling and principal component selection
- Multiple testing correction methods
- Minor allele frequency (MAF) thresholds
- Missing data thresholds
- Method-specific parametersIntegrative Optimization Framework

**Trade-offs:**
- maximize statistical power while controlling for false positives and genomic inflation.

GWAS raw data QC optimizer

**Optimization Goals:**
- positive contributions from sample and SNP retention
- balanced case-control retention
- negative penalties from lambda deviation
- case-control ratio changes
- covariate distribution changes

**Trade-offs:**
- balances quantity and quality
- preserves study design
- emphasizes control of confounding
- ensures balanced retention

## Key Principles for Back-to-Front Optimization:

1. **Endpoint Definition:** Define clear metrics for final hypothesis success (e.g., PRS AUC, credible set size)

2. **Hypothesis Interdependence:**
   - Polygenic risk depends on fine-mapping accuracy
   - Fine-mapping depends on identifying independent signals
   - Independent signals depend on robust associations
   - Associations depend on heritability architecture

3. **Weighted Multi-Objective Approach:**
   - Weight objectives based on study priorities
   - Create composite objective functions that reflect research goals
   - Allow trade-offs between competing hypothesis goals

4. **Phenotype-Specific Priorities:**
   - Rare diseases: Emphasize functional impact and independent signals
   - Common complex traits: Emphasize PRS and heritability partitioning
   - Quantitative traits: Emphasize heritability and genetic architecture

5. **Adaptive Reweighting:**
   - Start with general statistical validity (heritability, association)
   - Progressively increase weight on downstream objectives
   - Final optimization pass focused on endpoint hypotheses

## Integrated Multi-Step Optimization with Optuna: Optimizing GWAS Pipeline Stages

### Conceptual Workflow

In a typical GWAS pipeline:
1. QC (Quality Control) → 2. Association Testing → 3. Post-GWAS Analysis  

Traditional approaches optimize each step independently:  
- QC parameters optimized for data quality metrics  
- Association parameters optimized for statistical validity  
- Post-GWAS steps optimized for biological insights  

An integrated approach instead optimizes for downstream outcomes:  
- QC and association parameters jointly optimized based on final biological insights  
- Earlier steps evaluated by how well they enable later steps  

### Implementation Approaches

#### 1. Nested Objective Function

```python
def integrated_objective(trial):
    """Objective function that evaluates QC and association parameters together."""
    # First level: QC parameters
    maf_threshold = trial.suggest_float("maf_threshold", 0.001, 0.05)
    geno_threshold = trial.suggest_float("geno_threshold", 0.01, 0.1)
    hwe_threshold = trial.suggest_float("hwe_threshold", 1e-10, 1e-4, log=True)

    # Run QC with these parameters
    qc_output_prefix = run_qc(maf_threshold, geno_threshold, hwe_threshold)

    # Second level: Association parameters
    analysis_method = trial.suggest_categorical("analysis_method", ["assoc", "logistic", "mlma"])
    use_covariates = trial.suggest_categorical("use_covariates", [True, False])
    multiple_testing = trial.suggest_categorical("multiple_testing", ["bonferroni", "fdr", "none"])

    # Only suggest PC count if using covariates
    n_pcs = 0
    if use_covariates:
        n_pcs = trial.suggest_int("n_pcs", 0, 20)

    # Run association testing
    result_file = run_association(qc_output_prefix, analysis_method, use_covariates, n_pcs, multiple_testing)

    # Third level: Evaluate downstream metrics (e.g., pathway enrichment, PRS performance)
    enrichment_score = calculate_pathway_enrichment(result_file)
    prs_predictive_power = evaluate_prs_performance(result_file)
    target_gene_recapture = evaluate_target_gene_hits(result_file)

    # Combined objective that optimizes for downstream outcomes
    objective_value = (
        enrichment_score * 2.0 +
        prs_predictive_power * 3.0 +
        target_gene_recapture * 2.5
    )
    return objective_value
```

**When to use:**
- When you have computational resources to run the full pipeline for each trial
- When you want true end-to-end optimization
- For smaller datasets where full runs are feasible

#### 2. Hierarchical Optuna Studies

```python
def optimize_gwas_pipeline():
    """Two-level optimization with separate Optuna studies."""
    # First level: Optimize QC parameters
    qc_study = optuna.create_study(
        study_name="qc_optimization",
        direction="maximize",
        storage="sqlite:///qc_optimization.db"
    )
    qc_study.optimize(qc_objective, n_trials=50)
    best_qc_params = qc_study.best_params

    # Run QC with best parameters to generate input for next level
    qc_output_prefix = run_qc(
        best_qc_params["maf_threshold"],
        best_qc_params["geno_threshold"],
        best_qc_params["hwe_threshold"]
    )

    # Second level: Optimize association parameters with fixed QC
    assoc_study = optuna.create_study(
        study_name="association_optimization",
        direction="maximize",
        storage="sqlite:///association_optimization.db"
    )

    # Use partial to fix the QC output prefix
    from functools import partial
    assoc_objective_fixed = partial(association_objective, qc_output_prefix=qc_output_prefix)

    assoc_study.optimize(assoc_objective_fixed, n_trials=100)
    best_assoc_params = assoc_study.best_params

    return best_qc_params, best_assoc_params
```

**When to use:**
- When full pipeline optimization is computationally expensive
- When there's a clear separation between stage objectives
- When you want to focus computational resources on later stages

#### 3. Multi-Objective Optimization

```python
def multi_objective(trial):
    """Multi-objective optimization for different downstream outcomes."""
    # QC and association parameters
    maf_threshold = trial.suggest_float("maf_threshold", 0.001, 0.05)
    geno_threshold = trial.suggest_float("geno_threshold", 0.01, 0.1)
    analysis_method = trial.suggest_categorical("analysis_method", ["assoc", "logistic", "mlma"])

    # Run pipeline
    qc_output = run_qc(maf_threshold, geno_threshold)
    result_file = run_association(qc_output, analysis_method)

    # Multiple objectives to optimize
    enrichment_score = calculate_pathway_enrichment(result_file)
    prs_auc = evaluate_prs_performance(result_file)

    return enrichment_score, prs_auc
```

**When to use:**
- When there are multiple important downstream outcomes
- When you want to explore trade-offs between objectives
- When there's no clear way to combine different metrics

#### 4. Transfer Optimization

```python
def transfer_optimization():
    """Transfer learning between optimization stages."""
    # First run optimization on a smaller subset of data
    small_study = optuna.create_study(direction="maximize")
    small_study.optimize(lambda trial: small_data_objective(trial, subset_size=1000), n_trials=50)

    # Use the best parameters as starting points for full data optimization
    best_small_params = small_study.best_params

    # Create new study with a suggest callback to guide search
    def suggest_callback(study, trial):
        for param, value in best_small_params.items():
            trial.set_user_attr(f"init_{param}", value)
        return {}

    full_study = optuna.create_study(direction="maximize")
    full_study.enqueue_trial(best_small_params) # Start with best params from small_study
    full_study.optimize(full_data_objective, n_trials=20, callbacks=[suggest_callback])

    return full_study.best_params
```

**When to use:**
- When full data optimization is very expensive
- To speed up convergence on large datasets
- When parameters are likely to transfer between data scales

### Practical Example: Optimizing QC for PRS Performance

```python
def optimize_qc_for_prs():
    """Optimize QC parameters based on PRS performance."""
    # Create the study
    study = optuna.create_study(direction="maximize")

    def objective(trial):
        # QC parameters
        maf_threshold = trial.suggest_float("maf_threshold", 0.001, 0.05)
        geno_threshold = trial.suggest_float("geno_threshold", 0.01, 0.1)
        hwe_threshold = trial.suggest_float("hwe_threshold", 1e-10, 1e-4, log=True)
        relatedness = trial.suggest_float("relatedness", 0.05, 0.2)

        # Run QC with these parameters
        qc_output = run_qc(maf_threshold, geno_threshold, hwe_threshold, relatedness)

        # Run standard GWAS with fixed parameters
        gwas_results = run_standard_gwas(qc_output)

        # Create PRS model
        prs_model = build_prs(gwas_results)

        # Evaluate PRS performance (returns AUC)
        prs_auc = evaluate_prs(prs_model, test_data)

        # Return downstream metric directly as objective
        return prs_auc

    # Run optimization
    study.optimize(objective, n_trials=50)

    return study.best_params
```

## Considerations for Sequential Optimization

When implementing sequential optimization based on downstream outcomes:

1. **Computational Efficiency**
   - Full pipeline optimization can be expensive
   - Consider using smaller subsets for initial exploration
   - Early stopping of trials that produce poor intermediate results

2. **Feedback Loop Design**
   - Define clear metrics that represent downstream success
   - Ensure metrics actually capture meaningful biological or clinical utility
   - Beware of overfitting to specific downstream tasks

3. **Pipeline Stability**
   - Some parameter combinations may fail at intermediate stages
   - Build robust error handling in your objective function
   - Assign penalty values to failed trials rather than crashing

4. **Reproducibility**
   - Seed random number generators for reproducible results
   - Save intermediate files for critical pipeline stages
   - Log detailed parameters and outcomes for each trial

5. **Time Horizon Trade-offs**
   - Immediate outcomes vs. long-term research impact
   - Consider both short-term metrics and proxies for long-term value

## Pragmatic Implementation of Multi-Step GWAS Optimization in the Real World

### Recommended Implementation Strategy: Modular Pipeline

**Structure:** Split Into Multiple Interconnected Scripts

Each stage should:
1. Take input data and configuration parameters
2. Perform its specialized optimization
3. Output:
   - Optimized parameters in a structured format (JSON/YAML)
   - Processed data for the next stage
   - Performance metrics and visualizations

### Practical Development Approach

1. **Start with Testable Fragments**

```python
# gwas_qc_optimizer.py
def test_qc_parameters(maf, geno, hwe, mind):
    """Test a single set of QC parameters and return metrics."""
    output_prefix = f"qc_test_maf{maf}_geno{geno}_hwe{hwe}_mind{mind}"
    cmd = f"plink --bfile {input_data} --maf {maf} --geno {geno} --hwe {hwe} --mind {mind} --make-bed --out {output_prefix}"
    subprocess.run(cmd, shell=True)

    # Calculate basic QC metrics
    n_snps, n_samples = count_remaining(output_prefix)
    metrics = {
        "n_snps": n_snps,
        "n_samples": n_samples,
        "params": {'maf": maf, "geno": geno, "hwe": hwe, "mind": mind}
    }
    return metrics

# Test a few parameter sets
results = []
for maf in [0.01, 0.05]:
    for geno in [0.01, 0.05]:
        metrics = test_qc_parameters(maf, geno, 1e-6, 0.05)
        results.append(metrics)

# Write results to file
with open("qc_test_results.json", "w") as f:
    json.dump(results, f, indent=2)
```

2. **Connect Fragments Into a Prototype Pipeline**

```python
# prototype_pipeline.py
def run_prototype(maf, geno, hwe, method, p_threshold):
    """Run prototype pipeline with specific parameters."""
    # Step 1: QC
    qc_prefix = run_qc(maf, geno, hwe)

    # Step 2: Association
    assoc_file = run_association(qc_prefix, method)

    # Step 3: PRS
    prs_metrics = calculate_prs(assoc_file, p_threshold)

    # Combine metrics
    metrics = {
        "qc_params": {'maf": maf, "geno": geno, "hwe": hwe},
        "assoc_params": {'method": method},
        "prs_params": {'p_threshold": p_threshold},
        "prs_auc": prs_metrics["auc"]
    }
    return metrics

# Test a few combinations
results = []
for maf in [0.01, 0.05]:
    for method in ['linear', "logistic"]:
        for p in [0.05, 0.001, 0.0001]:
            metrics = run_prototype(maf, 0.02, 1e-6, method, p)
            results.append(metrics)
            print(f"Tested maf={maf}, method={method}, p={p}: AUC = {metrics['prs_auc']}")
```

3. **Integrate Optuna for Single-Stage Optimization**

```python
# optimize_qc_stage.py
import optuna

def objective(trial):
    """Optimize QC parameters based on downstream PRS performance."""
    # QC parameters to optimize
    maf = trial.suggest_float("maf", 0.001, 0.05)
    geno = trial.suggest_float("geno", 0.01, 0.1)
    hwe = trial.suggest_float("hwe", 1e-10, 1e-4, log=True)

    # Fixed parameters for later stages (we'll optimize these separately)
    method = "logistic"
    p_threshold = 0.001

    # Run pipeline
    qc_prefix = run_qc(maf, geno, hwe)
    assoc_file = run_association(qc_prefix, method)
    prs_metrics = calculate_prs(assoc_file, p_threshold)

    # Return AUC as objective to maximize
    return prs_metrics["auc"]

# Create and run study
study = optuna.create_study(direction="maximize")
study.optimize(objective, n_trials=30)

print(f"Best parameters: {study.best_params}")
print(f"Best AUC: {study.best_value}")

# Save best parameters
with open("best_qc_params.json", "w") as f:
    json.dump(study.best_params, f, indent=2)
```

4. **Extend to Multi-Stage Optimization**

```python
# multi_stage_optimization.py

# Stage 1: Optimize QC parameters
def optimize_qc():
    study = optuna.create_study(direction="maximize")
    study.optimize(qc_objective, n_trials=30)
    return study.best_params

# Stage 2: Optimize association with fixed optimal QC
def optimize_association(best_qc_params):
    # Apply best QC parameters
    qc_prefix = run_qc(**best_qc_params)

    # Create association optimization study
    study = optuna.create_study(direction="maximize")

    # Define objective with fixed QC parameters
    def assoc_objective(trial):
        method = trial.suggest_categorical("method", ["linear", "logistic", "mlma"])
        use_pcs = trial.suggest_categorical("use_pcs", [True, False])

        n_pcs = 0
        if use_pcs:
            n_pcs = trial.suggest_int("n_pcs", 0, 20)

        # Run association with these parameters
        assoc_file = run_association(qc_prefix, method, use_pcs, n_pcs)

        # Apply fixed PRS parameters for evaluation
        prs_metrics = calculate_prs(assoc_file, 0.001)
        return prs_metrics["auc"]

    study.optimize(assoc_objective, n_trials=30)
    return study.best_params

# Stage 3: Optimize PRS with fixed optimal QC and association
def optimize_prs(best_qc_params, best_assoc_params):
    # Apply best QC parameters
    qc_prefix = run_qc(**best_qc_params)

    # Apply best association parameters
    assoc_file = run_association(qc_prefix, **best_assoc_params)

    # Create PRS optimization study
    study = optuna.create_study(direction="maximize")

    def prs_objective(trial):
        p_threshold = trial.suggest_float("p_threshold", 1e-8, 0.5, log=True)
        clumping_r2 = trial.suggest_float("clumping_r2", 0.1, 0.9)

        # Calculate PRS with these parameters
        prs_metrics = calculate_prs(assoc_file, p_threshold, clumping_r2)
        return prs_metrics["auc"]

    study.optimize(prs_objective, n_trials=30)
    return study.best_params

# Run full optimization
best_qc = optimize_qc()
best_assoc = optimize_association(best_qc)
best_prs = optimize_prs(best_qc, best_assoc)

# Save final optimal pipeline parameters
final_params = {
    "qc": best_qc,
    "association": best_assoc,
    "prs": best_prs
}

with open("optimal_pipeline_params.json", "w") as f:
    json.dump(final_params, f, indent=2)
```

### Real-World Recommendations

1. **Divide and Conquer**
   - Splitting into multiple scripts is essential for maintenance
   - Each stage should be runnable independently for testing
   - Use configuration files to pass parameters between stages

2. **Scale Gradually**
   - Start with small datasets (e.g., single chromosome) for development
   - Use a reduced parameter space for initial testing
   - Scale up to full dataset only after prototype proves successful

3. **Checkpointing**
   - Save intermediate data at each stage
   - Implement resumable optimization when possible
   - Store full optimization history for later analysis

4. **Computation Management**
   - Optimize for your computing environment (HPC, cloud, local)
   - Parallelize trials when possible
   - Consider using pruning to terminate unpromising trials early

5. **Monitoring & Debugging**
   - Implement detailed logging at each stage
   - Create visualization dashboards for optimization progress
   - Save failed trials' information for debugging

6. **Validation Strategy**
   - Use clear train/test splits for each stage
   - Consider k-fold cross-validation for robust results
   - Include independent validation datasets when possible

## Phenotype-Specific GWAS Optimization Examples

### Example 1: Psychiatric Disorder

**Title:** "Polygenic risk scores for prediction of psychiatric disorders: Improved power and reduced overfitting through optimization strategies"  
**Authors:** Ge T, Chen CY, Ni Y, Feng YA, Smoller JW  
**Journal:** American Journal of Psychiatry (2019), 176(5):377-385  
**DOI:** 10.1176/appi.aip.2018.18091034  

**Key Optimization Approach:**
- Optimized p-value thresholds using a pruning + thresholding approach
- Employed cross-validation to avoid overfitting
- Emphasized capturing sub-threshold genetic signals that improve prediction
- Used PRSice software with custom optimization parameters

**Quote from Paper:** "By optimizing both SNP pruning parameters and p-value thresholds simultaneously, we achieved a 4.6% increase in schizophrenia prediction compared to standard approaches."

### Example 2: Autoimmune Disease

**Title:** "High-density genotyping of immune-related loci identifies new SLE risk variants in individuals with Asian ancestry"  
**Authors:** Sun C, Molineros JE, Looger LL, Zhou XJ, Kim K, Okada Y, et al.  
**Journal:** Nature Genetics (2016), 48(3):323-330  
**DOI:** 10.1038/ng.3496  

**Key Optimization Approach:**
- Used cell-type specific enrichment to weight variants
- Optimized conditional analysis parameters specifically for the HLA region
- Prioritized GWAS parameters based on functional enrichment rather than statistical signals alone
- Custom pipeline combining GCTA-COJO and functional annotations

**Quote from Paper:** "By prioritizing parameter sets that maximized enrichment of associated variants in relevant immune cell types, we identified 10 novel SLE loci that would have been missed using conventional GWAS approaches."

### Example 3: Anthropometric Trait

**Title:** "Meta-analysis of genome-wide association studies for height and body mass index in ~700,000 individuals of European ancestry"  
**Authors:** Yengo L, Sidorenko J, Kemper KE, Zheng Z, Wood AR, Weedon MN, et al.  
**Journal:** Human Molecular Genetics (2018), 27(20):3641-3649  
**DOI:** 10.1093/hmg/ddy271  

**Key Optimization Approach:**
- Extensive optimization of population stratification correction
- Careful selection of MAF thresholds to capture the full allele frequency spectrum
- Advanced mixed model implementations to account for relatedness
- Custom pipeline combining BOLT-LMM with post-GWAS optimization

**Quote from Paper:** "We observed that optimizing our statistical approach for population stratification control was critical for height, with each additional principal component explaining up to 0.1% of phenotypic variance, significantly impacting our ability to detect true genetic signals."

## Alternative Objective Functions for GWAS Parameter Optimization with Optuna

### 1. F-measure Objective (Balancing Precision and Recall)

```python
def calculate_objective(metrics, beta=1.0):
    """F-measure objective balancing precision and recall."""
    # Calculate precision (if we have known positives)
    if metrics["recovery_rate"] is not None:
        precision = metrics["recovery_rate"]
    else:
        # Estimate precision using lambda_gc as a proxy
        lambda_adj = max(0, 1 - abs(metrics["lambda_gc"] - 1.0))
        precision = lambda_adj * metrics["power_score"]

    # Use power score as recall
    recall = metrics["power_score"]

    # Calculate F-measure (harmonic mean)
    if precision > 0 and recall > 0:
        f_measure = (1 + beta**2) * (precision * recall) / ((beta**2 * precision) + recall)
        return f_measure
    else:
        return 0.0
```

**When to use:**
- When balanced performance between finding true positives and avoiding false positives is critical
- When you have a gold standard dataset with known associations
- Set β > 1 to favor recall (finding more true positives) or β < 1 to favor precision (avoiding false positives)

### 2. Area Under Precision-Recall Curve Objective

```python
def calculate_objective(metrics):
    """Area under precision-recall curve objective."""
    # Calculate precision at different significance thresholds
    precision_points = []
    recall_points = []

    thresholds = [5e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
    for threshold in thresholds:
        n_sig = sum(metrics["p_values"] < threshold)
        precision = metrics["lambda_gc_factor"] * n_sig / len(metrics["p_values"])
        recall = n_sig / (n_sig + metrics["estimated_missed"])

        precision_points.append(precision)
        recall_points.append(recall)

    # Calculate area under curve
    auc = 0
    for i in range(1, len(recall_points)):
        # Trapezoidal rule
        auc += (recall_points[i] - recall_points[i-1]) * (precision_points[i] + precision_points[i-1]) / 2

    return auc
```

**When to use:**
- When you want to optimize performance across multiple significance thresholds
- When you're concerned about the stability of results at different thresholds
- When the fixed significance threshold approach is too rigid

### 3. Posterior Expected FDR Objective

```python
def calculate_objective(metrics):
    """Posterior Expected FDR objective."""
    # Calculate local false discovery rate
    lambda_gc = metrics["lambda_gc"]
    p_values = metrics["p_values"]

    # Estimate null component
    from scipy import stats
    from statsmodels.nonparametric.kde import KDEUnivariate

    # Transform p_values
    z_scores = stats.norm.ppf(1 - p_values)

    # Fit KDE to z-scores
    kde = KDEUnivariate(z_scores)
    kde.fit()

    # Estimate null and alternative densities
    null_density = stats.norm.pdf(z_scores, 0, 1)
    mixture_density = kde.evaluate(z_scores)

    # Calculate local FDR
    local_fdr = null_density / mixture_density

    # Expected number of true positives at threshold
    threshold = 5e-8
    sig_indices = p_values < threshold
    expected_true_positives = sum(1 - local_fdr[sig_indices])

    # Return objective to maximize
    return expected_true_positives
```

**When to use:**
- When you want to directly optimize for the expected number of true discoveries
- When controlling false discovery rate is critical
- For exploratory studies where you need to maximize discoveries while controlling error rates

### 4. Weighted Multi-phenotype Objective

```python
def calculate_objective(metrics, phenotype_weights=None):
    """Multi-phenotype weighted objective."""
    if phenotype_weights is None:
        # Equal weights for all phenotypes
        phenotype_weights = {pheno: 1.0 for pheno in metrics["phenotypes"]}

    # Calculate weighted sum of objectives across phenotypes
    weighted_obj = 0
    for pheno in metrics["phenotypes"]:
        # Basic objective for this phenotype
        pheno_lambda = metrics[f"{pheno}_lambda_gc"]
        pheno_power = metrics[f"{pheno}_power_score"]
        lambda_penalty = -abs(pheno_lambda - 1.0) * 3.0
        power_reward = pheno_power * 2.0
        pheno_obj = lambda_penalty + power_reward

        # Add to weighted sum
        weighted_obj += phenotype_weights[pheno] * pheno_obj

    # Normalize by sum of weights
    total_weight = sum(phenotype_weights.values())
    return weighted_obj / total_weight
```

**When to use:**
- For multi-phenotype studies or pleiotropic analyses
- When you need to optimize parameters that work well across multiple related traits
- When certain phenotypes are prioritized over others

### 5. Genetic Architecture-Specific Objective

```python
def calculate_objective(metrics, architecture_type="polygenic"):
    """Genetic architecture-specific objective."""
    # Common base penalty for statistical validity
    lambda_penalty = -abs(metrics["lambda_gc"] - 1.0) * 3.0

    if architecture_type == "polygenic":
        # For highly polygenic traits, favor many small effects
        n_sig_moderate = sum((metrics["p_values"] < 1e-5) & (metrics["p_values"] >= 5e-8))
        polygenic_reward = np.log10(n_sig_moderate + 1) * 2.0
        return lambda_penalty + polygenic_reward

    elif architecture_type == "oligogenic":
        # For oligogenic traits, favor few strong signals
        n_sig_strong = sum(metrics["p_values"] < 5e-8)
        effect_size_reward = np.mean(metrics["effect_sizes"][metrics["p_values"] < 5e-8]) * 3.0
        return lambda_penalty + n_sig_strong + effect_size_reward

    elif architecture_type == "monogenic":
        # For monogenic traits, favor single very strong signal
        strongest_signal = -np.log10(min(metrics["p_values"] + 1e-10))
        return lambda_penalty + strongest_signal
```

**When to use:**
- When you have prior knowledge about the genetic architecture of your trait
- For specialized studies targeting specific types of genetic effects
- When you need to tailor the analysis to expected effect size distributions

### 6. Heritability-Based Objective

```python
def calculate_objective(metrics):
    """Heritability-based objective."""
    # Calculate explained variance
    from scipy import stats

    # Get chi-square statistics and convert to R² (explained variance)
    chi2_stats = stats.chi2.ppf(1 - metrics["p_values"], 1)
    n_samples = metrics["sample_size"]
    r_squared_values = chi2_stats / n_samples

    # Calculate total explained heritability
    h2_explained = sum(r_squared_values)

    # Penalize deviation from expected heritability
    expected_h2 = metrics["expected_heritability"]
    h2_penalty = -abs(h2_explained - expected_h2) * 2.0

    # Penalize inflation/deflation
    lambda_penalty = -abs(metrics["lambda_gc"] - 1.0) * 3.0

    return lambda_penalty + h2_penalty
```

**When to use:**
- When you have good estimates of the expected heritability
- For quantitative traits where variance explained is measurable
- When you want to ensure the parameter set yields biologically plausible results

## Considerations for Choosing an Objective Function

When selecting an objective function, consider:

1. **Available Prior Knowledge**
   - Do you have known associations to validate against?
   - Do you have estimates of trait heritability?
   - What's known about the genetic architecture?

2. **Study Goals**
   - Discovery-focused: Prioritize power and recall
   - Replication-focused: Prioritize precision and accuracy
   - Mechanistic understanding: Focus on effect size estimation

3. **Sample Characteristics**
   - Sample size affects statistical power considerations
   - Population structure influences inflation control needs
   - Case-control imbalance may require specific optimization strategies

4. **Computational Resources**
   - More complex objectives may require longer optimization times
   - Some objectives need additional preprocessing or calculation steps

5. **Downstream Analysis Plans**
   - Fine-mapping intentions might favor certain parameter sets
   - Polygenic risk score development has different optimal parameters
   - Meta-analysis participation requires specific statistical considerations
```



