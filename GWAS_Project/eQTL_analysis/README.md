# SNP-Expression Interaction Analysis Roadmap

This document outlines two complementary approaches for testing SNP effects on expression that depend on covariates (e.g., cell type, disease status).

## Approach 1: Traditional eQTL Interaction Analysis

### Overview
Traditional eQTL interaction analysis uses established statistical methods to directly test for interactions between SNPs, gene expression, and covariates.

### Implementation Steps

1. **Data Preparation**
   - Load normalized VCF file (genotypes)
   - Load normalized expression count matrix
   - Load sample metadata with covariates
   - Ensure sample IDs match across datasets
   - Handle missing values

2. **Statistical Testing**
   - Use lme4 package in R for linear mixed models
   - For each gene-SNP pair, test models with and without interaction terms
   - Compare models to determine significance of interaction effects
   - Apply multiple testing correction (Benjamini-Hochberg)

3. **Result Interpretation**
   - Sort interactions by adjusted p-values
   - Examine top significant interactions
   - Visualize interaction effects with boxplots/violin plots
   - Interpret effect sizes and directions

4. **Advantages**
   - Well-established, interpretable statistical framework
   - Direct test of interaction hypotheses
   - Provides clear p-values and effect sizes
   - Computationally efficient for targeted analyses
   - Easily interpretable results

5. **Limitations**
   - May miss complex, non-linear interactions
   - Computationally intensive for genome-wide analysis
   - Tests each interaction independently without leveraging shared patterns
   - May have limited power for rare variants or small effect sizes

## Approach 2: Variational Autoencoder (VAE) with SHAP

### Overview
This deep learning approach uses VAEs to learn latent representations of the data and SHAP (SHapley Additive exPlanations) to interpret how covariates influence the relationships between SNPs and gene expression.

### Implementation Steps

1. **Data Preparation**
   - Same initial steps as traditional approach
   - Select top variable genes and SNPs to manage computational complexity
   - Standardize/normalize data for neural network training

2. **VAE Architecture Design**
   - Design a conditional VAE that incorporates covariates
   - Create shared latent space for genotype and expression data
   - Include encoding and decoding networks for both data types

3. **Hyperparameter Optimization**
   - Use Optuna for Bayesian optimization of:
     - Latent dimension size
     - Hidden layer dimensions
     - Learning rate
     - Batch size
     - Beta parameter (KL divergence weight)

4. **Model Training**
   - Train the VAE to reconstruct both expression and genotype data
   - Monitor loss curves to ensure proper convergence
   - Extract latent representations for downstream analysis

5. **Latent Space Analysis**
   - Test associations between latent dimensions and covariates
   - Identify latent dimensions capturing interaction effects
   - Visualize latent space clustering by covariates

6. **SHAP Interpretation**
   - Apply SHAP to understand how latent representations and covariates influence predictions
   - Identify gene-covariate interactions with highest SHAP values
   - Create SHAP summary and dependency plots to visualize interactions

7. **Validation**
   - Compare with traditional eQTL results
   - Calculate agreement between methods
   - Identify interactions found by both approaches

8. **Advantages**
   - Can capture complex, non-linear relationships
   - Learns from patterns across all genes and SNPs simultaneously
   - May discover novel interactions missed by traditional methods
   - Potentially higher power for detecting weak but consistent signals
   - Reduces dimensionality of the problem

9. **Limitations**
   - Interpretability is more complex
   - Requires more computational resources
   - Hyperparameter tuning is crucial for performance
   - May not provide standard statistical measures (p-values)
   - Needs careful validation against established methods

## Comprehensive Analysis Plan

For the most robust analysis, we recommend:

1. **Start with traditional analysis for key genes/SNPs of interest**
   - Establish baseline results with well-understood methods
   - Generate clear statistical measures for primary hypotheses

2. **Apply VAE+SHAP approach to discover additional patterns**
   - Use deep learning to identify complex relationships
   - Leverage SHAP to interpret the model's findings

3. **Cross-validate findings between methods**
   - Prioritize interactions found by both approaches
   - Investigate discrepancies to understand methodological differences

4. **Functional interpretation**
   - Pathway analysis of genes with significant interactions
   - Literature review for biological context
   - Consider experimental validation of top findings

## Implementation Resources

Both approaches have been implemented in the provided code:

1. **Traditional eQTL**: R script using lme4 package
2. **VAE+SHAP**: Python implementation with PyTorch, Optuna, and SHAP

## Expected Outcomes

The combined analysis will provide:

1. Statistically rigorous tests of specific interaction hypotheses
2. Discovery of complex patterns and novel interactions
3. Visualizations highlighting how genetic effects vary across contexts
4. A ranked list of candidate interactions for further investigation
5. Insights into the biological mechanisms underlying context-dependent genetic regulation

By leveraging both approaches, we can gain a more comprehensive understanding of how genetic effects on gene expression depend on biological context.
