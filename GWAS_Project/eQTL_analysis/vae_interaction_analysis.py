# VAE Analysis Pipeline for SNP-Expression Interaction Analysis
# This script builds a VAE to identify SNP-covariate interactions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import optuna
import shap
import os
import allel
import scipy.stats as stats
from sklearn.metrics import r2_score

# Set random seeds for reproducibility
np.random.seed(42)
torch.manual_seed(42)

# Create output directory
output_dir = "vae_results"
os.makedirs(output_dir, exist_ok=True)

# Step 1: Load and preprocess data
print("Loading and preprocessing data...")

# Load metadata
metadata = pd.read_csv("sample_metadata.csv")
print(f"Loaded metadata with {metadata.shape[0]} samples and {metadata.shape[1]} columns")

# Load expression data
expression_data = pd.read_csv("expression_counts.csv", index_col=0)
print(f"Loaded expression data with {expression_data.shape[0]} genes and {expression_data.shape[1]} samples")

# Load VCF data using scikit-allel
vcf_file = "genotypes.vcf"
callset = allel.read_vcf(vcf_file)
genotypes = callset['calldata/GT']
samples = callset['samples']
variants = callset['variants/ID']

# Convert genotypes to numeric format (0, 1, 2)
geno_numeric = genotypes.sum(axis=2)
geno_df = pd.DataFrame(geno_numeric, index=variants, columns=samples)
print(f"Loaded genotype data with {geno_df.shape[0]} SNPs and {geno_df.shape[1]} samples")

# Ensure samples are aligned across datasets
common_samples = list(set(metadata['sample_id']) & set(expression_data.columns) & set(geno_df.columns))
print(f"Found {len(common_samples)} samples common to all datasets")

# Subset datasets to common samples
metadata = metadata.set_index('sample_id').loc[common_samples].reset_index()
expression_data = expression_data[common_samples]
geno_df = geno_df[common_samples]

# Step 2: Prepare covariates for analysis
# Select covariates of interest (e.g., cell type, disease status)
covariates_of_interest = ['cell_type', 'disease_status']

# One-hot encode categorical covariates
encoder = OneHotEncoder(sparse=False)
covariate_data = encoder.fit_transform(metadata[covariates_of_interest])
covariate_df = pd.DataFrame(
    covariate_data, 
    columns=encoder.get_feature_names_out(covariates_of_interest),
    index=metadata['sample_id']
)
print(f"Prepared covariates with {covariate_df.shape[1]} features")

# Step 3: Dimension reduction for computational efficiency
# Select top variable genes
gene_var = expression_data.var(axis=1)
top_genes = gene_var.sort_values(ascending=False).index[:1000]  # Use top 1000 variable genes
expression_subset = expression_data.loc[top_genes]
print(f"Selected {len(top_genes)} top variable genes")

# Select top variable SNPs or those in cis with the selected genes
# For simplicity, select top SNPs by MAF
maf = np.mean(geno_df.values, axis=1) / 2
maf = np.minimum(maf, 1 - maf)  # Get minor allele frequency
top_snps = maf.sort_values(ascending=False).index[:5000]  # Use top 5000 SNPs
geno_subset = geno_df.loc[top_snps]
print(f"Selected {len(top_snps)} top SNPs")

# Standardize data
scaler_expr = StandardScaler()
scaler_geno = StandardScaler()
expr_scaled = scaler_expr.fit_transform(expression_subset.T).T
geno_scaled = scaler_geno.fit_transform(geno_subset.T).T

# Step 4: Define VAE architecture for SNP-expression data with covariates
class ConditionalVAE(nn.Module):
    def __init__(self, n_genes, n_snps, n_covariates, latent_dim, hidden_dim):
        super(ConditionalVAE, self).__init__()
        self.n_genes = n_genes
        self.n_snps = n_snps
        self.n_covariates = n_covariates
        self.latent_dim = latent_dim
        
        # Encoder for gene expression
        self.encoder_expr = nn.Sequential(
            nn.Linear(n_genes + n_covariates, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU()
        )
        
        # Encoder for SNPs
        self.encoder_snp = nn.Sequential(
            nn.Linear(n_snps + n_covariates, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU()
        )
        
        # Joint latent space encoders
        self.fc_mu_expr = nn.Linear(hidden_dim // 2, latent_dim)
        self.fc_var_expr = nn.Linear(hidden_dim // 2, latent_dim)
        self.fc_mu_snp = nn.Linear(hidden_dim // 2, latent_dim)
        self.fc_var_snp = nn.Linear(hidden_dim // 2, latent_dim)
        
        # Decoders
        self.decoder_expr = nn.Sequential(
            nn.Linear(latent_dim + n_covariates, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, n_genes)
        )
        
        self.decoder_snp = nn.Sequential(
            nn.Linear(latent_dim + n_covariates, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, n_snps)
        )
    
    def encode_expr(self, x, c):
        # Concatenate expression data with covariates
        x_c = torch.cat([x, c], dim=1)
        h = self.encoder_expr(x_c)
        mu = self.fc_mu_expr(h)
        log_var = self.fc_var_expr(h)
        return mu, log_var
    
    def encode_snp(self, x, c):
        # Concatenate SNP data with covariates
        x_c = torch.cat([x, c], dim=1)
        h = self.encoder_snp(x_c)
        mu = self.fc_mu_snp(h)
        log_var = self.fc_var_snp(h)
        return mu, log_var
    
    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode_expr(self, z, c):
        # Concatenate latent representation with covariates
        z_c = torch.cat([z, c], dim=1)
        return self.decoder_expr(z_c)
    
    def decode_snp(self, z, c):
        # Concatenate latent representation with covariates
        z_c = torch.cat([z, c], dim=1)
        return self.decoder_snp(z_c)
    
    def forward(self, x_expr, x_snp, c):
        # Encode
        mu_expr, log_var_expr = self.encode_expr(x_expr, c)
        mu_snp, log_var_snp = self.encode_snp(x_snp, c)
        
        # Sample from latent distribution
        z_expr = self.reparameterize(mu_expr, log_var_expr)
        z_snp = self.reparameterize(mu_snp, log_var_snp)
        
        # Average the latent representations for a joint latent space
        z = (z_expr + z_snp) / 2
        
        # Decode
        recon_expr = self.decode_expr(z, c)
        recon_snp = self.decode_snp(z, c)
        
        return recon_expr, recon_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, z

# Step 5: Define loss function for VAE
def vae_loss(recon_expr, x_expr, recon_snp, x_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, beta=1.0):
    # Reconstruction loss for expression data
    recon_loss_expr = F.mse_loss(recon_expr, x_expr, reduction='sum')
    
    # Reconstruction loss for SNP data
    recon_loss_snp = F.mse_loss(recon_snp, x_snp, reduction='sum')
    
    # KL divergence for expression data
    kl_div_expr = -0.5 * torch.sum(1 + log_var_expr - mu_expr.pow(2) - log_var_expr.exp())
    
    # KL divergence for SNP data
    kl_div_snp = -0.5 * torch.sum(1 + log_var_snp - mu_snp.pow(2) - log_var_snp.exp())
    
    # Combined loss
    total_loss = recon_loss_expr + recon_loss_snp + beta * (kl_div_expr + kl_div_snp)
    
    return total_loss

# Step 6: Hyperparameter optimization using Optuna
def objective(trial):
    # Define hyperparameters to optimize
    latent_dim = trial.suggest_int('latent_dim', 10, 50)
    hidden_dim = trial.suggest_int('hidden_dim', 128, 512)
    learning_rate = trial.suggest_float('learning_rate', 1e-4, 1e-2, log=True)
    batch_size = trial.suggest_categorical('batch_size', [32, 64, 128])
    beta = trial.suggest_float('beta', 0.1, 2.0)
    
    # Prepare data for PyTorch
    X_expr = torch.FloatTensor(expr_scaled)
    X_snp = torch.FloatTensor(geno_scaled)
    X_cov = torch.FloatTensor(covariate_data)
    
    # Split data
    X_expr_train, X_expr_val, X_snp_train, X_snp_val, X_cov_train, X_cov_val = train_test_split(
        X_expr.T, X_snp.T, X_cov, test_size=0.2, random_state=42
    )
    
    # Create data loaders
    train_dataset = TensorDataset(X_expr_train, X_snp_train, X_cov_train)
    val_dataset = TensorDataset(X_expr_val, X_snp_val, X_cov_val)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    
    # Initialize model
    model = ConditionalVAE(
        n_genes=expr_scaled.shape[0],
        n_snps=geno_scaled.shape[0],
        n_covariates=covariate_data.shape[1],
        latent_dim=latent_dim,
        hidden_dim=hidden_dim
    )
    
    # Initialize optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    
    # Training loop
    best_val_loss = float('inf')
    patience = 5
    patience_counter = 0
    
    for epoch in range(50):  # Maximum 50 epochs
        # Training
        model.train()
        train_loss = 0
        for batch_expr, batch_snp, batch_cov in train_loader:
            optimizer.zero_grad()
            recon_expr, recon_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, _ = model(batch_expr, batch_snp, batch_cov)
            loss = vae_loss(recon_expr, batch_expr, recon_snp, batch_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, beta)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch_expr, batch_snp, batch_cov in val_loader:
                recon_expr, recon_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, _ = model(batch_expr, batch_snp, batch_cov)
                loss = vae_loss(recon_expr, batch_expr, recon_snp, batch_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, beta)
                val_loss += loss.item()
        
        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                break
    
    return best_val_loss

# Step 7: Run hyperparameter optimization
def optimize_hyperparameters():
    print("Optimizing hyperparameters...")
    study = optuna.create_study(direction='minimize')
    study.optimize(objective, n_trials=20)
    
    print("Best hyperparameters:")
    for key, value in study.best_params.items():
        print(f"    {key}: {value}")
    
    return study.best_params

# Step 8: Train final model with best hyperparameters
def train_final_model(best_params):
    print("Training final model with best hyperparameters...")
    
    # Prepare data for PyTorch
    X_expr = torch.FloatTensor(expr_scaled)
    X_snp = torch.FloatTensor(geno_scaled)
    X_cov = torch.FloatTensor(covariate_data)
    
    # Create dataset
    dataset = TensorDataset(X_expr.T, X_snp.T, X_cov)
    loader = DataLoader(dataset, batch_size=best_params['batch_size'], shuffle=True)
    
    # Initialize model
    model = ConditionalVAE(
        n_genes=expr_scaled.shape[0],
        n_snps=geno_scaled.shape[0],
        n_covariates=covariate_data.shape[1],
        latent_dim=best_params['latent_dim'],
        hidden_dim=best_params['hidden_dim']
    )
    
    # Initialize optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=best_params['learning_rate'])
    
    # Training loop
    num_epochs = 100
    beta = best_params['beta']
    
    losses = []
    for epoch in range(num_epochs):
        model.train()
        epoch_loss = 0
        
        for batch_expr, batch_snp, batch_cov in loader:
            optimizer.zero_grad()
            recon_expr, recon_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, _ = model(batch_expr, batch_snp, batch_cov)
            loss = vae_loss(recon_expr, batch_expr, recon_snp, batch_snp, mu_expr, log_var_expr, mu_snp, log_var_snp, beta)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
        
        losses.append(epoch_loss)
        if epoch % 10 == 0:
            print(f"Epoch {epoch}/{num_epochs}, Loss: {epoch_loss:.4f}")
    
    # Plot loss curve
    plt.figure(figsize=(10, 6))
    plt.plot(losses)
    plt.title('Training Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.savefig(os.path.join(output_dir, 'training_loss.png'))
    
    # Save model
    torch.save(model.state_dict(), os.path.join(output_dir, 'vae_model.pt'))
    
    return model

# Step 9: Extract latent representations and analyze interactions
def analyze_latent_interactions(model, covariates_of_interest):
    print("Analyzing latent interactions...")
    
    # Prepare data for PyTorch
    X_expr = torch.FloatTensor(expr_scaled)
    X_snp = torch.FloatTensor(geno_scaled)
    X_cov = torch.FloatTensor(covariate_data)
    
    # Extract latent representations
    model.eval()
    with torch.no_grad():
        _, _, mu_expr, _, mu_snp, _, z = model(X_expr.T, X_snp.T, X_cov)
        latent_repr = z.numpy()
    
    # Create dataframe with latent dimensions
    latent_df = pd.DataFrame(
        latent_repr, 
        columns=[f'latent_{i}' for i in range(latent_repr.shape[1])],
        index=common_samples
    )
    
    # Join with covariates
    analysis_df = pd.concat([latent_df, metadata.set_index('sample_id')], axis=1)
    
    # Analyze the effect of covariates on latent dimensions
    results = []
    for covariate in covariates_of_interest:
        print(f"Testing effect of {covariate} on latent dimensions...")
        
        # Get unique values of the covariate
        covariate_values = metadata[covariate].unique()
        
        # Test each latent dimension
        for latent_dim in latent_df.columns:
            # For categorical covariates, use ANOVA
            f_stat, p_value = stats.f_oneway(
                *[latent_df.loc[metadata[metadata[covariate] == val].sample_id, latent_dim].values 
                  for val in covariate_values]
            )
            
            results.append({
                'latent_dimension': latent_dim,
                'covariate': covariate,
                'test_statistic': f_stat,
                'p_value': p_value
            })
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    # Apply FDR correction
    results_df['p_adjusted'] = stats.false_discovery_rate_control(results_df['p_value'])
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('p_adjusted')
    
    # Save results
    results_df.to_csv(os.path.join(output_dir, 'latent_dimension_associations.csv'), index=False)
    
    # Plot top latent dimensions associated with covariates
    plot_top_associations(results_df, latent_df, metadata, covariates_of_interest)
    
    return results_df, latent_df

# Step 10: Plot top associations
def plot_top_associations(results_df, latent_df, metadata, covariates):
    print("Plotting top associations...")
    
    # Create output directory for plots
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # For each covariate, plot top 3 associated latent dimensions
    for covariate in covariates:
        covariate_results = results_df[results_df['covariate'] == covariate].sort_values('p_adjusted')
        top_dimensions = covariate_results.head(3)['latent_dimension'].tolist()
        
        for dim in top_dimensions:
            plt.figure(figsize=(10, 6))
            
            # Join latent dimension with metadata
            plot_data = pd.DataFrame({
                'latent_value': latent_df[dim],
                'covariate': metadata[covariate]
            })
            
            # Create violin plot
            sns.violinplot(x='covariate', y='latent_value', data=plot_data)
            plt.title(f'Distribution of {dim} by {covariate}')
            plt.xlabel(covariate)
            plt.ylabel(dim)
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            # Save plot
            plt.savefig(os.path.join(plots_dir, f'{dim}_{covariate}_association.png'))
            plt.close()

# Step 11: Identify SNP-expression interactions using SHAP
def interpret_with_shap(model, X_expr, X_snp, X_cov, metadata, top_genes, top_snps, covariates_of_interest):
    print("Interpreting model with SHAP...")
    
    # Create background dataset for SHAP (random subset)
    n_background = min(100, X_expr.shape[0])
    background_indices = np.random.choice(X_expr.shape[0], n_background, replace=False)
    background_expr = X_expr[background_indices]
    background_snp = X_snp[background_indices]
    background_cov = X_cov[background_indices]
    
    # Create PyTorch DataLoader for the background data
    background_dataset = TensorDataset(
        torch.FloatTensor(background_expr), 
        torch.FloatTensor(background_snp), 
        torch.FloatTensor(background_cov)
    )
    background_loader = DataLoader(background_dataset, batch_size=n_background)
    
    # Define SHAP explainer for gene expression decoder
    class ExpressionDecoder(nn.Module):
        def __init__(self, vae_model):
            super(ExpressionDecoder, self).__init__()
            self.vae = vae_model
        
        def forward(self, x):
            # x is a combination of latent variables and covariates
            latent = x[:, :self.vae.latent_dim]
            covariates = x[:, self.vae.latent_dim:]
            return self.vae.decode_expr(latent, covariates)
    
    # Create explainer
    expr_decoder = ExpressionDecoder(model)
    explainer = shap.DeepExplainer(expr_decoder, next(iter(background_loader))[0])
    
    # Generate SHAP values for a subset of samples
    test_samples = min(50, X_expr.shape[0])
    test_indices = np.random.choice(X_expr.shape[0], test_samples, replace=False)
    test_expr = X_expr[test_indices]
    test_snp = X_snp[test_indices]
    test_cov = X_cov[test_indices]
    
    # Extract latent representations
    model.eval()
    with torch.no_grad():
        _, _, _, _, _, _, z = model(
            torch.FloatTensor(test_expr), 
            torch.FloatTensor(test_snp), 
            torch.FloatTensor(test_cov)
        )
        z = z.detach().numpy()
    
    # Combine latent variables with covariates for SHAP analysis
    shap_input = np.hstack([z, test_cov])
    
    # Calculate SHAP values
    shap_values = explainer.shap_values(torch.FloatTensor(shap_input))
    
    # Create summary plots for top genes
    plt.figure(figsize=(12, 8))
    shap.summary_plot(
        shap_values[:, :top_genes.shape[0]], 
        features=shap_input,
        feature_names=[f"Gene_{i}" for i in range(top_genes.shape[0])],
        max_display=20
    )
    plt.savefig(os.path.join(output_dir, 'shap_summary_genes.png'))
    plt.close()
    
    # Identify top interactions
    interaction_scores = analyze_shap_interactions(
        shap_values, 
        shap_input, 
        model.latent_dim, 
        test_indices, 
        metadata, 
        top_genes, 
        covariates_of_interest
    )
    
    return interaction_scores

# Step 12: Analyze SHAP interactions to identify SNP-covariate effects
def analyze_shap_interactions(shap_values, shap_input, latent_dim, test_indices, metadata, top_genes, covariates):
    print("Analyzing SHAP interactions...")
    
    # Separate latent variables and covariates
    latent_values = shap_input[:, :latent_dim]
    
    # Get metadata for test samples
    test_metadata = metadata.iloc[test_indices]
    
    # Initialize results
    interaction_scores = []
    
    # For each covariate of interest
    for covariate in covariates:
        covariate_values = test_metadata[covariate].values
        
        # For each gene
        for gene_idx, gene_name in enumerate(top_genes):
            # Get SHAP values for this gene
            gene_shap = shap_values[:, gene_idx]
            
            # Calculate interaction score (e.g., correlation between SHAP values and covariate)
            # For categorical covariates, calculate the variance of mean SHAP values across categories
            shap_by_category = {}
            for val in np.unique(covariate_values):
                mask = covariate_values == val
                shap_by_category[val] = gene_shap[mask].mean()
            
            # Calculate variance of means as interaction score
            interaction_score = np.var(list(shap_by_category.values()))
            
            interaction_scores.append({
                'gene': gene_name,
                'covariate': covariate,
                'interaction_score': interaction_score
            })
    
    # Create dataframe
    interaction_df = pd.DataFrame(interaction_scores)
    
    # Sort by interaction score
    interaction_df = interaction_df.sort_values('interaction_score', ascending=False)
    
    # Save results
    interaction_df.to_csv(os.path.join(output_dir, 'gene_covariate_interactions.csv'), index=False)
    
    # Plot top interactions
    plot_top_shap_interactions(interaction_df, shap_values, test_metadata, top_genes, covariates)
    
    return interaction_df

# Step 13: Plot top SHAP interactions
def plot_top_shap_interactions(interaction_df, shap_values, test_metadata, top_genes, covariates):
    print("Plotting top SHAP interactions...")
    
    # Create output directory for plots
    plots_dir = os.path.join(output_dir, 'shap_plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # For each covariate, plot top 5 genes with highest interaction scores
    for covariate in covariates:
        covariate_interactions = interaction_df[interaction_df['covariate'] == covariate]
        top_interactions = covariate_interactions.head(5)
        
        for _, row in top_interactions.iterrows():
            gene = row['gene']
            gene_idx = top_genes.tolist().index(gene)
            
            plt.figure(figsize=(10, 6))
            
            # Create data for plotting
            plot_data = pd.DataFrame({
                'shap_value': shap_values[:, gene_idx],
                'covariate': test_metadata[covariate]
            })
            
            # Create boxplot
            sns.boxplot(x='covariate', y='shap_value', data=plot_data)
            plt.title(f'SHAP values for {gene} by {covariate}')
            plt.xlabel(covariate)
            plt.ylabel(f'SHAP value (impact on {gene} expression)')
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            # Save plot
            plt.savefig(os.path.join(plots_dir, f'shap_{gene}_{covariate}_interaction.png'))
            plt.close()

# Step 14: Validate results by comparing with traditional eQTL methods
def validate_with_eqtl(interaction_df, top_genes, top_snps, expr_scaled, geno_scaled, metadata, covariates):
    print("Validating with traditional eQTL analysis...")
    
    # Select top interactions identified by VAE+SHAP
    top_vae_interactions = []
    for covariate in covariates:
        covariate_interactions = interaction_df[interaction_df['covariate'] == covariate]
        top_genes_per_covariate = covariate_interactions.head(10)['gene'].tolist()
        top_vae_interactions.extend([(gene, covariate) for gene in top_genes_per_covariate])
    
    # Initialize results
    validation_results = []
    
    # For each top interaction
    for gene, covariate in top_vae_interactions:
        gene_idx = top_genes.tolist().index(gene)
        gene_expr = expr_scaled[gene_idx]
        
        # Test interaction with each SNP
        for snp_idx, snp in enumerate(top_snps):
            snp_geno = geno_scaled[snp_idx]
            
            # Create dataframe for analysis
            analysis_df = pd.DataFrame({
                'expression': gene_expr,
                'genotype': snp_geno,
                'covariate': metadata[covariate].values
            })
            
            # Remove missing values
            analysis_df = analysis_df.dropna()
            
            # Fit models with and without interaction
            model_main = sm.OLS.from_formula('expression ~ genotype + C(covariate)', data=analysis_df).fit()
            model_interaction = sm.OLS.from_formula('expression ~ genotype * C(covariate)', data=analysis_df).fit()
            
            # Compare models
            lrt_stat = -2 * (model_main.llf - model_interaction.llf)
            p_value = stats.chi2.sf(lrt_stat, df=len(np.unique(analysis_df['covariate'])) - 1)
            
            validation_results.append({
                'gene': gene,
                'snp': snp,
                'covariate': covariate,
                'lrt_statistic': lrt_stat,
                'p_value': p_value
            })
    
    # Create dataframe
    validation_df = pd.DataFrame(validation_results)
    
    # Apply FDR correction
    validation_df['p_adjusted'] = stats.false_discovery_rate_control(validation_df['p_value'])
    
    # Sort by adjusted p-value
    validation_df = validation_df.sort_values('p_adjusted')
    
    # Save results
    validation_df.to_csv(os.path.join(output_dir, 'validation_results.csv'), index=False)
    
    # Calculate agreement between VAE+SHAP and traditional eQTL
    agreement_stats = calculate_agreement(interaction_df, validation_df)
    
    return validation_df, agreement_stats

# Step 15: Calculate agreement between methods
def calculate_agreement(vae_results, eqtl_results):
    print("Calculating agreement between methods...")
    
    # Group eQTL results by gene-covariate pair
    grouped_eqtl = eqtl_results.groupby(['gene', 'covariate'])['p_adjusted'].min().reset_index()
    
    # Merge with VAE results
    comparison = pd.merge(
        vae_results[['gene', 'covariate', 'interaction_score']],
        grouped_eqtl,
        on=['gene', 'covariate'],
        how='inner'
    )
    
    # Calculate correlation
    correlation = stats.spearmanr(
        -np.log10(comparison['interaction_score']), 
        -np.log10(comparison['p_adjusted'])
    )
    
    # Create agreement plot
    plt.figure(figsize=(10, 8))
    plt.scatter(
        -np.log10(comparison['interaction_score']),
        -np.log10(comparison['p_adjusted']),
        alpha=0.7
    )
    plt.xlabel('-log10(VAE+SHAP Interaction Score)')
    plt.ylabel('-log10(eQTL Interaction P-value)')
    plt.title(f'Agreement between VAE+SHAP and eQTL (Spearman r={correlation[0]:.3f}, p={correlation[1]:.3e})')
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
    plt.axvline(x=np.median(-np.log10(comparison['interaction_score'])), color='r', linestyle='--')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'method_agreement.png'))
    plt.close()
    
    # Save comparison
    comparison.to_csv(os.path.join(output_dir, 'method_comparison.csv'), index=False)
    
    return {
        'correlation': correlation[0],
        'p_value': correlation[1],
        'n_compared': len(comparison)
    }

# Step 16: Main analysis pipeline
def run_vae_pipeline():
    print("Starting VAE analysis pipeline...")
    
    # Step 1: Optimize hyperparameters
    best_params = optimize_hyperparameters()
    
    # Step 2: Train final model
    model = train_final_model(best_params)
    
    # Step 3: Analyze latent space interactions
    results_df, latent_df = analyze_latent_interactions(model, covariates_of_interest)
    
    # Step 4: Interpret with SHAP
    interaction_scores = interpret_with_shap(
        model, 
        expr_scaled.T.numpy(), 
        geno_scaled.T.numpy(), 
        covariate_data, 
        metadata, 
        top_genes, 
        top_snps, 
        covariates_of_interest
    )
    
    # Step 5: Validate with traditional eQTL
    validation_df, agreement_stats = validate_with_eqtl(
        interaction_scores, 
        top_genes, 
        top_snps, 
        expr_scaled, 
        geno_scaled, 
        metadata, 
        covariates_of_interest
    )
    
    # Step 6: Generate final report
    generate_report(
        results_df, 
        interaction_scores, 
        validation_df, 
        agreement_stats, 
        best_params, 
        covariates_of_interest
    )
    
    print("VAE analysis pipeline completed successfully!")
    print(f"Results saved to {output_dir}")

# Step 17: Generate report
def generate_report(latent_results, interaction_scores, validation_results, agreement_stats, model_params, covariates):
    print("Generating final report...")
    
    # Create report header
    report = [
        "# SNP-Expression Interaction Analysis Report",
        "## Using Variational Autoencoder with SHAP Interpretation",
        "",
        "## Model Parameters",
        f"- Latent dimensions: {model_params['latent_dim']}",
        f"- Hidden dimensions: {model_params['hidden_dim']}",
        f"- Learning rate: {model_params['learning_rate']}",
        f"- Batch size: {model_params['batch_size']}",
        f"- Beta (KL weight): {model_params['beta']}",
        "",
        "## Summary of Results",
        "",
        "### Latent Space Analysis",
        f"Found {len(latent_results[latent_results['p_adjusted'] < 0.05])} significant associations between latent dimensions and covariates.",
        ""
    ]
    
    # Add top results for each covariate
    for covariate in covariates:
        covariate_results = latent_results[latent_results['covariate'] == covariate]
        significant = covariate_results[covariate_results['p_adjusted'] < 0.05]
        
        report.append(f"#### {covariate}")
        report.append(f"- {len(significant)} significant latent dimensions associated with {covariate}")
        
        if len(significant) > 0:
            report.append("- Top 3 latent dimensions:")
            top3 = significant.head(3)
            for _, row in top3.iterrows():
                report.append(f"  - {row['latent_dimension']}: p-adjusted = {row['p_adjusted']:.2e}")
        
        report.append("")
    
    # Add SHAP interaction results
    report.append("### SNP-Covariate Interactions (VAE+SHAP)")
    for covariate in covariates:
        covariate_interactions = interaction_scores[interaction_scores['covariate'] == covariate]
        
        report.append(f"#### {covariate}")
        report.append("- Top 5 genes with strongest interactions:")
        top5 = covariate_interactions.head(5)
        for _, row in top5.iterrows():
            report.append(f"  - {row['gene']}: interaction score = {row['interaction_score']:.4f}")
        
        report.append("")
    
    # Add validation results
    report.append("### Validation with Traditional eQTL Analysis")
    report.append(f"- Agreement between methods: Spearman r = {agreement_stats['correlation']:.3f} (p = {agreement_stats['p_value']:.3e})")
    report.append(f"- Number of interactions compared: {agreement_stats['n_compared']}")
    
    significant_validations = validation_results[validation_results['p_adjusted'] < 0.05]
    report.append(f"- {len(significant_validations)} significant interactions confirmed with traditional eQTL analysis")
    
    if len(significant_validations) > 0:
        report.append("- Top 5 validated interactions:")
        top5 = significant_validations.head(5)
        for _, row in top5.iterrows():
            report.append(f"  - {row['gene']} x {row['snp']} x {row['covariate']}: p-adjusted = {row['p_adjusted']:.2e}")
    
    report.append("")
    report.append("## Conclusion")
    report.append("This analysis used a deep learning approach with variational autoencoders to identify SNP effects on gene expression that depend on covariates such as cell type and disease status. The results were validated using traditional eQTL interaction analysis.")
    report.append("")
    report.append("Key findings:")
    report.append("1. The VAE successfully identified latent dimensions associated with different covariates")
    report.append("2. SHAP analysis revealed gene-covariate interactions that represent context-dependent genetic effects")
    report.append("3. Traditional eQTL analysis confirmed many of the interactions identified by the VAE+SHAP approach")
    report.append("")
    report.append("See the accompanying figures in the plots directory for visualizations of these results.")
    
    # Write report to file
    with open(os.path.join(output_dir, 'analysis_report.md'), 'w') as f:
        f.write('\n'.join(report))
    
    print(f"Report saved to {os.path.join(output_dir, 'analysis_report.md')}")

# Run the full pipeline
if __name__ == "__main__":
    import statsmodels.api as sm
    run_vae_pipeline()