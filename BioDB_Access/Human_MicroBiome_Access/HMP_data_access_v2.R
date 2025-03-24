# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install the main package and dependencies
BiocManager::install(c("curatedMetagenomicData", "phyloseq", "DESeq2", 
                      "ANCOMBC", "ggplot2", "dplyr", "tidyr", "vegan", 
                      "ComplexHeatmap", "ggrepel"))

# Load libraries
library(curatedMetagenomicData)
library(phyloseq)
library(DESeq2)
library(ANCOMBC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ComplexHeatmap)
library(ggrepel)
library(patchwork)  # For combining plots

#-----------------------------------------------
# 1. Data Retrieval and Exploration
#-----------------------------------------------

# List available HMP datasets
hmp_datasets <- availableDatasets()[grep("HMP", availableDatasets())]
print(hmp_datasets)

# Get the HMP stool and oral microbiome data
hmp_stool <- curatedMetagenomicData("HMP_2012.metaphlan_bugs_list.stool", 
                                    dryrun = FALSE)
hmp_oral <- curatedMetagenomicData("HMP_2012.metaphlan_bugs_list.oralcavity", 
                                  dryrun = FALSE)

# Extract the actual SummarizedExperiment objects
stool_se <- hmp_stool[[1]]
oral_se <- hmp_oral[[1]]

# Explore data structure
dim(stool_se)  # Number of taxa and samples
dim(oral_se)   # Number of taxa and samples

# Look at sample metadata
stool_meta <- colData(stool_se) %>% as.data.frame()
head(stool_meta)

# Combine the datasets for comparison
combined_se <- mergeData(hmp_stool, hmp_oral, sampleNames = TRUE)
combined_se <- combined_se[[1]]

# Add body site to sample metadata
colData(combined_se)$body_site <- ifelse(
  colData(combined_se)$body_subsite == "stool", 
  "Gut", 
  "Oral"
)

#-----------------------------------------------
# 2. Data Preprocessing
#-----------------------------------------------

# Filter low abundance taxa (keep taxa present in at least 10% of samples with rel. abundance > 0.01%)
keep_taxa <- rowSums(assay(combined_se) > 0.0001) > (ncol(combined_se) * 0.1)
filtered_se <- combined_se[keep_taxa,]

# Log transform abundances
log_transform <- function(x) {
  log_x <- log10(x + 1e-6)  # Add small pseudo-count to handle zeros
  return(log_x)
}

assay(filtered_se, "log10") <- log_transform(assay(filtered_se))

# Convert to phyloseq object for additional analyses
# First create the necessary components
otu_table <- otu_table(assay(filtered_se), taxa_are_rows = TRUE)
tax_table <- tax_table(as.matrix(rowData(filtered_se)))
sample_data <- sample_data(colData(filtered_se))

# Create phyloseq object
ps <- phyloseq(otu_table, tax_table, sample_data)

#-----------------------------------------------
# 3. Alpha Diversity Analysis
#-----------------------------------------------

# Calculate alpha diversity metrics
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$sample_id <- rownames(alpha_div)

# Add metadata
alpha_div <- left_join(
  alpha_div,
  as.data.frame(sample_data(ps)) %>% 
    select(sample_id, body_site, body_subsite, gender, age) %>%
    mutate(sample_id = rownames(sample_data(ps))),
  by = "sample_id"
)

# Plot alpha diversity by body site
alpha_plot <- ggplot(alpha_div, aes(x = body_site, y = Shannon, fill = body_site)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Shannon Diversity by Body Site",
       x = "Body Site",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none")

# Statistical test
alpha_test <- wilcox.test(Shannon ~ body_site, data = alpha_div)
alpha_pvalue <- paste("Wilcoxon p-value:", round(alpha_test$p.value, 6))

#-----------------------------------------------
# 4. Beta Diversity Analysis
#-----------------------------------------------

# Calculate Bray-Curtis distances
bray_dist <- phyloseq::distance(ps, method = "bray")

# Perform PERMANOVA
adonis_result <- adonis2(bray_dist ~ body_site, data = as.data.frame(sample_data(ps)))

# Perform PCoA
pcoa <- ordinate(ps, method = "PCoA", distance = bray_dist)

# Extract the percent variance explained by the first two axes
pcoa_var_explained <- pcoa$values$Eigenvalues / sum(pcoa$values$Eigenvalues)
pcoa_var_explained <- pcoa_var_explained[1:2] * 100

# Plot PCoA
beta_plot <- plot_ordination(ps, pcoa, color = "body_site") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "PCoA of Bray-Curtis Distances",
       x = paste0("PCo1 (", round(pcoa_var_explained[1], 1), "%)"),
       y = paste0("PCo2 (", round(pcoa_var_explained[2], 1), "%)"),
       color = "Body Site") +
  annotate("text", x = max(pcoa$vectors[,1]), y = min(pcoa$vectors[,2]),
           label = paste("PERMANOVA R² =", round(adonis_result$R2[1], 3), 
                         "\np-value =", format(adonis_result$`Pr(>F)`[1], scientific = TRUE, digits = 3)),
           hjust = 1, vjust = 0)

#-----------------------------------------------
# 5. Differential Abundance Analysis with ANCOM-BC
#-----------------------------------------------

# Prepare data for ANCOM-BC
ancom_data <- phyloseq_to_ancombc(ps, taxa_rank = "Genus")

# Run ANCOM-BC
ancom_bc_out <- ancombc(
  phyloseq = ancom_data,
  formula = "body_site",
  p_adj_method = "BH",
  zero_cut = 0.90,  # Filter genera absent in more than 90% of samples
  lib_cut = 1000,   # Minimum library size
  group = "body_site",
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.05
)

# Extract results
ancom_res <- ancom_bc_out$res
sig_taxa <- which(ancom_res$diff_abn$body_siteOral)
sig_genera <- rownames(ancom_res$diff_abn)[sig_taxa]

# Create a dataframe for plotting
ancom_df <- data.frame(
  taxon = rownames(ancom_res$diff_abn),
  log_fc = ancom_res$beta$body_siteOral,
  se = ancom_res$se$body_siteOral,
  p_val = ancom_res$p_val$body_siteOral,
  q_val = ancom_res$q_val$body_siteOral,
  diff_abundant = ancom_res$diff_abn$body_siteOral
)

# Add taxonomic information
taxa_info <- data.frame(tax_table(ancom_data))
ancom_df <- cbind(ancom_df, taxa_info[match(ancom_df$taxon, rownames(taxa_info)), c("Family", "Genus")])

# Volcano plot
volcano_plot <- ggplot(ancom_df, aes(x = log_fc, y = -log10(q_val), color = diff_abundant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey50", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = subset(ancom_df, diff_abundant & abs(log_fc) > 1),
    aes(label = Genus),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(title = "Differential Abundance: Oral vs Gut",
       x = "Log Fold Change (Oral vs Gut)",
       y = "-log10(q-value)",
       color = "Significant")

#-----------------------------------------------
# 6. Relative Abundance Analysis
#-----------------------------------------------

# Get top 20 genera by mean abundance across all samples
genus_abund <- psmelt(tax_glom(ps, taxrank = "Genus"))
genus_abund_summary <- genus_abund %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 20)

# Extract abundance data for these top genera
top_genera <- genus_abund_summary$Genus
top_genera_abund <- genus_abund %>%
  filter(Genus %in% top_genera) %>%
  mutate(Genus = factor(Genus, levels = top_genera))

# Calculate mean abundance by body site
top_genera_mean <- top_genera_abund %>%
  group_by(Genus, body_site) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = body_site, values_from = mean_abund, values_fill = 0) %>%
  mutate(log2_ratio = log2((Oral + 1e-5) / (Gut + 1e-5))) %>%
  arrange(desc(log2_ratio))

# Prepare for plotting
top_genera_long <- top_genera_mean %>%
  pivot_longer(cols = c(Gut, Oral), names_to = "body_site", values_to = "mean_abund") %>%
  mutate(Genus = factor(Genus, levels = top_genera_mean$Genus[order(top_genera_mean$log2_ratio)]))

# Create plot
abund_plot <- ggplot(top_genera_long, aes(x = Genus, y = mean_abund, fill = body_site)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Mean Abundance of Top 20 Genera",
       x = "",
       y = "Mean Relative Abundance",
       fill = "Body Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-----------------------------------------------
# 7. Heatmap of Top Taxa
#-----------------------------------------------

# Prepare data for heatmap (log-transformed relative abundance)
genus_mat <- as.matrix(otu_table(tax_glom(ps, taxrank = "Genus")))
genus_mat <- genus_mat[top_genera, ]
genus_mat <- log10(genus_mat + 1e-5)

# Get metadata for annotation
heatmap_meta <- data.frame(sample_data(ps))
heatmap_meta <- heatmap_meta[colnames(genus_mat), ]

# Create annotation for samples
sample_anno <- HeatmapAnnotation(
  body_site = heatmap_meta$body_site,
  body_subsite = heatmap_meta$body_subsite,
  col = list(
    body_site = c("Gut" = "#66C2A5", "Oral" = "#FC8D62"),
    body_subsite = c("stool" = "#66C2A5", 
                     "tongue_dorsum" = "#FC8D62", 
                     "supragingival_plaque" = "#8DA0CB", 
                     "buccal_mucosa" = "#E78AC3")
  )
)

# Create heatmap
heatmap <- Heatmap(
  genus_mat,
  name = "Log10 Abundance",
  col = colorRamp2(c(min(genus_mat), 0, max(genus_mat)), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  top_annotation = sample_anno,
  column_split = heatmap_meta$body_site,
  row_names_gp = gpar(fontsize = 8)
)

# Save heatmap
pdf("hmp_genus_heatmap.pdf", width = 12, height = 8)
draw(heatmap)
dev.off()

#-----------------------------------------------
# 8. Taxonomic Composition Barplots
#-----------------------------------------------

# Aggregate at phylum level
phylum_data <- tax_glom(ps, taxrank = "Phylum")
phylum_df <- psmelt(phylum_data)

# Calculate relative abundance by sample
phylum_df <- phylum_df %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance) * 100)

# Get top 9 phyla, group others
top_phyla <- phylum_df %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(rel_abund)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 9) %>%
  pull(Phylum)

phylum_df <- phylum_df %>%
  mutate(Phylum_grouped = ifelse(Phylum %in% top_phyla, Phylum, "Other"))

# Create stacked barplot
taxa_barplot <- ggplot(phylum_df, aes(x = Sample, y = rel_abund, fill = Phylum_grouped)) +
  geom_bar(stat = "identity") +
  facet_grid(~body_site, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Taxonomic Composition at Phylum Level",
       x = "Samples",
       y = "Relative Abundance (%)",
       fill = "Phylum")

#-----------------------------------------------
# 9. Combine Plots and Save Results
#-----------------------------------------------

# Combine alpha and beta diversity plots
diversity_plots <- alpha_plot + beta_plot +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Diversity Analysis of HMP Gut and Oral Microbiomes",
    caption = alpha_pvalue
  )

# Save main plots
ggsave("diversity_plots.pdf", diversity_plots, width = 12, height = 6)
ggsave("volcano_plot.pdf", volcano_plot, width = 10, height = 8)
ggsave("abundance_plot.pdf", abund_plot, width = 12, height = 6)
ggsave("taxa_barplot.pdf", taxa_barplot, width = 12, height = 6)

# Create a summary report
sink("hmp_analysis_summary.txt")
cat("# HMP Microbiome Analysis Summary\n\n")

cat("## Dataset Overview\n")
cat("- Stool samples:", nrow(stool_meta), "\n")
cat("- Oral samples:", nrow(colData(oral_se)), "\n")
cat("- Total taxa:", nrow(combined_se), "\n")
cat("- Taxa after filtering:", nrow(filtered_se), "\n\n")

cat("## Alpha Diversity\n")
cat(paste("- Wilcoxon test p-value:", round(alpha_test$p.value, 6), "\n"))
cat("- Mean Shannon Diversity in Gut:", mean(alpha_div$Shannon[alpha_div$body_site == "Gut"]), "\n")
cat("- Mean Shannon Diversity in Oral:", mean(alpha_div$Shannon[alpha_div$body_site == "Oral"]), "\n\n")

cat("## Beta Diversity\n")
cat("- PERMANOVA results:\n")
print(adonis_result)
cat("\n")

cat("## Differential Abundance\n")
cat("- Number of differentially abundant genera:", sum(ancom_df$diff_abundant), "\n")
cat("- Top genera enriched in Oral sites:\n")
oral_enriched <- ancom_df %>% 
  filter(diff_abundant & log_fc > 0) %>% 
  arrange(desc(log_fc)) %>% 
  head(10)
print(oral_enriched[, c("Genus", "log_fc", "q_val")])
cat("\n")

cat("- Top genera enriched in Gut:\n")
gut_enriched <- ancom_df %>% 
  filter(diff_abundant & log_fc < 0) %>% 
  arrange(log_fc) %>% 
  head(10)
print(gut_enriched[, c("Genus", "log_fc", "q_val")])
cat("\n")

sink()

# Return a named list of results
results <- list(
  phyloseq_object = ps,
  alpha_diversity = alpha_div,
  beta_diversity = pcoa,
  differential_abundance = ancom_df
)

# Print completion message
cat("Analysis complete. Results saved to current directory.\n")




