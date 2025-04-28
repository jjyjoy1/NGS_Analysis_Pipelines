# Snakefile for GWAS pipeline using GCTA

import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Define input and output directories
RAW_DATA = config["raw_data"]
OUTPUT_DIR = config["output_dir"]
PHENOTYPE_FILE = config["phenotype_file"]
COVARIATE_FILE = config["covariate_file"]
PHENOTYPE_NAME = config["phenotype_name"]
COVARIATE_NAMES = config["covariate_names"]
CHROMOSOMES = list(range(1, 23))

# Final target rule
rule all:
    input:
        qc_complete = os.path.join(OUTPUT_DIR, "qc_complete.bed"),
        h2_results = os.path.join(OUTPUT_DIR, "h2_results.hsq"),
        gwas_results = os.path.join(OUTPUT_DIR, "gwas_results.loco.mlma"),
        independent_signals = os.path.join(OUTPUT_DIR, "independent_signals.jma.cojo"),
        credible_sets = os.path.join(OUTPUT_DIR, "credible_sets.txt"),
        coloc_results = os.path.join(OUTPUT_DIR, "coloc_results.txt"),
        prs_results = os.path.join(OUTPUT_DIR, "prs_results.profile")

# 1. Sample Quality Control
rule sample_qc:
    input:
        bfile = [f"{RAW_DATA}.bed", f"{RAW_DATA}.bim", f"{RAW_DATA}.fam"]
    output:
        bfile = [os.path.join(OUTPUT_DIR, "sample_qc.bed"), 
                os.path.join(OUTPUT_DIR, "sample_qc.bim"), 
                os.path.join(OUTPUT_DIR, "sample_qc.fam")],
        sexcheck = os.path.join(OUTPUT_DIR, "sexcheck.sexcheck"),
        het = os.path.join(OUTPUT_DIR, "het.het"),
        het_outliers = os.path.join(OUTPUT_DIR, "het_outliers.txt"),
        ibd = os.path.join(OUTPUT_DIR, "ibd.genome"),
        related_to_remove = os.path.join(OUTPUT_DIR, "related_to_remove.txt"),
        pca = os.path.join(OUTPUT_DIR, "pca.eigenvec")
    params:
        out_prefix = os.path.join(OUTPUT_DIR, "sample_qc"),
        mind = config["qc"]["mind"],
        genome_min = config["qc"]["genome_min"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "sample_qc.log")
    shell:
        """
        # Check sample call rates
        plink --bfile {RAW_DATA} \
            --mind {params.mind} \
            --make-bed \
            --out {params.out_prefix}_mind 2>&1 | tee {log}
            
        # Sex check
        plink --bfile {params.out_prefix}_mind \
            --check-sex \
            --out {params.out_prefix}_sexcheck 2>&1 | tee -a {log}
            
        # Remove samples with sex discrepancies
        awk '{{ if ($5 == "PROBLEM") print $1,$2 }}' {params.out_prefix}_sexcheck.sexcheck > {output.sexcheck}
        
        plink --bfile {params.out_prefix}_mind \
            --remove {output.sexcheck} \
            --make-bed \
            --out {params.out_prefix}_sexchecked 2>&1 | tee -a {log}
        
        # Check for excessive heterozygosity
        plink --bfile {params.out_prefix}_sexchecked \
            --het \
            --out {params.out_prefix}_het 2>&1 | tee -a {log}
        
        # Identify samples with outlier heterozygosity
        Rscript -e '
            het <- read.table("{params.out_prefix}_het.het", header=TRUE);
            het$F <- (het$O.HOM - het$E.HOM)/het$E.HOM;
            mean_f <- mean(het$F);
            sd_f <- sd(het$F);
            outliers <- subset(het, F < mean_f-3*sd_f | F > mean_f+3*sd_f);
            write.table(outliers[,c(1,2)], "{output.het_outliers}", 
                        row.names=FALSE, col.names=FALSE, quote=FALSE);' 2>&1 | tee -a {log}
        
        # Remove heterozygosity outliers
        plink --bfile {params.out_prefix}_sexchecked \
            --remove {output.het_outliers} \
            --make-bed \
            --out {params.out_prefix}_het_filtered 2>&1 | tee -a {log}
        
        # Check for relatedness
        plink --bfile {params.out_prefix}_het_filtered \
            --indep-pairwise 50 5 0.2 \
            --out {params.out_prefix}_pruned 2>&1 | tee -a {log}
        
        plink --bfile {params.out_prefix}_het_filtered \
            --extract {params.out_prefix}_pruned.prune.in \
            --genome \
            --min {params.genome_min} \
            --out {params.out_prefix}_ibd 2>&1 | tee -a {log}
            
        cp {params.out_prefix}_ibd.genome {output.ibd}
        
        # Select one from each related cluster
        awk '{{ if ($9 > {params.genome_min}) print $1,$2,$3,$4,$9,$10 }}' {output.ibd} > {params.out_prefix}_related_pairs.txt
        
        Rscript -e '
            related <- read.table("{params.out_prefix}_related_pairs.txt", header=FALSE);
            samples_to_remove <- c();
            for(i in 1:nrow(related)) {{
                id1 <- paste(related[i,1], related[i,2], sep=":");
                id2 <- paste(related[i,3], related[i,4], sep=":");
                if(!(id1 %in% samples_to_remove) && !(id2 %in% samples_to_remove)) {{
                    samples_to_remove <- c(samples_to_remove, id2);
                }}
            }}
            ids <- do.call(rbind, strsplit(samples_to_remove, ":"));
            write.table(ids, "{output.related_to_remove}", 
                        row.names=FALSE, col.names=FALSE, quote=FALSE);' 2>&1 | tee -a {log}
        
        # Remove related individuals
        plink --bfile {params.out_prefix}_het_filtered \
            --remove {output.related_to_remove} \
            --make-bed \
            --out {params.out_prefix}_unrelated 2>&1 | tee -a {log}
        
        # Run PCA for population stratification
        plink --bfile {params.out_prefix}_unrelated \
            --pca 20 \
            --out {params.out_prefix}_pca 2>&1 | tee -a {log}
            
        cp {params.out_prefix}_pca.eigenvec {output.pca}
            
        # Create final sample QC files
        plink --bfile {params.out_prefix}_unrelated \
            --make-bed \
            --out {params.out_prefix} 2>&1 | tee -a {log}
        """

# 2. Variant Quality Control
rule variant_qc:
    input:
        bfile = rules.sample_qc.output.bfile
    output:
        bfile = [os.path.join(OUTPUT_DIR, "qc_complete.bed"), 
                os.path.join(OUTPUT_DIR, "qc_complete.bim"), 
                os.path.join(OUTPUT_DIR, "qc_complete.fam")],
        diffmiss = os.path.join(OUTPUT_DIR, "diffmiss.missing")
    params:
        in_prefix = os.path.join(OUTPUT_DIR, "sample_qc"),
        out_prefix = os.path.join(OUTPUT_DIR, "qc_complete"),
        geno = config["qc"]["geno"],
        maf = config["qc"]["maf"],
        hwe = config["qc"]["hwe"],
        diffmiss_pval = config["qc"]["diffmiss_pval"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "variant_qc.log")
    shell:
        """
        # Check variant call rates
        plink --bfile {params.in_prefix} \
            --geno {params.geno} \
            --make-bed \
            --out {params.out_prefix}_geno 2>&1 | tee {log}
        
        # Filter by Minor Allele Frequency
        plink --bfile {params.out_prefix}_geno \
            --maf {params.maf} \
            --make-bed \
            --out {params.out_prefix}_maf 2>&1 | tee -a {log}
        
        # Filter by Hardy-Weinberg Equilibrium
        plink --bfile {params.out_prefix}_maf \
            --hwe {params.hwe} \
            --make-bed \
            --out {params.out_prefix}_hwe 2>&1 | tee -a {log}
        
        # Check for SNPs with differential missingness
        plink --bfile {params.out_prefix}_hwe \
            --test-missing \
            --out {params.out_prefix}_diffmiss 2>&1 | tee -a {log}
            
        cp {params.out_prefix}_diffmiss.missing {output.diffmiss}
        
        # Remove SNPs with differential missingness
        awk '{{ if ($5 < {params.diffmiss_pval}) print $2 }}' {output.diffmiss} > {params.out_prefix}_diffmiss_snps.txt
        
        plink --bfile {params.out_prefix}_hwe \
            --exclude {params.out_prefix}_diffmiss_snps.txt \
            --make-bed \
            --out {params.out_prefix} 2>&1 | tee -a {log}
        """

# 3. Imputation Preparation
rule imputation_prep:
    input:
        bfile = rules.variant_qc.output.bfile
    output:
        vcfs = expand(os.path.join(OUTPUT_DIR, "imputation", "chr{chrom}.vcf.gz"), chrom=CHROMOSOMES)
    params:
        in_prefix = os.path.join(OUTPUT_DIR, "qc_complete"),
        out_dir = os.path.join(OUTPUT_DIR, "imputation")
    log:
        os.path.join(OUTPUT_DIR, "logs", "imputation_prep.log")
    shell:
        """
        mkdir -p {params.out_dir}
        
        for CHR in {{1..22}}; do
            plink --bfile {params.in_prefix} \
                --chr $CHR \
                --recode vcf-4.2 \
                --out {params.out_dir}/chr$CHR
                
            bgzip {params.out_dir}/chr$CHR.vcf
            tabix -p vcf {params.out_dir}/chr$CHR.vcf.gz
        done 2>&1 | tee {log}
        """

# 4. Process Imputed Data (assuming imputation has been done externally)
rule process_imputed:
    input:
        vcfs = rules.imputation_prep.output.vcfs
    output:
        bfile = [os.path.join(OUTPUT_DIR, "imputed_merged.bed"), 
                os.path.join(OUTPUT_DIR, "imputed_merged.bim"), 
                os.path.join(OUTPUT_DIR, "imputed_merged.fam")]
    params:
        in_dir = os.path.join(OUTPUT_DIR, "imputation"),
        out_prefix = os.path.join(OUTPUT_DIR, "imputed_merged"),
        tmp_dir = os.path.join(OUTPUT_DIR, "tmp")
    log:
        os.path.join(OUTPUT_DIR, "logs", "process_imputed.log")
    shell:
        """
        mkdir -p {params.tmp_dir}
        
        # Convert imputed VCFs to PLINK format
        for CHR in {{1..22}}; do
            plink --vcf {params.in_dir}/chr$CHR.vcf.gz \
                --make-bed \
                --out {params.tmp_dir}/imputed_chr$CHR
        done
        
        # Create merge list
        echo "{params.tmp_dir}/imputed_chr1" > {params.tmp_dir}/merge_list.txt
        for CHR in {{2..22}}; do
            echo "{params.tmp_dir}/imputed_chr$CHR" >> {params.tmp_dir}/merge_list.txt
        done
        
        # Merge all chromosomes
        plink --merge-list {params.tmp_dir}/merge_list.txt \
            --make-bed \
            --out {params.out_prefix} 2>&1 | tee {log}
        """

# 5. Create Genomic Relationship Matrix
rule make_grm:
    input:
        bfile = rules.process_imputed.output.bfile
    output:
        grm = [os.path.join(OUTPUT_DIR, "grm.grm.bin"), 
              os.path.join(OUTPUT_DIR, "grm.grm.id"), 
              os.path.join(OUTPUT_DIR, "grm.grm.N.bin")]
    params:
        in_prefix = os.path.join(OUTPUT_DIR, "imputed_merged"),
        out_prefix = os.path.join(OUTPUT_DIR, "grm"),
        threads = config["resources"]["gcta_threads"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "make_grm.log")
    shell:
        """
        # Prune SNPs for GRM
        plink --bfile {params.in_prefix} \
            --indep-pairwise 50 5 0.2 \
            --out {params.in_prefix}_prune 2>&1 | tee {log}
        
        # Create GRM
        gcta64 --bfile {params.in_prefix} \
            --extract {params.in_prefix}_prune.prune.in \
            --make-grm \
            --out {params.out_prefix} \
            --thread-num {params.threads} 2>&1 | tee -a {log}
        """

# 6. Estimate Heritability
rule estimate_heritability:
    input:
        grm = rules.make_grm.output.grm,
        phenotype = PHENOTYPE_FILE,
        covariate = COVARIATE_FILE
    output:
        hsq = os.path.join(OUTPUT_DIR, "h2_results.hsq")
    params:
        grm_prefix = os.path.join(OUTPUT_DIR, "grm"),
        out_prefix = os.path.join(OUTPUT_DIR, "h2_results"),
        phenotype_name = PHENOTYPE_NAME,
        covariate_names = COVARIATE_NAMES,
        threads = config["resources"]["gcta_threads"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "estimate_heritability.log")
    shell:
        """
        gcta64 --grm {params.grm_prefix} \
            --pheno {input.phenotype} \
            --mpheno {params.phenotype_name} \
            --covar {input.covariate} \
            --qcovar {input.covariate} \
            --reml \
            --out {params.out_prefix} \
            --thread-num {params.threads} 2>&1 | tee {log}
        """

# 7. Run GWAS with GCTA-MLMA-LOCO
rule gwas_mlma_loco:
    input:
        bfile = rules.process_imputed.output.bfile,
        grm = rules.make_grm.output.grm,
        phenotype = PHENOTYPE_FILE,
        covariate = COVARIATE_FILE
    output:
        mlma = os.path.join(OUTPUT_DIR, "gwas_results.loco.mlma")
    params:
        bfile_prefix = os.path.join(OUTPUT_DIR, "imputed_merged"),
        grm_prefix = os.path.join(OUTPUT_DIR, "grm"),
        out_prefix = os.path.join(OUTPUT_DIR, "gwas_results"),
        phenotype_name = PHENOTYPE_NAME,
        covariate_names = COVARIATE_NAMES,
        threads = config["resources"]["gcta_threads"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "gwas_mlma_loco.log")
    shell:
        """
        gcta64 --mlma-loco \
            --bfile {params.bfile_prefix} \
            --grm {params.grm_prefix} \
            --pheno {input.phenotype} \
            --mpheno {params.phenotype_name} \
            --covar {input.covariate} \
            --qcovar {input.covariate} \
            --out {params.out_prefix} \
            --thread-num {params.threads} 2>&1 | tee {log}
        """

# 8. Identify Independent Signals with COJO
rule cojo_analysis:
    input:
        bfile = rules.process_imputed.output.bfile,
        gwas = rules.gwas_mlma_loco.output.mlma
    output:
        cojo = os.path.join(OUTPUT_DIR, "independent_signals.jma.cojo")
    params:
        bfile_prefix = os.path.join(OUTPUT_DIR, "imputed_merged"),
        out_prefix = os.path.join(OUTPUT_DIR, "independent_signals"),
        p_threshold = config["analysis"]["cojo_p_threshold"],
        threads = config["resources"]["gcta_threads"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "cojo_analysis.log")
    shell:
        """
        gcta64 --bfile {params.bfile_prefix} \
            --cojo-file {input.gwas} \
            --cojo-slct \
            --cojo-p {params.p_threshold} \
            --out {params.out_prefix} \
            --thread-num {params.threads} 2>&1 | tee {log}
        """

# 9. Create Credible Sets for Each Locus
rule credible_sets:
    input:
        bfile = rules.process_imputed.output.bfile,
        gwas = rules.gwas_mlma_loco.output.mlma,
        indep = rules.cojo_analysis.output.cojo
    output:
        credible_sets = os.path.join(OUTPUT_DIR, "credible_sets.txt")
    params:
        bfile_prefix = os.path.join(OUTPUT_DIR, "imputed_merged"),
        out_prefix = os.path.join(OUTPUT_DIR, "credible_sets"),
        window_size = config["analysis"]["credible_set_window"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "credible_sets.log")
    script:
        "scripts/credible_sets.R"

# 10. Colocalization Analysis with eQTL Data
rule colocalization:
    input:
        gwas = rules.gwas_mlma_loco.output.mlma,
        indep = rules.cojo_analysis.output.cojo,
        eqtl = config["eqtl_file"]
    output:
        coloc = os.path.join(OUTPUT_DIR, "coloc_results.txt")
    params:
        out_prefix = os.path.join(OUTPUT_DIR, "coloc"),
        window_size = config["analysis"]["coloc_window"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "colocalization.log")
    script:
        "scripts/colocalization.R"

# 11. Generate Polygenic Risk Score
rule create_prs:
    input:
        bfile = rules.process_imputed.output.bfile,
        gwas = rules.gwas_mlma_loco.output.mlma
    output:
        profile = os.path.join(OUTPUT_DIR, "prs_results.profile")
    params:
        bfile_prefix = os.path.join(OUTPUT_DIR, "imputed_merged"),
        out_prefix = os.path.join(OUTPUT_DIR, "prs_results"),
        p_threshold = config["analysis"]["prs_p_threshold"]
    log:
        os.path.join(OUTPUT_DIR, "logs", "create_prs.log")
    shell:
        """
        # Extract effect sizes from GWAS
        awk '{{if ($13 < {params.p_threshold}) print $2,$8}}' {input.gwas} > {params.out_prefix}_effects.txt
        
        # Create PRS
        plink --bfile {params.bfile_prefix} \
            --score {params.out_prefix}_effects.txt 1 2 header \
            --out {params.out_prefix} 2>&1 | tee {log}
        """
