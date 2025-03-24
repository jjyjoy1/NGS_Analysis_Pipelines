#!/usr/bin/env Rscript
# TCGA Data Download Script using TCGAbiolinks
# This script provides functions to download TCGA data using TCGAbiolinks package
# Install required packages if not already installed

# Check and install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "readr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load required libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
})

# Function to list available TCGA projects
list_tcga_projects <- function() {
  projects <- TCGAbiolinks::getGDCprojects()
  tcga_projects <- projects %>% 
    filter(grepl("^TCGA", project_id)) %>%
    select(project_id, name, tumor, tissue)
  
  return(tcga_projects)
}

# Function to search for files based on criteria
search_files <- function(project_id = NULL, 
                         data_category = NULL, 
                         data_type = NULL, 
                         workflow_type = NULL,
                         experimental_strategy = NULL,
                         limit = 100) {
  
  # Build query parameters
  query_params <- list()
  
  # Add project filter if provided
  if (!is.null(project_id)) {
    query_params$project <- project_id
  } else {
    # Default to TCGA program if no specific project
    query_params$program <- "TCGA"
  }
  
  # Add other filters if provided
  if (!is.null(data_category)) {
    query_params$data.category <- data_category
  }
  
  if (!is.null(data_type)) {
    query_params$data.type <- data_type
  }
  
  if (!is.null(workflow_type)) {
    query_params$workflow.type <- workflow_type
  }
  
  if (!is.null(experimental_strategy)) {
    query_params$experimental.strategy <- experimental_strategy
  }
  
  # Execute the query
  tryCatch({
    query <- do.call(TCGAbiolinks::GDCquery, query_params)
    result <- getResults(query)
    
    # Limit results if specified
    if (nrow(result) > limit) {
      result <- result[1:limit, ]
    }
    
    return(result)
  }, error = function(e) {
    message("Error in search: ", e$message)
    return(NULL)
  })
}

# Function to download gene expression data
download_expression_data <- function(project_id, 
                                    workflow_type = "HTSeq - Counts",
                                    data_type = "Gene Expression Quantification",
                                    save_path = "downloads") {
  
  # Create directory if it doesn't exist
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  # Set up query for RNA-Seq data
  message(paste0("Querying gene expression data for project ", project_id, "..."))
  query <- TCGAbiolinks::GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = data_type,
    workflow.type = workflow_type,
    legacy = FALSE
  )
  
  # Download the data
  message("Downloading data...")
  TCGAbiolinks::GDCdownload(
    query = query,
    method = "api",
    files.per.chunk = 10,
    directory = save_path
  )
  
  # Prepare the data
  message("Preparing gene expression data...")
  expr_data <- TCGAbiolinks::GDCprepare(
    query = query,
    directory = save_path
  )
  
  # Save as TSV or RDS
  output_file_rds <- file.path(save_path, paste0(project_id, "_", gsub(" - ", "_", workflow_type), ".rds"))
  saveRDS(expr_data, file = output_file_rds)
  message(paste0("Data saved as RDS: ", output_file_rds))
  
  # Extract counts matrix and sample metadata
  if (is(expr_data, "SummarizedExperiment")) {
    # Extract counts matrix
    counts_matrix <- assay(expr_data)
    output_file_counts <- file.path(save_path, paste0(project_id, "_", gsub(" - ", "_", workflow_type), "_counts.tsv"))
    write.table(counts_matrix, file = output_file_counts, sep = "\t", quote = FALSE, col.names = NA)
    message(paste0("Counts matrix saved as TSV: ", output_file_counts))
    
    # Extract sample metadata
    sample_info <- colData(expr_data) %>% as.data.frame()
    output_file_samples <- file.path(save_path, paste0(project_id, "_", gsub(" - ", "_", workflow_type), "_samples.tsv"))
    write.table(sample_info, file = output_file_samples, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste0("Sample metadata saved as TSV: ", output_file_samples))
  }
  
  return(expr_data)
}

# Function to download clinical data
download_clinical_data <- function(project_id, save_path = "downloads") {
  # Create directory if it doesn't exist
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  message(paste0("Querying clinical data for project ", project_id, "..."))
  
  # Query clinical data
  clinical_query <- TCGAbiolinks::GDCquery(
    project = project_id,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    legacy = FALSE
  )
  
  # Download the data
  message("Downloading clinical data...")
  TCGAbiolinks::GDCdownload(
    query = clinical_query,
    method = "api",
    directory = save_path
  )
  
  # Get clinical data in a processed format
  message("Processing clinical data...")
  clinical <- TCGAbiolinks::GDCquery_clinic(project = project_id, type = "clinical")
  
  # Save the clinical data
  output_file <- file.path(save_path, paste0(project_id, "_clinical.tsv"))
  write.table(clinical, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message(paste0("Clinical data saved as TSV: ", output_file))
  
  return(clinical)
}

# Function to download mutation data
download_mutation_data <- function(project_id, 
                                  workflow_type = "MuSE Variant Aggregation and Masking",
                                  save_path = "downloads") {
  
  # Create directory if it doesn't exist
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  message(paste0("Querying mutation data for project ", project_id, "..."))
  
  # Query mutation data
  mutation_query <- TCGAbiolinks::GDCquery(
    project = project_id,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    workflow.type = workflow_type,
    legacy = FALSE
  )
  
  # Download the data
  message("Downloading mutation data...")
  TCGAbiolinks::GDCdownload(
    query = mutation_query,
    method = "api",
    directory = save_path
  )
  
  # Prepare the mutation data
  message("Preparing mutation data...")
  mutations <- TCGAbiolinks::GDCprepare(
    query = mutation_query,
    directory = save_path
  )
  
  # Save mutation data
  output_file <- file.path(save_path, paste0(project_id, "_mutations.tsv"))
  write.table(mutations, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message(paste0("Mutation data saved as TSV: ", output_file))
  
  return(mutations)
}

# Function to download methylation data
download_methylation_data <- function(project_id, save_path = "downloads") {
  # Create directory if it doesn't exist
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  message(paste0("Querying DNA methylation data for project ", project_id, "..."))
  
  # Query methylation data
  methyl_query <- TCGAbiolinks::GDCquery(
    project = project_id,
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    legacy = FALSE
  )
  
  # Download the data
  message("Downloading methylation data...")
  TCGAbiolinks::GDCdownload(
    query = methyl_query,
    method = "api",
    directory = save_path
  )
  
  # Prepare the methylation data
  message("Preparing methylation data...")
  methylation <- TCGAbiolinks::GDCprepare(
    query = methyl_query,
    directory = save_path
  )
  
  # Save methylation data if it's a SummarizedExperiment
  if (is(methylation, "SummarizedExperiment")) {
    beta_values <- assay(methylation)
    output_file <- file.path(save_path, paste0(project_id, "_methylation_beta.tsv"))
    write.table(beta_values, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
    message(paste0("Methylation beta values saved as TSV: ", output_file))
    
    # Save annotation
    probe_info <- rowData(methylation) %>% as.data.frame()
    output_file_probes <- file.path(save_path, paste0(project_id, "_methylation_probes.tsv"))
    write.table(probe_info, file = output_file_probes, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste0("Probe annotation saved as TSV: ", output_file_probes))
    
    # Save sample info
    sample_info <- colData(methylation) %>% as.data.frame()
    output_file_samples <- file.path(save_path, paste0(project_id, "_methylation_samples.tsv"))
    write.table(sample_info, file = output_file_samples, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste0("Sample information saved as TSV: ", output_file_samples))
  } else {
    # Just save whatever we got
    output_file <- file.path(save_path, paste0(project_id, "_methylation.rds"))
    saveRDS(methylation, file = output_file)
    message(paste0("Methylation data saved as RDS: ", output_file))
  }
  
  return(methylation)
}

# Main function to handle command line arguments
main <- function() {
  # Parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
    cat("TCGA Data Download Script\n")
    cat("Usage:\n")
    cat("  Rscript tcga_download.R list                                  - List all TCGA projects\n")
    cat("  Rscript tcga_download.R expression TCGA-BRCA [output_dir]     - Download gene expression data\n")
    cat("  Rscript tcga_download.R clinical TCGA-BRCA [output_dir]       - Download clinical data\n")
    cat("  Rscript tcga_download.R mutation TCGA-BRCA [output_dir]       - Download mutation data\n")
    cat("  Rscript tcga_download.R methylation TCGA-BRCA [output_dir]    - Download methylation data\n")
    cat("  Rscript tcga_download.R all TCGA-BRCA [output_dir]            - Download all data types\n")
    return(invisible(NULL))
  }
  
  # Process commands
  command <- args[1]
  
  if (command == "list") {
    message("Listing TCGA projects...")
    projects <- list_tcga_projects()
    print(projects)
    
    # Save projects to CSV
    write.csv(projects, file = "tcga_projects.csv", row.names = FALSE)
    message("Projects saved to tcga_projects.csv")
    
  } else {
    # For data download commands, we need at least a project ID
    if (length(args) < 2) {
      stop("Project ID is required for download commands")
    }
    
    project_id <- args[2]
    
    # Set output directory (default or user-specified)
    output_dir <- if (length(args) >= 3) args[3] else "downloads"
    
    # Create a project-specific subdirectory
    project_dir <- file.path(output_dir, project_id)
    dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Execute the requested command
    if (command == "expression") {
      message("Downloading gene expression data for ", project_id)
      download_expression_data(project_id, save_path = project_dir)
      
    } else if (command == "clinical") {
      message("Downloading clinical data for ", project_id)
      download_clinical_data(project_id, save_path = project_dir)
      
    } else if (command == "mutation") {
      message("Downloading mutation data for ", project_id)
      download_mutation_data(project_id, save_path = project_dir)
      
    } else if (command == "methylation") {
      message("Downloading methylation data for ", project_id)
      download_methylation_data(project_id, save_path = project_dir)
      
    } else if (command == "all") {
      message("Downloading all data types for ", project_id)
      
      # Download each data type
      download_clinical_data(project_id, save_path = project_dir)
      download_expression_data(project_id, save_path = project_dir)
      download_mutation_data(project_id, save_path = project_dir)
      download_methylation_data(project_id, save_path = project_dir)
      
      message("All data downloaded for ", project_id)
    } else {
      stop("Unknown command: ", command)
    }
  }
}

# Run the main function if script is executed
if (!interactive()) {
  main()
}

