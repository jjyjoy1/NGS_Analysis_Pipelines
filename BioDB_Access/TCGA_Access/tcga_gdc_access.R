#!/usr/bin/env Rscript
# TCGA Data Download using GenomicDataCommons (GDC) R package
# This offers an alternative approach to TCGAbiolinks

# Check and install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("GenomicDataCommons", "dplyr", "readr", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load libraries
suppressPackageStartupMessages({
  library(GenomicDataCommons)
  library(dplyr)
  library(readr)
  library(data.table)
})

# Check GDC status
gdc_status <- function() {
  status <- GenomicDataCommons::status()
  message("GDC Status:")
  message("- API Status: ", ifelse(status$status == "OK", "Available", "Unavailable"))
  message("- API Version: ", status$`tag:github.com,2008:Repository/35925266@`$version)
  message("- Data Release: ", status$data_release)
  message("- Projects: ", status$`tag:github.com,2008:Repository/35925266@`$projects$count)
  message("- Cases: ", status$`tag:github.com,2008:Repository/35925266@`$cases$count)
  message("- Files: ", status$`tag:github.com,2008:Repository/35925266@`$files$count)
  message("- Annotations: ", status$`tag:github.com,2008:Repository/35925266@`$annotations$count)
  
  return(status)
}

# List all TCGA projects
list_tcga_projects <- function() {
  message("Fetching TCGA projects from GDC...")
  
  projects <- projects() %>%
    GenomicDataCommons::filter(program.name == "TCGA") %>%
    results_all()
  
  # Extract relevant information
  project_data <- data.frame(
    project_id = projects$project_id,
    name = projects$name,
    disease_type = projects$disease_type,
    primary_site = projects$primary_site,
    case_count = sapply(projects$summary, function(x) x$case_count),
    stringsAsFactors = FALSE
  )
  
  return(project_data)
}

# Get available data categories for a project
get_data_categories <- function(project_id) {
  message("Fetching available data categories for ", project_id, "...")
  
  files <- files() %>%
    GenomicDataCommons::filter(cases.project.project_id == project_id) %>%
    facet("data_category") %>%
    aggregations()
  
  categories <- data.frame(
    data_category = files$data_category$key,
    file_count = files$data_category$doc_count,
    stringsAsFactors = FALSE
  )
  
  return(categories)
}

# Get data types for a specific data category in a project
get_data_types <- function(project_id, data_category) {
  message("Fetching data types for ", project_id, " in category ", data_category, "...")
  
  files <- files() %>%
    GenomicDataCommons::filter(cases.project.project_id == project_id) %>%
    GenomicDataCommons::filter(data_category == data_category) %>%
    facet("data_type") %>%
    aggregations()
  
  types <- data.frame(
    data_type = files$data_type$key,
    file_count = files$data_type$doc_count,
    stringsAsFactors = FALSE
  )
  
  return(types)
}

# Get available workflow types for a project, data category, and data type
get_workflow_types <- function(project_id, data_category, data_type) {
  message("Fetching workflow types for ", project_id, " in category ", data_category, 
          ", type ", data_type, "...")
  
  files <- files() %>%
    GenomicDataCommons::filter(cases.project.project_id == project_id) %>%
    GenomicDataCommons::filter(data_category == data_category) %>%
    GenomicDataCommons::filter(data_type == data_type) %>%
    facet("analysis.workflow_type") %>%
    aggregations()
  
  if (!is.null(files$`analysis.workflow_type`)) {
    types <- data.frame(
      workflow_type = files$`analysis.workflow_type`$key,
      file_count = files$`analysis.workflow_type`$doc_count,
      stringsAsFactors = FALSE
    )
  } else {
    types <- data.frame(
      workflow_type = character(0),
      file_count = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(types)
}

# Search for files based on criteria
search_files <- function(project_id = NULL, 
                         data_category = NULL, 
                         data_type = NULL, 
                         workflow_type = NULL,
                         limit = 1000) {
  
  message("Searching for files with specified criteria...")
  
  # Start building the query
  query <- files()
  
  # Add filters based on provided parameters
  if (!is.null(project_id)) {
    query <- query %>% GenomicDataCommons::filter(cases.project.project_id == project_id)
  }
  
  if (!is.null(data_category)) {
    query <- query %>% GenomicDataCommons::filter(data_category == data_category)
  }
  
  if (!is.null(data_type)) {
    query <- query %>% GenomicDataCommons::filter(data_type == data_type)
  }
  
  if (!is.null(workflow_type)) {
    query <- query %>% GenomicDataCommons::filter(`analysis.workflow_type` == workflow_type)
  }
  
  # Execute the query with size limit
  file_results <- query %>% 
    GenomicDataCommons::select(c("file_id", "file_name", "cases.case_id", "cases.submitter_id", 
             "cases.project.project_id", "data_category", "data_type", 
             "analysis.workflow_type", "file_size", "access")) %>%
    GenomicDataCommons::results(size = limit)
  
  return(file_results)
}

# Download files based on file IDs
download_files <- function(file_ids, save_path = "downloads") {
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  message("Downloading ", length(file_ids), " files to ", save_path, "...")
  
  # Download files
  downloaded_files <- GenomicDataCommons::gdcdata(file_ids, directory = save_path, progress = TRUE)
  
  message("Download complete. Files saved to ", save_path)
  return(downloaded_files)
}

# Function to download gene expression data
download_expression_data <- function(project_id, 
                                    data_type = "Gene Expression Quantification",
                                    workflow_type = "HTSeq - Counts",
                                    save_path = "downloads",
                                    max_files = 1000) {
  
  # Create directory
  project_dir <- file.path(save_path, project_id, "Expression")
  dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Searching for gene expression data (", workflow_type, ") for ", project_id, "...")
  
  # Find expression files
  expr_files <- search_files(
    project_id = project_id,
    data_category = "Transcriptome Profiling",
    data_type = data_type,
    workflow_type = workflow_type,
    limit = max_files
  )
  
  if (length(expr_files) == 0) {
    message("No expression files found for the specified criteria.")
    return(NULL)
  }
  
  message("Found ", length(expr_files), " expression files.")
  
  # Download files
  file_ids <- expr_files$file_id
  downloaded_files <- download_files(file_ids, save_path = project_dir)
  
  # Create a manifest file
  manifest <- data.frame(
    file_id = expr_files$file_id,
    file_name = expr_files$file_name,
    case_id = sapply(expr_files$cases, function(x) ifelse(length(x$submitter_id) > 0, x$submitter_id[1], NA)),
    file_size = expr_files$file_size,
    local_path = downloaded_files,
    stringsAsFactors = FALSE
  )
  
  manifest_file <- file.path(project_dir, paste0(project_id, "_expression_manifest.tsv"))
  write_tsv(manifest, manifest_file)
  message("Manifest file saved to: ", manifest_file)
  
  return(manifest)
}

# Function to download clinical data
download_clinical_data <- function(project_id, save_path = "downloads") {
  # Create directory
  project_dir <- file.path(save_path, project_id, "Clinical")
  dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Retrieving clinical data for ", project_id, "...")
  
  # Get all cases for the project
  cases <- cases() %>%
    GenomicDataCommons::filter(project.project_id == project_id) %>%
    GenomicDataCommons::expand("demographic") %>%
    GenomicDataCommons::expand("diagnoses") %>%
    GenomicDataCommons::expand("exposures") %>%
    results_all()
  
  if (length(cases) == 0) {
    message("No cases found for project ", project_id)
    return(NULL)
  }
  
  message("Found ", length(cases), " cases. Processing clinical data...")
  
  # Extract demographics
  demographics <- do.call(rbind, lapply(1:length(cases), function(i) {
    case <- cases[i]
    demo <- case$demographic
    if (is.null(demo)) return(NULL)
    
    # Convert to data frame with case_id
    demo_df <- as.data.frame(demo)
    demo_df$case_id <- case$case_id
    demo_df$submitter_id <- case$submitter_id
    return(demo_df)
  }))
  
  # Extract diagnoses
  diagnoses <- do.call(rbind, lapply(1:length(cases), function(i) {
    case <- cases[i]
    diags <- case$diagnoses
    if (length(diags) == 0) return(NULL)
    
    # Process all diagnoses for this case
    diag_list <- lapply(diags, function(diag) {
      diag_df <- as.data.frame(diag)
      diag_df$case_id <- case$case_id
      diag_df$submitter_id <- case$submitter_id
      return(diag_df)
    })
    
    return(do.call(rbind, diag_list))
  }))
  
  # Extract exposures
  exposures <- do.call(rbind, lapply(1:length(cases), function(i) {
    case <- cases[i]
    exps <- case$exposures
    if (length(exps) == 0) return(NULL)
    
    # Process all exposures for this case
    exp_list <- lapply(exps, function(exp) {
      exp_df <- as.data.frame(exp)
      exp_df$case_id <- case$case_id
      exp_df$submitter_id <- case$submitter_id
      return(exp_df)
    })
    
    return(do.call(rbind, exp_list))
  }))
  
  # Save demographic data
  if (!is.null(demographics) && nrow(demographics) > 0) {
    demographics_file <- file.path(project_dir, paste0(project_id, "_demographics.tsv"))
    write_tsv(demographics, demographics_file)
    message("Demographics saved to: ", demographics_file)
  } else {
    message("No demographic data available.")
  }
  
  # Save diagnoses data
  if (!is.null(diagnoses) && nrow(diagnoses) > 0) {
    diagnoses_file <- file.path(project_dir, paste0(project_id, "_diagnoses.tsv"))
    write_tsv(diagnoses, diagnoses_file)
    message("Diagnoses saved to: ", diagnoses_file)
  } else {
    message("No diagnoses data available.")
  }
  
  # Save exposures data
  if (!is.null(exposures) && nrow(exposures) > 0) {
    exposures_file <- file.path(project_dir, paste0(project_id, "_exposures.tsv"))
    write_tsv(exposures, exposures_file)
    message("Exposures saved to: ", exposures_file)
  } else {
    message("No exposures data available.")
  }
  
  # Also download clinical supplement files
  message("Searching for clinical supplement files...")
  clinical_files <- search_files(
    project_id = project_id,
    data_category = "Clinical",
    data_type = "Clinical Supplement",
    limit = 1000
  )
  
  if (length(clinical_files) > 0) {
    message("Found ", length(clinical_files), " clinical supplement files.")
    file_ids <- clinical_files$file_id
    download_files(file_ids, save_path = project_dir)
  } else {
    message("No clinical supplement files found.")
  }
  
  return(list(
    demographics = demographics,
    diagnoses = diagnoses,
    exposures = exposures,
    clinical_files = clinical_files
  ))
}

# Function to download mutation data
download_mutation_data <- function(project_id, save_path = "downloads") {
  # Create directory
  project_dir <- file.path(save_path, project_id, "Mutations")
  dir.create(project_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("Searching for mutation data for ", project_id, "...")
  
  # Get all available workflow types for mutations
  workflow_types <- get_workflow_types(
    project_id = project_id,
    data_category = "Simple Nucleotide Variation",
    data_type = "Masked Somatic Mutation"
  )
  
  if (nrow(workflow_types) == 0) {
    message("No mutation data found for project ", project_id)
    return(NULL)
  }
  
  message("Found ", nrow(workflow_types), " mutation workflow types. Downloading files...")
  
  # Download files for each workflow type
  results <- list()
  
  for (i in 1:nrow(workflow_types)) {
    workflow <- workflow_types$workflow_type[i]
    message("Processing workflow: ", workflow)
    
    mutation_files <- search_files(
      project_id = project_id,
      data_category = "Simple Nucleotide Variation",
      data_type = "Masked Somatic Mutation",
      workflow_type = workflow,
      limit = 1000
    )
    
    if (length(mutation_files) > 0) {
      file_ids <- mutation_files$file_id
      workflow_dir <- file.path(project_dir, make.names(workflow))
      dir.create(workflow_dir, showWarnings = FALSE, recursive = TRUE)
      
      downloaded_files <- download_files(file_ids, save_path = workflow_dir)
      
      # Create a manifest file
      manifest <- data.frame(
        file_id = mutation_files$file_id,
        file_name = mutation_files$file_name,
        case_id = sapply(mutation_files$cases, function(x) ifelse(length(x$submitter_id) > 0, x$submitter_id[1], NA)),
        file_size = mutation_files$file_size,
        local_path = downloaded_files,
        stringsAsFactors = FALSE
      )
      
      manifest_file <- file.path(workflow_dir, paste0(project_id, "_", make.names(workflow), "_manifest.tsv"))
      write_tsv(manifest, manifest_file)
      message("Manifest file saved to: ", manifest_file)
      
      results[[workflow]] <- manifest
    }
  }
  
  return(results)
}

# Main function to handle command line arguments
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || args[1] == "-h" || args[1] == "--help") {
    cat("TCGA Data Download Script using GenomicDataCommons\n")
    cat("Usage:\n")
    cat("  Rscript gdc_download.R status                               - Check GDC API status\n")
    cat("  Rscript gdc_download.R list                                 - List all TCGA projects\n")
    cat("  Rscript gdc_download.R categories TCGA-BRCA                 - List data categories for a project\n")
    cat("  Rscript gdc_download.R types TCGA-BRCA 'Data Category'      - List data types for a category\n")
    cat("  Rscript gdc_download.R workflow TCGA-BRCA 'Data Category' 'Data Type' - List workflows\n")
    cat("  Rscript gdc_download.R search TCGA-BRCA 'Data Category' 'Data Type' ['Workflow'] - Search files\n")
    cat("  Rscript gdc_download.R expression TCGA-BRCA [output_dir]    - Download gene expression data\n")
    cat("  Rscript gdc_download.R clinical TCGA-BRCA [output_dir]      - Download clinical data\n")
    cat("  Rscript gdc_download.R mutation TCGA-BRCA [output_dir]      - Download mutation data\n")
    return(invisible(NULL))
  }
  
  # Process commands
  command <- args[1]
  
  if (command == "status") {
    status <- gdc_status()
    print(status)
    
  } else if (command == "list") {
    projects <- list_tcga_projects()
    print(projects)
    
    # Save projects to CSV
    output_file <- "tcga_projects.csv"
    write_csv(projects, output_file)
    message("Projects saved to ", output_file)
    
  } else if (command == "categories") {
    if (length(args) < 2) stop("Project ID is required")
    project_id <- args[2]
    
    categories <- get_data_categories(project_id)
    print(categories)
    
  } else if (command == "types") {
    if (length(args) < 3) stop("Project ID and data category are required")
    project_id <- args[2]
    data_category <- args[3]
    
    types <- get_data_types(project_id, data_category)
    print(types)
    
  } else if (command == "workflow") {
    if (length(args) < 4) stop("Project ID, data category, and data type are required")
    project_id <- args[2]
    data_category <- args[3]
    data_type <- args[4]
    
    workflows <- get_workflow_types(project_id, data_category, data_type)
    print(workflows)
    
  } else if (command == "search") {
    if (length(args) < 2) stop("Project ID is required")
    project_id <- args[2]
    data_category <- if (length(args) >= 3) args[3] else NULL
    data_type <- if (length(args) >= 4) args[4] else NULL
    workflow_type <- if (length(args) >= 5) args[5] else NULL
    
    files <- search_files(project_id, data_category, data_type, workflow_type)
    print(head(files))
    message("Found ", length(files), " files.")
    
    # Save search results
    output_file <- "search_results.csv"
    write_csv(as.data.frame(files), output_file)
    message("Search results saved to ", output_file)
    
  } else if (command == "expression") {
    if (length(args) < 2) stop("Project ID is required")
    project_id <- args[2]
    output_dir <- if (length(args) >= 3) args[3] else "downloads"
    
    download_expression_data(project_id, save_path = output_dir)
    
  } else if (command == "clinical") {
    if (length(args) < 2) stop("Project ID is required")
    project_id <- args[2]
    output_dir <- if (length(args) >= 3) args[3] else "downloads"
    
    download_clinical_data(project_id, save_path = output_dir)
    
  } else if (command == "mutation") {
    if (length(args) < 2) stop("Project ID is required")
    project_id <- args[2]
    output_dir <- if (length(args) >= 3) args[3] else "downloads"
    
    download_mutation_data(project_id, save_path = output_dir)
    
  } else {
    stop("Unknown command: ", command)
  }
}

# Run the main function if script is executed directly
if (!interactive()) {
  main()
}


