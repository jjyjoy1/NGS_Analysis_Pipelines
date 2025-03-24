#!/usr/bin/env python3
"""
TCGA Data Download Script
This script provides functions to download TCGA data directly using the GDC API.
Requirements:
- requests
- pandas
- tqdm

Install with: pip install requests pandas tqdm
"""

import os
import json
import requests
import pandas as pd
from tqdm import tqdm
import argparse
import time
import tarfile
import shutil
import gzip

# GDC API endpoints
GDC_API_BASE = "https://api.gdc.cancer.gov/"
GDC_DATA_ENDPOINT = GDC_API_BASE + "data/"
GDC_FILES_ENDPOINT = GDC_API_BASE + "files/"
GDC_CASES_ENDPOINT = GDC_API_BASE + "cases/"
GDC_PROJECTS_ENDPOINT = GDC_API_BASE + "projects/"

def query_gdc_api(endpoint, params=None, method="GET", data=None):
    """Query the GDC API with proper error handling"""
    try:
        if method == "GET":
            response = requests.get(endpoint, params=params, timeout=60)
        else:  # POST
            headers = {"Content-Type": "application/json"}
            response = requests.post(endpoint, headers=headers, data=json.dumps(data), timeout=60)
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error querying GDC API: {e}")
        return None

def list_tcga_projects():
    """List all available TCGA projects"""
    params = {
        "filters": json.dumps({
            "op": "and",
            "content": [{
                "op": "=",
                "content": {
                    "field": "program.name",
                    "value": "TCGA"
                }
            }]
        }),
        "size": 100,
        "from": 0
    }
    
    response = query_gdc_api(GDC_PROJECTS_ENDPOINT, params=params)
    if not response:
        return []
    
    projects = []
    for project in response.get('data', {}).get('hits', []):
        projects.append({
            'project_id': project.get('project_id'),
            'name': project.get('name'),
            'disease_type': project.get('disease_type'),
            'primary_site': project.get('primary_site'),
            'case_count': project.get('summary', {}).get('case_count', 0)
        })
    
    return pd.DataFrame(projects)

def search_files(project_id=None, data_category=None, data_type=None, workflow_type=None, limit=100):
    """
    Search for files in the GDC database with optional filters
    
    Parameters:
    - project_id: TCGA project identifier (e.g., 'TCGA-BRCA')
    - data_category: e.g., 'Transcriptome Profiling', 'Clinical', 'DNA Methylation'
    - data_type: e.g., 'Gene Expression Quantification', 'Methylation Beta Value'
    - workflow_type: e.g., 'HTSeq - Counts', 'MuSE Variant Aggregation and Masking'
    - limit: Maximum number of results to return
    
    Returns:
    - DataFrame with file information
    """
    filters = {"op": "and", "content": []}
    
    # Add program filter for TCGA
    filters["content"].append({
        "op": "=",
        "content": {
            "field": "cases.project.program.name",
            "value": "TCGA"
        }
    })
    
    # Add optional filters
    if project_id:
        filters["content"].append({
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": project_id
            }
        })
    
    if data_category:
        filters["content"].append({
            "op": "=",
            "content": {
                "field": "data_category",
                "value": data_category
            }
        })
    
    if data_type:
        filters["content"].append({
            "op": "=",
            "content": {
                "field": "data_type",
                "value": data_type
            }
        })
    
    if workflow_type:
        filters["content"].append({
            "op": "=",
            "content": {
                "field": "analysis.workflow_type",
                "value": workflow_type
            }
        })
    
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,"
                 "analysis.workflow_type,file_size,access,cases.project.project_id",
        "format": "JSON",
        "size": limit,
        "from": 0
    }
    
    response = query_gdc_api(GDC_FILES_ENDPOINT, params=params)
    if not response:
        return pd.DataFrame()
    
    files = []
    for file_info in response.get('data', {}).get('hits', []):
        cases = file_info.get('cases', [])
        case_ids = [case.get('submitter_id') for case in cases]
        project_ids = [case.get('project', {}).get('project_id') for case in cases]
        
        files.append({
            'file_id': file_info.get('file_id'),
            'file_name': file_info.get('file_name'),
            'case_ids': ','.join(filter(None, case_ids)),
            'project_id': ','.join(filter(None, project_ids)),
            'data_category': file_info.get('data_category'),
            'data_type': file_info.get('data_type'),
            'workflow_type': file_info.get('analysis', {}).get('workflow_type'),
            'file_size': file_info.get('file_size'),
            'access': file_info.get('access')
        })
    
    return pd.DataFrame(files)

def download_files(file_ids, output_dir="downloads", unpack=False):
    """
    Download files from GDC using their file IDs
    
    Parameters:
    - file_ids: List of file IDs to download
    - output_dir: Directory where files will be saved
    - unpack: Whether to unpack tar.gz and gz files after download
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare data for the POST request
    data = {"ids": file_ids}
    
    print(f"Requesting download for {len(file_ids)} files...")
    response = requests.post(GDC_DATA_ENDPOINT, data=json.dumps(data), headers={"Content-Type": "application/json"})
    
    if response.status_code != 200:
        print(f"Error: {response.status_code} - {response.text}")
        return
    
    # Handle chunked download with progress bar
    response_head = response.headers
    content_disposition = response_head.get("Content-Disposition")
    
    if "attachment" in content_disposition:
        # Single file download
        filename = content_disposition.split("filename=")[1].strip("\"'")
        filepath = os.path.join(output_dir, filename)
        
        with open(filepath, "wb") as f, tqdm(
            desc=f"Downloading {filename}",
            total=int(response_head.get("Content-Length", 0)),
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=1024*1024):  # 1MB chunks
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
        
        print(f"Downloaded: {filepath}")
        
        # Unpack if requested
        if unpack and filepath.endswith((".tar.gz", ".gz")):
            _unpack_file(filepath, output_dir)
            
    else:
        # Multiple files as a tarball
        with open(os.path.join(output_dir, "gdc_download.tar.gz"), "wb") as f, tqdm(
            desc="Downloading files",
            total=int(response_head.get("Content-Length", 0)),
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=1024*1024):  # 1MB chunks
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
        
        archive_path = os.path.join(output_dir, "gdc_download.tar.gz")
        print(f"Downloaded archive: {archive_path}")
        
        # Unpack the archive
        if unpack:
            print(f"Unpacking {archive_path}...")
            with tarfile.open(archive_path) as tar:
                tar.extractall(path=output_dir)
            print(f"Files extracted to {output_dir}")

def _unpack_file(filepath, output_dir):
    """Helper function to unpack a compressed file"""
    if filepath.endswith(".tar.gz"):
        print(f"Unpacking {filepath}...")
        with tarfile.open(filepath) as tar:
            tar.extractall(path=output_dir)
        print(f"Extracted to {output_dir}")
        
    elif filepath.endswith(".gz"):
        print(f"Unpacking {filepath}...")
        output_file = filepath[:-3]  # Remove .gz extension
        with gzip.open(filepath, 'rb') as f_in, open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Extracted to {output_file}")

def download_clinical_data(project_id, output_dir="downloads"):
    """Download clinical data for a specific TCGA project"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Search for clinical files
    files_df = search_files(
        project_id=project_id,
        data_category="Clinical",
        data_type="Clinical Supplement",
        limit=100
    )
    
    if files_df.empty:
        print(f"No clinical data found for project {project_id}")
        return
    
    # Download the clinical files
    download_files(files_df['file_id'].tolist(), output_dir=output_dir, unpack=True)
    
    print(f"Clinical data for {project_id} downloaded to {output_dir}")

def download_expression_data(project_id, workflow_type="HTSeq - Counts", output_dir="downloads"):
    """Download gene expression data for a specific TCGA project"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Search for gene expression files
    files_df = search_files(
        project_id=project_id,
        data_category="Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type=workflow_type,
        limit=1000
    )
    
    if files_df.empty:
        print(f"No gene expression data found for project {project_id} with workflow {workflow_type}")
        return
    
    # Download the expression files
    download_files(files_df['file_id'].tolist(), output_dir=output_dir, unpack=True)
    
    print(f"Gene expression data for {project_id} downloaded to {output_dir}")

def download_methylation_data(project_id, output_dir="downloads"):
    """Download DNA methylation data for a specific TCGA project"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Search for methylation files
    files_df = search_files(
        project_id=project_id,
        data_category="DNA Methylation",
        data_type="Methylation Beta Value",
        limit=1000
    )
    
    if files_df.empty:
        print(f"No methylation data found for project {project_id}")
        return
    
    # Download the methylation files
    download_files(files_df['file_id'].tolist(), output_dir=output_dir, unpack=True)
    
    print(f"Methylation data for {project_id} downloaded to {output_dir}")

def download_mutation_data(project_id, output_dir="downloads"):
    """Download somatic mutation data for a specific TCGA project"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Search for mutation files
    files_df = search_files(
        project_id=project_id,
        data_category="Simple Nucleotide Variation",
        data_type="Masked Somatic Mutation",
        workflow_type="MuSE Variant Aggregation and Masking",
        limit=1000
    )
    
    if files_df.empty:
        print(f"No mutation data found for project {project_id}")
        return
    
    # Download the mutation files
    download_files(files_df['file_id'].tolist(), output_dir=output_dir, unpack=True)
    
    print(f"Mutation data for {project_id} downloaded to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Download TCGA data from GDC")
    
    # Command arguments
    parser.add_argument("--list-projects", action="store_true", help="List all TCGA projects")
    parser.add_argument("--search", action="store_true", help="Search for files")
    parser.add_argument("--download", action="store_true", help="Download files by ID")
    parser.add_argument("--clinical", action="store_true", help="Download clinical data")
    parser.add_argument("--expression", action="store_true", help="Download gene expression data")
    parser.add_argument("--methylation", action="store_true", help="Download DNA methylation data")
    parser.add_argument("--mutation", action="store_true", help="Download mutation data")
    
    # Options
    parser.add_argument("--project", type=str, help="TCGA project ID (e.g., TCGA-BRCA)")
    parser.add_argument("--category", type=str, help="Data category")
    parser.add_argument("--data-type", type=str, help="Data type")
    parser.add_argument("--workflow", type=str, help="Workflow type")
    parser.add_argument("--file-ids", type=str, help="Comma-separated list of file IDs to download")
    parser.add_argument("--limit", type=int, default=100, help="Limit for search results")
    parser.add_argument("--output-dir", type=str, default="downloads", help="Output directory for downloads")
    parser.add_argument("--unpack", action="store_true", help="Unpack downloaded files")
    
    args = parser.parse_args()
    
    if args.list_projects:
        print("Listing TCGA projects...")
        projects_df = list_tcga_projects()
        if not projects_df.empty:
            print(projects_df)
        else:
            print("No projects found or error occurred.")
    
    elif args.search:
        print("Searching for files...")
        files_df = search_files(
            project_id=args.project,
            data_category=args.category,
            data_type=args.data_type,
            workflow_type=args.workflow,
            limit=args.limit
        )
        
        if not files_df.empty:
            print(f"Found {len(files_df)} files:")
            print(files_df)
            # Save search results
            csv_file = os.path.join(args.output_dir, "search_results.csv")
            os.makedirs(args.output_dir, exist_ok=True)
            files_df.to_csv(csv_file, index=False)
            print(f"Search results saved to {csv_file}")
        else:
            print("No files found matching the criteria.")
    
    elif args.download and args.file_ids:
        file_ids = [fid.strip() for fid in args.file_ids.split(",")]
        download_files(file_ids, output_dir=args.output_dir, unpack=args.unpack)
    
    elif args.clinical and args.project:
        download_clinical_data(args.project, output_dir=args.output_dir)
    
    elif args.expression and args.project:
        download_expression_data(
            args.project, 
            workflow_type=args.workflow or "HTSeq - Counts", 
            output_dir=args.output_dir
        )
    
    elif args.methylation and args.project:
        download_methylation_data(args.project, output_dir=args.output_dir)
    
    elif args.mutation and args.project:
        download_mutation_data(args.project, output_dir=args.output_dir)
    
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

