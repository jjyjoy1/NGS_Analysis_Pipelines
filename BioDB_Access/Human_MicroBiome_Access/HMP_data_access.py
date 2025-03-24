#!/usr/bin/env python3
"""
HMP Data Download Tool
Downloads data from the Human Microbiome Project using a manifest file.
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
from pathlib import Path

def install_hmp_client():
    """Install the HMP client if not already installed"""
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "hmp_client"])
        print("HMP client installed successfully")
    except subprocess.CalledProcessError:
        print("Failed to install HMP client. Please install manually with 'pip install hmp_client'")
        sys.exit(1)

def check_hmp_client():
    """Check if HMP client is installed"""
    try:
        subprocess.check_call(["hmp_client", "--help"], 
                             stdout=subprocess.DEVNULL, 
                             stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def download_from_manifest(manifest_file, output_dir, max_threads=4):
    """Download files listed in the manifest file"""
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Check if manifest file exists
    if not os.path.exists(manifest_file):
        print(f"Error: Manifest file '{manifest_file}' not found")
        return False
    
    # Build command
    cmd = ["hmp_client", "-m", manifest_file, "--destination", output_dir, 
           "-p", str(max_threads), "download"]
    
    # Execute command
    print(f"Downloading files from manifest '{manifest_file}'...")
    try:
        subprocess.check_call(cmd)
        print("Download completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error downloading files: {e}")
        return False

def list_manifest_contents(manifest_file):
    """List contents of a manifest file"""
    try:
        # Try to read as TSV first
        manifest = pd.read_csv(manifest_file, sep='\t')
        print(f"Manifest contains {len(manifest)} files")
        
        # Print summary by data type if available
        if 'sample_type' in manifest.columns:
            print("\nFiles by sample type:")
            print(manifest['sample_type'].value_counts())
        
        # Print total size if available
        if 'size' in manifest.columns:
            total_size = manifest['size'].sum() / (1024**3)  # Convert to GB
            print(f"\nTotal download size: {total_size:.2f} GB")
            
        return True
    except Exception as e:
        print(f"Error reading manifest file: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Download data from the Human Microbiome Project")
    parser.add_argument("-m", "--manifest", required=True, help="Path to the manifest file")
    parser.add_argument("-o", "--output", default="hmp_data", help="Output directory for downloaded files")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Maximum number of download threads")
    parser.add_argument("-l", "--list", action="store_true", help="List contents of manifest without downloading")
    parser.add_argument("--install-client", action="store_true", help="Install HMP client if not already installed")
    
    args = parser.parse_args()
    
    # Check if HMP client is installed
    if not check_hmp_client():
        if args.install_client:
            install_hmp_client()
        else:
            print("HMP client not found. Run with --install-client to install automatically.")
            sys.exit(1)
    
    # List manifest contents if requested
    if args.list:
        list_manifest_contents(args.manifest)
    else:
        # Download files
        download_from_manifest(args.manifest, args.output, args.threads)

if __name__ == "__main__":
    main()

