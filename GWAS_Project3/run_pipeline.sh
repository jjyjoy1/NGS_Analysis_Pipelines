#!/bin/bash

# Run GWAS Pipeline using Snakemake
# This script provides a convenient way to run the entire workflow

# Default values
CONFIG="config.yaml"
PROFILE="slurm"
CORES=1
DRY_RUN=false
REASON=false
UNLOCK=false
RERUN_INCOMPLETE=true
FORCE=false

# Print usage information
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -h, --help               Show this help message"
    echo "  -c, --config FILE        Config file (default: config.yaml)"
    echo "  -p, --profile PROFILE    Profile to use: 'slurm' or 'local' (default: slurm)"
    echo "  --cores N                Number of cores for local execution (default: 1)"
    echo "  -n, --dry-run            Dry run (don't execute commands)"
    echo "  -r, --reason             Print the reason for rule execution"
    echo "  -u, --unlock             Unlock the working directory"
    echo "  -f, --force              Force the execution of rules"
    echo "  --no-rerun-incomplete    Do not rerun incomplete jobs"
    echo ""
    echo "Example:"
    echo "  $0 --config my_config.yaml --profile slurm"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            usage
            exit 0
            ;;
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -r|--reason)
            REASON=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        --no-rerun-incomplete)
            RERUN_INCOMPLETE=false
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Check if config file exists
if [ ! -f "$CONFIG" ]; then
    echo "Error: Config file '$CONFIG' not found."
    exit 1
fi

# Create logs directory
mkdir -p logs
mkdir -p logs/slurm

# Build snakemake command
SNAKEMAKE_CMD="snakemake --configfile $CONFIG"

# Add profile if specified
if [ "$PROFILE" == "slurm" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --profile profile/slurm"
elif [ "$PROFILE" == "local" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --cores $CORES"
else
    echo "Error: Profile '$PROFILE' not recognized. Use 'slurm' or 'local'."
    exit 1
fi

# Add optional flags
if [ "$DRY_RUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD -n"
fi

if [ "$REASON" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD -r"
fi

if [ "$UNLOCK" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --unlock"
fi

if [ "$RERUN_INCOMPLETE" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --rerun-incomplete"
fi

if [ "$FORCE" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD -F"
fi

# Print command
echo "Running: $SNAKEMAKE_CMD"

# Execute snakemake
eval $SNAKEMAKE_CMD

# Check exit status
STATUS=$?
if [ $STATUS -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "Pipeline failed with exit code $STATUS"
    exit $STATUS
fi


