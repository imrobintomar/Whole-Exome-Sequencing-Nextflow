#!/bin/bash

# ==========================================
# WES Pipeline Runner with Optimized JVM Settings
# ==========================================
# This script runs the Whole Exome Sequencing pipeline
# with proper Java memory configuration to prevent StackOverflowError
# ==========================================

set -e  # Exit on error

# Color codes for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored messages
print_header() {
    echo -e "${BLUE}╔════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║  WES Pipeline Runner - v1.0            ║${NC}"
    echo -e "${BLUE}╚════════════════════════════════════════╝${NC}"
}

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# JVM Settings for Nextflow
# These settings prevent StackOverflowError and optimize memory usage
export NXF_OPTS="-Xms2g -Xmx8g -Xss4m"
# -Xms2g  : Initial heap size (2GB)
# -Xmx8g  : Maximum heap size (8GB) - adjust based on your system
# -Xss4m  : Thread stack size (4MB) - fixes StackOverflowError

print_header
echo ""

# Display JVM settings
print_info "JVM Settings configured:"
echo "  - Initial Heap: 2GB"
echo "  - Max Heap: 8GB"
echo "  - Stack Size: 4MB"
echo ""

# Check if Nextflow is available
if ! command -v nextflow &> /dev/null; then
    print_error "Nextflow is not installed or not in PATH"
    print_info "Please install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html"
    exit 1
fi

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

print_info "Working directory: $SCRIPT_DIR"
echo ""

# Parse command line arguments
RESUME_FLAG=""
PROFILE=""
ADDITIONAL_ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -resume|--resume)
            RESUME_FLAG="-resume"
            print_info "Resume mode enabled - will skip completed tasks"
            shift
            ;;
        -profile|--profile)
            PROFILE="-profile $2"
            print_info "Using profile: $2"
            shift 2
            ;;
        -clean|--clean)
            print_warn "Cleaning work directory..."
            rm -rf work/
            rm -rf .nextflow/
            print_info "Work directory cleaned"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -resume, --resume       Resume from last successful step"
            echo "  -profile PROFILE        Use specific profile (docker, singularity, etc.)"
            echo "  -clean, --clean         Clean work directory before running"
            echo "  -h, --help             Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                     # Run pipeline from scratch"
            echo "  $0 -resume             # Resume failed pipeline"
            echo "  $0 -profile docker     # Run with Docker profile"
            echo "  $0 -clean -resume      # Clean and resume"
            exit 0
            ;;
        *)
            ADDITIONAL_ARGS="$ADDITIONAL_ARGS $1"
            shift
            ;;
    esac
done

# Check if main.nf exists
if [ ! -f "main.nf" ]; then
    print_error "main.nf not found in current directory"
    exit 1
fi

# Check if nextflow.config exists
if [ ! -f "nextflow.config" ]; then
    print_warn "nextflow.config not found - using default configuration"
fi

# Build the command
CMD="nextflow run main.nf $RESUME_FLAG $PROFILE $ADDITIONAL_ARGS"

print_info "Executing pipeline..."
echo -e "${BLUE}Command:${NC} $CMD"
echo ""
print_info "Pipeline started at: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Run the pipeline
$CMD
EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    print_info "✓ Pipeline completed successfully at: $(date '+%Y-%m-%d %H:%M:%S')"
    print_info "Check results in the output directory"
else
    print_error "✗ Pipeline failed with exit code: $EXIT_CODE"
    print_info "To resume from last successful step, run:"
    echo "  $0 -resume"
fi

exit $EXIT_CODE
