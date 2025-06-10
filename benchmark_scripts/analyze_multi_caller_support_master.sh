#!/bin/bash

# Master script to analyze multi-caller support rates across multiple samples.
# Generates violin plots and heatmap tables for SV support analysis across N callers.

set -e  # Exit on error

# Function to show usage
show_usage() {
    echo "Usage: $0 --mapping MAPPING_FILE --input-dir INPUT_DIR [--n-matches N_MATCHES] [--svtypes SVTYPES] [--output-dir OUTPUT_DIR] [--use-regular-tsv]"
    echo ""
    echo "Analyze multi-caller support rates across all available samples."
    echo ""
    echo "Arguments:"
    echo "  --mapping MAPPING_FILE    Path to mapping file with SR_ID and LR_ID columns"
    echo "  --input-dir INPUT_DIR     Input directory containing multi-caller results"
    echo "  --n-matches N_MATCHES     Maximum number of caller matches to analyze (generates 1+, 2+, ..., N+ subfigures) (default: 2)"
    echo "  --svtypes SVTYPES         Comma-separated list of SV types to include (e.g., DEL,DUP,INS) (default: all types)"
    echo "  --output-dir OUTPUT_DIR   Output directory for plots and tables (default: multi_caller_support_analysis)"
    echo "  --use-regular-tsv         Use regular .coalesced.tsv files instead of .combined.coalesced.tsv files"
    echo ""
    echo "Output files:"
    echo "  - violin_combined_support_N_callers.png         - Combined violin plot with 1+ through N+ subfigures"
    echo "  - heatmap_us_rm_regions_N_callers.png          - Heatmap for US/RM regions with 1+ through N+ subfigures"
    echo "  - heatmap_sd_sr_regions_N_callers.png          - Heatmap for SD/SR regions with 1+ through N+ subfigures"
    echo "  - support_rates_us_rm_regions_X_callers.csv    - Support rate data for US/RM regions (X = 1, 2, ..., N)"
    echo "  - support_rates_sd_sr_regions_X_callers.csv    - Support rate data for SD/SR regions (X = 1, 2, ..., N)"
    echo "  - summary_stats_X_callers.txt                  - Summary statistics (X = 1, 2, ..., N)"
    echo ""
    echo "Examples:"
    echo "  $0 --mapping Mapping.tsv --input-dir cutesv_sniffles_results --n-matches 2"
    echo "  $0 --mapping Mapping.tsv --input-dir cutesv_sniffles_results --n-matches 3 --svtypes DEL,DUP,INS"
}

# Check if required tools are available
check_requirements() {
    if ! command -v python3 &> /dev/null; then
        echo "Error: python3 is required but not found in PATH"
        exit 1
    fi
    
    # Check if required Python packages are available
    python3 -c "import pandas, numpy, matplotlib, seaborn" 2>/dev/null || {
        echo "Error: Required Python packages not found (pandas, numpy, matplotlib, seaborn)"
        echo "Please install them using: pip install pandas numpy matplotlib seaborn"
        exit 1
    }
}

# Verify file exists
verify_file_exists() {
    local file_path="$1"
    local description="$2"
    
    if [ ! -f "$file_path" ]; then
        echo "Error: $description not found: $file_path"
        return 1
    fi
    
    if [ ! -r "$file_path" ]; then
        echo "Error: Cannot read $description: $file_path"
        return 1
    fi
    
    return 0
}

# Main script
if [ "$#" -lt 2 ]; then
    show_usage
    exit 1
fi

# Parse arguments
MAPPING_FILE=""
INPUT_DIR=""
N_MATCHES=2
SVTYPES=""
OUTPUT_DIR="multi_caller_support_analysis"
USE_REGULAR_TSV=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --mapping)
            MAPPING_FILE="$2"
            shift 2
            ;;
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --n-matches)
            N_MATCHES="$2"
            shift 2
            ;;
        --svtypes)
            SVTYPES="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --use-regular-tsv)
            USE_REGULAR_TSV=true
            shift
            ;;
        --help|-h)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$MAPPING_FILE" ] || [ -z "$INPUT_DIR" ]; then
    echo "Error: --mapping and --input-dir are required arguments"
    show_usage
    exit 1
fi

# Check requirements
check_requirements

# Verify mapping file exists
verify_file_exists "$MAPPING_FILE" "Mapping file" || exit 1

# Verify input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Get the directory where this script is located (should be benchmark_scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_SCRIPT="${SCRIPT_DIR}/analyze_multi_caller_support.py"

# Verify the analysis script exists
verify_file_exists "$ANALYSIS_SCRIPT" "Analysis script" || exit 1

echo ""
echo "========================================================"
echo "MULTI-CALLER SUPPORT ANALYSIS PIPELINE"
echo "========================================================"
echo "Mapping file: $MAPPING_FILE"
echo "Input directory: $INPUT_DIR"
echo "N-matches (max): $N_MATCHES"
echo "SV types filter: ${SVTYPES:-all types}"
echo "Output directory: $OUTPUT_DIR"
echo "Use regular TSV files: $USE_REGULAR_TSV"
echo ""

# Prepare arguments for Python script
PYTHON_ARGS="--mapping $MAPPING_FILE --input-dir $INPUT_DIR --n-matches $N_MATCHES --output-dir $OUTPUT_DIR"

if [ -n "$SVTYPES" ]; then
    PYTHON_ARGS="$PYTHON_ARGS --svtypes $SVTYPES"
fi

if [ "$USE_REGULAR_TSV" = false ]; then
    PYTHON_ARGS="$PYTHON_ARGS --use-combined"
fi

# Run the analysis
echo "Running multi-caller support analysis..."
python3 "$ANALYSIS_SCRIPT" $PYTHON_ARGS

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================"
    echo "ANALYSIS COMPLETE"
    echo "========================================================"
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    echo "  - violin_combined_support_${N_MATCHES}_callers.png"
    echo "  - heatmap_us_rm_regions_${N_MATCHES}_callers.png"
    echo "  - heatmap_sd_sr_regions_${N_MATCHES}_callers.png"
    for ((i=1; i<=N_MATCHES; i++)); do
        echo "  - support_rates_us_rm_regions_${i}_callers.csv"
        echo "  - support_rates_sd_sr_regions_${i}_callers.csv"
        echo "  - summary_stats_${i}_callers.txt"
    done
    echo "========================================================"
else
    echo ""
    echo "Error: Analysis failed. Check the output above for details."
    exit 1
fi 