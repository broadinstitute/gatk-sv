#!/bin/bash

# Process multiple samples through the SV callset evaluation pipeline.
# Takes a mapping file of SR_ID to LR_ID and processes a specified number of samples.

set -e  # Exit on error

# Get the directory where this script is located (should be benchmark_scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GATK_SV_ROOT="$(dirname "$SCRIPT_DIR")"

# Check if required tools are available
check_requirements() {
    local required_tools=("bcftools" "tabix" "svtk" "bgzip")
    local missing_tools=()
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        echo "Error: The following required tools are not installed or not in PATH: ${missing_tools[*]}"
        exit 1
    fi
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

# Wait for file to exist
wait_for_file() {
    local file_path="$1"
    local timeout="$2"
    local interval=5
    local elapsed=0
    
    while [ ! -f "$file_path" ] && [ $elapsed -lt $timeout ]; do
        sleep $interval
        elapsed=$((elapsed + interval))
        echo "Waiting for file $file_path... ($elapsed/$timeout seconds)"
    done
    
    if [ ! -f "$file_path" ]; then
        echo "Error: File $file_path not found after $timeout seconds"
        return 1
    fi
    
    return 0
}

# Process short-read data for a sample
process_sr_sample() {
    local sr_id="$1"
    local lr_id="$2"
    local base_dir="$3"
    
    echo ""
    echo "--- SHORT-READ PROCESSING ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    mkdir -p "$sample_dir"
    
    local sr_vcf="${sample_dir}/${sr_id}.sr.vcf.gz"
    local sr_filtered_vcf="${sample_dir}/${sr_id}.sr.filtered.vcf.gz"
    local sr_filtered_revised_vcf="${sample_dir}/${sr_id}.sr.filtered.revised.vcf.gz"
    local sr_bed="${sample_dir}/${sr_id}.sr.bed"
    local sr_annotated_bed="${sample_dir}/${sr_id}.sr.annotated.bed"
    
    # Use configurable input VCF path
    local input_vcf="${base_dir}/${SR_VCF_PATH}"
    
    local ref_dir="${SCRIPT_DIR}/input"
    local rm_ref="${ref_dir}/hg38.RM.sorted.merged.bed.gz"
    local sd_ref="${ref_dir}/hg38.SD.sorted.merged.bed.gz"
    local sr_ref="${ref_dir}/hg38.SR.sorted.merged.bed.gz"
    
    # Verify required files exist
    verify_file_exists "$rm_ref" "Reference file (RM)" || return 1
    verify_file_exists "$sd_ref" "Reference file (SD)" || return 1
    verify_file_exists "$sr_ref" "Reference file (SR)" || return 1
    verify_file_exists "$input_vcf" "Input VCF file" || return 1
    
    echo "Extracting sample $sr_id from joint VCF..."
    bcftools view -s "$sr_id" -Oz -o "$sr_vcf" "$input_vcf" || return 1
    tabix "$sr_vcf" || return 1
    
    echo "Filtering short-read variants..."
    bcftools view -i 'GT!="0/0" && GT!="./."' -Oz -o "$sr_filtered_vcf" "$sr_vcf" || return 1
    tabix "$sr_filtered_vcf" || return 1
    
    echo "Annotating SV types..."
    "${SCRIPT_DIR}/annotate_svtype.sh" --in "$sr_filtered_vcf" --out "$sr_filtered_revised_vcf" || return 1
    
    echo "Converting to BED format..."
    svtk vcf2bed "$sr_filtered_revised_vcf" "$sr_bed" --info SVLEN --info SVTYPE --info AF --no-samples || return 1
    
    echo "Extracting filter information..."
    python3 "${SCRIPT_DIR}/extract_filter_info.py" --vcf "$sr_filtered_revised_vcf" --output "${sample_dir}/${sr_id}.sr.filter_info.tsv" || return 1
    
    echo "Annotating genomic context..."
    "${GATK_SV_ROOT}/module18_annotate_genomic_context/annotate_genomic_context.sh" "$sr_bed" "$sr_annotated_bed" "$rm_ref" "$sd_ref" "$sr_ref" || return 1
    
    return 0
}

# Process long-read data for a sample
process_lr_sample() {
    local sr_id="$1"
    local lr_id="$2"
    local base_dir="$3"
    local lr_caller="$4"
    
    echo ""
    echo "--- LONG-READ PROCESSING ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    mkdir -p "$sample_dir"
    
    # Use VCF from specified LR caller directory with appropriate naming
    # Different callers have different filename patterns
    local lr_source_file
    case "$lr_caller" in
        "sniffles")
            lr_source_file="${base_dir}/${LR_VCF_PATH}/${lr_id}.hg38.${lr_caller}.sv.phased.vcf.gz"
            ;;
        "cutesv")
            lr_source_file="${base_dir}/${LR_VCF_PATH}/${lr_id}.hg38.${lr_caller}.sv.vcf.gz"
            ;;
        *)
            # Use configurable filename pattern
            lr_source_file="${base_dir}/${LR_VCF_PATH}/${LR_VCF_PATTERN}"
            # Replace placeholders
            lr_source_file="${lr_source_file//\{LR_ID\}/$lr_id}"
            lr_source_file="${lr_source_file//\{LR_CALLER\}/$lr_caller}"
            ;;
    esac
    
    echo "Source LR VCF: $lr_source_file"
    
    local lr_vcf="${sample_dir}/${sr_id}.lr.vcf.gz"
    local lr_filtered_vcf="${sample_dir}/${sr_id}.lr.filtered.vcf.gz"
    local lr_bed="${sample_dir}/${sr_id}.lr.bed"
    local lr_annotated_bed="${sample_dir}/${sr_id}.lr.annotated.bed"
    
    local ref_dir="${SCRIPT_DIR}/input"
    local rm_ref="${ref_dir}/hg38.RM.sorted.merged.bed.gz"
    local sd_ref="${ref_dir}/hg38.SD.sorted.merged.bed.gz"
    local sr_ref="${ref_dir}/hg38.SR.sorted.merged.bed.gz"
    
    # Verify required files exist
    verify_file_exists "$rm_ref" "Reference file (RM)" || return 1
    verify_file_exists "$sd_ref" "Reference file (SD)" || return 1
    verify_file_exists "$sr_ref" "Reference file (SR)" || return 1
    verify_file_exists "$lr_source_file" "Long-read VCF file" || return 1
    
    echo "Processing long-read data for $lr_id..."
    
    echo "Filtering long-read variants..."
    bcftools view -i 'GT!="0/0" && GT!="./."' -Oz -o "$lr_vcf" "$lr_source_file" || return 1
    tabix "$lr_vcf" || return 1
    
    echo "Annotating SV types..."
    "${SCRIPT_DIR}/annotate_svtype.sh" --in "$lr_vcf" --out "$lr_filtered_vcf" || return 1
    
    echo "Converting to BED format..."
    svtk vcf2bed "$lr_filtered_vcf" "$lr_bed" --info SVLEN --info SVTYPE --no-samples || return 1
    
    echo "Extracting filter information..."
    python3 "${SCRIPT_DIR}/extract_filter_info.py" --vcf "$lr_filtered_vcf" --output "${sample_dir}/${sr_id}.lr.filter_info.tsv" || return 1
    
    echo "Annotating genomic context..."
    "${GATK_SV_ROOT}/module18_annotate_genomic_context/annotate_genomic_context.sh" "$lr_bed" "$lr_annotated_bed" "$rm_ref" "$sd_ref" "$sr_ref" || return 1
    
    return 0
}

# Run callset comparisons
compare_callsets_sample() {
    local sr_id="$1"
    local lr_id="$2" 
    local base_dir="$3"
    
    echo ""
    echo "--- CALLSET COMPARISONS ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    
    # Check if required input files exist
    verify_file_exists "${sample_dir}/${sr_id}.sr.bed" "SR BED file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.bed" "LR BED file" || return 1
    
    echo "Processing BED files for comparison..."
    
    # Build the command with optional limit arguments
    cmd="python3 \"${SCRIPT_DIR}/process_bed_files.py\" --sample \"$sr_id\" --sr \"${sample_dir}/${sr_id}.sr.bed\" --lr \"${sample_dir}/${sr_id}.lr.bed\""
    
    if [ -n "$SR_LIMIT" ]; then
        cmd="$cmd --sr-limit \"$SR_LIMIT\""
    fi
    
    if [ -n "$LR_LIMIT" ]; then
        cmd="$cmd --lr-limit \"$LR_LIMIT\""
    fi
    
    eval "$cmd" || return 1
    
    # Verify output files were created
    verify_file_exists "${sample_dir}/${sr_id}.sr.query.bed" "SR query BED" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.sr.ref.bed" "SR reference BED" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.query.bed" "LR query BED" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.ref.bed" "LR reference BED" || return 1
    
    echo "Running callset comparisons..."
    
    # Determine which comparisons to run based on COMPARISON_DIRECTIONS
    if [[ "$COMPARISON_DIRECTIONS" == *"SR"* ]]; then
        echo "Running LR query vs SR reference comparison..."
        "${SCRIPT_DIR}/compare_callsets_V2.sh" -O "${sample_dir}/${sr_id}.lr_query.bed" -p lr_query -d "$DISTANCE" "${sample_dir}/${sr_id}.lr.query.bed" "${sample_dir}/${sr_id}.sr.ref.bed" || return 1
    fi
    
    if [[ "$COMPARISON_DIRECTIONS" == *"LR"* ]]; then
        echo "Running SR query vs LR reference comparison..."
        "${SCRIPT_DIR}/compare_callsets_V2.sh" -O "${sample_dir}/${sr_id}.sr_query.bed" -p sr_query -d "$DISTANCE" "${sample_dir}/${sr_id}.sr.query.bed" "${sample_dir}/${sr_id}.lr.ref.bed" || return 1
    fi
    
    return 0
}

# Analyze match rates using analyze_match_rates.py
analyze_rates_sample() {
    local sr_id="$1"
    local lr_id="$2" 
    local base_dir="$3"
    
    echo ""
    echo "--- MATCH RATE ANALYSIS ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    
    # Check if required input files exist
    verify_file_exists "${sample_dir}/${sr_id}.sr.annotated.bed" "SR annotated file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.annotated.bed" "LR annotated file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.sr.filter_info.tsv" "SR filter info file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.filter_info.tsv" "LR filter info file" || return 1
    
    # Check for comparison files based on COMPARISON_DIRECTIONS and analyze each separately
    local analysis_run=false
    
    if [[ "$COMPARISON_DIRECTIONS" == *"SR"* ]] && [[ -f "${sample_dir}/${sr_id}.lr_query.bed" ]]; then
        echo "Analyzing LR query match rates..."
        # LR query uses SR annotated bed and SR filter info (since LR query contains SR variant IDs)
        python3 "${SCRIPT_DIR}/analyze_match_rates.py" \
            --query-bed "${sample_dir}/${sr_id}.lr_query.bed" \
            --annotated-bed "${sample_dir}/${sr_id}.sr.annotated.bed" \
            --filter-info "${sample_dir}/${sr_id}.sr.filter_info.tsv" \
            --output-prefix "${sample_dir}/lr_query" \
            --query-type lr_query || return 1
        analysis_run=true
    fi
    
    if [[ "$COMPARISON_DIRECTIONS" == *"LR"* ]] && [[ -f "${sample_dir}/${sr_id}.sr_query.bed" ]]; then
        echo "Analyzing SR query match rates..."
        # SR query uses LR annotated bed and LR filter info (since SR query contains LR variant IDs)
        python3 "${SCRIPT_DIR}/analyze_match_rates.py" \
            --query-bed "${sample_dir}/${sr_id}.sr_query.bed" \
            --annotated-bed "${sample_dir}/${sr_id}.lr.annotated.bed" \
            --filter-info "${sample_dir}/${sr_id}.lr.filter_info.tsv" \
            --output-prefix "${sample_dir}/sr_query" \
            --query-type sr_query || return 1
        analysis_run=true
    fi
    
    if [[ "$analysis_run" == "false" ]]; then
        echo "Error: No comparison files found. Please run COMPARE_CALLSETS process first."
        return 1
    fi
    
    echo "Match rate analysis completed!"
    return 0
}

# Process a single sample with specified process types
process_sample() {
    local sr_id="$1"
    local lr_id="$2"
    local base_dir="$3"
    local processes="$4"
    
    echo ""
    echo "=========================================="
    echo "PROCESSING SAMPLE: $sr_id (LR: $lr_id)"
    echo "Processes: $processes"
    echo "=========================================="
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    mkdir -p "$sample_dir"
    
    # Run specified processes
    if [[ "$processes" == *"PROCESS_SR"* ]]; then
        if ! process_sr_sample "$sr_id" "$lr_id" "$base_dir"; then
            return 1
        fi
    fi
    
    if [[ "$processes" == *"PROCESS_LR"* ]]; then
        if ! process_lr_sample "$sr_id" "$lr_id" "$base_dir" "$LR_CALLER"; then
            return 1
        fi
    fi
    
    if [[ "$processes" == *"COMPARE_CALLSETS"* ]]; then
        if ! compare_callsets_sample "$sr_id" "$lr_id" "$base_dir"; then
            return 1
        fi
    fi
    
    if [[ "$processes" == *"ANALYZE_RATES"* ]]; then
        if ! analyze_rates_sample "$sr_id" "$lr_id" "$base_dir"; then
            return 1
        fi
    fi
    
    echo ""
    echo "✓ Sample $sr_id processing completed successfully"
    echo "=========================================="
    return 0
}

# Show usage information
show_usage() {
    echo "Usage: $0 --mapping MAPPING_FILE --base-dir BASE_DIR [OPTIONS]"
    echo ""
    echo "Long-read SV benchmarking pipeline for processing multiple samples."
    echo ""
    echo "Required arguments:"
    echo "  --mapping MAPPING_FILE          Path to mapping file with SR_ID and LR_ID columns"
    echo "  --base-dir BASE_DIR             Base directory for input data"
    echo ""
    echo "Optional arguments:"
    echo "  --num-samples NUM_SAMPLES       Number of samples to process (default: 3)"
    echo "  --lr-caller LR_CALLER           Long-read caller: sniffles, cutesv, svim, etc. (default: sniffles)"
    echo "  --comparison-directions DIRS    Comparison directions: SR, LR, or SR,LR (default: SR,LR)"
    echo "  --sr-limit SR_LIMIT             SR variant limit for testing (default: no limit)"
    echo "  --lr-limit LR_LIMIT             LR variant limit for testing (default: no limit)"
    echo "  --distance DISTANCE             Maximum distance between breakpoints (default: 50)"
    echo "  --sr-vcf-path SR_VCF_PATH       Path to SR VCF relative to base-dir (default: all_batches.filter_genotypes.sanitized.135.overlaps_tagged.vcf.gz)"
    echo "  --lr-vcf-path LR_VCF_PATH       Directory for LR VCFs relative to base-dir (default: {lr-caller}_vcfs)"
    echo "  --lr-vcf-pattern PATTERN        Filename pattern for LR VCFs (default: auto-detect, use {LR_ID} and {LR_CALLER} as placeholders)"
    echo "  --output-dir OUTPUT_DIR         Output directory (default: {lr-caller}_results)"
    echo ""
    echo "Process types (comma-separated):"
    echo "  PROCESS_SR    - Process short-read data (VCF filtering, annotation, BED conversion)"
    echo "  PROCESS_LR    - Process long-read data (VCF filtering, annotation, BED conversion)"
    echo "  COMPARE_CALLSETS - Run callset comparison"
    echo "  ANALYZE_RATES - Run match rate analysis (requires comparison results)"
    echo ""
    echo "Examples:"
    echo "  # Basic usage with defaults"
    echo "  $0 --mapping Mapping.tsv --base-dir ."
    echo ""
    echo "  # Custom SR VCF and distance parameter"
    echo "  $0 --mapping Mapping.tsv --base-dir . --sr-vcf-path custom_sr_callset.vcf.gz --distance 100"
    echo ""
    echo "  # Custom LR caller with pattern"
    echo "  $0 --mapping Mapping.tsv --base-dir . --lr-caller pbsv --lr-vcf-pattern '{LR_ID}.pbsv.vcf.gz'"
    echo ""
    echo "  # Process only 5 samples, run only analysis"
    echo "  $0 --mapping Mapping.tsv --base-dir . --num-samples 5 --processes ANALYZE_RATES"
}

# Main script
if [ "$#" -eq 0 ]; then
    show_usage
    exit 1
fi

# Parse arguments
MAPPING_FILE=""
BASE_DIR=""
NUM_SAMPLES=3
PROCESSES="PROCESS_SR,PROCESS_LR,COMPARE_CALLSETS,ANALYZE_RATES"  # Default: do everything
LR_CALLER=""
COMPARISON_DIRECTIONS=""
SR_LIMIT=""
LR_LIMIT=""
DISTANCE=""
SR_VCF_PATH=""
LR_VCF_PATH=""
LR_VCF_PATTERN=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --mapping)
            MAPPING_FILE="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --num-samples)
            NUM_SAMPLES="$2"
            shift 2
            ;;
        --processes)
            PROCESSES="$2"
            shift 2
            ;;
        --lr-caller)
            LR_CALLER="$2"
            shift 2
            ;;
        --comparison-directions)
            COMPARISON_DIRECTIONS="$2"
            shift 2
            ;;
        --sr-limit)
            SR_LIMIT="$2"
            shift 2
            ;;
        --lr-limit)
            LR_LIMIT="$2"
            shift 2
            ;;
        --distance)
            DISTANCE="$2"
            shift 2
            ;;
        --sr-vcf-path)
            SR_VCF_PATH="$2"
            shift 2
            ;;
        --lr-vcf-path)
            LR_VCF_PATH="$2"
            shift 2
            ;;
        --lr-vcf-pattern)
            LR_VCF_PATTERN="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --help|-h)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$MAPPING_FILE" ] || [ -z "$BASE_DIR" ]; then
    echo "Error: --mapping and --base-dir are required arguments"
    show_usage
    exit 1
fi

# Check requirements only if we need processing tools
if [[ "$PROCESSES" == *"PROCESS_"* ]]; then
    check_requirements
fi

# Verify mapping file exists
verify_file_exists "$MAPPING_FILE" "Mapping file" || exit 1

# Set default values for new arguments
if [ -z "$LR_CALLER" ]; then
    LR_CALLER="sniffles"
fi

if [ -z "$COMPARISON_DIRECTIONS" ]; then
    COMPARISON_DIRECTIONS="SR,LR"
fi

if [ -z "$DISTANCE" ]; then
    DISTANCE="50"
fi

if [ -z "$SR_VCF_PATH" ]; then
    SR_VCF_PATH="all_batches.filter_genotypes.sanitized.135.overlaps_tagged.vcf.gz"
fi

# Automatically set LR_VCF_PATH and OUTPUT_DIR based on LR_CALLER if not specified
if [ -z "$LR_VCF_PATH" ]; then
    LR_VCF_PATH="${LR_CALLER}_vcfs"
fi

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${LR_CALLER}_results"
fi

# Set default LR VCF pattern if not specified and using custom caller
if [ -z "$LR_VCF_PATTERN" ] && [[ "$LR_CALLER" != "sniffles" && "$LR_CALLER" != "cutesv" ]]; then
    LR_VCF_PATTERN="{LR_ID}.hg38.{LR_CALLER}.sv.vcf.gz"
fi

# Validate comparison directions
case "$COMPARISON_DIRECTIONS" in
    "SR"|"LR"|"SR,LR"|"LR,SR")
        ;;
    *)
        echo "Error: Invalid comparison directions '$COMPARISON_DIRECTIONS'. Use SR, LR, or SR,LR"
        exit 1
        ;;
esac

echo ""
echo "========================================================"
echo "LR BENCHMARKING PIPELINE"
echo "========================================================"
echo "Mapping file: $MAPPING_FILE"
echo "Base directory: $BASE_DIR"
echo "Number of samples to process: $NUM_SAMPLES"
echo "Processes to run: $PROCESSES"
echo "LR caller: $LR_CALLER"
echo "SR VCF path: $SR_VCF_PATH"
echo "LR VCF path: $LR_VCF_PATH"
if [ -n "$LR_VCF_PATTERN" ]; then
    echo "LR VCF pattern: $LR_VCF_PATTERN"
fi
echo "Output directory: $OUTPUT_DIR"
echo "Comparison directions: $COMPARISON_DIRECTIONS"
echo "Distance parameter: $DISTANCE"
echo "SR limit: ${SR_LIMIT:-none}"
echo "LR limit: ${LR_LIMIT:-none}"
echo ""

# Process samples
successful_samples=0
failed_samples=0
total_processed=0

# Read mapping file and process samples
while IFS=$'\t' read -r sr_id lr_id || [ -n "$sr_id" ]; do
    if [ "$total_processed" -ge "$NUM_SAMPLES" ]; then
        break
    fi
    
    if [ "$sr_id" = "SR_ID" ]; then
        continue  # Skip header
    fi
    
    ((total_processed++))
    echo ""
    echo "SAMPLE $total_processed/$NUM_SAMPLES"
    
    if process_sample "$sr_id" "$lr_id" "$BASE_DIR" "$PROCESSES"; then
        ((successful_samples++))
        echo ""
        echo "✓ SUCCESS: Sample $sr_id processed successfully"
    else
        ((failed_samples++))
        echo ""
        echo "✗ FAILED: Sample $sr_id processing failed"
    fi
done < "$MAPPING_FILE"

echo ""
echo "========================================================"
echo "PROCESSING COMPLETE"
echo "========================================================"
echo "Successfully processed: $successful_samples samples"
echo "Failed to process: $failed_samples samples"
echo "========================================================" 