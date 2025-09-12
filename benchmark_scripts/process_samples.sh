#!/bin/bash

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
    local input_vcf="${base_dir}/annotated.manual_review.vcf.gz"
    
    local ref_dir="../../gatk-sv/module18_annotate_genomic_context/references"
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
    bcftools view -i 'FILTER="PASS" && GT!="0/0" && GT!="./."' -Oz -o "$sr_filtered_vcf" "$sr_vcf" || return 1
    tabix "$sr_filtered_vcf" || return 1
    
    echo "Annotating SV types..."
    "./annotate_svtype.sh" --in "$sr_filtered_vcf" --out "$sr_filtered_revised_vcf" || return 1
    
    echo "Converting to BED format..."
    svtk vcf2bed "$sr_filtered_revised_vcf" "$sr_bed" --info SVLEN --info SVTYPE --info AF --no-samples || return 1
    
    echo "Extracting filter information..."
    /Users/kjaising/miniconda3/envs/venv/bin/python "./extract_filter_info.py" --vcf "$sr_filtered_revised_vcf" --output "${sample_dir}/${sr_id}.sr.filter_info.tsv" || return 1
    
    echo "Annotating genomic context..."
    "${ref_dir}/../annotate_genomic_context.sh" "$sr_bed" "$sr_annotated_bed" --rm "$rm_ref" --sd "$sd_ref" --sr "$sr_ref" || return 1
    
    return 0
}

process_lr_sample() {
    local sr_id="$1"
    local lr_id="$2"
    local base_dir="$3"
    local lr_caller="$4"
    
    echo ""
    echo "--- LONG-READ PROCESSING ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    mkdir -p "$sample_dir"

    local lr_source_file="${LR_VCF_PATH}/${lr_id}.hg38.${lr_caller}.sv.vcf.gz"    
    local lr_vcf="${sample_dir}/${sr_id}.lr.vcf.gz"
    local lr_filtered_vcf="${sample_dir}/${sr_id}.lr.filtered.vcf.gz"
    local lr_bed="${sample_dir}/${sr_id}.lr.bed"
    local lr_annotated_bed="${sample_dir}/${sr_id}.lr.annotated.bed"
    
    local ref_dir="../../gatk-sv/module18_annotate_genomic_context/references"
    local rm_ref="${ref_dir}/hg38.RM.sorted.merged.bed.gz"
    local sd_ref="${ref_dir}/hg38.SD.sorted.merged.bed.gz"
    local sr_ref="${ref_dir}/hg38.SR.sorted.merged.bed.gz"
    
    verify_file_exists "$rm_ref" "Reference file (RM)" || return 1
    verify_file_exists "$sd_ref" "Reference file (SD)" || return 1
    verify_file_exists "$sr_ref" "Reference file (SR)" || return 1
    verify_file_exists "$lr_source_file" "Long-read VCF file" || return 1
    
    echo "Processing long-read data for $lr_id..."
    
    echo "Filtering long-read variants..."
    bcftools view -i 'GT!="0/0" && GT!="./."' -Oz -o "$lr_vcf" "$lr_source_file" || return 1
    tabix "$lr_vcf" || return 1
    
    echo "Annotating SV types..."
    "./annotate_svtype.sh" --in "$lr_vcf" --out "$lr_filtered_vcf" || return 1
    
    echo "Converting to BED format..."
    svtk vcf2bed "$lr_filtered_vcf" "$lr_bed" --info SVLEN --info SVTYPE --no-samples || return 1
    
    echo "Extracting filter information..."
    /Users/kjaising/miniconda3/envs/venv/bin/python "./extract_filter_info.py" --vcf "$lr_filtered_vcf" --output "${sample_dir}/${sr_id}.lr.filter_info.tsv" || return 1
    
    echo "Annotating genomic context..."
    "${ref_dir}/../annotate_genomic_context.sh" "$lr_bed" "$lr_annotated_bed" --rm "$rm_ref" --sd "$sd_ref" --sr "$sr_ref" || return 1
    
    return 0
}

compare_callsets_sample() {
    local sr_id="$1"
    local lr_id="$2" 
    local base_dir="$3"
    
    echo ""
    echo "--- CALLSET COMPARISONS ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    
    verify_file_exists "${sample_dir}/${sr_id}.sr.bed" "SR BED file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.bed" "LR BED file" || return 1
    
    echo "Processing BED files for comparison..."
    
    cmd="python \"./process_bed_files.py\" --sample \"$sr_id\" --sr \"${sample_dir}/${sr_id}.sr.bed\" --lr \"${sample_dir}/${sr_id}.lr.bed\""
    
    if [ -n "$SR_LIMIT" ]; then
        cmd="$cmd --sr-limit \"$SR_LIMIT\""
    fi
    
    if [ -n "$LR_LIMIT" ]; then
        cmd="$cmd --lr-limit \"$LR_LIMIT\""
    fi
    
    eval "$cmd" || return 1
    
    verify_file_exists "${sample_dir}/${sr_id}.sr.query.bed" "SR query BED" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.sr.ref.bed" "SR reference BED" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.query.bed" "LR query BED" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.ref.bed" "LR reference BED" || return 1
    
    echo "Running callset comparisons..."
    
    if [[ "$COMPARISON_DIRECTIONS" == *"SR"* ]]; then
        echo "Running LR query vs SR reference comparison..."
        "../../gatk-sv/benchmark_scripts/compare_callsets_V2.sh" -O "${sample_dir}/${sr_id}.lr_query.bed" -p lr_query "${sample_dir}/${sr_id}.lr.query.bed" "${sample_dir}/${sr_id}.sr.ref.bed" || return 1
    fi
    
    if [[ "$COMPARISON_DIRECTIONS" == *"LR"* ]]; then
        echo "Running SR query vs LR reference comparison..."
        "../../gatk-sv/benchmark_scripts/compare_callsets_V2.sh" -O "${sample_dir}/${sr_id}.sr_query.bed" -p sr_query "${sample_dir}/${sr_id}.sr.query.bed" "${sample_dir}/${sr_id}.lr.ref.bed" || return 1
    fi
    
    return 0
}

analyze_rates_sample() {
    local sr_id="$1"
    local lr_id="$2" 
    local base_dir="$3"
    
    echo ""
    echo "--- MATCH RATE ANALYSIS ---"
    
    local sample_dir="${OUTPUT_DIR}/${sr_id}"
    
    verify_file_exists "${sample_dir}/${sr_id}.sr.annotated.bed" "SR annotated file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.annotated.bed" "LR annotated file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.sr.filter_info.tsv" "SR filter info file" || return 1
    verify_file_exists "${sample_dir}/${sr_id}.lr.filter_info.tsv" "LR filter info file" || return 1
    
    local analysis_run=false
    
    if [[ "$COMPARISON_DIRECTIONS" == *"SR"* ]] && [[ -f "${sample_dir}/${sr_id}.lr_query.bed" ]]; then
        echo "Analyzing LR query match rates..."
        python "./analyze_match_rates.py" \
            --query-bed "${sample_dir}/${sr_id}.lr_query.bed" \
            --annotated-bed "${sample_dir}/${sr_id}.sr.annotated.bed" \
            --filter-info "${sample_dir}/${sr_id}.sr.filter_info.tsv" \
            --output-prefix "${sample_dir}/lr_query" \
            --query-type lr_query || return 1
        analysis_run=true
    fi
    
    if [[ "$COMPARISON_DIRECTIONS" == *"LR"* ]] && [[ -f "${sample_dir}/${sr_id}.sr_query.bed" ]]; then
        echo "Analyzing SR query match rates..."
        python "./analyze_match_rates.py" \
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

# Main script
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 --mapping MAPPING_FILE --base-dir BASE_DIR [--num-samples NUM_SAMPLES] [--processes PROCESS_LIST] [--lr-caller LR_CALLER] [--comparison-directions COMPARISON_DIRECTIONS] [--sr-limit SR_LIMIT] [--lr-limit LR_LIMIT]"
    echo ""
    echo "Required arguments:"
    echo "  --mapping MAPPING_FILE          - Path to mapping file with SR_ID and LR_ID columns"
    echo "  --base-dir BASE_DIR             - Base directory for input data"
    echo ""
    echo "Optional arguments:"
    echo "  --num-samples NUM_SAMPLES       - Number of samples to process (default: 3)"
    echo "  --lr-caller LR_CALLER           - Long-read caller: sniffles, cutesv, svim, etc. (default: sniffles)"
    echo "  --comparison-directions DIRS    - Comparison directions: SR, LR, or SR,LR (default: SR,LR)"
    echo "  --sr-limit SR_LIMIT            - SR limit for processing (default: no limit)"
    echo "  --lr-limit LR_LIMIT            - LR limit for processing (default: no limit)"
    echo ""
    echo "Process types (comma-separated):"
    echo "  PROCESS_SR    - Process short-read data (VCF filtering, annotation, BED conversion)"
    echo "  PROCESS_LR    - Process long-read data (VCF filtering, annotation, BED conversion)"
    echo "  COMPARE_CALLSETS - Run callset comparison"
    echo "  ANALYZE_RATES - Run match rate analysis (requires comparison results)"
    echo ""
    echo "Note: Output directory and LR VCF path are automatically determined from lr-caller:"
    echo "  For --lr-caller cutesv: output-dir=cutesv_results, vcf-path=cutesv_vcfs"
    echo "  For --lr-caller sniffles: output-dir=sniffles_results, vcf-path=sniffles_vcfs"
    echo ""
    echo "Example: $0 --mapping Mapping.tsv --base-dir . --processes ANALYZE_RATES --lr-caller cutesv"
    exit 1
fi

MAPPING_FILE=""
BASE_DIR=""
NUM_SAMPLES=3
PROCESSES="PROCESS_SR,PROCESS_LR,COMPARE_CALLSETS,ANALYZE_RATES"
LR_CALLER=""
COMPARISON_DIRECTIONS=""
SR_LIMIT=""
LR_LIMIT=""

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
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ "$PROCESSES" == *"PROCESS_"* ]]; then
    check_requirements
fi

verify_file_exists "$MAPPING_FILE" "Mapping file" || exit 1

if [ -z "$LR_CALLER" ]; then
    LR_CALLER="sniffles"
fi

if [ -z "$COMPARISON_DIRECTIONS" ]; then
    COMPARISON_DIRECTIONS="SR,LR"
fi

LR_VCF_PATH="${BASE_DIR}/${LR_CALLER}_vcfs"
OUTPUT_DIR="${BASE_DIR}/${LR_CALLER}_results"

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
echo "LR VCF path: $LR_VCF_PATH"
echo "Output directory: $OUTPUT_DIR"
echo "Comparison directions: $COMPARISON_DIRECTIONS"
echo "SR limit: $SR_LIMIT"
echo "LR limit: $LR_LIMIT"
echo ""

successful_samples=0
failed_samples=0
total_processed=0
while IFS=$'\t' read -r sr_id lr_id || [ -n "$sr_id" ]; do
    if [ "$total_processed" -ge "$NUM_SAMPLES" ]; then
        break
    fi
    
    if [ "$sr_id" = "SR_ID" ]; then
        continue
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