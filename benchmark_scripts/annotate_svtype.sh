#!/bin/bash

# Script to replace ALT field with SVTYPE from INFO field
# Usage: ./annotate_svtype.sh --in input.vcf.gz --out output.vcf.gz

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --in)
            INPUT_VCF="$2"
            shift 2
            ;;
        --out)
            OUTPUT_VCF="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --in input.vcf.gz --out output.vcf.gz"
            echo "  --in   Input VCF file (can be compressed)"
            echo "  --out  Output VCF file (will be compressed)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$INPUT_VCF" || -z "$OUTPUT_VCF" ]]; then
    echo "Error: Both --in and --out arguments are required"
    echo "Usage: $0 --in input.vcf.gz --out output.vcf.gz"
    exit 1
fi

# Check if input file exists
if [[ ! -f "$INPUT_VCF" ]]; then
    echo "Error: Input file '$INPUT_VCF' not found"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
mkdir -p "$OUTPUT_DIR"

# Create temporary directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Annotating SV types in VCF: $(basename "$INPUT_VCF")"

# Process the VCF file in one step
{
    # Extract and output header
    bcftools view -h "$INPUT_VCF"
    
    # Process records: replace ALT field with SVTYPE from INFO
    bcftools view -H "$INPUT_VCF" | awk 'BEGIN{FS=OFS="\t"} {
        # Parse INFO field to extract SVTYPE
        n = split($8, info_fields, ";");
        svtype = "";
        for(i=1; i<=n; i++){
            if(info_fields[i] ~ /^SVTYPE=/){
                svtype = substr(info_fields[i], 8);
                break;
            }
        }
        
        # If SVTYPE not found, keep original ALT or use "."
        if(svtype == ""){
            svtype = ($5 == "" || $5 == ".") ? "." : $5;
        }
        
        # Replace ALT field (column 5) with SVTYPE
        $5 = svtype;
        print;
    }'
} | bgzip -c > "$OUTPUT_VCF"

# Index the output VCF
tabix -p vcf "$OUTPUT_VCF"

echo "SV type annotation completed: $(basename "$OUTPUT_VCF")"