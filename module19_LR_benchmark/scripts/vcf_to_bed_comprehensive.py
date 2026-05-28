#!/usr/bin/env python3
"""
Convert VCF to BED format with specific columns and INFO fields.

Columns included:
- chr, pos, end, id, filter (standard BED columns)
- Selected INFO fields: allele_length, allele_type, SOURCE, AF, AC, AN,
    gnomAD_V4_match_type, gnomAD_V4_match_ID, gnomAD_V4_match_source,
    dbGaP_ID, MERGE_TYPE, TRID, REGION
- All INFO fields containing "PREDICTED_"
- VEP fields: 2nd column as vep and 4th column as vep_gene
"""

import sys
import gzip
import argparse
from typing import Dict, Tuple


def open_vcf(path: str):
    """Open VCF file (handles .gz compression)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info_field(info_str: str) -> Dict[str, str]:
    """Parse INFO field from VCF into a dictionary."""
    info_dict = {}
    if info_str == ".":
        return info_dict
    
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def extract_vep_columns(vep_str: str, col2_idx: int = 1, col4_idx: int = 3) -> Tuple[str, str]:
    """
    Extract specific columns from VEP field (pipe-delimited).
    
    Args:
        vep_str: VEP field value (comma-separated entries, each pipe-delimited)
        col2_idx: Index for 2nd column (0-based, default 1 for column 2)
        col4_idx: Index for 4th column (0-based, default 3 for column 4)
    
    Returns:
        Two comma-separated strings containing the selected VEP columns.
    """
    if not vep_str or vep_str == ".":
        return ".", "."

    vep_values = []
    vep_gene_values = []
    # VEP can have multiple entries separated by commas
    entries = vep_str.split(",")

    for entry in entries:
        parts = entry.split("|")

        # Extract 2nd column (index 1)
        if len(parts) > col2_idx:
            vep_values.append(parts[col2_idx])
        else:
            vep_values.append(".")

        # Extract 4th column (index 3)
        if len(parts) > col4_idx:
            vep_gene_values.append(parts[col4_idx])
        else:
            vep_gene_values.append(".")

    return ",".join(vep_values), ",".join(vep_gene_values)


def vcf_to_bed(vcf_file: str, bed_file: str, vep_col2: int = 2, 
               vep_col4: int = 4, skip_missing: bool = False) -> None:
    """
    Convert VCF to BED format.
    
    Args:
        vcf_file: Input VCF file path
        bed_file: Output BED file path
        vep_col2: Column number to extract from VEP (2nd column default, 1-indexed)
        vep_col4: Column number to extract from VEP (4th column default, 1-indexed)
        skip_missing: If True, skip lines where critical fields are missing
    """
    # Convert to 0-based indices
    vep_col2_idx = vep_col2 - 1
    vep_col4_idx = vep_col4 - 1
    
    required_info_fields = [
        "allele_length",
        "allele_type",
        "SOURCE",
        "AF",
        "AC",
        "AN",
        "gnomAD_V4_match_type",
        "gnomAD_V4_match_ID",
        "gnomAD_V4_match_source",
        "dbGaP_ID",
        "MERGE_TYPE",
        "TRID",
        "REGION"
    ]
    
    # Header for output BED file
    header_cols = [
        "chr", "pos", "end", "id", "filter",
        *required_info_fields,
        "PREDICTED_fields",
        "vep",
        "vep_gene"
    ]
    
    try:
        with open_vcf(vcf_file) as infile, open(bed_file, 'w') as outfile:
            # Write header
            outfile.write("\t".join(header_cols) + "\n")
            
            line_num = 0
            variants_written = 0
            
            for line in infile:
                line = line.rstrip("\n")
                
                # Skip empty lines and comments (except header)
                if not line or line.startswith("##"):
                    continue
                
                # Skip header line
                if line.startswith("#CHROM"):
                    continue
                
                line_num += 1
                
                # Parse VCF line
                fields = line.split("\t")
                
                if len(fields) < 8:
                    print(f"Warning: Line {line_num} has fewer than 8 fields, skipping.", 
                          file=sys.stderr)
                    continue
                
                chrom = fields[0]
                pos = fields[1]
                vcf_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filt = fields[6]
                info_str = fields[7]
                
                # Parse INFO field
                info_dict = parse_info_field(info_str)
                
                # Determine END position
                # Try to get from INFO[END], otherwise use pos + ref length
                if "END" in info_dict:
                    end_pos = info_dict["END"]
                else:
                    try:
                        end_pos = str(int(pos) + len(ref) - 1)
                    except (ValueError, TypeError):
                        end_pos = pos
                
                # Extract required INFO fields
                info_values = []
                for field in required_info_fields:
                    info_values.append(info_dict.get(field, "."))
                
                # Extract PREDICTED_ fields
                predicted_fields = []
                for key in sorted(info_dict.keys()):
                    if "PREDICTED_" in key:
                        predicted_fields.append(f"{key}={info_dict[key]}")
                predicted_str = ";".join(predicted_fields) if predicted_fields else "."
                
                # Extract VEP field
                vep_str = info_dict.get("VEP", ".")
                vep_value, vep_gene_value = extract_vep_columns(vep_str, vep_col2_idx, vep_col4_idx)
                
                # Build output line
                output_fields = [
                    chrom,
                    pos,
                    end_pos,
                    vcf_id,
                    filt,
                    *info_values,
                    predicted_str,
                    vep_value,
                    vep_gene_value
                ]
                
                outfile.write("\t".join(output_fields) + "\n")
                variants_written += 1
        
        print(f"Successfully converted {variants_written} variants from {vcf_file} to {bed_file}", 
              file=sys.stderr)
    
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF to BED format with specific INFO fields and VEP annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.vcf output.bed
  %(prog)s input.vcf.gz output.bed
  %(prog)s input.vcf output.bed --vep-col2 1 --vep-col4 3

VEP columns:
  Default extraction uses columns 2 and 4 from pipe-delimited VEP entries.
  Use --vep-col2 and --vep-col4 to specify different columns (1-indexed).
        """
    )
    
    parser.add_argument("vcf_file", help="Input VCF file (supports .gz compression)")
    parser.add_argument("bed_file", help="Output BED file")
    parser.add_argument("--vep-col2", type=int, default=2,
                        help="Column number for 2nd VEP extraction (1-indexed, default: 2)")
    parser.add_argument("--vep-col4", type=int, default=4,
                        help="Column number for 4th VEP extraction (1-indexed, default: 4)")
    
    args = parser.parse_args()
    
    vcf_to_bed(args.vcf_file, args.bed_file, args.vep_col2, args.vep_col4)


if __name__ == "__main__":
    main()
