import pysam
import re
import argparse
import os

def extract_end2_from_alt(alt):
    match = re.search(r"[\[\]]chr[0-9XYM]+:([0-9]+)[\[\]]", alt)
    return int(match.group(1)) if match else None

def process_vcf(input_vcf, output_vcf):
    # Open input VCF
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Ensure INFO field includes END2
    if "END2" not in vcf_in.header.info:
        vcf_in.header.add_meta("INFO", items=[
            ("ID", "END2"), 
            ("Number", "1"), 
            ("Type", "Integer"), 
            ("Description", "Extracted END2 from ALT field")
        ])

    # Open output VCF
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    # Process each variant record
    for record in vcf_in:
        if record.info.get("SVTYPE") == "BND" and record.info.get("CHR2") != record.chrom:
            end2_value = extract_end2_from_alt(str(record.alts[0]))
            if end2_value:
                record.info["END2"] = end2_value  # Add END2 field
        vcf_out.write(record)

    # Close VCFs
    vcf_in.close()
    vcf_out.close()

    # Index the output VCF if it's compressed (.vcf.gz)
    if output_vcf.endswith(".gz"):
        pysam.tabix_index(output_vcf, preset="vcf")

def main():
    parser = argparse.ArgumentParser(description="Extract END2 from ALT for BND variants in a VCF file.")
    parser.add_argument("-v", "--vcf", required=True, help="Path to input VCF (.vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True, help="Path to output VCF (.vcf or .vcf.gz)")

    args = parser.parse_args()

    # Ensure input file exists
    if not os.path.exists(args.vcf):
        print(f"Error: Input file '{args.vcf}' not found.")
        exit(1)

    # Run the processing function
    process_vcf(args.vcf, args.output)
    print(f"Updated VCF written to: {args.output}")

if __name__ == "__main__":
    main()
