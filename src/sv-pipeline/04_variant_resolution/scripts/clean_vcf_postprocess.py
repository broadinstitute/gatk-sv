#!/bin/python

import argparse
import pysam

# Constants
EV = 'EV'
SVTYPE = 'SVTYPE'
ME = 'ME'
UNR = 'UNR'
FILTER_VCF_INFO_LINES = {'BND_DEPTH', 'BND_MATEID', 'SPLIT_READS', 'PAIRED_END_READS', 'CLUSTER_MEMBER_IDS'}
FILTER_VCF_LINES = {'ID=UNR', 'ID=BND_DEPTH', 'ID=BND_MATEID', 'ID=CLUSTER_MEMBER_IDS', 'ID=PAIRED_END_READS', 'ID=SPLIT_READS'}

def modify_header(header):
    new_header = pysam.VariantHeader()

    # Copy over header lines, excluding some
    for line in header.records:
        include_line = True
        if line.type == 'INFO' and line.get('ID') in FILTER_VCF_INFO_LINES:
            include_line = False
        elif line.type == 'FORMAT' and line.get('ID') == EV:
            include_line = False
        elif line.type == 'ALT' and line.get('ID') == UNR:
            include_line = False
        elif any(fv_line in str(line) for fv_line in FILTER_VCF_LINES):
            include_line = False
        if include_line:
            new_header.add_line(str(line))

    # Add new header line for EV
    new_header.add_line('##FORMAT=<ID=EV,Number=1,Type=String,Description="Classes of evidence supporting final genotype">')

    # Add samples to header
    for sample in header.samples:
        new_header.add_sample(sample)

    return new_header

def process_record(record):
    record = cleanse_info_fields(record)
    record = process_svtype(record)
    return record

def cleanse_info_fields(record):
    for field in FILTER_VCF_INFO_LINES:
        if field in record.info:
            del record.info[field]
    return record

def process_svtype(record):
    svtype = record.info.get(SVTYPE, None)

    # Check for mobile element in alleles
    has_mobile_element = False
    if record.alts:
        for allele in record.alts:
            if allele.startswith('<') and allele.endswith('>'):
                symbol = allele[1:-1]
                if symbol == ME:
                    has_mobile_element = True
                    break

    # If SVTYPE is missing or variant has mobile element, skip processing
    if svtype is None or has_mobile_element:
        return record

    # Update alleles
    ref_allele = record.ref
    alt_allele = f'<{svtype}>'
    record.alleles = (ref_allele, alt_allele)

    # Update genotypes
    for sample in record.samples:
        genotype = record.samples[sample]
        gt = genotype.get('GT', (None, None))

        # Count number of alt alleles
        alt_count = sum(1 for allele_index in gt if allele_index is not None and allele_index > 0)

        # Update GT accordingly
        if alt_count == 1:
            genotype['GT'] = (0, 1)
        elif alt_count == 2:
            genotype['GT'] = (1, 1)

    return record

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Clean VCF post-processing.')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF name')
    parser.add_argument('input_vcf', help='Input VCF file')
    args = parser.parse_args()

    # Open input VCF
    vcf_in = pysam.VariantFile(args.input_vcf)

    # Modify header
    new_header = modify_header(vcf_in.header)

    # Open output VCF
    vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=new_header)

    # Process and write variants
    for record in vcf_in:
        record = process_record(record)
        vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()
