#!/bin/python

import argparse
import pysam
import gzip


EV = 'EV'
SVTYPE = 'SVTYPE'
ME = 'ME'
UNR = 'UNR'
FILTER_VCF_INFO_LINES = {
    'BND_DEPTH', 'BND_MATEID', 'SPLIT_READS', 'PAIRED_END_READS', 
    'CLUSTER_MEMBER_IDS', 'MULTIALLELIC', 'UNRESOLVED', 'VARGQ', 
    'EVENT', 'REVISED_EVENT', 'MULTI_CNV'
}
FILTER_VCF_TEXT_LINES = {
    'CIPOS', 'CIEND', 'RMSSTD', 'source', 'bcftools', 'GATKCommandLine', 'fileformat'
}

# TODO: Remove INFO fields in advance of script: 'MULTI_CNV', 'VARGQ', 'REVISED_EVENT'

def cleanse_header(header):
    new_header = pysam.VariantHeader()

    for line in header.records:
        include_line = True
        if line.type == 'INFO' and line.get('ID') in FILTER_VCF_INFO_LINES:
            include_line = False
        elif any(fv_line in str(line) for fv_line in FILTER_VCF_TEXT_LINES):
            include_line = False
        elif line.type == 'FORMAT' and line.get('ID') == EV:
            include_line = False
        elif line.type == 'ALT' and line.get('ID') == UNR:
            include_line = False
        if include_line:
            new_header.add_line(str(line))
    
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

    # Skip if variant has mobile element
    has_mobile_element = False
    if record.alts:
        for allele in record.alts:
            if allele.startswith('<') and allele.endswith('>'):
                symbol = allele[1:-1]
                if symbol == ME:
                    has_mobile_element = True
                    break
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

        alt_count = sum(1 for allele_index in gt if allele_index is not None and allele_index > 0)
        if alt_count == 1:
            genotype['GT'] = (0, 1)
        elif alt_count == 2:
            genotype['GT'] = (1, 1)

    return record


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='CleanVcf postprocessing.')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    args = parser.parse_args()

    # Open input VCF
    if args.input_vcf.endswith('.gz'):
        vcf_in = pysam.VariantFile(gzip.open(args.input_vcf, 'rt'))
    else:
        vcf_in = pysam.VariantFile(args.input_vcf)
    new_header = cleanse_header(vcf_in.header)
    
    # Open output file
    if args.output_vcf.endswith('.gz'):
        vcf_out = pysam.VariantFile(args.output_vcf, 'wz', header=new_header)
    else:
        vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=new_header)
    
    # Process records
    for record in vcf_in:
        record = process_record(record)
        vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()
