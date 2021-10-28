#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
"""

import argparse
import pysam


def merge_annotations(vcf, protein_coding_vcf, lincRNA_vcf, promoter_vcf, noncoding_vcf):

    record_iter = zip(vcf, protein_coding_vcf, lincRNA_vcf,
                      promoter_vcf, noncoding_vcf)

    gene_classes = ['PROTEIN_CODING', 'LINCRNA']

    coding_classes = 'LOF DUP_LOF COPY_GAIN MSV_EXON_OVR DUP_PARTIAL INTRONIC INV_SPAN UTR'.split()

    for records in record_iter:
        # First check that all records have same ID
        ids = [r.id for r in records]
        if len(set(ids)) != 1:
            msg = 'Mismatched records: {0}'.format(','.join(set(ids)))
            raise Exception(msg)

        # Use record from unannotated VCF with updated header as base
        base = records[0]

        # Copy coding annotations over from each gene class
        for gene_class, record in zip(gene_classes, records[1:6]):
            for coding_class in coding_classes:
                if coding_class in record.info.keys():
                    info = '{0}__{1}'.format(gene_class, coding_class)
                    base.info[info] = record.info[coding_class]

            # Additionally tag nearest TSS and intergenic for protein coding
            if gene_class == 'PROTEIN_CODING':
                if 'INTERGENIC' in record.info.keys():
                    base.info['PROTEIN_CODING__NEAREST_TSS'] = record.info['NEAREST_TSS']
                    base.info['PROTEIN_CODING__INTERGENIC'] = True

        # Annotate promoters of any genes if there wasn't a coding hit
        promoter_record = records[3]
        promoters = promoter_record.info.get(
            'NONCODING_SPAN', ()) + promoter_record.info.get('NONCODING_BREAKPOINT', ())
        kept_promoters = []
        for gene in promoters:
            coding_hit = False
            for coding_class in coding_classes:
                info = 'PROTEIN_CODING__{0}'.format(coding_class)
                if info in base.info.keys() and gene in base.info[info]:
                    coding_hit = True
            if not coding_hit:
                kept_promoters.append(gene)

        if len(kept_promoters) > 0:
            base.info['PROTEIN_CODING__PROMOTER'] = kept_promoters

        noncoding_record = records[4]
        if 'NONCODING_SPAN' in noncoding_record.info.keys():
            base.info['NONCODING_SPAN'] = noncoding_record.info['NONCODING_SPAN']
        if 'NONCODING_BREAKPOINT' in noncoding_record.info.keys():
            base.info['NONCODING_BREAKPOINT'] = noncoding_record.info['NONCODING_BREAKPOINT']

        yield base


def update_header(vcf):
    GENCODE_INFO = [
        '##INFO=<ID={0}__LOF,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a loss-of-function effect.">',
        '##INFO=<ID={0}__DUP_LOF,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a loss-of-function effect via intragenic exonic duplication.">',
        '##INFO=<ID={0}__COPY_GAIN,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a copy-gain effect.">',
        '##INFO=<ID={0}__DUP_PARTIAL,Number=.,Type=String,Description="Gene(s) which are partially overlapped by an SV\'s duplication, such that an unaltered copy is preserved.">',
        '##INFO=<ID={0}__MSV_EXON_OVR,Number=.,Type=String,Description="Gene(s) on which the multiallelic SV would be predicted to have a LOF, DUP_LOF, COPY_GAIN, or DUP_PARTIAL annotation if the SV were biallelic.">',
        '##INFO=<ID={0}__INTRONIC,Number=.,Type=String,Description="Gene(s) where the SV was found to lie entirely within an intron.">',
        '##INFO=<ID={0}__INV_SPAN,Number=.,Type=String,Description="Gene(s) which are entirely spanned by an SV\'s inversion.">',
        '##INFO=<ID={0}__UTR,Number=.,Type=String,Description="Gene(s) for which the SV is predicted to disrupt a UTR.">'
    ]

    OTHER_INFO = [
        '##INFO=<ID=NONCODING_SPAN,Number=.,Type=String,Description="Classes of noncoding elements spanned by SV.">',
        '##INFO=<ID=NONCODING_BREAKPOINT,Number=.,Type=String,Description="Classes of noncoding elements disrupted by SV breakpoint.">',
        '##INFO=<ID=PROTEIN_CODING__NEAREST_TSS,Number=.,Type=String,Description="Nearest transcription start site to intragenic variants.">',
        '##INFO=<ID=PROTEIN_CODING__INTERGENIC,Number=0,Type=Flag,Description="SV does not overlap coding sequence.">',
        '##INFO=<ID=PROTEIN_CODING__PROMOTER,Number=.,Type=String,Description="Genes whose promoter sequence (1 kb) was disrupted by SV.">'
    ]

    classes = ['PROTEIN_CODING', 'LINCRNA']

    for info in GENCODE_INFO:
        for c in classes:
            vcf.header.add_line(info.format(c))

    for info in OTHER_INFO:
        vcf.header.add_line(info)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('protein_coding_vcf')
    # parser.add_argument('antisense_vcf')
    parser.add_argument('lincRNA_vcf')
    # parser.add_argument('processed_transcript_vcf')
    # parser.add_argument('pseudogene_vcf')
    parser.add_argument('promoter_vcf')
    parser.add_argument('noncoding_vcf')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    protein_coding_vcf = pysam.VariantFile(args.protein_coding_vcf)
    # antisense_vcf = pysam.VariantFile(args.antisense_vcf)
    lincRNA_vcf = pysam.VariantFile(args.lincRNA_vcf)
    # processed_transcript_vcf = pysam.VariantFile(args.processed_transcript_vcf)
    # pseudogene_vcf = pysam.VariantFile(args.pseudogene_vcf)
    promoter_vcf = pysam.VariantFile(args.promoter_vcf)
    noncoding_vcf = pysam.VariantFile(args.noncoding_vcf)

    update_header(vcf)
    fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    for record in merge_annotations(vcf, protein_coding_vcf, lincRNA_vcf,
                                    promoter_vcf, noncoding_vcf):
        fout.write(record)


if __name__ == '__main__':
    main()
