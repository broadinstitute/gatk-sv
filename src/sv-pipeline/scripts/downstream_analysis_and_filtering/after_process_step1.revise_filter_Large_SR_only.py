#!modify filter label of SR only PASS SVs:

import argparse
import pysam


def revise_filter_Large_SR_only(vcf, fout, size_cutoff=500, sample_proportion_cutoff=0.1):
    records = []
    for record in vcf:
        if 'PASS' in record.filter.keys():
            if record.info['SVLEN'] > size_cutoff:
                SR_only_list = []
                gt_list = []
                for sample in record.samples.values():
                    if sample['RD_CN'] == 2 and sample['PE_GT'] == 0 and sample['PE_GQ'] == 99 and sample['SR_GT'] != 0:
                        SR_only_list.append(1)
                        gt_list.append(sample['GT'])
                    else:
                        SR_only_list.append(0)
                        gt_list.append(sample['GT'])
                if sum(SR_only_list) / len(SR_only_list) > sample_proportion_cutoff:
                    record.filter.clear()
                    record.filter.add('LARGE_SR_ONLY')
                    # print([str(j) for j in [record.chrom, record.pos,record.id, sum(SR_only_list)/len(SR_only_list)]])
        records.append(record)
    for record in records:
        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    # parser.add_argument('PCRPLUS_samples', help='List of PCRPLUS sample IDs.')
    # parser.add_argument('male_samples', help='List of male sample IDs.')
    parser.add_argument('fout', help='Output file (supports "stdout").')
    args = parser.parse_args()

    # Add new FILTER lines to VCF header
    vcf = pysam.VariantFile(args.vcf)

    NEW_FILTERS = ['##FILTER=<ID=LARGE_SR_ONLY,Description="Site genotyped as non-reference ' +
                   'but has only support from SR evidences.">']

    header = vcf.header
    for filt in NEW_FILTERS:
        header.add_line(filt)

    fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)
    revise_filter_Large_SR_only(vcf, fout)
    fout.close()


if __name__ == '__main__':
    main()
