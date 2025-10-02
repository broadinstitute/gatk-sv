#!/usr/bin/env python3

import argparse
import pysam
import os

def reclassify_cnv_genotypes(input_vcf, output_vcf):
    vcf_in = pysam.VariantFile(input_vcf)
    
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
    
    for record in vcf_in:
        if 'CNV' not in record.id:
            vcf_out.write(record)
            continue
        
        svtype = record.info.get('SVTYPE')
        if svtype not in ['DEL', 'DUP']:
            vcf_out.write(record)
            continue
        
        for sample in record.samples:
            sample_data = record.samples[sample]
            
            if 'ECN' not in sample_data or 'RD_CN' not in sample_data:
                continue
            
            ecn = sample_data['ECN']
            rd_cn = sample_data['RD_CN']
            
            if isinstance(ecn, tuple):
                ecn = ecn[0]
            if isinstance(rd_cn, tuple):
                rd_cn = rd_cn[0]

            if ecn is None or rd_cn is None:
                continue
            
            if svtype == 'DEL':
                if ecn <= rd_cn:
                    sample_data['GT'] = (0, 0)
                elif rd_cn == ecn - 1:
                    sample_data['GT'] = (0, 1)
                elif rd_cn < ecn - 1:
                    sample_data['GT'] = (1, 1)
            
            elif svtype == 'DUP':
                if ecn >= rd_cn:
                    sample_data['GT'] = (0, 0)
                elif rd_cn == ecn + 1:
                    sample_data['GT'] = (0, 1)
                elif rd_cn > ecn + 1:
                    sample_data['GT'] = (1, 1)
        
        vcf_out.write(record)
    
    vcf_in.close()
    vcf_out.close()

def main():
    parser = argparse.ArgumentParser(description='Reclassify CNV genotypes based on ECN and RD_CN values')
    parser.add_argument('input_vcf', help='Input VCF file (.vcf or .vcf.gz)')
    parser.add_argument('output_vcf', help='Output VCF file (.vcf or .vcf.gz)')
    
    args = parser.parse_args()
    
    reclassify_cnv_genotypes(args.input_vcf, args.output_vcf)

if __name__ == '__main__':
    main()