#!python
#script to extract subset of SVs that are of specific type and fell within certain genomic regions from vcf
# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import os
import pysam
import argparse

def SVID_GC_readin(SVID_GC,GenomicContext):
	fin=open(SVID_GC)
	out=[]
	for line in fin:
		pin=line.strip().split()
		if pin[1] == GenomicContext:
			out.append(pin[0])
	fin.close()
	return(out)

def vcf_split(vcf_file, out_file, SVID_list, svtype):
	vcf = pysam.VariantFile(vcf_file)
	header  = vcf.header
	header.add_line('##source=cleanvcf')
	out = pysam.VariantFile(out_file, 'w', header = header)
	for record in vcf:
		if record.id in SVID_list and record.info['SVTYPE'] == svtype:
			record.info["STRANDS"] = "+-"
			out.write(record)
	vcf.close()
	out.close()

def main():
	parser = argparse.ArgumentParser(description='Split vcf by genomic context and sv type')
	parser.add_argument('vcf', metavar='', type=str,help='name of vcf file to process')
	parser.add_argument('output', metavar='', type=str,help='name of output file')
	parser.add_argument('svid_gc', metavar='', type=str,help='2 column sample with SVID and genomic context')
	parser.add_argument('svtype', metavar='', type=str,help='type of SVs to retain')
	parser.add_argument('genomic_context', metavar='', type=str,help='genomic context to retain')

	args = parser.parse_args()
	SVID_list = SVID_GC_readin(args.svid_gc, args.genomic_context)
	vcf_split(args.vcf, args.output, SVID_list, args.svtype)


if __name__ == '__main__':
	main()

