#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Contact: Xuefang Zhao <xzhao@broadinstitute.org> 

"""
Annotate sample-level Boost scores directly into an input VCF
"""

import csv
import numpy as np
import argparse
import pysam
import sys

def SVID_info_readin(SVID_info):
	fin=open(SVID_info)
	out = {}
	for line in fin:
		pin=line.strip().split()
		out[pin[0]] = pin[1]
	fin.close()
	return(out)

def samp_batch_readin(sample_batch):
	fin=open(sample_batch)
	out = {}
	for line in fin:
		pin=line.strip().split()
		if not pin[1] in out.keys():
			out[pin[1]] = []
		out[pin[1]].append(pin[0])
	fin.close()
	return(out)

def main():
	"""
	Command-line main block
	"""

	# Parse command line arguments and options
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--vcf', required=True, help='VCF to be annotated. Required.')
	parser.add_argument('--SVID_info', required=True, help='two-column .tsv of ' +'file including the SVID that failed batch effect' +'Required.')
	parser.add_argument('--SVID_filter', required=True, help='two-column .tsv of ' +'file including the SVID to have filter column changed' +'Required.')
	parser.add_argument('--SVID_format', required=True, help='two-column .tsv of ' +'file including the SVID to have format column changed' +'Required.')
	parser.add_argument('--sample_batch', required=True, help='two-column .tsv of ' +'file including the sample ID and corresponding batches' +'Required.')
	parser.add_argument('-o', '--outfile', default='stdout', help='Path to output ' + 'VCF. Default: stdout.')
	args = parser.parse_args()


	#readin batch effect results
	SVID_info = SVID_info_readin(args.SVID_info)
	SVID_filter = SVID_info_readin(args.SVID_filter)
	SVID_format = SVID_info_readin(args.SVID_format)
	samp_batch = samp_batch_readin(args.sample_batch)

	#Open connection to input VCF
	if args.vcf in '- stdin'.split():
		vcf = pysam.VariantFile(sys.stdin)
	else:
		vcf = pysam.VariantFile(args.vcf)

	# Add new INFO line to VCF header
	vcf.header.add_meta('INFO',
						items=[('ID', "BATCH_EFFECT"), ('Number', "1"), ('Type', "String"),
								('Description', "Batch effect results")])

	vcf.header.add_meta('FILTER',items=[('ID', "PCR_AF_BIAS"), 
								('Description', "SVs showed significant bias between PCR+ and PCR- batches")])

	vcf.header.add_meta('FILTER',items=[('ID', "UNSTABLE_AF_ESTIMATE_PCRMINUS"), 
								('Description', "SVs showed significant bias between pairs of PCR- batches")])


	fout = pysam.VariantFile(args.outfile, 'w', header = vcf.header)
	for record in vcf:
		if record.id in SVID_info.keys():
				record.info['BATCH_EFFECT'] = SVID_info[record.id]
		if record.id in SVID_filter.keys():
				record.filter.add(SVID_filter[record.id])
		if record.id in SVID_format.keys():
				batch_key = SVID_format[record.id].split(',')
				sample_list = []
				for batch in batch_key:
					sample_list += samp_batch['gnomAD-SV_v3_' + batch]
				for sample in sample_list:
					if sample in record.samples.keys():
						record.samples[sample]['GT'] = (None, None)
		fout.write(record)
	vcf.close()
	fout.close()

if __name__ == '__main__':
	main()


