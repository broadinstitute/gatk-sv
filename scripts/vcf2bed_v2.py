#!/bin/python

import sys
import pysam

MINGQ = 20
GQ_SUM = 2

vcf = pysam.VariantFile(sys.argv[1], 'r')

for record in vcf:
    if record.info['SVTYPE'] in ['DEL', 'INS', 'BND']:
        carriers = [s for s in record.samples.keys() if sum(record.samples[s]['GT']) == GQ_SUM and record.samples[s]['GQ'] > MINGQ]
    elif record.info['SVTYPE'] in ['DUP']:
        carriers = [s for s in record.samples.keys() if record.samples[s]['CN'] - record.samples[s]['NCN'] == GQ_SUM and record.samples[s]['GQ'] > MINGQ]
    print('\t'.join([str(x) for x in [record.chrom, record.start, record.stop, record.info['SVTYPE'], ",".join(carriers)]]))
