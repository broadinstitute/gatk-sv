#!/bin/python

import sys
import pysam

GQ_MIN = 0
P_MIN = 0.9

vcf = pysam.VariantFile(sys.argv[1], 'r')
sys.stdout.write(str(vcf.header))

for record in vcf:
    gq_list = sorted([record.samples[s]['GQ'] for s in record.samples.keys()])
    gq = gq_list[int(len(gq_list) / 2)]
    has_p_min = ('PPE' in record.info and record.info['PPE'] > P_MIN) or \
                ('PSR1' in record.info and 'PSR2' in record.info and record.info['PSR1'] > P_MIN and record.info['PSR2'] > P_MIN)
    if gq > GQ_MIN and has_p_min:
        sys.stdout.write(str(record))
