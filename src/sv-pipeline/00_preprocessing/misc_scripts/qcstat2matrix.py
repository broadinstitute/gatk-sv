#!/bin/python

import sys
import collections

[_, fil, batch, ev, chroms_file] = sys.argv
dct = collections.defaultdict(dict)

with open(fil, 'r') as f:
    f.readline()
    for line in f:
        dat = line.rstrip().split()
        sample = dat[0]
        chrom = dat[1]
        number = dat[2]
        dct[sample][chrom] = number

with open(chroms_file, 'r') as f:
    chroms = [line.strip() for line in f]

with open(batch + '.' + ev + '.QC_matrix.txt', 'w') as g:
    g.write("#batch\tID\t" + '\t'.join(chroms) + '\n')
    for sample in dct.keys():
        g.write(batch + '\t' + sample)
        g.write('\t' + '\t'.join([dct[sample].get(c, '0') for c in chroms]))
        g.write('\n')
