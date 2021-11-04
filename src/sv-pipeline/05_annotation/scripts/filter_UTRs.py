#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
For every gene in the provided GTF, subtract UTR sequence from the exonic
features. Requires the input GTF be sorted by gene_id (`sort -k10,10`).
"""

import argparse
import sys
import re
import pybedtools as pbt


def parse_attribute(attribute):
    """
    Convert GTF entry attribute to dictionary
    """
    data = attribute.strip(';').split(';')
    data = [d.strip() for d in data]
    pairs = [d.split() for d in data]
    attribute_dict = {k: v.strip('"') for k, v in pairs}

    return attribute_dict


class GTFEntry:
    def __init__(self, seqname, source, feature, start, end, score, strand,
                 frame, attribute, line):
        """
        Simple representation of an entry in a GTF file
        """
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = parse_attribute(attribute)
        self.line = line


def GTFParser(gtf):
    """
    Simple GTF parser - for each line, convert to GTFEntry
    """
    for line in gtf:
        yield GTFEntry(*line.strip().split('\t'), line)


def reformat_bedtool_GTF(entry):
    """
    pybedtools doesn't preserve the GTF format - replaces spaces in the
    attribute with tabs and eliminates the space after 'gene_id'
    """
    data = str(entry).split()
    metadata = '\t'.join(data[:8])
    attribute = ' '.join(data[9:])
    gene = re.sub('"', ' "', data[8], count=1)

    line = metadata + '\t' + gene + ' ' + attribute + '\n'

    return line


def update_exons(exons, UTRs):
    """
    Subtract UTR sequence from exonic features
    """
    exon_bt = pbt.BedTool(''.join(e.line for e in exons), from_string=True)
    UTR_bt = pbt.BedTool(''.join(e.line for e in UTRs), from_string=True)

    update = exon_bt.subtract(UTR_bt).saveas()

    for f in str(update).strip().split('\n'):
        line = reformat_bedtool_GTF(f)
        yield GTFEntry(*str(line).strip().split('\t'), str(line))


def filter_UTRs(gtf, fout):
    """
    For every gene in a GTF, subtract UTR sequence from exonic features. Leave
    other features untouched.
    """
    curr_gene = None
    exons = []
    UTRs = []

    for entry in gtf:
        if curr_gene is None:
            curr_gene = entry.attribute['gene_id']

        if curr_gene != entry.attribute['gene_id']:
            for UTR in UTRs:
                fout.write(UTR.line)
            for exon in update_exons(exons, UTRs):
                fout.write(exon.line)

            curr_gene = entry.attribute['gene_id']
            exons = []
            UTRs = []

        if entry.feature == 'exon':
            exons.append(entry)
        elif entry.feature == 'UTR':
            UTRs.append(entry)
        else:
            fout.write(entry.line)

    for UTR in UTRs:
        fout.write(UTR.line)
    for exon in update_exons(exons, UTRs):
        fout.write(exon.line)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gtf')
    parser.add_argument('fout')
    args = parser.parse_args()

    if args.gtf in '- stdin'.split():
        gtf = GTFParser(sys.stdin)
    else:
        gtf = GTFParser(open(args.gtf))

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    filter_UTRs(gtf, fout)


if __name__ == '__main__':
    main()
