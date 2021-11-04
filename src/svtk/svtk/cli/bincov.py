#!/usr/bin/env python3


"""
Calculates non-duplicate primary-aligned binned coverage
of a chromosome from an input BAM file
"""

# Import libraries
import argparse
import sys
from subprocess import call
import pysam
import pybedtools
import gzip
import shutil
import os


def countable_read(read):
    """
    Requirements to include a read when computing coverage.

    Require non-duplicate primary read alignments.
    """
    ok = (not read.is_duplicate and
          not read.is_unmapped and
          not read.is_secondary and
          not read.is_supplementary and
          read.reference_start > 0 and
          read.next_reference_start)  # TODO: check if this should be >0
    return ok


# Function to return read or fragment intervals from pysam.AlignmentFile
def filter_mappings(bam, mode='nucleotide'):
    """
    Generates bed intervals from a bam for a specific chromosome corresponding
    either to read coordinates or fragment coordinates

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Input bam
    mode : str
        'physical' or 'nucleotide' (default: 'physical')

    Returns
    ------
    mappings : BedTool
        Read or fragment intervals (depending on mode)
    """

    # Sanity check mode
    if mode not in 'nucleotide physical'.split():
        msg = 'Invalid mode: {0} (options: nucleotide, physical)'
        raise ValueError(msg.format(mode))

    # For nucleotide mode, return non-duplicate primary read mappings
    for read in bam:
        if not countable_read(read):
            continue

        if mode == 'nucleotide':
            yield '\t'.join([read.reference_name,
                             str(read.reference_start),
                             str(read.reference_end)]) + '\n'
        else:
            if read.is_read1 and read.is_proper_pair:
                fstart, fend = sorted([int(read.reference_start),
                                       int(read.next_reference_start)])
                if fstart < fend:
                    yield '\t'.join([read.reference_name,
                                     str(fstart), str(fend)]) + '\n'


# Function to evaluate nucleotide or physical coverage
def binCov(bam, chr, binsize, mode='nucleotide', overlap=0.05,
           blacklist=None, presubbed=False, oldBT=False):
    """
    Generates non-duplicate, primary-aligned nucleotide or physical coverage
    in regular bin sizes on a specified chromosome from a coordinate-sorted
    bamfile

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Input bam
    chr : string
        Chromosome to evaluate
    binsize : int
        Size of bins in bp
    mode : str
        Evaluate 'nucleotide' or 'physical' coverage
    overlap : float
        Maximum tolerated blacklist overlap before excluding bin
    blacklist : string
        Path to blacklist BED file
    presubbed : boolean
        Has the bam already been subsetted to the desired chromosome?
    oldBT : boolean
        Are you using a version of bedtools pre-2.24.0?

    Returns
    ------
    coverage : pybedtools.BedTool
        chr, start, end, coverage
    """

    # Create coverage bins and convert to BedTool
    maxchrpos = {d['SN']: d['LN'] for d in bam.header['SQ']}[chr]
    bin_starts = range(0, maxchrpos - binsize, binsize)
    bin_stops = range(binsize, maxchrpos, binsize)
    bins = []
    for i in range(0, len(bin_starts) - 1):
        bins.append([chr, bin_starts[i], bin_stops[i]])
    bins = pybedtools.BedTool(bins)

    # Remove bins that have at least 5% overlap with blacklist by size
    if blacklist is not None:
        blist = pybedtools.BedTool(blacklist)
        bins_filtered = bins.intersect(blist, v=True, f=overlap)
    else:
        bins_filtered = bins

    # Filter bam
    if presubbed:
        mappings = filter_mappings(bam, mode)
    else:
        mappings = filter_mappings(bam.fetch(chr), mode)
    bambed = pybedtools.BedTool(mappings)

    # Generate & return coverage
    if oldBT:
        coverage = bambed.coverage(bins_filtered, counts=True)
    else:
        coverage = bins_filtered.coverage(bambed, counts=True, sorted=True)
    return coverage


# Main function
def main(argv):
    # Add arguments
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk bincov')
    parser.add_argument('bam', type=str,
                        help='Input bam')
    parser.add_argument('chr', help='Contig to evaluate')
    parser.add_argument('cov_out', help='Output bed file of raw coverage')
    parser.add_argument('-n', '--norm_out', type=str,
                        help='Output bed file of normalized coverage')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size (bp) [1000]')
    parser.add_argument('-m', '--mode', default='nucleotide',
                        choices=['nucleotide', 'physical'],
                        help='Type of coverage to calculate '
                             '[nucleotide]')
    parser.add_argument('-x', '--blacklist', type=str, default=None,
                        help='BED file of regions to exclude')
    parser.add_argument('-z', '--gzip', default=False, action='store_true',
                        help='Compress output bed files')
    parser.add_argument('-p', '--presubsetted', dest='presubbed',
                        action='store_true', default=False,
                        help='Input bam is already subsetted to desired chr')
    parser.add_argument('-v', '--overlap', type=float, default=0.05,
                        help='Maximum fraction of each bin permitted to '
                        'overlap with blacklist regions. [0.05]')
    parser.add_argument('--oldBT', dest='oldBT', default=False,
                        action='store_true',
                        help='Using a bedtools version pre-2.24.0')
    parser.set_defaults(presubbed=False)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Correct filename for py3/py2 string inconsistency
    if args.bam.endswith("'") and args.bam.startswith("b'"):
        filename = args.bam[2:-1]
    else:
        filename = args.bam

    # Stores bam input as pysam.AlignmentFile
    bamfile = pysam.AlignmentFile(filename, 'rb')

    # Get coverage & write out
    coverage = binCov(bamfile, args.chr, args.binsize,
                      args.mode, args.overlap, args.blacklist,
                      args.presubbed, args.oldBT)
    coverage.saveas(args.cov_out)
    call('sort -Vk1,1 -k2,2n -o ' + args.cov_out + ' ' + args.cov_out,
         shell=True)

    # Gzip if optioned
    if args.gzip:
        f_in = open(args.cov_out, 'rb')
        f_out = gzip.open(args.cov_out + '.gz', 'wb')
        shutil.copyfileobj(f_in, f_out)
        os.remove(args.cov_out)

    # Normalize coverage (if optioned) & write out
    if args.norm_out is not None:
        ncoverage = coverage.to_dataframe(names='chr start end cov'.split())
        medcov = ncoverage.loc[ncoverage['cov'] > 0, 'cov'].median()
        ncoverage['cov'] = ncoverage['cov'] / medcov
        ncoverage.to_csv(args.norm_out, sep='\t', index=False, header=False)
        call(' '.join(['sort -Vk1,1 -k2,2n -o', args.norm_out,
                       args.norm_out]), shell=True)

    # Gzip if optioned
    if args.gzip is True:
        f_in = open(args.norm_out, 'rb')
        f_out = gzip.open(args.norm_out + '.gz', 'wb')
        shutil.copyfileobj(f_in, f_out)
        os.remove(args.norm_out)


if __name__ == '__main__':
    main(sys.argv[1:])
