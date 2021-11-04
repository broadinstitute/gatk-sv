#!/usr/bin/env python

"""
Calculates non-duplicate primary-aligned binned coverage
of a chromosome from an input BAM file
"""

# Import libraries
import argparse
from subprocess import call
import pysam
import pybedtools
import gzip
import shutil
import os

# Define exception class for invalid coverage modes


class InvalidModeError(Exception):
    """Invalid coverage mode"""

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
        raise InvalidModeError('Invalid mode: ' + mode +
                               ' (options: nucleotide, physical)')

    # For nucleotide mode, return non-duplicate primary read mappings
    for read in bam:
        if (not any([read.is_duplicate, read.is_unmapped,
                     read.is_secondary, read.is_supplementary]) and
                all([read.reference_start > 0, read.next_reference_start])):
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
    if presubbed is True:
        mappings = filter_mappings(bam, mode)
    else:
        mappings = filter_mappings(bam.fetch(chr), mode)
    bambed = pybedtools.BedTool(mappings)

    # Generate & return coverage
    if oldBT is True:
        coverage = bambed.coverage(bins_filtered, counts=True)
    else:
        coverage = bins_filtered.coverage(bambed, counts=True, sorted=True)
    return coverage


# Main function
def main():
    # Add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam', type=str,
                        help='Input bam (or sam/cram, but requires appropriate flag)')
    parser.add_argument('chr', help='Contig to evaluate')
    parser.add_argument('cov_out', help='Output bed file of raw coverage')
    parser.add_argument('-S', '--SAM', default=False, action='store_true',
                        help='Input file is in sam format')
    parser.add_argument('-C', '--CRAM', default=False, action='store_true',
                        help='Input file is in cram format')
    parser.add_argument('-I', '--index_path', type=str, default=None,
                        help='Bam/cram index file')
    parser.add_argument('-z', '--gzip', dest='gzip', default=False,
                        action='store_true', help='Gzip output files'
                        ' bed files')
    parser.add_argument('-n', '--norm_out', type=str,
                        help='Output normalized coverage')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size, in bp (default: 1000)')
    parser.add_argument('-m', '--mode', default='nucleotide',
                        choices=['nucleotide', 'physical'],
                        help='Evaluate nucleotide or physical coverage '
                             '(default: nucleotide)')
    parser.add_argument('-x', '--blacklist', type=str, default=None,
                        help='BED file of regions to ignore')
    # parser.add_argument('-p', '--presubsetted', dest='presubbed',
    #                     action='store_true', help='Boolean flag to indicate'
    #                     ' if input bam is already subsetted to desired chr',
    #                     default=False)
    parser.add_argument('-v', '--overlap', nargs=1, type=float, default=0.05,
                        help='Maximum tolerated blacklist overlap before '
                        'excluding bin')
    parser.add_argument('--oldBT', dest='oldBT', default=False,
                        action='store_true', help='Flag to indicate if you are'
                        ' using a BEDTools version pre-2.24.0')
    parser.set_defaults(presubbed=False)
    args = parser.parse_args()

    # Correct filename for py3/py2 string inconsistency
    if args.bam.endswith("'") and args.bam.startswith("b'"):
        filename = args.bam[2:-1]
    else:
        filename = args.bam

    # Set appropriate read mode for bam/sam/cram input
    if args.SAM:
        read_mode = 'r'
    elif args.CRAM:
        read_mode = 'rc'
    else:
        read_mode = 'rb'

    # Read bamfile as pysam.AlignmentFile
    bamfile = pysam.AlignmentFile(
        filename, read_mode, index_filename=args.index_path)

    # Get coverage & write out
    coverage = binCov(bamfile, args.chr, args.binsize,
                      args.mode, args.overlap, args.blacklist,
                      args.presubbed, args.oldBT)
    coverage.saveas(args.cov_out)
    call('sort -Vk1,1 -k2,2n -o ' + args.cov_out + ' ' + args.cov_out,
         shell=True)
    # Gzip if optioned
    if args.gzip:
        with open(args.cov_out, 'rb') as f_in, gzip.open(args.cov_out + '.gz',
                                                         'wb') as f_out:
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
        if args.gzip:
            with open(args.norm_out, 'rb') as f_in, gzip.open(args.norm_out + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(args.norm_out)


# Main block
if __name__ == '__main__':
    main()
