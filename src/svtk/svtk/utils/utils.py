# -*- coding: utf-8 -*-
#
"""
utils.py

Helper functions for svtk.
"""

from collections import deque
import pysam
import pybedtools as pbt


NULL_GT = [(0, 0), (None, None), (0, ), (None, )]


def is_smaller_chrom(chrA, chrB):
    """
    Test if chrA is naturally less than chrB

    Returns True if chrA == chrB so comparison will default to position
    """

    if chrA.startswith('chr'):
        chrA = chrA[3:]
    if chrB.startswith('chr'):
        chrB = chrB[3:]

    # Numeric comparison, if possible
    if chrA.isdigit() and chrB.isdigit():
        return int(chrA) <= int(chrB)

    # String comparison for X/Y
    elif not chrA.isdigit() and not chrB.isdigit():
        return chrA <= chrB

    # Numeric is always less than X/Y
    else:
        return chrA.isdigit()


def is_biallelic(record):
    """
    Check if record is biallelic
    """
    if 'MULTIALLELIC' not in record.filter \
            and len(record.alleles) <= 2 \
            and record.info['SVTYPE'] not in 'CNV MCNV'.split():
        return True
    else:
        return False


def recip(startA, endA, startB, endB, frac):
    """
    Test if two intervals share a specified reciprocal overlap.
    """

    if frac == 0:
        return True

    start = max(startA, startB)
    end = min(endA, endB)
    olen = end - start
    lenA = endA - startA
    lenB = endB - startB

    try:
        lapA = olen / float(lenA)
        lapB = olen / float(lenB)
    except ZeroDivisionError:
        return False

    return (olen > 0) and (lapA >= frac) and (lapB >= frac)


def make_bnd_alt(chrom, pos, strands, ref_base='N'):
    """
    Make ALT for BND record in accordance with VCF specification.
    """

    p = '{0}:{1}'.format(chrom, pos)
    t = ref_base

    if strands == '++':
        fmt = '{1}]{0}]'
    elif strands == '+-':
        fmt = '{1}[{0}['
    elif strands == '-+':
        fmt = ']{0}]{1}'
    elif strands == '--':
        fmt = '[{0}[{1}'

    return fmt.format(p, t)


def get_called_samples(record, include_null=False):
    """
    Return list of samples with variant call

    Parameters
    ----------
    record : pysam.VariantRecord
    include_null : bool
        Include samples without an explicit reference (0/0) call (i.e. ./.)

    Returns
    -------
    samples : list of str
        Sorted list of sample IDs with a variant call
    """

    samples = deque()
    for sample in record.samples.keys():
        if record.samples[sample]['GT'] not in NULL_GT:
            samples.append(sample)

    if record.info.get('SVTYPE', None) == 'CNV':
        for sample in record.samples.keys():
            if record.samples[sample]['CN'] != 2:
                samples.append(sample)

    return sorted(samples)


# TODO: check if record is CPX and make entry per complex interval
def vcf2bedtool(vcf, split_bnd=True, include_samples=False,
                include_strands=True, split_cpx=False, include_infos=None,
                annotate_ins=True, report_alt=False, svtypes=None,
                no_sort_coords=False, simple_sinks=False,
                include_unresolved=True, include_filters=False):
    """
    Wrap VCF as a bedtool. Necessary as pybedtools does not support SV in VCF.

    Parameters
    ----------
    vcf : str or pysam.VariantFile
    split_bnd : bool, optional
        Provide two records for each BND, one per breakend
    include_samples : bool, optional
        Provide comma-delimited list of called samples
    include_strands : bool, optional
        Provide breakpoint strandedness
    include_infos : list of str, optional
        INFO fields to add as columns in output. If "ALL" is present in the
        list, all INFO fields will be reported.
    annotate_ins : bool, optional
        Rename SVTYPE of insertion records to DEL to annotate sink.
    report_alt : bool, optional
        Report record's ALT as SVTYPE in bed
    svtypes : list of str, optional
        Whitelist of SV types to restrict generated bed to
    no_sort_coords : bool, optional
        Do not sort start & end coordinates
    simple_sinks : bool, optional
        Treat all insertion sinks as single-bp windows
    include_unresolved : bool, optional
        Output unresolved variants
    include_filters : bool, optional
        Output FILTER field after INFO fields

    Returns
    -------
    bt : pybedtools.BedTool
        SV converted to Bedtool. Starts of BND records are assigned as pos - 1.
        Note that BED is 0-based but VCF is 1-based, so 1bp is subtracted from
        all pos values when converted from VCF to BED
        Included columns: chrom, start, end, name, svtype, strands
    """

    cpx_ins_classes = 'dDUP dDUP_iDEL INS_iDEL'.split()

    if not isinstance(vcf, pysam.VariantFile):
        vcf = pysam.VariantFile(vcf)

    entry = '{chrom}\t{start}\t{end}\t{name}\t{svtype}'
    if include_strands:
        entry += '\t{strands}'
    if include_samples:
        entry += '\t{samples}'
    if include_infos:
        if 'ALL' in include_infos:
            include_infos = vcf.header.info.keys()
        entry += '\t{infos}'
    if include_filters:
        entry += '\t{filters}'
    entry += '\n'

    def _format_info(info):
        if info is None:
            return 'NA'
        elif isinstance(info, tuple) or isinstance(info, list):
            return ','.join([str(x) for x in info])
        else:
            return str(info)

    # Convert each record in vcf to bed entry
    def _converter():
        for record in vcf:
            if svtypes is not None and record.info['SVTYPE'] not in svtypes:
                continue
            if not include_unresolved:
                if 'UNRESOLVED' in record.info.keys() \
                        or 'UNRESOLVED_TYPE' in record.info.keys() \
                        or 'UNRESOLVED' in record.filter:
                    continue

            chrom = record.chrom
            name = record.id

            # Set start & end coordinates to appropriate sorted order
            # for all records (to not break bedtools)
            if no_sort_coords:
                start = int(record.pos)
                end = int(record.stop)
            else:
                start, end = sorted([int(record.pos), int(record.stop)])

            # Subtract 1bp from pos to convert to 0-based BED vs 1-based VCF
            start = max([0, int(start) - 1])

            if report_alt:
                svtype = record.alts[0].strip('<>')
            else:
                svtype = record.info['SVTYPE']

            if include_strands:
                try:
                    strands = record.info.get('STRANDS', None)
                except ValueError:
                    strands = None
                strands = '.' if strands is None else strands
            if include_samples:
                samples = ','.join(get_called_samples(record))
            if include_infos:
                infos = deque()
                for key in include_infos:
                    # Can't access END through info
                    if key == 'END':
                        infos.append(record.stop)
                    else:
                        infos.append(record.info.get(key))

                # reformat for tabular output
                infos = [_format_info(v) for v in infos]
                infos = '\t'.join(infos)
            if include_filters:
                filters = [f for f in record.filter]
                filters = ','.join(filters)

            if record.info.get('SVTYPE', None) == 'BND':
                # First end of breakpoint
                start = max([0, int(record.pos) - 1])
                end = record.pos
                yield entry.format(**locals())

                # Second end of breakpoint
                if split_bnd:
                    chrom = record.info['CHR2']
                    start = max([0, int(record.stop) - 1])
                    end = record.stop
                    yield entry.format(**locals())

            elif record.info.get('SVTYPE', None) == 'INS':
                # Only yield insertion sinks for now
                # Treat them as deletions
                # TODO: rename CPX_INTERVALS to SOURCE for insertions
                if annotate_ins:
                    svtype = 'DEL'
                # if not no_sort_coords:
                #     start, end = sorted([start, end])
                # Reduce insertion sinks to single-bp intervals if optioned
                if simple_sinks:
                    start = max([0, int(record.pos)])
                    end = start + 1
                yield entry.format(**locals())

            elif record.info.get('SVTYPE', None) == 'CTX':
                start = max([0, int(record.pos) - 1])
                end = start + 1
                yield entry.format(**locals())

                # Second end of breakpoint
                if split_bnd:
                    chrom = record.info['CHR2']
                    start = max([0, int(record.stop) - 1])
                    end = record.stop
                    yield entry.format(**locals())

            # Deconstruct complex intervals, if optioned
            elif 'CPX_INTERVALS' in record.info and split_cpx:
                # If complex, all constituent intervals are in CPX_INTERVALS
                for interval in record.info['CPX_INTERVALS']:
                    svtype, region = interval.split('_')
                    chrom, coords = region.split(':')
                    start, end = coords.split('-')
                    start = max([0, int(start) - 1])
                    yield entry.format(**locals())
                # If complex insertion, return insertion point as 1bp DEL
                if record.info.get('CPX_TYPE', None) in cpx_ins_classes:
                    svtype = 'DEL'
                    chrom = record.chrom
                    start = max([0, int(record.pos) - 1])
                    end = record.pos
                    yield entry.format(**locals())

            # elif (record.info.get('SVTYPE', None) == 'CPX' and
            #       'CPX_TYPE' in record.info.keys()):
            #     if (record.info.get('CPX_TYPE', None) in cpx_ins_classes):
            #         if annotate_ins:
            #             svtype = 'DEL'
            #         yield entry.format(**locals())

            #     if split_cpx:
            #         if 'dDUP' in record.info.get('CPX_TYPE', None):
            #             svtype = 'DUP'
            #         else:
            #             svtype = 'INS'
            #         source = record.info.get('SOURCE')
            #         region = source.split('_')[1]
            #         chrom, coords = region.split(':')
            #         start, end = coords.split('-')
            #         yield entry.format(**locals())

            else:
                if not no_sort_coords:
                    start, end = sorted([start, end])
                    if start == end:
                        end += 1
                yield entry.format(**locals())

    return pbt.BedTool(_converter()).saveas()


def set_null(record, sample):
    """Remove sample call from a VariantRecord"""
    dat = record.samples[sample].items()

    # Set genotype to no-call
    record.samples[sample]['GT'] = (0, 0)

    for fmt, value in dat:
        if fmt == 'GT':
            continue

        # Get type and count of FORMAT
        n = record.format[fmt].number
        dtype = record.format[fmt].type

        # Set null value based on FORMAT type
        if dtype in 'Integer Float'.split():
            null_val = 0
        elif dtype in 'String Character'.split():
            null_val = ''
        else:
            raise ValueError('Invalid VCF FORMAT type: {0}'.format(dtype))

        # Set the appropriate count of values
        if n == 1:
            record.samples[sample][fmt] = null_val
        elif n == '.':
            record.samples[sample][fmt] = (null_val, )
        else:
            record.samples[sample][fmt] = tuple(null_val for i in range(n))


def samples_overlap(samplesA, samplesB, upper_thresh=0.5, lower_thresh=0.5):
    """
    Test if two sets of samples are sufficiently similar.

    The two sets are tested to determine which fraction of their respective
    samples appears in the other set. The set with greater overlap must have at
    least upper_thresh% of its samples appear in the other set, while the set with lesser
    overlap must have at least lower_thresh% of its samples appear in the other set.

    Arguments
    ---------
    samplesA : list of str OR pysam.VariantRecord
    samplesB : list of str OR pysam.VariantRecord
    """

    if isinstance(samplesA, pysam.VariantRecord):
        samplesA = get_called_samples(samplesA)
        samplesB = get_called_samples(samplesB)

    # Get lists of called samples for each record
    if samplesA is not set:
        samplesA = set(samplesA)
    if samplesB is not set:
        samplesB = set(samplesB)

    # Compute fraction of each record's samples which are shared
    if len(samplesA) > 0 and len(samplesB) > 0:
        shared = samplesA & samplesB
        fracA = len(shared) / len(samplesA)
        fracB = len(shared) / len(samplesB)
    else:
        fracA = 0
        fracB = 0

    min_frac, max_frac = sorted([fracA, fracB])

    return min_frac >= lower_thresh and max_frac >= upper_thresh
