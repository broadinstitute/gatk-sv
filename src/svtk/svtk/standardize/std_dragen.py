# -*- coding: utf-8 -*-
#
"""
std_dragen.py

Standardize a DRAGEN-SV record.
"""


from collections import deque
from svtk.utils import is_smaller_chrom, parse_bnd_pos, parse_bnd_strands
from .standardize import VCFStandardizer


def checkInversion(record):
    isInv3 = False
    isInv5 = False
    mateChr = ""
    matePos = -1

    def getMateInfo(splitChar):
        items = record.alts[0].split(splitChar)
        mateInfo = items[1].split(':')
        matePos = mateInfo[-1]
        mateChr = ":".join(mateInfo[:-1])
        matePos = int(matePos)
        return mateChr, matePos

    if record.alts[0].startswith('['):
        mateChr, matePos = getMateInfo('[')
        if mateChr == record.chrom:
            isInv5 = True
    elif record.alts[0].endswith(']'):
        mateChr, matePos = getMateInfo(']')
        if mateChr == record.chrom:
            isInv3 = True

    return isInv3, isInv5, matePos


@VCFStandardizer.register('dragen')
class DragenStandardizer(VCFStandardizer):
    def standardize_records(self):
        """
        Filter Dragen VCF.

        Skip mated events that are not marked with SECONDARY tag.
        """

        # Track IDs of observed records
        mate_IDs = deque()

        # Iterate over records
        for record in self.filter_raw_vcf():
            # Filter unmarked SECONDARY on same chromosome
            if 'MATEID' in record.info:
                mate_ID = record.info['MATEID'][0]

                # Skip records with an observed mate
                if mate_ID in mate_IDs:
                    continue

                mate_IDs.append(record.id)

            yield self.standardize_record(record)

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Dragen record.

        1. Define SVTYPE.
        2. Define ID.
        3. Define CHR2 and END.
        3. Define STRANDS.
        4. Define SVLEN.
        5. Define ALGORITHMS.
        """

        # Retrieve record data
        svtype = raw_rec.info['SVTYPE']
        isInv3, isInv5, _ = checkInversion(raw_rec)
        if isInv3 or isInv5:
            svtype = 'INV'

        # Convert INS to DUP for small DUP:TANDEM events
        if 'DUP:TANDEM' in raw_rec.id:
            svtype = 'DUP'

        # Define SVTYPE
        std_rec.info['SVTYPE'] = svtype

        # Define ID
        std_rec.id = std_rec.id.replace(':', '_')

        # Define CHR2 and END
        if svtype == 'BND' or svtype == 'INV':
            chrA, posA = raw_rec.chrom, raw_rec.pos
            chrB, posB = parse_bnd_pos(raw_rec.alts[0])
            if not is_smaller_chrom(chrA, chrB):
                posA, posB = posB, posA
                chrA, chrB = chrB, chrA
                std_rec.chrom = chrA
                std_rec.pos = posA
        elif svtype == 'INS':
            chrB = raw_rec.chrom
            posB = raw_rec.pos + 1
        elif svtype == 'DUP':
            chrB = raw_rec.chrom
            if isinstance(raw_rec.info.get('SVLEN', 0), tuple):
                svlen = raw_rec.info['SVLEN'][0]
            else:
                svlen = raw_rec.info.get('SVLEN', 0)
            posB = raw_rec.pos + svlen
        else:
            chrB = raw_rec.chrom
            posB = raw_rec.stop
        std_rec.info['CHR2'] = chrB
        std_rec.stop = posB

        # Define STRANDS
        if svtype == 'INV':
            if isInv3:
                strands = '++'
            else:
                strands = '--'
        elif svtype == 'BND':
            strands = parse_bnd_strands(raw_rec.alts[0])
        elif svtype == 'DEL':
            strands = '+-'
        elif svtype == 'DUP':
            strands = '-+'
        elif svtype == 'INS':
            strands = '+-'
        if not is_smaller_chrom(std_rec.chrom, std_rec.info['CHR2']):
            strands = strands[::-1]
        std_rec.info['STRANDS'] = strands

        # Define SVLEN
        if svtype == 'BND' and std_rec.chrom != std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = -1
        elif svtype == 'INS':
            std_rec.info['SVLEN'] = raw_rec.info.get('SVLEN', -1)
        else:
            std_rec.info['SVLEN'] = std_rec.stop - std_rec.pos

        # Define ALGORITHMS
        std_rec.info['ALGORITHMS'] = ['dragen']

        # Retain QUAL score
        std_rec.qual = raw_rec.qual

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Retain GT and GQ.
        """

        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            gt = raw_rec.samples[sample]['GT']
            gq = raw_rec.samples[sample]['GQ']
            if self.call_null_sites:
                if gt == (None, None):
                    gt = (0, 1)
                if gt == (None,):
                    gt = (1,)

            std_rec.samples[std_sample]['GT'] = gt
            std_rec.samples[std_sample]['GQ'] = gq

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize Dragen ALT field.

        When the full ref/alt sequence is specified for deletions or
        insertions, replace with N and <SVTYPE>.
        """

        # Format BND ALT
        std_rec = super().standardize_alts(std_rec, raw_rec)

        # Replace ALT sequence with svtype tag
        svtype = std_rec.info['SVTYPE']
        simple_alt = f'<{svtype}>'
        if svtype != 'BND':
            std_rec.alts = (simple_alt, )

        # Set reference to null for non-BND variants
        stop = std_rec.stop
        std_rec.ref = 'N'
        std_rec.stop = stop

        return std_rec
