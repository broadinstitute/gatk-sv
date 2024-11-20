# -*- coding: utf-8 -*-
#
"""
std_dragen.py

Standardize a Dragen record.
"""


from collections import deque
from svtk.utils import is_smaller_chrom, parse_bnd_pos, parse_bnd_strands
from .standardize import VCFStandardizer


@VCFStandardizer.register('dragen')
class DragenStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Dragen record.

        1. Replace colons in ID with underscores
        2. Define CHR2 and END
        3. Add strandedness
        4. Add SVLEN
        """

        # Replace colons in the ID
        std_rec.id = std_rec.id.replace(':', '_')

        svtype = raw_rec.info['SVTYPE']
        std_rec.info['SVTYPE'] = svtype

        # Define CHR2 and END
        if 'END' in raw_rec.info:
            end = raw_rec.info['END']
        else:
            end = raw_rec.stop

        std_rec.stop = end
        chr2 = raw_rec.chrom
        std_rec.info['CHR2'] = chr2

        # Strand parsing
        if svtype == 'BND':
            # Need to parse BND ALT
            strands = parse_bnd_strands(raw_rec.alts[0])
            chr2, end = parse_bnd_pos(raw_rec.alts[0])
            std_rec.info['CHR2'] = chr2
            std_rec.stop = end
            if not is_smaller_chrom(std_rec.chrom, chr2):
                # Swap chrom and pos
                std_rec.chrom, std_rec.info['CHR2'] = chr2, std_rec.chrom
                std_rec.pos, std_rec.stop = end, std_rec.pos
        elif svtype == 'DEL':
            strands = '+-'
        elif svtype == 'DUP':
            strands = '-+'
        elif svtype == 'INS':
            # Treat 'DUPSVLEN' as DUP rather than INS
            if 'DUPSVLEN' in raw_rec.info:
                svtype = 'DUP'
                std_rec.info['SVTYPE'] = svtype
                strands = '-+'
            else:
                strands = '+-'
        else:
            # Default strands
            strands = '+-'

        if not is_smaller_chrom(std_rec.chrom, std_rec.info['CHR2']):
            strands = strands[::-1]
        std_rec.info['STRANDS'] = strands

        # SVLEN
        if 'SVLEN' in raw_rec.info:
            svlen = raw_rec.info['SVLEN']
            if isinstance(svlen, list):
                svlen = svlen[0]
            std_rec.info['SVLEN'] = svlen
        else:
            if svtype == 'BND' and std_rec.chrom != std_rec.info['CHR2']:
                std_rec.info['SVLEN'] = -1
            elif svtype == 'INS':
                std_rec.info['SVLEN'] = -1
            else:
                std_rec.info['SVLEN'] = std_rec.stop - std_rec.pos

        std_rec.info['ALGORITHMS'] = ['dragen']

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize ALT.

        When the full ref/alt sequence is specified for deletions or
        insertions, replace with N and <SVTYPE>
        """

        # Format BND ALT
        std_rec = super().standardize_alts(std_rec, raw_rec)

        # Replace ALT sequence with svtype tag
        svtype = std_rec.info['SVTYPE']
        simple_alt = f'<{svtype}>'
        if svtype != 'BND':
            std_rec.alts = (simple_alt, )

        # Set reference to 'N' for non-BND variants
        # BNDs aren't symbolic so setting `ref` will reset `stop`
        stop = std_rec.stop
        std_rec.ref = 'N'
        std_rec.stop = stop

        return std_rec
