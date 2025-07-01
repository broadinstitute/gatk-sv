# -*- coding: utf-8 -*-
#
"""
std_manta.py

Standardize a Manta record.
"""


from collections import deque
from svtk.utils import is_smaller_chrom, parse_bnd_pos, parse_bnd_strands
from .standardize import VCFStandardizer


@VCFStandardizer.register('manta')
class MantaStandardizer(VCFStandardizer):
    def standardize_records(self):
        """
        Filter Manta VCF.

        Mated events are not marked with SECONDARY tag; must filter manually.
        """

        mate_IDs = deque()
        for record in self.filter_raw_vcf():
            # Filter unmarked SECONDARY on same chromosome
            # TODO: Check if necessary to filter diff chromosomes
            if 'MATEID' in record.info:
                mate_ID = record.info['MATEID'][0]

                # Skip records with an observed mate
                if mate_ID in mate_IDs:
                    continue

                # Track IDs of observed records
                mate_IDs.append(record.id)

            yield self.standardize_record(record)

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Manta record.

        1) Replace colons in ID with underscores (otherwise breaks VCF parsing)
        2) Define CHR2 and END
        3) Add strandedness
        4) Add SVLEN
        """

        # Colons in the ID can break parsing
        std_rec.id = std_rec.id.replace(':', '_')

        svtype = raw_rec.info['SVTYPE']
        std_rec.info['SVTYPE'] = svtype

        # Define CHR2 and END
        if svtype == 'BND':
            chr2, end = parse_bnd_pos(raw_rec.alts[0])
            chrom, pos = raw_rec.chrom, raw_rec.pos
            if not is_smaller_chrom(chrom, chr2):
                pos, end = end, pos
                chrom, chr2 = chr2, chrom
                std_rec.pos = pos
                std_rec.chrom = chrom
            std_rec.info['CHR2'] = chr2
            std_rec.stop = std_rec.pos
            std_rec.info['END2'] = end
        elif svtype == 'INS':
            std_rec.info['CHR2'] = raw_rec.chrom
            std_rec.stop = raw_rec.pos + 1
        else:
            std_rec.info['CHR2'] = raw_rec.chrom
            std_rec.stop = raw_rec.stop

        # Strand parsing
        if svtype == 'INV':
            if 'INV3' in raw_rec.info.keys():
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

        if svtype == 'BND' and std_rec.chrom != std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = -1
        elif svtype == 'INS':
            std_rec.info['SVLEN'] = raw_rec.info.get('SVLEN', -1)
        else:
            std_rec.info['SVLEN'] = std_rec.stop - std_rec.pos

        std_rec.info['ALGORITHMS'] = ['manta']

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
        simple_alt = '<{0}>'.format(svtype)
        if svtype != 'BND':
            std_rec.alts = (simple_alt, )

        # Set reference to null
        # BNDs aren't symbolic so setting `ref` will reset `stop`
        stop = std_rec.stop
        std_rec.ref = 'N'
        std_rec.stop = stop

        return std_rec
