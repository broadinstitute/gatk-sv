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
    def standardize_records(self):
        """
        Filter Dragen VCF.

        Skip mated events that are not marked with SECONDARY tag.
        """

        mate_IDs = deque()
        for record in self.filter_raw_vcf():
            # Filter unmarked SECONDARY on same chromosome
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
        Standardize Dragen record.

        1. Replace colons in ID with underscores.
        2. Define CHR2 and END.
        3. Define strandedness.
        4. Define SVLEN.
        5. Define ALGORITHMS.
        """

        # Update ID
        std_rec.id = std_rec.id.replace(':', '_')

        # Update SVTYPE
        svtype = raw_rec.info['SVTYPE']
        std_rec.info['SVTYPE'] = svtype

        # Update CHR2 and END
        if svtype == 'BND':
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
        else:
            chrB = raw_rec.chrom
            posB = raw_rec.stop
        std_rec.info['CHR2'] = chrB
        std_rec.stop = posB

        # Update STRANDS
        if svtype == 'BND':
            strands = parse_bnd_strands(raw_rec.alts[0])
        elif svtype == 'DEL':
            strands = '+-'
        elif svtype == 'DUP':
            strands = '-+'
        elif svtype == 'INS':  # Treat DUPSVLEN as DUP
            if 'DUPSVLEN' in raw_rec.info:
                svtype = 'DUP'
                std_rec.info['SVTYPE'] = svtype
                strands = '-+'
            else:
                strands = '+-'
        else:  # Default
            strands = '+-'
        if not is_smaller_chrom(std_rec.chrom, std_rec.info['CHR2']):
            strands = strands[::-1]
        std_rec.info['STRANDS'] = strands

        # Update SVLEN
        if svtype == 'BND' and std_rec.chrom != std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = -1
        elif svtype == 'INS':
            std_rec.info['SVLEN'] = raw_rec.info.get('SVLEN', -1)
        else:
            std_rec.info['SVLEN'] = std_rec.stop - std_rec.pos

        # Update ALGORITHMS
        std_rec.info['ALGORITHMS'] = ['dragen']

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
