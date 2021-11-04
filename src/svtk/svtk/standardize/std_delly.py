# -*- coding: utf-8 -*-
#
"""
std_delly.py

Standardize a Delly record.
"""

from .standardize import VCFStandardizer
from svtk.utils import is_smaller_chrom


@VCFStandardizer.register('delly')
class DellyStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Delly record.

        1) Rename 'TRA' to 'BND'.
        2) Swap CHROM and CHR2 in translocations.
        3) Add END.
        4) Rename 'CT' to 'STRANDS' and convert notation.
        5) Compute SVLEN.
        6) Add ALGORITHMS.
        7) Standardize ALT to VCF spec.
        """

        # Rename TRA to BND
        svtype = raw_rec.info['SVTYPE']
        if svtype == 'TRA':
            svtype = 'BND'
        std_rec.info['SVTYPE'] = svtype
        END = std_rec.stop
        if "<" not in std_rec.alts[0]:
            std_rec.alts = ('<' + svtype + '>', )
        std_rec.ref = "N"
        std_rec.stop = END
        # Convert strandedness notation
        raw_strands = raw_rec.info['CT']
        if raw_strands == '5to3':
            strands = '-+'
        elif raw_strands == '3to3':
            strands = '++'
        elif raw_strands == '5to5':
            strands = '--'
        elif raw_strands == '3to5':
            strands = '+-'
        elif raw_strands == 'NtoN' and svtype == 'INS':
            strands = '+-'
        else:
            msg = 'Improper strands ({0}) in record {1}'
            raise Exception(msg.format(raw_strands, raw_rec.id))

        std_rec.info['STRANDS'] = strands

        pos, end = raw_rec.pos, raw_rec.stop
        if pos == 0:
            pos = 1

        # Swap CHR2/CHROM if necessary and update ALT
        if svtype == 'BND':
            chrom, chr2 = raw_rec.chrom, raw_rec.info['CHR2']

            # swap chr2/chrom, pos/end, and reverse strandedness
            if not is_smaller_chrom(chrom, chr2):
                if end == 0:
                    end = 1
                std_rec.pos, end = end, pos
                std_rec.chrom, chr2 = chr2, chrom
                std_rec.info['STRANDS'] = strands[::-1]

        else:
            chr2 = raw_rec.chrom

        # Add CHR2 and END
        std_rec.stop = end
        std_rec.info['CHR2'] = chr2

        # Add SVLEN
        if std_rec.chrom == std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = std_rec.stop - std_rec.pos
        else:
            std_rec.info['SVLEN'] = -1

        std_rec.info['ALGORITHMS'] = ['delly']

        return std_rec
