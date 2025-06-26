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

        # Add placeholder ALT and REF to avoid pysam errors
        if "<" not in std_rec.alts[0]:
            std_rec.alts = ('<' + svtype + '>', )
        std_rec.ref = "N"

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
                pos, end = end, pos
                chrom, chr2 = chr2, chrom
                std_rec.info['STRANDS'] = strands[::-1]
            
            std_rec.chrom = chrom
            std_rec.pos = pos
            std_rec.info['CHR2'] = chr2
            std_rec.info['END2'] = end
            std_rec.stop = pos

        else:
            chr2 = raw_rec.chrom
            std_rec.info['CHR2'] = chr2
            std_rec.stop = end

        # Add SVLEN
        if std_rec.chrom == std_rec.info['CHR2']:
            sv_end = std_rec.info.get('END2') if svtype == 'BND' else std_rec.stop
            std_rec.info['SVLEN'] = sv_end - std_rec.pos
        else:
            std_rec.info['SVLEN'] = -1

        std_rec.info['ALGORITHMS'] = ['delly']

        return std_rec
