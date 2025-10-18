# -*- coding: utf-8 -*-
#
"""
std_depth.py

Standardize a depth record from Dragen CNV.
"""


from svtk.utils import is_smaller_chrom
from .standardize import VCFStandardizer, NULL_GT


@VCFStandardizer.register('depth')
class DragenCnvStandardizer(VCFStandardizer):
    def standardize_records(self):
        """
        Filter Dragen CNV VCF.

        Skip REF events.
        """

        # Iterate over records
        for record in self.filter_raw_vcf():
            # Skip records that are reference
            if 'SVTYPE' not in record.info:
                continue

            # Reset filters
            record.filter.clear()

            yield self.standardize_record(record)

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Dragen CNV record.

        1. Define SVTYPE.
        2. Define ID.
        3. Define CHR2 and END.
        3. Define STRANDS.
        4. Define SVLEN.
        5. Define ALGORITHMS.
        """

        # Define SVTYPE
        svtype = raw_rec.alts[0].replace('<', '').replace('>', '')
        std_rec.info['SVTYPE'] = svtype

        # Define ID
        std_rec.id = std_rec.id.replace(':', '_').replace('-', '_')

        # Define CHR2 and END
        std_rec.info['CHR2'] = raw_rec.chrom
        std_rec.stop = raw_rec.stop

        # Define STRANDS
        if svtype == 'DEL':
            strands = '+-'
        elif svtype == 'DUP':
            strands = '-+'
        if not is_smaller_chrom(std_rec.chrom, std_rec.info['CHR2']):
            strands = strands[::-1]
        std_rec.info['STRANDS'] = strands

        # Define SVLEN
        std_rec.info['SVLEN'] = raw_rec.stop - raw_rec.pos

        # Define ALGORITHMS
        std_rec.info['ALGORITHMS'] = ['depth']

        # Retain QUAL score
        std_rec.qual = raw_rec.qual

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Parse ./1 GTs for DUPs into 1/1.
        """

        source = std_rec.info['ALGORITHMS'][0]

        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            gt = raw_rec.samples[sample]['GT']
            if self.call_null_sites:
                if gt == (None, None):
                    gt = (0, 1)
                if gt == (None,):
                    gt = (1,)

            if None in gt or any(allele == "." for allele in gt):
                gt = (1, 1)
            std_rec.samples[std_sample]['GT'] = gt

            if gt not in NULL_GT:
                std_rec.samples[std_sample][source] = 1
            else:
                std_rec.samples[std_sample][source] = 0

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize Dragen CNV ALT field.

        Replace with N and SVTYPE.
        """

        # Set ALT to SVTYPE tag
        svtype = std_rec.info['SVTYPE']
        simple_alt = f'<{svtype}>'
        std_rec.alts = (simple_alt, )

        # Set REF to null
        stop = std_rec.stop
        std_rec.ref = 'N'
        std_rec.stop = stop

        return std_rec
