# -*- coding: utf-8 -*-
#
"""
std_wham.py

Standardize WHAM record.
"""


from .standardize import VCFStandardizer


@VCFStandardizer.register('wham')
class WhamStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Wham record.
        """

        # Add placeholder ID if no ID is present
        if std_rec.id is None:
            std_rec.id = 'SV'
        else:
            # colons in the ID break VCF parsing
            std_rec.id = std_rec.id.replace(':', '_')

            # commas conflict with VCF format field with num='.'
            std_rec.id = std_rec.id.replace(',', ';')

        # SV types already standard
        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']

        # No CTX events
        std_rec.info['CHR2'] = std_rec.chrom
        std_rec.stop = raw_rec.stop

        # Strand not provided for inv/tloc
        if std_rec.info['SVTYPE'] == 'DEL':
            strands = '+-'
        elif std_rec.info['SVTYPE'] == 'DUP':
            strands = '-+'
        else:
            strands = '.'
        std_rec.info['STRANDS'] = strands

        # SVLEN is a list and can be negative
        std_rec.info['SVLEN'] = abs(raw_rec.info['SVLEN'][0])

        std_rec.info['ALGORITHMS'] = ['wham']

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Parse called samples from TAGS field
        """

        source = std_rec.info['ALGORITHMS'][0]

        # Any sample in TAGS field is considered to be called
        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            if sample in raw_rec.info['TAGS']:
                std_rec.samples[std_sample]['GT'] = (0, 1)
                std_rec.samples[std_sample][source] = 1
            else:
                std_rec.samples[std_sample]['GT'] = (0, 0)
                std_rec.samples[std_sample][source] = 0

        return std_rec
