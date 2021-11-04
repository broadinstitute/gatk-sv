# -*- coding: utf-8 -*-
#
"""
std_s.py

Standardize Lumpy records.
"""


from .standardize import VCFStandardizer, parse_bnd_pos
from svtk.utils import is_smaller_chrom


@VCFStandardizer.register('smoove')
class SmooveStandardizer(VCFStandardizer):
    def standardize_records(self):
        for record in self.filter_raw_vcf():
            # Split inversion events into their constituent breakpoints
            # For each strandedness in record, make a corresponding std record
            strands = record.info['STRANDS']
            for i, strand in enumerate(strands):
                record.info['STRANDS'] = (strand, )
                std_rec = self.std_vcf.new_record()
                std_rec = self.standardize_record(std_rec, record)

                # Some variants have stranded pairs that don't match their
                # SV type
                if std_rec.info['SVTYPE'] == 'DEL':
                    if std_rec.info['STRANDS'] != '+-':
                        continue
                if std_rec.info['SVTYPE'] == 'DUP':
                    if std_rec.info['STRANDS'] != '-+':
                        continue
                if std_rec.info['SVTYPE'] == 'INV':
                    if std_rec.info['STRANDS'] not in '++ --'.split():
                        continue

                # Tag split record
                std_rec.id += 'abcd'[i]

                yield std_rec

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Lumpy record.

        1) Add CHR2, END
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']

        # Strip per-strand counts
        std_rec.info['STRANDS'] = raw_rec.info['STRANDS'][0].split(':')[0]

        # Parse CHR2 and END
        if std_rec.info['SVTYPE'] == 'BND':
            chr2, end = parse_bnd_pos(std_rec.alts[0])

            # swap chr2/chrom, pos/end, and reverse strandedness
            if not is_smaller_chrom(std_rec.chrom, chr2):
                std_rec.pos, end = end, std_rec.pos
                std_rec.chrom, chr2 = chr2, std_rec.chrom
                std_rec.info['STRANDS'] = std_rec.info['STRANDS'][::-1]
        else:
            chr2, end = raw_rec.chrom, raw_rec.stop

        std_rec.info['CHR2'] = chr2
        std_rec.stop = end

        # Add SVLEN
        if std_rec.chrom == std_rec.info['CHR2']:
            std_rec.info['SVLEN'] = end - std_rec.pos
        else:
            std_rec.info['SVLEN'] = -1

        std_rec.info['ALGORITHMS'] = ['smoove']

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Parse called samples from TAGS field
        """

        source = std_rec.info['ALGORITHMS'][0]

        # Any sample in TAGS field is considered to be called
        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            if raw_rec.samples[sample]['GT'] == (0, 1):
                std_rec.samples[std_sample]['GT'] = (0, 1)
                std_rec.samples[std_sample][source] = 1
            else:
                std_rec.samples[std_sample]['GT'] = (0, 0)
                std_rec.samples[std_sample][source] = 0

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize ALT.

        When a BND is intrachromosomal, convert to appropriate SVTYPE
        """

        # Format BND ALT
        std_rec = super().standardize_alts(std_rec, raw_rec)

        stop = std_rec.stop
        if (std_rec.info['SVTYPE'] == 'BND' and
                std_rec.chrom == std_rec.info['CHR2']):
            if std_rec.info['STRANDS'] == '+-':
                std_rec.info['SVTYPE'] = 'DEL'
                std_rec.alts = ('<DEL>', )
            elif std_rec.info['STRANDS'] == '-+':
                std_rec.info['SVTYPE'] = 'DUP'
                std_rec.alts = ('<DUP>', )
            else:
                std_rec.info['SVTYPE'] = 'INV'
                std_rec.alts = ('<INV>', )
        std_rec.stop = stop
        return std_rec
