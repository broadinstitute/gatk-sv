# -*- coding: utf-8 -*-
#
"""
standardize.py

Standardize a VCF of SV calls.

Each record corresponds to a single SV breakpoint and will have the following
INFO fields, with specified constraints:
  SVTYPE:  SV type [DEL,DUP,INV,BND]
  CHR2:    Secondary chromosome [Must be lexicographically greater than CHROM]
  END:     SV end position (or position on CHR2 in translocations)
  STRANDS: Breakpoint strandedness [++,+-,-+,--]
  SVLEN:   SV length (-1 if translocation)
"""

import os
import tempfile
import pysam
from svtk.utils import make_bnd_alt, NULL_GT, parse_bnd_pos


def any_called(record):
    null_GTs = [(0, 0), (None, None), (0, ), (None, )]

    def _is_called(sample):
        return record.samples[sample]['GT'] not in null_GTs

    return any([_is_called(sample) for sample in record.samples])


class VCFStandardizer:
    subclasses = {}

    def __init__(self, raw_vcf, std_vcf, sample_names, prefix=None, min_size=50,
                 include_reference_sites=False, call_null_sites=False):
        """
        Standardize a VCF.

        Parameters
        ----------
        raw_vcf : pysam.VariantFile
            Input VCF.
        std_vcf : pysam.VariantFile
            Standardized VCF (containing no samples in the header).
        sample_names : list of str
            Sample names to use in header (must match order in the raw VCF)
        prefix : str, optional
            ID prefix to assign to standardized records
        min_size : int, optional
            Minimum SV size to report
        include_reference_sites : bool, optional
            Keep sites where no sample has a variant call (i.e. all 0/0 or ./.)
        call_null_sites : bool, optional
            Convert null genotypes (./.) to variant calls (0/1)
        """
        self.raw_vcf = raw_vcf
        self.std_vcf = std_vcf
        self.prefix = prefix
        self.min_size = min_size
        self.include_reference_sites = include_reference_sites
        self.call_null_sites = call_null_sites

        self.std_sample_names = list(sample_names)
        if len(std_vcf.header.samples) > 0 and list(std_vcf.header.samples) != self.std_sample_names:
            raise ValueError("Input std_vcf must have empty samples header, or the samples must be identical to "
                             "provided sample_names")

        num_out_samples = len(sample_names)
        num_raw_samples = len(self.raw_vcf.header.samples)
        if num_out_samples != num_raw_samples:
            raise ValueError("There are %d standardized sample names but the "
                             "raw vcf contains %d samples." % (num_out_samples,
                                                               num_raw_samples))

    @staticmethod
    def get_header_from_template(template_file, samples_list):
        # pysam can no longer open up a VCF header with FORMAT but no samples, so copy template to temporary file
        # and add samples, then open and return header
        append_samples = '\t'.join(samples_list)
        with tempfile.NamedTemporaryFile(suffix=".vcf") as named_temp_file:
            # use sed to add the sample ids to the end of the main header line, and redirect to temporary file
            os.system(f"sed 's/FORMAT$/FORMAT\t{append_samples}/' {template_file} > {named_temp_file.name}")
            # use pysam to open the temporary file and get the header
            return pysam.VariantFile(named_temp_file.name).header

    @classmethod
    def register(cls, source):
        def decorator(subclass):
            cls.subclasses[source] = subclass
            return subclass
        return decorator

    @classmethod
    def create(cls, source, *args):
        if source not in cls.subclasses:
            msg = 'No standardizer defined for {0}'.format(source)
            raise ValueError(msg)
        return cls.subclasses[source](*args)

    def filter_raw_vcf(self):
        """
        Filter records.

        Records are excluded if:
        1) They exist on contigs excluded from new file
        2) They are tagged as SECONDARY
        """
        for record in self.raw_vcf:
            if record.chrom not in self.std_vcf.header.contigs:
                continue

            if 'CHR2' in record.info:
                chr2 = record.info.get('CHR2', None)
                if chr2 not in self.std_vcf.header.contigs:
                    continue

            # Filter on chr2 if a breakend
            if '[' in record.alts[0] or ']' in record.alts[0]:
                chr2, end = parse_bnd_pos(record.alts[0])
                if chr2 not in self.std_vcf.header.contigs:
                    continue

            # Skip SECONDARY events to avoid double counting
            if 'SECONDARY' in record.info.keys():
                continue

            yield record

    def standardize_records(self):
        """
        Standardize every record in a VCF.

        Yields
        ------
        std_rec : pysam.VariantRecord
            Standardized records
        """
        for record in self.filter_raw_vcf():
            yield self.standardize_record(record)

    def standardize_vcf(self):
        """
        Standardize a VCF of SV records.

        Any filtering of records should be implemented in this method.

        Yields
        ------
        std_rec : pysam.VariantRecord
            Standardized records
        """

        if len(self.std_vcf.header.samples) == 0:
            # Add provided sample names to std header
            for sample in self.std_sample_names:
                self.std_vcf.header.add_sample(sample)

        idx = 1
        for std_rec in self.standardize_records():
            # Apply size filter (but keep breakends (SVLEN=-1))
            if 0 < std_rec.info['SVLEN'] < self.min_size:
                continue

            # Exclude insertions of unknown SVLEN
            if std_rec.info['SVTYPE'] == 'INS' and std_rec.info['SVLEN'] == -1:
                continue

            # Exclude sites with no called samples unless requested otherwise
            if not any_called(std_rec) and not self.include_reference_sites:
                continue

            # Filter unstranded breakpoints
            if std_rec.info['STRANDS'] not in '++ +- -+ --'.split():
                continue

            # Assign new variant IDs
            if self.prefix is not None:
                std_rec.id = '{0}_{1}'.format(self.prefix, idx)
                idx += 1

            yield std_rec

    def standardize_record(self, raw_rec):
        """
        Create a standardized copy of a VCF record.

        Parameters
        ----------
        raw_rec : pysam.VariantRecord

        Returns
        -------
        std_rec : pysam.VariantRecord
        """

        # Construct a new record and copy basic VCF fields
        std_rec = self.std_vcf.new_record()
        std_rec.chrom = raw_rec.chrom

        # pysam/htslib require non-negative pos
        if raw_rec.pos == 0:
            std_rec.pos = 1
        else:
            std_rec.pos = raw_rec.pos
        std_rec.id = raw_rec.id
        std_rec.ref = raw_rec.ref
        std_rec.alts = raw_rec.alts

        # Strip filters
        std_rec.filter.add('PASS')

        # Standardize the required INFO fields
        std_rec = self.standardize_info(std_rec, raw_rec)
        std_rec = self.standardize_alts(std_rec, raw_rec)
        std_rec = self.standardize_format(std_rec, raw_rec)

        return std_rec

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize VCF record INFO.

        When implementing this function, assume the basic data fields (CHROM,
        POS, ID, REF, ALTS, and FILTER) have been copied directly from the
        original record.

        The default implementation is a placeholder for testing and should be
        overridden with the logic for each algorithm's formatting.

        Parameters
        ----------
        std_rec : pysam.VariantRecord
            Standardized record with populated basic data.
        raw_rec : pysam.VariantRecord
            Raw record to be standardized.
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']
        std_rec.info['CHR2'] = raw_rec.chrom
        std_rec.stop = raw_rec.pos + 1
        std_rec.info['SVLEN'] = 0
        std_rec.info['ALGORITHMS'] = ['source']

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Copy desired FORMAT fields to new record.

        By default, copy GT and tag source FORMAT appropriately.

        Note: self.std_sample_names order must match raw_rec.samples
        """

        source = std_rec.info['ALGORITHMS'][0]

        # Add per-sample genotypes (ignoring other FORMAT fields)
        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            gt = raw_rec.samples[sample]['GT']
            if self.call_null_sites:
                if gt == (None, None):
                    gt = (0, 1)
                if gt == (None,):
                    gt = (1,)
            std_rec.samples[std_sample]['GT'] = gt

            if gt not in NULL_GT:
                std_rec.samples[std_sample][source] = 1
            else:
                std_rec.samples[std_sample][source] = 0

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize ALT field.

        Default behavior is to standardize BND alt to VCF spec and leave
        other SVTYPE alts untouched.
        """

        # Standardize tloc ALT after SVTYPE and CHR2/END are standardized
        if std_rec.info['SVTYPE'] == 'BND':
            alt = make_bnd_alt(std_rec.info['CHR2'], std_rec.info['END2'],
                               std_rec.info['STRANDS'])
            stop = std_rec.stop
            std_rec.alts = (alt, )
            std_rec.stop = stop

        return std_rec
