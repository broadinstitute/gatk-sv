# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
"""
svfile.py

Wrap the pysam API to permit clustering of standardized SV VCF records.
"""

import numpy as np
from collections import defaultdict
from .utils import recip, make_bnd_alt, update_best_genotypes, get_called_samples
from .genomeslink import GSNode


class SVFile(object):
    def __init__(self, vcf):
        """
        Wrapper for standardized VCF files.

        Parameters
        ----------
        vcf : pysam.VariantFile
        """
        self.reader = vcf
        self.filename = vcf.filename.decode('utf-8')
        self.samples = list(self.reader.header.samples)

        # Confirm all standard INFO fields are present
        required_info = 'SVTYPE CHR2 END STRANDS SVLEN ALGORITHMS'.split()
        for info in required_info:
            if info not in self.reader.header.info.keys():
                msg = "Required INFO field {0} not found in file {1}"
                msg = msg.format(info, self.filename)
                raise KeyError(msg)

        # Unfortunately no way to index into "source" metadata record
        # via pysam API, must manually check all header records
        self.sources = None
        for hrec in self.reader.header.records:
            if hrec.key == 'source':
                self.sources = hrec.value.split(',')

        if self.sources is None:
            msg = "Source not specified in header of {0}"
            msg = msg.format(self.filename)
            raise KeyError(msg)

    def fetch(self, chrom, start=None, end=None):
        """
        Fetch only calls from specified region of SVFile.

        Requires Tabix index. Updates SVFile in place.

        Parameters
        ----------
        chrom : str
        start : int, optional
        end : int, optional
            Required if start specified
        """
        if start is not None and end is None:
            msg = 'Start {}:{} specified but no end coordinate provided'
            msg = msg.format(chrom, start)
            raise ValueError(msg)

        # First check if VCF is empty
        try:
            pos = self.reader.tell()
            next(self.reader)
            self.reader.seek(pos)
        except StopIteration:
            self.reader = iter(())
            return

        # Then check if index is present
        try:
            self.reader = self.reader.fetch(chrom, start, end)
        except ValueError:
            msg = 'No index found for {0}'.format(self.filename)
            raise FileNotFoundError(msg)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        record = next(self.reader)
        return SVRecord(record)


class SVRecord(GSNode):
    """
    Clusterable VCF record.
    """

    def __init__(self, record):
        """
        record : pysam.VariantRecord
            Must specify 'CHR2' and 'END' in INFO
        """

        self.record = record
        self.sources = record.info['ALGORITHMS']
        self.called_samples = None

        chrA = record.chrom
        posA = record.pos
        chrB = record.info['CHR2']
        posB = record.stop
        name = record.id

        super().__init__(chrA, posA, chrB, posB, name)

    def get_called_samples_set(self):
        if self.called_samples is not None:
            return self.called_samples
        self.called_samples = set(get_called_samples(self.record))
        return self.called_samples

    def overlaps(self, other, frac=0.0):
        """
        Check if two records meet minimum reciprocal overlap.

        1) DEL, DUP, INV - reciprocal overlap calculated between SV intervals
        2) INS - intervals are calculated as insertion site plus SVLEN
        3) BND - always true
        """

        if self.is_tloc:
            return True
        if self.svtype == 'INS':
            # If either record's insertion length is unknown, consider
            # overlap met
            svlens = self.record.info['SVLEN'], other.record.info['SVLEN']
            if svlens[0] == -1 or svlens[1] == -1:
                return True
            # Otherwise, model insertion region as (coord + ins length)
            else:
                posBs = (self.posA + svlens[0], other.posA + svlens[1])
        else:
            posBs = (self.posB, other.posB)

        return recip(self.posA, posBs[0], other.posA, posBs[1], frac)

    @property
    def svtype(self):
        """
        Returns
        -------
        svtype : str
            One of {DEL, DUP, INV, BND}
        """
        return self.record.info['SVTYPE']

    @property
    def is_tloc(self):
        return self.chrA != self.chrB

    def __hash__(self):
        return id(self)


class SVRecordCluster:
    def __init__(self, records):
        self.records = records

    def sources(self):
        """
        Return list of source algorithms in clustered records.

        By default, searches for an ALGORITHMS INFO field

        Returns
        -------
        sources : list of str
        """
        call_sources = set()
        for record in self.records:
            sources = record.record.info.get('ALGORITHMS')
            if sources:
                call_sources = call_sources.union(sources)

        return sorted(call_sources)

    def merge_record_data(self, new_record):
        """
        Aggregate metadata (coordinates, alts, and INFO) of clustered records.

        * Original record IDs are preserved in a new INFO field
        * svtype is reported for base record

        Parameters
        ----------
        new_record : pysam.VariantRecord
            Blank record to fill with aggregated data

        Returns
        -------
        new_record : pysam.VariantRecord
            VCF record populated with aggregate cluster data
        """

        # Secondary records have duplicate POS/END/INFO
        # TODO: move secondary filtering elsewhere

        if len(self.records) == 0:
            return None

        base_record = self.records[0]

        new_record.chrom = base_record.chrA
        new_record.ref = base_record.record.ref
        new_record.info['SVTYPE'] = base_record.svtype
        new_record.info['CHR2'] = base_record.chrB

        # TODO: Check if strands should be merged
        new_record.info['STRANDS'] = base_record.record.info['STRANDS']

        # Merge coordinates
        POS, END, CIPOS, CIEND = self.merge_pos()
        new_record.pos = POS
        new_record.stop = END
        #  new_record.info['CIPOS'] = CIPOS
        #  new_record.info['CIEND'] = CIEND

        # Assign alts, updating translocation alt based on merged coordinates
        if new_record.info['SVTYPE'] == 'BND':
            strands = new_record.info['STRANDS']
            alt = make_bnd_alt(base_record.chrB, END, strands)
            new_record.alts = (alt, )
            new_record.stop = END
        if new_record.info['SVTYPE'] == 'INS':
            alts = set()
            for record in self.records:
                alts = alts.union(record.record.alts)
            if len(alts) > 1:
                alts = tuple(a for a in alts if a != '<INS>')
            new_record.alts = alts
            new_record.stop = END
        else:
            new_record.alts = base_record.record.alts
            new_record.stop = END

        # SVLEN for intra-chromosomal is -1
        if base_record.is_tloc:
            new_record.info['SVLEN'] = -1
        # SVLEN for insertions is median of observed insertions
        elif base_record.record.info['SVTYPE'] == 'INS':
            svlens = [r.record.info['SVLEN'] for r in self.records]
            svlens = [svlen for svlen in svlens if svlen > -1]
            if len(svlens) == 0:
                svlen = -1
            else:
                svlen = int(np.around(np.median(svlens)))
            new_record.info['SVLEN'] = svlen
        else:
            new_record.info['SVLEN'] = END - POS

        # QUAL, FILTER currently unused
        new_record.filter.add('PASS')

        # Report cluster RMSSTD
        #  new_record.info['RMSSTD'] = self.rmsstd

        # List of aggregate sources
        new_record.info['ALGORITHMS'] = self.sources()

        # If merging, all will have same CLUSTER ID
        if 'CLUSTER' in base_record.record.info:
            new_record.info['CLUSTER'] = base_record.record.info['CLUSTER']

        return new_record

    def merge_record_infos(self, new_record, header):
        """
        Aggregate INFO fields in child records
        """
        PROTECTED_INFOS = ('SVTYPE CHR2 END STRANDS SVLEN ALGORITHMS CIPOS '
                           'CIEND RMSSTD MEMBERS').split()
        records = [r.record for r in self.records]
        infos = defaultdict(list)

        for record in records:
            for info, value in record.info.items():
                if info not in PROTECTED_INFOS:
                    infos[info].append(value)

        for info, values in infos.items():
            if header.info[info].type == 'Flag':
                new_record.info[info] = True
            elif header.info[info].type == 'String':
                if header.info[info].number == '.':
                    new_record.info[info] = sorted(
                        set([v for vlist in values for v in vlist]))
                elif header.info[info].number == 1:
                    new_record.info[info] = ','.join(
                        sorted(set([v for vlist in values for v in vlist])))
                else:
                    new_record.info[info] = [
                        ','.join(vlist) for vlist in zip(values)]

            # TODO merge numeric INFO
            elif info == 'varGQ':
                new_record.info[info] = max(values)
            else:
                pass

        return new_record

    def merge_record_formats(self, new_record, sourcelist,
                             preserve_genotypes=False):
        """
        Aggregate sample genotype data across records.

        1) Set GT to 0/1 for samples called in any record, 0/0 otherwise.
        2) For each provided source, set a corresponding FORMAT field to 1
           for samples called by the source. If the FORMAT field is available
           in the record being merged, use per-sample data. If the FORMAT field
           is not available, use the record's ALGORITHMS to determine
           support for all samples called in that record.

        Parameters
        ----------
        new_record : pysam.VariantRecord
            Record to populate with FORMAT data
        sourcelist : list of str
            List of all sources to add a FORMAT field for

        Returns
        -------
        new_record : pysam.VariantRecord
            Populated record
        """

        # Seed with null values
        for sample in new_record.samples:
            new_record.samples[sample]['GT'] = (0, 0)

            #  for source in sourcelist:
            #  new_record.samples[sample][source] = 0

        # Update with called samples
        if preserve_genotypes:
            for fmt in new_record.format.keys():
                if fmt != 'GT':
                    del new_record.format[fmt]

            # then overwrite genotypes of non-multiallelic sites as necessary
            records = [r.record for r in self.records]
            update_best_genotypes(new_record, records,
                                  preserve_multiallelic=True)

        # TODO: optionally permit ./. instead of rejecting
        # I think that was an issue with one caller, maybe handle in preproc
        else:
            null_GTs = [(0, 0), (None, None), (0, ), (None, )]
            for record in self.records:
                for sample in record.record.samples:
                    gt = record.record.samples[sample]['GT']

                    # Skip samples without a call
                    if gt in null_GTs:
                        continue

                    # Otherwise call the sample in the new record
                    new_record.samples[sample]['GT'] = (0, 1)

        return new_record

    def merge_pos(self):
        """
        Compute aggregate POS/END of clustered SVRecords.

        Defaults to computing median of POS and END over constituent
        records. CIPOS/CIEND are calculated as (MIN - MEDIAN, MAX - MEDIAN)
        over POS and END.

        Returns
        -------
        POS : int
        END : int
        CIPOS : list of int
        CIEND : list of int
        """

        # Position bounds
        MIN_POS = min(rec.posA for rec in self.records)
        MAX_POS = max(rec.posA for rec in self.records)
        MIN_END = min(rec.posB for rec in self.records)
        MAX_END = max(rec.posB for rec in self.records)

        POS = int(np.median([rec.posA for rec in self.records]))
        END = int(np.median([rec.posB for rec in self.records]))
        CIPOS = [MIN_POS - POS, MAX_POS - POS]
        CIEND = [MIN_END - END, MAX_END - END]

        return POS, END, CIPOS, CIEND

    @property
    def rmsstd(self):
        """
        Root-mean-square standard deviation of cluster coordinates
        """

        if hasattr(self, '_rmsstd'):
            return self._rmmstd

        starts = np.array([record.posA for record in self.records])
        ends = np.array([record.posB for record in self.records])

        def _meanSS(X):
            mu = np.mean(X)
            return np.sum((X - mu) ** 2) / len(X)

        SS = _meanSS(starts) + _meanSS(ends)
        self._rmmstd = np.sqrt(SS)

        return self._rmmstd
