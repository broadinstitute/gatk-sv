#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Resolve clustered records into a complex SV
"""

import numpy as np
import pybedtools as pbt
import svtk.utils as svu
from .cpx_inv import classify_complex_inversion
from .cpx_tloc import classify_simple_translocation, classify_insertion


class ComplexSV:
    def __init__(self, records, cytobands, mei_bed, SR_only_cutoff):
        """
        Parameters
        ----------
        records : list of pysam.VariantRecord
            Clustered records to resolve
        cytobands : pysam.TabixFile
            Cytoband bed file (to classify interchromosomal)
        mei_bed : pybedtools.BedTool
        """

        self.records = records
        self.cytobands = cytobands
        self.mei_bed = mei_bed
        #  self.rdtest = rdtest

        self.organize_records()

        self.make_record()
        self.resolve(SR_only_cutoff=SR_only_cutoff)
        self.clean_record()

    def organize_records(self):
        self.inversions = [
            r for r in self.records if r.info['SVTYPE'] == 'INV']
        self.tlocs = [r for r in self.records if r.chrom != r.info['CHR2']]
        self.breakends = [r for r in self.records if (r.chrom == r.info['CHR2']) and
                          (r.info['SVTYPE'] == 'BND')]
        self.insertions = [
            r for r in self.records if r.info['SVTYPE'] == 'INS']

        cnvtypes = 'DEL DUP'.split()
        self.cnvs = [r for r in self.records if r.info['SVTYPE'] in cnvtypes]
        self.dels = [r for r in self.records if r.info['SVTYPE'] == 'DEL']
        self.dups = [r for r in self.records if r.info['SVTYPE'] == 'DUP']

    def remove_SR_only_breakpoints(self):
        def _is_SR_only(record):
            return record.info.get('EVIDENCE', None) == ('SR', )

        if 'EVIDENCE' in self.records[0].header.info.keys():
            clean_records = [r for r in self.records if not _is_SR_only(r)]
            if len(clean_records) > 0:
                self.records = clean_records

    def resolve(self, SR_only_cutoff):
        self.set_cluster_type()
        if self.cluster_type == 'CANDIDATE_INVERSION':
            self.resolve_inversion(SR_only_cutoff=SR_only_cutoff)
            if self.svtype == 'UNR':
                self.set_unresolved()
                self.vcf_record.info['SVTYPE'] = 'BND'
                self.vcf_record.info['UNRESOLVED_TYPE'] = 'SR_ONLY_LARGE_INVERSION'
                self.cluster_type = 'SR_ONLY_LARGE_INVERSION'
                for r in self.records:
                    r.info['UNRESOLVED'] = True
                    r.info['UNRESOLVED_TYPE'] = 'SR_ONLY_LARGE_INVERSION'
                    r.info['SVTYPE'] = 'BND'
                    if 'CPX_TYPE' in r.info.keys():
                        r.info.pop('CPX_TYPE')
                    if 'CPX_INTERVALS' in r.info.keys():
                        r.info.pop('CPX_INTERVALS')
        elif self.cluster_type == 'CANDIDATE_TRANSLOCATION':
            self.resolve_translocation()
        elif self.cluster_type == 'CANDIDATE_INSERTION':
            self.resolve_insertion()
        elif self.cluster_type == 'CANDIDATE_DISDUP':
            self.resolve_disdup()
        elif self.cluster_type == 'RESOLVED_INSERTION':
            self.report_simple_insertion()
        elif self.cluster_type == 'SINGLE_INSERTION_STRIP_CNVS':
            self.report_insertion_strip_CNVs()
        elif self.cluster_type == 'TANDEM_DUPLICATION_FLANKING_INSERTIONS':
            self.report_manta_tandem_dup()

        # Second pass through the record after excluding SR-only breakpoints
        # if first pass is unsuccessful
        else:
            # Collect SR-only records that are not remaining in existing cluster
            # Note: if the cluster only had a single breakpoint to begin with,
            # it will not be removed by .remove_SR_only_breakpoints(), even if
            # it is SR only. The below code is necessary to prevent duplicating
            # SR-only single-enders
            def _is_SR_only(record):
                return record.info.get('EVIDENCE', None) == ('SR', )
            non_sr_only_rec_ids = [
                r.id for r in self.records if not _is_SR_only(r)]
            if len(non_sr_only_rec_ids) > 0:
                sr_only_recs = [r for r in self.records if _is_SR_only(
                    r) and r.id not in non_sr_only_rec_ids]
            else:
                sr_only_recs = []

            self.remove_SR_only_breakpoints()
            self.organize_records()
            self.set_cluster_type()

            if self.cluster_type in 'CANDIDATE_INVERSION CANDIDATE_TRANSLOCATION '.split() + \
                'CANDIDATE_INSERTION RESOLVED_INSERTION'.split() \
                    and len(self.records) > 0:
                if self.cluster_type == 'CANDIDATE_INVERSION':
                    self.resolve_inversion(SR_only_cutoff=SR_only_cutoff)
                    if self.svtype == 'UNR':
                        self.set_unresolved()
                        self.vcf_record.info['SVTYPE'] = 'BND'
                        self.vcf_record.info['UNRESOLVED_TYPE'] = 'SR_ONLY_LARGE_INVERSION'
                        self.cluster_type = 'SR_ONLY_LARGE_INVERSION'
                        for r in self.records:
                            r.info['UNRESOLVED'] = True
                            r.info['UNRESOLVED_TYPE'] = 'SR_ONLY_LARGE_INVERSION'
                            r.info['SVTYPE'] = 'BND'
                            if 'CPX_TYPE' in r.info.keys():
                                r.info.pop('CPX_TYPE')
                            if 'CPX_INTERVALS' in r.info.keys():
                                r.info.pop('CPX_INTERVALS')
                elif self.cluster_type == 'CANDIDATE_TRANSLOCATION':
                    self.resolve_translocation()
                elif self.cluster_type == 'CANDIDATE_INSERTION':
                    self.resolve_insertion()
                elif self.cluster_type == 'RESOLVED_INSERTION':
                    self.report_simple_insertion()
                elif self.cluster_type == 'CANDIDATE_DISDUP':
                    self.resolve_disdup()

            else:
                # Add back SR-only records if cluster is still unresolved,
                # and if SR-only record doesnt match only existing cluster
                if len(sr_only_recs) > 0:
                    for r in sr_only_recs:
                        self.records.append(r)
                self.organize_records()
                self.set_cluster_type()
                self.set_unresolved()

                if self.cluster_type == 'SINGLE_ENDER':
                    self.vcf_record.info['SVTYPE'] = 'BND'
                if self.cluster_type == 'INVERSION_SINGLE_ENDER':
                    self.vcf_record.info['SVTYPE'] = 'BND'
                for r in self.records:
                    r.info['UNRESOLVED'] = True
                    r.info['UNRESOLVED_TYPE'] = self.cpx_type

    def clean_record(self):
        """
        Merge and clean metadata
        """
        sources = set(s for r in self.records for s in r.info['ALGORITHMS'])

        # some variants throw an error if you try to overwrite ALGORITHMS info
        # without removing it via `pop` first. don't remove with `del`, it
        # will break the ability to set ALGORITHMS at all (Invalid INFO field)
        # bcf_update_info: Assertion `!inf->vptr_free' failed.
        self.vcf_record.info.pop('ALGORITHMS')
        self.vcf_record.info['ALGORITHMS'] = ','.join(sorted(sources))

        self.vcf_record.info['MEMBERS'] = tuple(
            sorted(r.id for r in self.records))

        varGQs = []
        for record in self.records:
            if 'varGQ' in record.info.keys():
                varGQs.append(record.info['varGQ'])
        if len(varGQs) > 0 and 'varGQ' in self.vcf_record.header.info.keys():
            self.vcf_record.info['varGQ'] = max(varGQs)

        # if resolved as a CNV, ensure RD_CN and RD_GQ are set
        if len(self.records) > 1 and self.vcf_record.info['SVTYPE'] in ['DEL', 'DUP', 'CNV'] and len(self.cnvs) > 0:
            cnv_record = self.cnvs[0]
            if 'RD_CN' in cnv_record.format.keys() and 'RD_GQ' in cnv_record.format.keys():
                for sample in self.vcf_record.samples:
                    self.vcf_record.samples[sample]['RD_CN'] = cnv_record.samples[sample]['RD_CN']
                    self.vcf_record.samples[sample]['RD_GQ'] = cnv_record.samples[sample]['RD_GQ']

    @property
    def record_ids(self):
        return [r.id for r in self.records]

    def set_unresolved(self):
        self.svtype = 'UNR'
        self.cpx_type = self.cluster_type

    def resolve_inversion(self, SR_only_cutoff):
        if self.inversions[0].info['STRANDS'] == '++':
            FF, RR = self.inversions
        else:
            RR, FF = self.inversions

        self.cpx_type, cnvs = classify_complex_inversion(FF, RR, self.cnvs)
        self.records = [FF, RR] + cnvs

        if self.cpx_type == 'INV':
            self.svtype = 'INV'
        elif self.cpx_type == 'UNK':
            self.svtype = 'UNR'
        elif self.cpx_type == 'COMPLEX_INS':
            self.svtype = 'UNR'
        elif 'INS' in self.cpx_type:
            self.svtype = 'INS'
        else:
            self.svtype = 'CPX'

        # If event is >1kb and both breakpoints are SR-only,
        # leave variant as unresolved
        if 'EVIDENCE' in FF.info.keys() and 'EVIDENCE' in RR.info.keys():
            if FF.info['EVIDENCE'] == ('SR',) and RR.info['EVIDENCE'] == ('SR',):
                if max(FF.stop, RR.stop) - min(FF.pos, RR.pos) > SR_only_cutoff:
                    self.cpx_type = 'UNK'
                    self.svtype = 'UNR'

        # Overall variant start/end
        if self.svtype in ['INV', 'CPX', 'UNR']:
            self.vcf_record.pos = min(FF.pos, RR.pos)
            self.vcf_record.stop = max(FF.stop, RR.stop)

            self.vcf_record.info['SVLEN'] = abs(self.vcf_record.stop -
                                                self.vcf_record.pos)

            if self.svtype != 'UNR':
                cpx_intervals = make_inversion_intervals(FF, RR, self.cnvs,
                                                         self.cpx_type)
                self.vcf_record.info['CPX_INTERVALS'] = cpx_intervals

        elif self.svtype == 'INS':
            #   C B   A D
            # -->|<----|-->
            if self.cpx_type == 'DUP5/INS3':
                source_start, source_end = RR.pos, FF.pos
                sink_start, sink_end = FF.stop, RR.stop

            #   A D   C B
            # -->|<----|-->
            elif self.cpx_type == 'DUP3/INS5':
                source_start, source_end = RR.stop, FF.stop
                sink_start, sink_end = FF.pos, RR.pos

            # first check for overlap with MEI
            is_mei = check_mei_overlap(self.vcf_record.chrom, source_start,
                                       source_end, self.mei_bed)

            if is_mei:
                self.cpx_type = 'MEI_' + self.cpx_type.split('/')[1]

            self.vcf_record.pos = sink_start
            self.vcf_record.stop = sink_end

            # As in MELT, use length of inserted sequence as SVLEN
            self.vcf_record.info['SVLEN'] = abs(source_end - source_start)

            interval = 'INV_{0}:{1}-{2}'
            source = interval.format(self.vcf_record.chrom,
                                     source_start, source_end)
            self.vcf_record.info['SOURCE'] = source

        # Setting alts removes END, so do it up front
        self.vcf_record.alts = ('<{0}>'.format(self.svtype), )
        self.vcf_record.info['SVTYPE'] = self.svtype
        self.vcf_record.info['CPX_TYPE'] = self.cpx_type

    def resolve_translocation(self):
        # Force to ++/-- or +-/-+ ordering
        plus, minus = sorted(self.tlocs, key=lambda t: t.info['STRANDS'])
        armA, armB = get_arms(plus, self.cytobands)

        self.cpx_type = classify_simple_translocation(plus, minus)

        if 'INS' in self.cpx_type:
            self.svtype = 'INS'
        elif self.cpx_type in ['TLOC_MISMATCH_CHROM', 'CTX_UNR']:
            self.svtype = 'UNR'
        elif self.cpx_type == 'CTX_PP/QQ':
            # Don't report sites where posA/posB are identical at each bkpt
            if plus.pos == minus.pos and plus.stop == minus.stop:
                self.svtype = 'UNR'
                self.cpx_type += '_DUPLICATE_COORDS'
            elif armA == armB:
                self.svtype = 'CTX'
            else:
                self.svtype = 'UNR'
                self.cpx_type += '_MISMATCH'
        elif self.cpx_type == 'CTX_PQ/QP':
            if plus.pos == minus.pos and plus.stop == minus.stop:
                self.svtype = 'UNR'
                self.cpx_type += '_DUPLICATE_COORDS'
            elif armA != armB:
                self.svtype = 'CTX'
            else:
                self.svtype = 'UNR'
                self.cpx_type += '_MISMATCH'
        else:
            raise Exception('Invalid cpx type: ' + self.cpx_type)

        # Setting alts removes END, so do it up front
        self.vcf_record.alts = ('<{0}>'.format(self.svtype), )
        self.vcf_record.info['SVTYPE'] = self.svtype
        self.vcf_record.info['CPX_TYPE'] = self.cpx_type

        if self.svtype == 'CTX':
            self.vcf_record.chrom = plus.chrom
            self.vcf_record.pos = plus.pos
            self.vcf_record.info['CHR2'] = plus.info['CHR2']
            self.vcf_record.info['END2'] = plus.stop
            self.vcf_record.stop = plus.stop
            self.vcf_record.info['SVLEN'] = -1

        elif self.svtype == 'INS':
            if 'B2A' in self.cpx_type:
                sink_chrom, source_chrom = plus.chrom, plus.info['CHR2']
            else:
                sink_chrom, source_chrom = plus.info['CHR2'], plus.chrom

            if self.cpx_type == 'CTX_INS_B2A':
                sink_start = plus.pos
                sink_end = minus.pos
                source_start = plus.stop
                source_end = minus.stop
            elif self.cpx_type == 'CTX_INV_INS_B2A':
                sink_start = plus.pos
                sink_end = minus.pos
                source_start = minus.stop
                source_end = plus.stop
            elif self.cpx_type == 'CTX_INS_A2B':
                sink_start = minus.stop
                sink_end = plus.stop
                source_start = minus.pos
                source_end = plus.pos
            elif self.cpx_type == 'CTX_INV_INS_A2B':
                sink_start = plus.stop
                sink_end = minus.stop
                source_start = minus.pos
                source_end = plus.pos

            self.vcf_record.chrom = sink_chrom
            self.vcf_record.pos = sink_start
            self.vcf_record.stop = sink_end

            self.vcf_record.info['CHR2'] = source_chrom
            self.vcf_record.info['END2'] = sink_end
            self.vcf_record.info['SVLEN'] = abs(source_end - source_start)

            interval = '{0}_{1}:{2}-{3}'
            if 'INV' in self.cpx_type:
                interval_type = 'INV'
            else:
                interval_type = 'INS'

            source = interval.format(interval_type, source_chrom,
                                     source_start, source_end)
            self.vcf_record.info['SOURCE'] = source

    def resolve_insertion(self):
        plus, minus = sorted(self.breakends, key=lambda t: t.info['STRANDS'])
        self.cpx_type = classify_insertion(plus, minus)
        if self.cpx_type == 'INS_UNCLASSIFIED':
            self.svtype = 'UNR'
            return
        else:
            self.svtype = 'INS'

        # Setting alts removes END, so do it up front
        self.vcf_record.alts = ('<{0}>'.format(self.svtype), )
        self.vcf_record.info['SVTYPE'] = self.svtype
        self.vcf_record.info['CPX_TYPE'] = self.cpx_type

        sink_chrom = plus.chrom
        source_chrom = plus.chrom

        if self.cpx_type == 'INS_B2A':
            sink_start = plus.pos
            sink_end = minus.pos
            source_start = plus.stop
            source_end = minus.stop
        elif self.cpx_type == 'INS_A2B':
            sink_start = minus.stop
            sink_end = plus.stop
            source_start = minus.pos
            source_end = plus.pos

        self.vcf_record.chrom = sink_chrom
        self.vcf_record.pos = sink_start
        self.vcf_record.stop = sink_end

        self.vcf_record.info['CHR2'] = source_chrom
        self.vcf_record.info['SVLEN'] = abs(source_end - source_start)

        interval = 'INS_{0}:{1}-{2}'
        source = interval.format(source_chrom, source_start, source_end)
        self.vcf_record.info['SOURCE'] = source

        if len(self.insertions) == 1:
            b_algs = self.vcf_record.info['ALGORITHMS']
            mei_algs = self.insertions[0].info['ALGORITHMS']
            algs = tuple(sorted(set(b_algs).union(mei_algs)))
            self.vcf_record = self.insertions[0]
            self.vcf_record.info['ALGORITHMS'] = algs

    def resolve_disdup(self):
        plus, minus = sorted(self.cnvs, key=lambda t: t.info['STRANDS'])
        self.cpx_type = classify_insertion(plus, minus)

        if self.cpx_type == 'INS_UNCLASSIFIED':
            self.svtype = 'UNR'
            return
        else:
            self.svtype = 'INS'

        # Setting alts removes END, so do it up front
        self.vcf_record.alts = ('<{0}>'.format(self.svtype), )
        self.vcf_record.info['SVTYPE'] = self.svtype
        self.vcf_record.info['CPX_TYPE'] = self.cpx_type

        sink_chrom = plus.chrom
        source_chrom = plus.chrom

        if self.cpx_type == 'INS_B2A':
            sink_start = plus.pos
            sink_end = minus.pos
            source_start = plus.stop
            source_end = minus.stop
        elif self.cpx_type == 'INS_A2B':
            sink_start = minus.stop
            sink_end = plus.stop
            source_start = minus.pos
            source_end = plus.pos

        # RLC Note: no longer need this code, as this will now be handled in
        # the complex regenotyping WDL
        # # Don't report insertions with large deletions at insertion site
        # if sink_end - sink_start >= 100:
        #     self.set_unresolved()
        #     return

        self.vcf_record.chrom = sink_chrom
        self.vcf_record.pos = sink_start
        self.vcf_record.stop = sink_end

        self.vcf_record.info['CHR2'] = source_chrom
        self.vcf_record.info['SVLEN'] = abs(source_end - source_start)

        interval = 'INS_{0}:{1}-{2}'
        source = interval.format(source_chrom, source_start, source_end)
        self.vcf_record.info['SOURCE'] = source

        if len(self.insertions) == 1:
            b_algs = self.vcf_record.info['ALGORITHMS']
            mei_algs = self.insertions[0].info['ALGORITHMS']
            algs = tuple(sorted(set(b_algs).union(mei_algs)))
            self.vcf_record = self.insertions[0]
            self.vcf_record.info['ALGORITHMS'] = algs

    def report_simple_insertion(self):
        if len(self.cnvs) == 1:
            if self.cnvs[0].info['SVTYPE'] == 'DUP':
                record = self.cnvs[0]
                self.svtype = 'DUP'
                self.vcf_record.pos = record.pos
                self.vcf_record.stop = record.stop
                self.vcf_record.id = record.id
                self.vcf_record.alts = record.alts
                self.vcf_record.info['SVTYPE'] = self.svtype
                self.vcf_record.info['CHR2'] = record.info['CHR2']
                self.vcf_record.info['SVLEN'] = record.info['SVLEN']
            else:
                record = self.insertions[0]
                self.cpx_type = record.alts[0].strip('<>')
                self.svtype = 'INS'
                self.vcf_record.pos = record.pos
                self.vcf_record.stop = record.stop
                self.vcf_record.alts = record.alts
                self.vcf_record.id = record.id
                self.vcf_record.info['SVTYPE'] = self.svtype
                self.vcf_record.info['CPX_TYPE'] = self.cpx_type
                self.vcf_record.info['CHR2'] = record.info['CHR2']
                self.vcf_record.info['SVLEN'] = record.info['SVLEN']
        else:
            record = self.insertions[0]
            self.cpx_type = record.alts[0].strip('<>')
            self.svtype = 'INS'
            self.vcf_record.pos = record.pos
            self.vcf_record.stop = record.stop
            self.vcf_record.id = record.id
            self.vcf_record.alts = record.alts
            self.vcf_record.info['SVTYPE'] = self.svtype
            self.vcf_record.info['CPX_TYPE'] = self.cpx_type
            self.vcf_record.info['CHR2'] = record.info['CHR2']
            self.vcf_record.info['SVLEN'] = record.info['SVLEN']
        self.vcf_record.info['MEMBERS'] = tuple(r.id for r in self.records)

    # Handles situations where a single insertion clusters with one or more CNVs
    def report_insertion_strip_CNVs(self):
        # Desired behavior: if cluster contains a single duplication, report it
        # and merge INS into the MEMBERS tag for the duplication.
        # Otherwise, mark cluster as svtype = 'SPLIT' and pass entire cluster
        # to resolve.py, where they will be parsed and reported appropriately
        if len(self.cnvs) == 1 and len(self.dups) == 1 \
                and self.cnvs[0].info['SVTYPE'] == 'DUP':
            record = self.cnvs[0]
            self.svtype = 'DUP'
            self.vcf_record.pos = record.pos
            self.vcf_record.stop = record.stop
            self.vcf_record.id = record.id
            self.vcf_record.alts = record.alts
            self.vcf_record.info['SVTYPE'] = self.svtype
            self.vcf_record.info['CHR2'] = record.info['CHR2']
            self.vcf_record.info['SVLEN'] = record.info['SVLEN']
            self.vcf_record.info['MEMBERS'] = tuple(r.id for r in self.records)
        elif len(self.cnvs) > 0:
            self.svtype = 'SPLIT'
        else:
            self.svtype = 'INS'

    # Where Manta calls two insertions flanking a duplication, report just the dup
    def report_manta_tandem_dup(self):
        record = self.dups[0]
        self.cpx_type = record.alts[0].strip('<>')
        self.svtype = 'DUP'
        self.vcf_record.pos = record.pos
        self.vcf_record.stop = record.stop
        self.vcf_record.alts = record.alts
        self.vcf_record.info['SVTYPE'] = self.svtype
        self.vcf_record.info['CPX_TYPE'] = self.cpx_type
        self.vcf_record.info['CHR2'] = record.info['CHR2']
        self.vcf_record.info['SVLEN'] = record.info['SVLEN']

    def set_cluster_type(self):
        # Restrict to double-ended inversion events with appropriate
        # strand pairing
        #  n_invs = len(self.inversions)
        #  n_tlocs = len(self.tlocs)
        #  n_bnds = len(self.breakends)

        class_counts = [len(records) for records in
                        [self.inversions, self.tlocs, self.breakends, self.insertions]]

        paired = np.array([count == 2 for count in class_counts])
        absent = np.array([count == 0 for count in class_counts])

        del_counts = [len(records) for records in [self.dels]]
        dup_counts = [len(records) for records in [self.dups]]

        # If one class is paired and rest are absent
        if all(paired ^ absent) and len(np.where(paired)[0]) == 1 and len(self.cnvs) == 0:
            idx = np.where(paired)[0][0]
            if idx == 0:
                if (self.inversions[0].info['STRANDS'] ==
                        self.inversions[1].info['STRANDS']):
                    self.cluster_type = 'CONGRUENT_INV_STRANDS'
                else:
                    self.cluster_type = 'CANDIDATE_INVERSION'
            elif idx == 1:
                if len(self.cnvs) > 0:
                    self.cluster_type = 'TLOC_WITH_CNV'
                elif ok_tloc_strands(*self.tlocs):
                    self.cluster_type = 'CANDIDATE_TRANSLOCATION'
                else:
                    self.cluster_type = 'STRAND_MISMATCH_TLOC'
            elif idx == 2:
                if len(self.cnvs) > 0:
                    self.cluster_type = 'INS_WITH_CNV'
                elif ok_ins_strands(*self.breakends):
                    self.cluster_type = 'CANDIDATE_INSERTION'
                else:
                    self.cluster_type = 'STRAND_MISMATCH_INS'
            elif idx == 3:
                # Catch case where Manta emits two insertions and a duplication
                if len(self.dups) == 1 and len(self.insertions) == 2:
                    self.cluster_type = 'TANDEM_DUPLICATION_FLANKING_INSERTIONS'
                else:
                    self.cluster_type = 'MULTIPLE_RESOLVED_INSERTIONS'

        elif sum(class_counts) == 0:
            self.cluster_type = 'ERROR_CNV_ONLY'
            if sum(del_counts) == 1 and sum(dup_counts) == 1:
                self.cluster_type = 'CANDIDATE_DISDUP'
        elif sum(class_counts) == 1:
            if len(self.insertions) == 1:
                if len(self.cnvs) > 0:
                    self.cluster_type = 'SINGLE_INSERTION_STRIP_CNVS'
                else:
                    self.cluster_type = 'RESOLVED_INSERTION'
            else:
                if len(self.cnvs) > 0:
                    self.cluster_type = 'MIXED_BREAKENDS'
                elif len(self.inversions) > 0:
                    self.cluster_type = 'INVERSION_SINGLE_ENDER_' + \
                        self.inversions[0].info['STRANDS']
                else:
                    self.cluster_type = 'SINGLE_ENDER'
        elif sum(class_counts) >= 2:
            if len(self.insertions) >= 1 and len(self.breakends) >= 1:
                self.cluster_type = 'RESOLVED_INSERTION'
            else:
                self.cluster_type = 'MIXED_BREAKENDS'
        else:
            self.cluster_type = 'ERROR_UNCLASSIFIED'

    def make_record(self):
        self.vcf_record = self.records[0].copy()
        if len(self.records) > 1:
            svu.update_best_genotypes(
                self.vcf_record, self.records, preserve_multiallelic=False)


def ok_tloc_strands(tloc1, tloc2):
    strand1, strand2 = sorted([t.info['STRANDS'] for t in (tloc1, tloc2)])

    return ((strand1 == '++' and strand2 == '--') or
            (strand1 == '+-' and strand2 == '-+'))


def ok_ins_strands(bnd1, bnd2):
    strand1, strand2 = sorted([t.info['STRANDS'] for t in (bnd1, bnd2)])

    return (strand1 == '+-') and (strand2 == '-+')


def get_arms(record, cytobands):
    regionA = '{0}:{1}-{1}'.format(record.chrom, record.pos)
    regionB = '{0}:{1}-{1}'.format(record.info['CHR2'], record.stop)

    def _get_arm(region):
        arm = next(cytobands.fetch(region))
        return arm.split()[3][0]

    return _get_arm(regionA), _get_arm(regionB)


def make_inversion_intervals(FF, RR, cnvs, cpx_type):
    intervals = []
    chrom = FF.chrom

    interval = '{svtype}_{chrom}:{start}-{end}'

    # First add 5' CNV
    if cpx_type.startswith('del'):
        svtype = 'DEL'
        start = FF.pos
        end = RR.pos
        intervals.append(interval.format(**locals()))

    if cpx_type.startswith('dup'):
        svtype = 'DUP'
        start = RR.pos
        end = FF.pos
        intervals.append(interval.format(**locals()))

    # Then add inversion
    svtype = 'INV'
    start = RR.pos
    end = FF.stop
    intervals.append(interval.format(**locals()))

    # Finally add 3' CNV
    if cpx_type.endswith('del'):
        svtype = 'DEL'
        start = FF.stop
        end = RR.stop
        intervals.append(interval.format(**locals()))

    if cpx_type.endswith('dup'):
        svtype = 'DUP'
        start = RR.stop
        end = FF.stop
        intervals.append(interval.format(**locals()))

    return intervals


def check_mei_overlap(chrom, start, end, mei_bed):
    """
    Check if putative insertion is covered by MEIs
    """

    bedline = '{0}\t{1}\t{2}\n'.format(chrom, start, end)
    bed = pbt.BedTool(bedline, from_string=True)

    cov_bed = bed.coverage(mei_bed).saveas()
    i = next(cov_bed.intervals)
    cov = float(i.fields[6])

    return cov >= 0.5
