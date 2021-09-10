#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


def choose_best_genotype(sample, records):
    """
    Return record where sample has best non-reference genotype, if available

    Parameters
    ----------
    sample : str
    records : list of pysam.VariantRecord

    Returns
    -------
    best_record : pysam.VariantRecord
    """

    best_GT = (0, 0)
    best_GQ = 0
    best_record = None

    if not any(['GQ' in record.samples[sample].keys() for record in records]):
        for record in records:
            if record.samples[sample]['GT'] != (0, 0):
                return record
        return records[0]

    # Pick best non-reference genotype
    # for record in records:
    #    if record.samples[sample]['GQ'] >= best_GQ:
        # if record is non-reference , use it
        # or if it's a higher GQ for a reference call, use it
    #        if record.samples[sample]['GT'] != (0, 0) or best_GT == (0, 0):
    #            best_GT = record.samples[sample]['GT']
    #            best_GQ = record.samples[sample]['GQ']
    #            best_record = record

    for record in records:
        # if found non-ref GT, replace GT and GQ together
        if record.samples[sample]['GT'] != (0, 0) and best_GT == (0, 0):
            best_GT = record.samples[sample]['GT']
            best_GQ = record.samples[sample]['GQ']
            best_record = record
        elif record.samples[sample]['GT'] == (0, 0) and best_GT != (0, 0):
            continue
        # if new GT  = best_GT, while found a higher GQ, replace GQ
        else:
            if record.samples[sample]['GQ'] >= best_GQ:
                best_GQ = record.samples[sample]['GQ']
                best_record = record

    return best_record


def check_multiallelic(records):
    """
    Returns true if any record is multiallelic

    Parameters
    ----------
    records : list of pysam.VariantRecord
    """
    for record in records:
        if record.alts[0] == '<CN0>':
            return True
        #  for sample in record.samples:
            #  GT = record.samples[sample]['GT']
            #  if GT[0] > 1 or GT[1] > 1:
            #  return True

    return False


def make_multiallelic_alts(records):
    """
    Make list of alts corresponding to record with highest observed copy number

    Parameters
    ----------
    records : list of pysam.VariantRecord

    Returns
    -------
    alts : tuple of str
    """

    max_CN = 2

    is_bca = records[0].info['SVTYPE'] not in 'DEL DUP'.split()

    for record in records:
        if record.alts[0] == '<CN0>':
            CN = int(record.alts[-1].strip('<CN>'))
            if CN > max_CN:
                max_CN = CN

    if is_bca:
        return tuple(['<CN1>'] + ['<CN%d>' % i for i in range(1, max_CN + 1)])
    else:
        return tuple(['<CN0>'] + ['<CN%d>' % i for i in range(1, max_CN + 1)])


def update_best_genotypes(new_record, records, preserve_multiallelic=False):
    """
    For each sample in record, update GT and other formats with best genotype

    Parameters
    ----------
    new_record : pysam.VariantRecord
    records : list of SVRecord
    preserve_multiallelic : bool
        If any record is multiallelic, make all genotypes multiallelic
    """

    if preserve_multiallelic:
        is_multiallelic = check_multiallelic(records)
    else:
        is_multiallelic = False

    if is_multiallelic:
        new_record.alts = make_multiallelic_alts(records)

    new_record_formats = new_record.header.formats.keys()
    for sample in new_record.samples:
        best_record = choose_best_genotype(sample, records)
        for key in best_record.format.keys():
            if key in new_record_formats:
                # If any record is multiallelic, replace non-multiallelic
                # genotypes with multiallelic equivalents
                if key == 'GT':
                    GT = best_record.samples[sample][key]
                    if is_multiallelic:
                        if GT == (0, 1):
                            new_record.samples[sample][key] = (0, 2)
                        elif GT == (1, 1):
                            new_record.samples[sample][key] = (2, 2)
                        else:
                            new_record.samples[sample][key] = GT
                    else:
                        GT = tuple([min(x, 1) for x in GT])
                        new_record.samples[sample][key] = GT
                elif key == 'EV':
                    if 'EV' in new_record.samples[sample].keys():
                        curr_ev = new_record.samples[sample][key]
                    else:
                        curr_ev = 0

                    if curr_ev is None:
                        curr_ev = 0

                    for record in records:
                        if 'EV' in record.samples[sample].keys():
                            ev = record.samples[sample][key]
                            if ev is None:
                                ev = 0
                            curr_ev = curr_ev | ev
                    new_record.samples[sample][key] = curr_ev
                else:
                    new_record.samples[sample][key] = best_record.samples[sample][key]
