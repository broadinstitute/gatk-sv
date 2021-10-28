#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Annotate a VCF of structural variants with a list of genomic elements.
"""

import pandas as pd


def split_gencode_fields(field_str):
    fields = field_str.strip(';').split('; ')
    fields = [x.split(' ') for x in fields]
    fields = dict(map(lambda s: s.strip('"'), f) for f in fields)
    return fields


def intersection_type(variant, element, filetype='bed'):
    """
    list of str
        chrom, start, end, ID
    """

    variant_start, variant_end = [int(x) for x in variant[1:3]]

    if filetype == 'bed':
        element_start, element_end = [int(x) for x in element[1:3]]
    elif filetype == 'gtf':
        element_start, element_end = [int(x) for x in element[3:5]]

    if variant_start > element_start and variant_end < element_end:
        return 'BOTH-INSIDE'
    elif variant_start > element_start and variant_start < element_end:
        return 'ONE-INSIDE'
    elif variant_end > element_start and variant_end < element_end:
        return 'ONE-INSIDE'
    else:
        return 'SPAN'


def disruption_type(hit_type, svtype):
    disruptions = {
        'DEL': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
        'DUP': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'COPY_GAIN'},
        'CNV': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'COPY_GAIN'},
        'INV': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'SPAN'},
        'BND': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
        'CTX': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
        'ANEUPLOIDY': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
        'CPX': {
            'BOTH-INSIDE': 'DISRUPTING',
            'ONE-INSIDE': 'DISRUPTING',
            'SPAN': 'DISRUPTING'},
    }

    return disruptions[svtype][hit_type]


def annotate_intersection(sv, elements, filetype='gtf'):
    """
    Parameters
    ----------
    sv : pbt.BedTool
        SV breakpoints and CNV intervals
    gencode : pbt.BedTool
        Gencode annotations
    """

    # Number of fields in SV bedtool
    N_BED_FIELDS = 6

    # Check intersection with gene boundaries
    sect = sv.intersect(elements, wa=True, wb=True)

    def _annotate():
        for hit in sect.intervals:
            variant = hit.fields[:N_BED_FIELDS]
            variant_ID = variant[3]
            svtype = variant[4]

            # Noncoding data
            element = hit.fields[N_BED_FIELDS:]

            # Gencode data
            if filetype == 'gtf':
                fields = split_gencode_fields(element[-1])
                gene_name = fields['gene_name']
                element_type = element[2]
            # Noncoding elements
            else:
                gene_name = element[3]
                element_type = 'noncoding'

            hit_type = intersection_type(variant, element, filetype)
            disrupt_type = disruption_type(hit_type, svtype)

            yield (variant_ID, svtype, gene_name, element_type, hit_type,
                   disrupt_type)

    columns = 'name svtype gene_name element_type hit_type disrupt_type'
    columns = columns.split()
    hits = pd.DataFrame.from_records(_annotate(), columns=columns)

    return hits
