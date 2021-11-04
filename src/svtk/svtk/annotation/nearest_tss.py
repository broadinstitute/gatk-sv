#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import warnings
from .annotate_intersection import split_gencode_fields


def annotate_nearest_tss(sv, gencode):
    """
    Annotate each variant record with its nearest TSS.

    Parameters
    ----------
    sv : pbt.BedTool
        columns = (chrom, start, end, name, svtype, strands)
    gencode : pbt.BedTool
        Gencode gene annotations (GTF) of 9 columns

    Returns
    -------
    nearest_tss : pd.DataFrame
        Columns = (name, svtype, gene_name, effect)
    """

    def _make_tss(feature):
        strand = feature.fields[6]
        if strand == '-':
            feature.start = feature.end

        feature.end = feature.start + 1
        return feature

    transcripts = gencode.filter(lambda r: r.fields[2] == 'transcript')
    tss = transcripts.each(_make_tss).saveas()

    nearest_tss = sv.sort().closest(tss.sort(), D='b', id=True).saveas()

    sv_names = ['chrom', 'start', 'end', 'name', 'svtype', 'strand']
    # A fix by SHuang on March 1st 2019:
    #   I really don't know how to provide sensible names to tss
    tss_names = ['chr2', 'id2', 'class', 'start2', 'end2', 'dot', 'strand2',
                 'another_dot', 'gencode_fields']
    df_names = sv_names + tss_names + ['distance']
    # Suppress warning about column names
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        nearest_tss = nearest_tss.to_dataframe(disable_auto_names=True,
                                               names=df_names)

    def _get_gene_name(field_str):
        fields = split_gencode_fields(field_str)
        return fields['gene_name']

    nearest_tss['gene_name'] = nearest_tss['gencode_fields'].apply(
        _get_gene_name)
    nearest_tss = nearest_tss[['name', 'svtype', 'gene_name']].copy()
    nearest_tss['effect'] = 'NEAREST_TSS'

    return nearest_tss
