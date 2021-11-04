# -*- coding: utf-8 -*-
#

"""

"""

import pandas as pd
import pysam
from .annotate_intersection import annotate_intersection
from .classify_effect import classify_effect
from .nearest_tss import annotate_nearest_tss
import svtk.utils as svu


def annotate_gencode(sv, gencode):
    # Annotate Gencode hits and predicted effects
    hits = annotate_intersection(sv, gencode, filetype='gtf')
    effects = classify_effect(hits)
    effects = effects.loc[effects.effect != 'GENE_OTHER'].copy()

    # Annotate nearest TSS
    tss = annotate_nearest_tss(sv, gencode)

    # Only include TSS if no genic hit observed
    if len(effects) > 0:
        tss = tss.loc[~tss.name.isin(effects.name)].copy()

    # Merge annotations
    effects = pd.concat([effects, tss])
    return effects


def annotate_noncoding(sv, noncoding):
    # Annotate noncoding hits.
    # For now, any overlap gets annotated
    noncoding_hits = annotate_intersection(sv, noncoding, filetype='bed')
    noncoding_hits = noncoding_hits.drop_duplicates()

    if noncoding_hits.shape[0] > 0:
        if len(noncoding_hits.hit_type == 'SPAN') != 0:
            noncoding_hits.loc[noncoding_hits.hit_type ==
                               'SPAN', 'effect'] = 'NONCODING_SPAN'
        if len(noncoding_hits.hit_type != 'SPAN') != 0:
            noncoding_hits.loc[noncoding_hits.hit_type !=
                               'SPAN', 'effect'] = 'NONCODING_BREAKPOINT'
        noncoding_cols = 'name svtype gene_name effect'.split()
        effects = noncoding_hits[noncoding_cols].drop_duplicates()
    else:
        # no hits found, add an empty 'effect' column
        noncoding_hits['effect'] = ''
        effects = noncoding_hits

    return effects


def annotate(sv, gencode, noncoding):
    """
    Gene disrupting if:
       1) Breakpoint hits exon
       2) Inverted segment lies inside gene AND hits exon
       3) One end of inverted segment lies inside gene AND other does not

    Notes:
        * Intragenic variants are annotated with the nearest TSS.
        * Inversions which entirely span a gene are considered "no_effect" and
        excluded from annotation
        * Duplications which overlap the start OR end of a gene boundary
        without fully spanning the boundary are considered "no_effect" and
        excluded from annotation

    Parameters
    ----------
    sv : pbt.BedTool
        SV breakpoints and CNV intervals
    gencode : pbt.BedTool
        Gencode annotations
    """

    coding_anno = None if gencode is None else annotate_gencode(sv, gencode)
    noncoding_anno = None if noncoding is None else annotate_noncoding(
        sv, noncoding)

    effects = pd.concat([coding_anno, noncoding_anno])

    # Aggregate genic effects by variant ID
    effects = effects.pivot_table(index='name',
                                  columns='effect',
                                  values='gene_name',
                                  aggfunc=lambda s: ','.join(sorted(set(s))),
                                  fill_value='NA')

    return effects


GENCODE_INFO = [
    '##INFO=<ID=LOF,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a loss-of-function effect.">',
    '##INFO=<ID=DUP_LOF,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a loss-of-function effect via intragenic exonic duplication.">',
    '##INFO=<ID=COPY_GAIN,Number=.,Type=String,Description="Gene(s) on which the SV is predicted to have a copy-gain effect.">',
    '##INFO=<ID=DUP_PARTIAL,Number=.,Type=String,Description="Gene(s) which are partially overlapped by an SV\'s duplication, such that an unaltered copy is preserved.">',
    '##INFO=<ID=MSV_EXON_OVR,Number=.,Type=String,Description="Gene(s) on which the multiallelic SV would be predicted to have a LOF, DUP_LOF, COPY_GAIN, or DUP_PARTIAL annotation if the SV were biallelic.">',
    '##INFO=<ID=INTRONIC,Number=.,Type=String,Description="Gene(s) where the SV was found to lie entirely within an intron.">',
    '##INFO=<ID=INV_SPAN,Number=.,Type=String,Description="Gene(s) which are entirely spanned by an SV\'s inversion.">',
    '##INFO=<ID=UTR,Number=.,Type=String,Description="Gene(s) for which the SV is predicted to disrupt a UTR.">',
    '##INFO=<ID=NEAREST_TSS,Number=.,Type=String,Description="Nearest transcription start site to intragenic variants.">',
    '##INFO=<ID=promoter,Number=0,Type=Flag,Description="Nearest promoter.">',
    '##INFO=<ID=INTERGENIC,Number=0,Type=Flag,Description="SV does not overlap coding sequence.">'
]

NONCODING_INFO = [
    '##INFO=<ID=NONCODING_SPAN,Number=.,Type=String,Description="Classes of noncoding elements spanned by SV.">',
    '##INFO=<ID=NONCODING_BREAKPOINT,Number=.,Type=String,Description="Classes of noncoding elements disrupted by SV breakpoint.">',
]


def annotate_vcf(vcf, gencode, noncoding, annotated_vcf):
    """
    Parameters
    ----------
    vcf : pysam.VariantFile
    gencode : pbt.BedTool
        Gencode gene annotations
    noncoding : pbt.BedTool
        Noncoding elements
    annotated_vcf : str
        Path to output VCF
    """

    # Add metadata lines for annotations
    header = vcf.header

    if gencode is not None:
        for line in GENCODE_INFO:
            header.add_line(line)
    if noncoding is not None:
        for line in NONCODING_INFO:
            header.add_line(line)

    # Open output file
    fout = pysam.VariantFile(annotated_vcf, 'w', header=header)

    # Annotate genic hits
    if isinstance(vcf.filename, bytes):
        fname = vcf.filename.decode()
    else:
        fname = vcf.filename
    sv = svu.vcf2bedtool(fname, split_bnd=True, split_cpx=True,
                         simple_sinks=True, include_unresolved=False)

    effects = annotate(sv, gencode, noncoding)
    effects = effects.to_dict(orient='index')
    # Add results to variant records and save
    for record in vcf:
        anno = effects.get(record.id)
        if anno is None:
            fout.write(record)
            continue

        # Handle general catch-all intersection for MULTIALLELIC variants
        if 'MULTIALLELIC' in record.filter:
            multi_ovr = []
            for info, genelist in anno.items():
                if info in 'LOF DUP_LOF COPY_GAIN DUP_PARTIAL'.split():
                    if genelist != 'NA':
                        for gene in genelist.split(','):
                            if gene not in multi_ovr:
                                multi_ovr.append(gene)
                else:
                    if genelist != 'NA':
                        record.info[info] = genelist
            if len(multi_ovr) > 0:
                record.info['MSV_EXON_OVR'] = ','.join(multi_ovr)
        else:
            for info, genelist in anno.items():
                if genelist != 'NA':
                    record.info[info] = genelist

        if 'NEAREST_TSS' in record.info:
            record.info['INTERGENIC'] = True

        fout.write(record)

    fout.close()
