#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Identify the ENSEMBL canonical transcript for each protein-coding Gencode gene.

ENSEMBL determines the canonical transcript based on the following hierarchy
(definition obtained from http://www.ensembl.org/Help/Glossary?id=346):
    1) Longest CCDS translation with no stop codons.
    2) If no (1), choose the longest Ensembl/Havana merged translation with no
       stop codons.
    3) If no (2), choose the longest translation with no stop codons.
    4) If no translation, choose the longest non-protein-coding transcript.

Known to work on Gencode v19 files.
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/
    gencode.v19.annotation.gtf.gz
    gencode.v19.pc_translations.fa.gz
    gencode.v19.pc_transcripts.fa.gz
    gencode.v19.metadata.Transcript_source.gz
"""

import argparse
import gzip
from contextlib import ExitStack
import pandas as pd


def select_canonical_transcript(txs):
    """
    Select ENSEMBL-defined canonical transcript.

    Parameters
    ----------
    txs : pd.DataFrame
        Group of transcripts associated with a single gene ID.

    Returns
    -------
    canonical_tx : pd.Series
        Canonical transcript.
    """

    # Choose longer translation, with transcript length as tiebreaker
    sizes = 'translation_length transcript_length'.split()

    # 4) Longest non-protein-coding transcript
    translated = txs.loc[txs.transcript_type == 'protein_coding']
    if translated.shape[0] == 0:
        return txs.sort_values('seq_length', ascending=False).iloc[0]

    # 1) Longest CCDS translation with no stop codons
    CCDS = translated.loc[translated.tags.str.contains('CCDS')]
    if CCDS.shape[0] > 0:
        return CCDS.sort_values(sizes, ascending=False).iloc[0]

    # 2) Longest Ensembl/Havana merged translation
    ensembl_havana = translated.loc[
        (translated.transcript_id == 'ensembl_havana_transcript')]
    if ensembl_havana.shape[0] > 0:
        return ensembl_havana.sort_values(sizes, ascending=False).iloc[0]

    # 2b) Longest Ensembl/Havana merged translation - "proj" prefix
    ensembl_havana = translated.loc[
        (translated.transcript_id == 'proj_ensembl_havana_transcript')]
    if ensembl_havana.shape[0] > 0:
        return ensembl_havana.sort_values(sizes, ascending=False).iloc[0]

    # 3) Longest translation
    return translated.sort_values(sizes, ascending=False).iloc[0]


def parse_gencode_transcripts(gtffile):
    """
    Parse metadata for all Gencode transcripts of protein-coding genes.

    Parameters
    ----------
    gtffile : file
        Gencode annotation GTF

    Returns
    -------
    metadata : pd.DataFrame
        Columns = (gene_id, transcript_id, gene_type, gene_status, gene_name,
                   transcript_type, transcript_status, transcript_name)
    """

    def _parse():
        for line in gtffile:
            if line.startswith('#'):
                continue

            data = line.strip().split('\t')

            # Report end - start as length of element sequence
            # GTF is closed interval, so add 1
            seq_length = int(data[4]) - int(data[3]) + 1

            # gene/transcript/exon/CDS/UTR/start_codon/stop_codon
            element = data[2]

            # Annotation fields
            fields = data[8].strip(';').split('; ')

            # First 8 fields are consistent metadata across all records
            field_tuples = [f.split() for f in fields[:8]]
            field_dict = {f[0]: f[1].strip('"') for f in field_tuples}

            # Parse tags from remaining fields
            tags = []
            for f in fields[8:]:
                if '=' in f:
                    continue
                field, val = f.split()
                if field == 'tag':
                    tags.append(val.strip('"'))
            tags = ','.join(tags) if len(tags) > 0 else '.'

            # Restrict to transcripts of protein-coding genes
            if (element == 'transcript' and
                    field_dict['gene_type'] == 'protein_coding'):

                entries = [f[1].strip('"') for f in field_tuples]
                entries = entries + [seq_length, tags]

                yield entries

    columns = ('gene_id transcript_id gene_type gene_status gene_name '
               'transcript_type transcript_status transcript_name '
               'seq_length tags'.split())

    return pd.DataFrame.from_records(_parse(), columns=columns)


def parse_gencode_fasta(fastafile, len_type=None):
    """
    Read metadata from record descriptions in Gencode fasta.

    Expects following Gencode format for protein-coding translation and
    protein-coding transcript sequence fastas.

    >tx_id|gene_id|havana_gene|havana_tx|tx_name|gene_name|seq_length

    Parameters
    ----------
    fastafile : file
    len_type : str, optional
        Prefix to add to "length" column

    Returns
    -------
    metadata : pd.DataFrame
        Columns = (transcript_id, gene_id, havana_gene, havana_transcript,
                   transcript_name, gene_name, length)
    """

    # Parse first 7 "|"-delimited entries from each fasta record description
    def _parse():
        for line in fastafile:
            if not line.startswith('>'):
                continue

            desc = line.strip().strip('>').split('|')

            yield tuple(desc[:7])

    # Add transcript/translation prefix to length column title if desired
    columns = ('transcript_id gene_id havana_gene havana_transcript '
               'transcript_name gene_name').split()
    if len_type is None:
        columns = columns + ['length']
    else:
        columns = columns + [len_type + '_length']

    return pd.DataFrame.from_records(_parse(), columns=columns)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('annotation_gtf', help='Gencode gene annotations.')
    parser.add_argument('pc_translations_fa', help='Gencode protein-coding '
                        'transcript translation sequences. Parsed for '
                        'translation lengths.')
    parser.add_argument('pc_transcripts_fa', help='Gencode protein-coding '
                        'transcript sequences. Parsed for transcript lengths.')
    parser.add_argument('transcript_source', help='Gencode transcript source '
                        'metadata.')
    parser.add_argument('canonical_transcripts', help='Tab-delimited table of '
                        'canonical transcripts.')
    args = parser.parse_args()

    with ExitStack() as cm:
        # Parse translation metadata and lengths
        tl_fasta = cm.enter_context(gzip.open(args.pc_translations_fa, 'rt'))
        translations = parse_gencode_fasta(tl_fasta, len_type='translation')

        # Parse transcript metadata and lengths
        tx_fasta = cm.enter_context(gzip.open(args.pc_transcripts_fa, 'rt'))
        transcripts = parse_gencode_fasta(tx_fasta, len_type='transcript')

        # Parse transcript sources
        names = 'transcript_id transcript_source'.split()
        sources = pd.read_table(args.transcript_source, names=names)

        # Merge translation/transcript metadata and sources into one dataframe
        tx_cols = 'transcript_id transcript_length'.split()
        metadata = pd.merge(translations, transcripts[tx_cols],
                            on='transcript_id', how='left')
        metadata = pd.merge(metadata, sources,
                            on='transcript_id', how='left')

        # Parse all annotated transcripts of protein-coding genes
        anno_gtf = cm.enter_context(gzip.open(args.annotation_gtf, 'rt'))
        transcripts = parse_gencode_transcripts(anno_gtf)

        # Merge length and source info into set of transcripts
        metadata_cols = ('transcript_id translation_length transcript_length '
                         'transcript_source').split()
        transcripts = pd.merge(transcripts, metadata[metadata_cols],
                               on='transcript_id', how='left')

    #  transcripts.to_csv('transcript_data2.csv', index=False)

    # Filter to canonical transcripts and save
    canon = transcripts.groupby('gene_id')\
                       .apply(select_canonical_transcript)
    canon.to_csv(args.canonical_transcripts, index=False, sep='\t',
                 na_rep='NA', float_format='%0.f')


if __name__ == '__main__':
    main()
