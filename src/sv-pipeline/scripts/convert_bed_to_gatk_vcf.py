#!/bin/env python

import argparse
import gzip
import pysam
import sys

from operator import attrgetter
from typing import Optional, List, Text, Set, Dict, IO, Iterable

# Sort batches with matching contig/start
BATCH_SORT_SPEC = [
    ('contig', False),
    ('start', False),
    ('end', False),
    ('svtype', False)
]


# Sorts list xs by specified attributes
def multisort(xs, specs):
    for key, reverse in reversed(specs):
        xs.sort(key=attrgetter(key), reverse=reverse)
    return xs


def create_header(contigs: Set[Text],
                  samples: Iterable[Text]) -> pysam.VariantHeader:
    """
    Creates header for GATK CNVs froms scratch.

    Parameters
    ----------
    contigs: Set[Text]
        set of reference contigs
    samples: Iterable[Text]
        set of sample names

    Returns
    -------
    header: pysam.VariantHeader
        gatk-style header
    """
    header = pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    for sample in samples:
        header.add_sample(sample)
    for contig in contigs:
        header.add_line(f"##contig=<ID={contig}>")
    header.add_line('##ALT=<ID=DEL,Description="Deletion">')
    header.add_line('##ALT=<ID=DUP,Description="Duplication">')
    header.add_line('##ALT=<ID=CNV,Description="Copy number variant">')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">')
    header.add_line('##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of affected segment on the reference">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">')
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    header.add_line('##source=depth')
    return header


class CNVRecord:

    def __init__(self,
                 bed_record: Text):
        """
        Initializes from a bed record

        Parameters
        ----------
        record: Text
            BED record
        ploidy_dict: Dict[Text, Dict[Text, int]]
            map from sample to contig to ploidy

        Returns
        -------
        header: pysam.VariantRecord
            gatk-style record
        """
        cols = bed_record.strip().split('\t')
        self.contig = cols[0]
        self.start = int(cols[1])
        self.end = int(cols[2]) + 1
        self.name = cols[3]
        self.sample = cols[4]
        self.svtype = cols[5]


class CNVWriter:

    def __init__(self,
                 bed_path: Text,
                 vcf_path: Text,
                 contigs: Set[Text],
                 samples: Iterable[Text],
                 vid_prefix: Text,
                 ploidy_dict: Dict[Text, Dict[Text, int]],
                 contig_len_dict: Dict[Text, int]):
        self.bed = self.__open_bed(bed_path)
        header = create_header(contigs=contigs, samples=samples)
        self.vcf_out = pysam.VariantFile(vcf_path, mode='w', header=header)
        self.vid_index = 0
        self.vid_prefix = vid_prefix
        self.ploidy_dict = ploidy_dict
        self.contig_len_dict = contig_len_dict

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.bed.close()
        self.vcf_out.close()

    def __open_bed(self, path: Text) -> IO[Text]:
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        else:
            return open(path, 'r')

    def collapse_batch(self, batch: List[CNVRecord]):
        contig = batch[0].contig
        start = max(1, batch[0].start)
        end = batch[0].end
        # Clip interval if enabled
        if self.contig_len_dict:
            end = min(self.contig_len_dict[contig], end)
        svtype = batch[0].svtype
        if svtype != 'DEL' and svtype != 'DUP' and svtype != 'CNV':
            raise ValueError(f"Unsupported SV type: {svtype}")
        alleles = ('N', '<' + svtype + '>')
        vcf_record = self.vcf_out.new_record(contig=contig, start=start, stop=end, alleles=alleles)
        vcf_record.info['SVTYPE'] = svtype
        vcf_record.info['ALGORITHMS'] = "depth"
        carriers = set(r.sample for r in batch)
        for sample in self.vcf_out.header.samples:
            ploidy = self.ploidy_dict[sample][contig]
            if ploidy > 2:
                raise ValueError(f"Unsupported {contig} ploidy for sample {sample}: {ploidy}")
            genotype = vcf_record.samples[sample]
            genotype['ECN'] = ploidy
            if ploidy == 0:
                genotype['GT'] = ()
                genotype['CN'] = 0
                continue
            if sample in carriers:
                if ploidy == 1:
                    genotype['GT'] = (1,)
                    if svtype == 'DEL':
                        genotype['CN'] = 0
                    else:
                        # Treat CNV and DUP the same
                        genotype['CN'] = 2
                elif ploidy == 2:
                    # Mark as het for carrier status, but we don't know the actual genotype from the bed file
                    genotype['GT'] = (0, 1)
                    if svtype == 'DEL':
                        genotype['CN'] = 1
                    else:
                        # Treat CNV and DUP the same
                        genotype['CN'] = 3
            else:
                genotype['CN'] = ploidy
                if ploidy == 1:
                    genotype['GT'] = (0,)
                elif ploidy == 2:
                    genotype['GT'] = (0, 0)
        return vcf_record

    def collapse_and_write_batch(self, batch: List[CNVRecord]):
        sorted_batch = multisort(xs=batch, specs=BATCH_SORT_SPEC)
        current_group_spec = None
        current_group = []
        for r in sorted_batch:
            group_spec = (r.contig, r.start, r.end, r.svtype)
            if current_group_spec != group_spec:
                if len(current_group) > 0:
                    record = self.collapse_batch(batch=current_group)
                    record.id = f"{self.vid_prefix}{self.vid_index}"
                    self.vid_index += 1
                    self.vcf_out.write(record)
                current_group_spec = group_spec
                current_group = [r]
            else:
                current_group.append(r)
        if len(current_group) > 0:
            record = self.collapse_batch(batch=current_group)
            record.id = f"{self.vid_prefix}{self.vid_index}"
            self.vid_index += 1
            self.vcf_out.write(record)

    def convert_bed_to_vcf(self):
        last_contig = None
        last_start = None
        batch = []  # batch of identical sites (across multiple samples)
        for line in self.bed:
            if line.startswith("#"):
                continue
            record = CNVRecord(line)
            if last_contig == record.contig and last_start == record.start:
                batch.append(record)
            else:
                self.collapse_and_write_batch(batch=batch)
                last_contig = record.contig
                last_start = record.start
                batch = [record]
        self.collapse_and_write_batch(batch=batch)


def __read_list(path: Text) -> List[Text]:
    with open(path, 'r') as f:
        return [line.strip().split('\t')[0] for line in f]


def __parse_ploidy_table(path: Text) -> Dict[Text, Dict[Text, int]]:
    """
    Parses tsv of sample ploidy values.

    Parameters
    ----------
    path: Text
        table path

    Returns
    -------
    header: Dict[Text, Dict[Text, int]]
        map of sample to contig to ploidy, i.e. Dict[sample][contig] = ploidy
    """
    ploidy_dict = dict()
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            tokens = line.strip().split('\t')
            sample = tokens[0]
            ploidy_dict[sample] = {header[i]: int(tokens[i]) for i in range(1, len(header))}
    return ploidy_dict


def __parse_contig_lengths(path: Text) -> Dict[Text, int]:
    """
    Parses fai file for contig lengths.

    Parameters
    ----------
    path: Text
        fai path

    Returns
    -------
    header: Dict[Text, int]
        map of contig name to length
    """
    contig_len_dict = dict()
    with open(path, 'r') as f:
        for line in f:
            tokens = line.strip().split('\t')
            contig = tokens[0]
            length = int(tokens[1])
            contig_len_dict[contig] = length
    return contig_len_dict


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Convert a BED file of CNVs to a GATK-style VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bed", type=str, required=True,
                        help="BED file containing CNVs with columns: chr, start, end, name, sample, svtype, sources")
    parser.add_argument("--contigs", type=str, required=True,
                        help="List of contigs")
    parser.add_argument("--samples", type=str, required=True,
                        help="List of samples")
    parser.add_argument("--vid-prefix", type=str, required=True,
                        help="Variant ID prefix")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--ploidy-table", type=str, required=True,
                        help="Tab-delimited table of sample ploidies. The table should have a header row where the "
                             "first column is SAMPLE, and the remaining columns are contig names. For each row "
                             "thereafter, the first column is the sample name, and remaining columns are the contig "
                             "ploidy values for that sample.")
    parser.add_argument("--fai", type=str, help="If provided, clips coordinates to [1, contig_len]")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    if parsed_arguments.fai and not parsed_arguments.fai.endswith(".fai"):
        raise ValueError("--fai file must end with .fai")
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    contigs = __read_list(arguments.contigs)
    samples = __read_list(arguments.samples)
    ploidy_dict = __parse_ploidy_table(arguments.ploidy_table)
    if arguments.fai:
        contig_len_dict = __parse_contig_lengths(arguments.fai)
    else:
        contig_len_dict = None
    with CNVWriter(bed_path=arguments.bed,
                   vcf_path=arguments.out,
                   samples=samples,
                   vid_prefix=arguments.vid_prefix,
                   contigs=contigs,
                   ploidy_dict=ploidy_dict,
                   contig_len_dict=contig_len_dict) as writer:
        writer.convert_bed_to_vcf()


if __name__ == "__main__":
    main()
