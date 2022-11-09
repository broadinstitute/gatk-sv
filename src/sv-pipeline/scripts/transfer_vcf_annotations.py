#!/bin/env python

import argparse
import re
import pysam
import sys
from typing import Optional, List, Text, Set, Dict, Tuple
from collections import deque


def create_header(header_in: pysam.VariantHeader,
                  header_ann: pysam.VariantHeader,
                  transfer_infos: Set[Text],
                  transfer_formats: Set[Text]) -> pysam.VariantHeader:
    """
    Adds annotation metadata lines to the input header. Also checks that the fields exist in the annotation vcf.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header
    header_ann: pysam.VariantHeader
        annotation vcf header
    transfer_infos: Set[Text]
        set of info fields to annotate
    transfer_formats: Set[Text]
        set of format fields to annotate

    Returns
    -------
    header: pysam.VariantHeader
        output header
    """
    header = pysam.VariantHeader()
    for sample in header_in.samples:
        header.add_sample(sample)
    for line in header_in.records:
        header.add_line(str(line))

    # Check annotation fields exist
    for k in transfer_infos:
        if k not in header_ann.info:
            raise ValueError(f"INFO field {k} not found in annotation vcf header")
    for k in transfer_formats:
        if k not in header_ann.formats:
            raise ValueError(f"FORMAT field {k} not found in annotation vcf header")

    # Add annotation headerlines
    for line in header_ann.records:
        if len(line.attrs) > 0 and 'ID' in line.keys() and \
                (line['ID'] in transfer_infos or line['ID'] in transfer_formats):
            header.add_line(str(line))
    return header


def get_annotations(ann_vcf: pysam.VariantFile,
                    header_out: pysam.VariantHeader,
                    transfer_infos: Set[Text],
                    transfer_formats: Set[Text]) -> Tuple[Dict, Dict]:
    """
    Gathers info/format metadata from a source vcf

    Parameters
    ----------
    ann_vcf: pysam.VariantFile
        Source vcf with annotations
    header_out: pysam.VariantHeader
        Output vcf header
    transfer_infos: Set[Text]
        set of info fields to annotate
    transfer_formats: Set[Text]
        set of format fields to annotate

    Returns
    -------
    (dict, dict)
        pair of nested dicts with INFO and FORMAT fields: {vid: {info: value}}, {vid: {sample: {format: value}}}
    """
    sample_list = sorted(list(set(header_out.samples).intersection(set(ann_vcf.header.samples))))
    info_list = sorted(list(transfer_infos))
    format_list = sorted(list(transfer_formats))
    if len(sample_list) == 0 and len(transfer_formats) > 0:
        raise ValueError(f"There are FORMAT fields to transfer but no samples in common between the vcfs")
    info_dict = dict()
    format_dict = dict()
    for record in ann_vcf:
        if record.id in info_dict:
            raise ValueError(f"Encountered duplicate variant id: {record.id}")
        info_dict[record.id] = [record.info.get(k, None) for k in info_list]
        sample_formats = [record.samples[s] for s in sample_list]
        format_dict[record.id] = \
            [[f.get(k, None) for k in format_list] for f in sample_formats]
    return info_dict, format_dict, sample_list, info_list, format_list


def annotate(record: pysam.VariantRecord,
             vcf_out: pysam.VariantFile,
             info_dict: Tuple[Dict, Dict],
             format_dict: Tuple[Dict, Dict],
             sample_list: List[Text],
             info_list: List[Text],
             format_list: List[Text]) -> pysam.VariantRecord:
    """
    Annotates record with INFO and FORMAT fields

    Parameters
    ----------
    record: pysam.VariantRecord
        input record
    vcf_out: pysam.VariantFile
        output vcf object
    info_dict: Set[Text]
        info field dict from get_annotations: {vid: [value]} i.e. vid -> info_index
    format_dict: Set[Text]
        format field dict from get_annotations: {vid: [[value]]} i.e. vid -> sample_index -> format_index
    sample_list: List[Text]
        list of sample IDs as ordered in format_dict
    info_list: List[Text]
        list of format keys as ordered in format_dict
    format_list: List[Text]
        list of format keys as ordered in format_dict

    Returns
    -------
    header: pysam.VariantRecord
        annotated record
    """
    new_record = vcf_out.new_record(alleles=record.alleles, contig=record.contig, filter=record.filter, id=record.id,
                                    info=record.info, start=record.start, qual=record.qual, stop=record.stop)
    for sample, dest in new_record.samples.items():
        source = record.samples[sample]
        for key, val in source.items():
            dest[key] = val
    if record.id in info_dict:
        for key, val in zip(info_list, info_dict[record.id]):
            if val is not None:
                new_record.info[key] = val
    if record.id in format_dict:
        for sample_id, sample_format in zip(sample_list, format_dict[record.id]):
            new_sample_format = new_record.samples[sample_id]
            for key, val in zip(format_list, sample_format):
                if val is not None:
                    new_sample_format[key] = val
    return new_record


def cmp_coord(x_chrom: Text,
              x_pos: int,
              y_chrom: Text,
              y_pos: int,
              contig_indices: Dict[Text, int]) -> int:
    x_chrom_i = contig_indices[x_chrom]
    y_chrom_i = contig_indices[y_chrom]
    if x_chrom_i < y_chrom_i:
        return -1
    elif x_chrom_i > y_chrom_i:
        return 1
    elif x_pos < y_pos:
        return -1
    elif x_pos > y_pos:
        return 1
    else:
        return 0


class ContigComparator(str):
    def __lt__(self, other):
        num_self = re.sub("[^0-9]", "", self)
        num_other = re.sub("[^0-9]", "", other)
        if num_self == "" and num_other == "":
            return str(self) < str(other)  # e.g. chrX, chrY
        elif num_self == "":
            return False
        elif num_other == "":
            return True
        else:
            return int(num_self) < int(num_other)


def pairs(xs, ys, header):
    xs = deque(xs)
    ys = deque(ys)

    # Sometimes contigs aren't sorted properly in the header
    sorted_contig_list = [c for c in header.contigs]
    sorted_contig_list.sort(key=ContigComparator)
    contig_indices = {c: sorted_contig_list.index(c) for c in sorted_contig_list}

    chrom = None
    pos = None
    x_buff = dict()
    y_buff = dict()
    while xs and ys:
        x = xs[0]
        y = ys[0]
        if chrom is None:
            cxy = cmp_coord(x.chrom, x.pos, y.chrom, y.pos, contig_indices)
            if cxy <= 0:
                chrom = x.chrom
                pos = x.pos
            else:
                chrom = y.chrom
                pos = y.pos
        else:
            cxy = cmp_coord(x.chrom, x.pos, y.chrom, y.pos, contig_indices)
            cx = cmp_coord(x.chrom, x.pos, chrom, pos, contig_indices)
            cy = cmp_coord(y.chrom, y.pos, chrom, pos, contig_indices)
            if cx < 0:
                raise ValueError('X pos less than current pos')
            elif cy < 0:
                raise ValueError('Y pos less than current pos')
            if cx == 0:
                x_buff[x.id] = xs.popleft()
            if cy == 0:
                y_buff[y.id] = ys.popleft()
            if cx > 0 and cy > 0:
                for key, val in x_buff.items():
                    yield val, y_buff.get(key, None)
                x_buff = dict()
                y_buff = dict()
                if cxy <= 0:
                    chrom = x.chrom
                    pos = x.pos
                else:
                    chrom = y.chrom
                    pos = y.pos
    for key, val in x_buff.items():
        yield val, y_buff.get(key, None)
    for x in xs:
        yield x, None
    for y in ys:
        yield None, y


def annotate_v2(in_record: pysam.VariantRecord,
                ann_record: pysam.VariantRecord,
                vcf_out: pysam.VariantFile,
                info_list: List[Text],
                format_list: List[Text]) -> pysam.VariantRecord:
    """
    Annotates record with INFO and FORMAT fields

    Parameters
    ----------
    record: pysam.VariantRecord
        input record
    vcf_out: pysam.VariantFile
        output vcf object
    info_dict: Set[Text]
        info field dict from get_annotations: {vid: [value]} i.e. vid -> info_index
    format_dict: Set[Text]
        format field dict from get_annotations: {vid: [[value]]} i.e. vid -> sample_index -> format_index
    sample_list: List[Text]
        list of sample IDs as ordered in format_dict
    info_list: List[Text]
        list of format keys as ordered in format_dict
    format_list: List[Text]
        list of format keys as ordered in format_dict

    Returns
    -------
    header: pysam.VariantRecord
        annotated record
    """
    new_record = vcf_out.new_record(alleles=in_record.alleles, contig=in_record.contig, filter=in_record.filter, id=in_record.id,
                                    info=in_record.info, start=in_record.start, qual=in_record.qual, stop=in_record.stop)
    sample_list = sorted(list(set(vcf_out.header.samples).intersection(set(ann_record.samples))))
    for sample, dest in new_record.samples.items():
        source = in_record.samples[sample]
        for key, val in source.items():
            dest[key] = val
    for key in info_list:
        val = ann_record.info.get(key, None)
        if val is not None:
            new_record.info[key] = val
    for sample_id in sample_list:
        ann_format = ann_record.samples[sample_id]
        new_sample_format = new_record.samples[sample_id]
        for key in format_list:
            val = ann_format.get(key, None)
            if val is not None:
                new_sample_format[key] = val
    return new_record


def run(vcf_in,
        vcf_ann,
        vcf_out,
        info_keys,
        format_keys):
    for in_rec, ann_rec in pairs(vcf_in, vcf_ann, vcf_out.header):
        if in_rec is None:
            continue
        if ann_rec is None:
            vcf_out.write(in_rec)
        else:
            if in_rec.id != ann_rec.id:
                raise ValueError(f"Record ids do not match: {in_rec.id} {ann_rec.id}")
            if in_rec.chrom != ann_rec.chrom or in_rec.pos != ann_rec.pos:
                raise ValueError(f"Record positions do not match for {in_rec.id}: "
                                 f"{in_rec.chrom}:{in_rec.pos} {ann_rec.chrom}:{ann_rec.pos}")
            vcf_out.write(annotate_v2(in_record=in_rec, ann_record=ann_rec, vcf_out=vcf_out, info_list=info_keys,
                                  format_list=format_keys))


def __parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return list()
    else:
        return arg.split(',')


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Annotate the given vcf with fields from another, matching by variant id",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf", type=str,
                        help="VCF to annotate")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--ann-vcf", type=str, required=True,
                        help="Vcf with source annotations")
    parser.add_argument("--formats", type=str,
                        help="Comma-delimited list of FORMAT fields to transfer (e.g. SR_GT,SR_GQ)")
    parser.add_argument("--infos", type=str,
                        help="Comma-delimited list of INFO fields to transfer (e.g. AC,AN,AF)")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    transfer_formats = set(__parse_arg_list(arguments.formats))
    transfer_infos = set(__parse_arg_list(arguments.infos))
    if len(transfer_infos) + len(transfer_formats) == 0:
        raise ValueError('No INFO or FORMAT fields were specified')

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf_in, pysam.VariantFile(arguments.ann_vcf) as vcf_ann:
        header = create_header(
            header_in=vcf_in.header,
            header_ann=vcf_ann.header,
            transfer_formats=transfer_formats,
            transfer_infos=transfer_infos
        )
        with pysam.VariantFile(arguments.out, mode='w', header=header) as vcf_out:
            print("Applying annotations...")
            run(vcf_in=vcf_in, vcf_ann=vcf_ann, vcf_out=vcf_out, info_keys=transfer_infos, format_keys=transfer_formats)


if __name__ == "__main__":
    main()
