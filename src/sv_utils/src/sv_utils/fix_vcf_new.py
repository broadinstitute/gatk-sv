#!/usr/bin/env python

import sys
import os
import warnings
import enum
import argparse
import pysam
import attrs
from typing import List, Text, Optional, Iterator


from sv_utils import common, genomics_io


# noinspection PyArgumentList
class FlippedIntervalStrategy(enum.Enum):
    end_to_pos = enum.auto()
    flip = enum.auto()


class Default:
    num_threads = os.cpu_count()
    encoding = "utf-8"
    index_output_vcf = True
    error_on_mismatched_chr2 = False
    flipped_interval_strategy = FlippedIntervalStrategy.end_to_pos
    gq_scale: Optional[float] = None


class VcfKeys:
    gt = genomics_io.VcfKeys.gt
    gq = genomics_io.VcfKeys.gq
    svtype = genomics_io.VcfKeys.svtype
    end = genomics_io.VcfKeys.end
    bnd_contig_2 = genomics_io.VcfKeys.bnd_contig_2
    bnd_end_2 = genomics_io.VcfKeys.bnd_end_2
    source = genomics_io.VcfKeys.source


BND_TYPES = frozenset({"BND", "CTX"})
GAIN_TYPES = frozenset({"INS", "CPX", "MEI", "DUP"})
VALID_EMPTY_INTERVAL_TYPES = frozenset({"INS", "BND", "CPX", "CTX"})


@attrs.frozen(slots=True)
class SortableRecordLine:
    record_line: str
    chrom: str
    pos: int
    end: int


def fix_vcf(
        input_vcf: str,
        output_vcf: str,
        index_output_vcf: bool = Default.index_output_vcf,
        error_on_mismatched_chr2: bool = Default.error_on_mismatched_chr2,
        flipped_interval_strategy: FlippedIntervalStrategy = Default.flipped_interval_strategy,
        encoding: str = Default.encoding
):
    f"""
    Fix issues with VCF header and record spec compliance that may prevent GATK from reading input VCF.
    This may alter correct record order, so sort corrected VCF records.
    Create tabix index for output VCF if requested.
    Args:
        input_vcf: str
            Path to input VCF
        output_vcf: str
            Path to save corrected VCF
        index_output_vcf: bool (default={Default.index_output_vcf})
            If True, create tabix index for output_vcf
        error_on_mismatched_chr2: bool (Default={Default.error_on_mismatched_chr2})
            If {VcfKeys.bnd_contig_2} does not match info in SOURCE, then throw an error if true, otherwise warn
        flipped_interval_strategy: FlippedIntervalStrategy (Default={Default.flipped_interval_strategy})
            If END < POS, fix the interval with the defined strategy
        encoding: str (default={Default.encoding})
            Format for encoding bytes for header processing
    """
    with pysam.BGZFile(input_vcf, "rb", index=None) as vcf_in:
        peekable_vcf_in = common.PeekableIter(
            vcf_line.decode(encoding) for vcf_line in vcf_in
        )
        header_lines = []
        while peekable_vcf_in.has_next() and peekable_vcf_in.peek_next().startswith("#"):
            header_lines.append(next(peekable_vcf_in))

        with pysam.BGZFile(output_vcf, "wb", index=None) as vcf_out:
            for fixed_header_line in _get_fixed_header(header_lines):
                vcf_out.write(f"{fixed_header_line}\n".encode(encoding))

            for fixed_record_line in _low_mem_cheat_sort(
                fix_record_line(
                    record_line,
                    error_on_mismatched_chr2=error_on_mismatched_chr2,
                    flipped_interval_strategy=flipped_interval_strategy,
                )
                for record_line in peekable_vcf_in
            ):
                vcf_out.write(f"{fixed_record_line.record_line}\n".encode(encoding))

    if index_output_vcf:
        pysam.tabix_index(output_vcf, preset="vcf", force=True)


def _get_fixed_header(header_lines: list[str]) -> list[str]:
    """Fix common problems with header lines"""
    # guarantee that the correct header lines are present. If they are, don't try to parse them,
    # just remove them. Either way, add them into the header so they are there correctly.
    mandatory_header_replacements = {
        f"##INFO=<ID={VcfKeys.bnd_end_2},":
            '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">',
        f"##INFO=<ID=CHR2":
            '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">',

    }
    # replace these header lines with the correct values if they are present, otherwise don't do
    # anything
    optional_header_replacements = {
        f"##FORMAT=<ID=EV,":
        '##FORMAT=<ID=EV,Number=.,Type=String,Description="Classes of evidence supporting final genotype">',
    }
    for line_index, line in enumerate(header_lines):
        for start, fix_line in list(mandatory_header_replacements.items()):
            if line.startswith(start):
                header_lines[line_index] = fix_line
                mandatory_header_replacements.pop(start)
                break
        for start, fix_line in list(optional_header_replacements.items()):
            if line.startswith(start):
                header_lines[line_index] = fix_line
                optional_header_replacements.pop(start)
                break
    for fix_line in mandatory_header_replacements.values():
        header_lines.append(fix_line)
    return header_lines


def fix_record_line(
        record_line: str,
        error_on_mismatched_chr2: bool = Default.error_on_mismatched_chr2,
        flipped_interval_strategy: FlippedIntervalStrategy = Default.flipped_interval_strategy,
) -> SortableRecordLine:
    f"""
    For supplied VCF record, detect and fix errors with pos/end, and allele designation
    Args:
        record_line: text of variant line in VCF
        error_on_mismatched_chr2: (Default={Default.error_on_mismatched_chr2})
            If {VcfKeys.bnd_contig_2} does not match info in SOURCE, then throw an error if true,
            otherwise warn
        flipped_interval_strategy: (Default={Default.flipped_interval_strategy})
            If END < POS, fix the interval with the defined strategy
    Returns:
        fixed_record: pysam.VariantRecord
            Corrected variant line
    """
    columns = record_line.split("\t", 7)
    info = {
        word[0]: word[1]
        for field in columns[7].split(";")
        for word in field.split("=", 1)
    }
    modified_info = False
    modified_pos = False
    pos = int(columns[1])
    if VcfKeys.bnd_contig_2 in info:
        if info[VcfKeys.svtype] in BND_TYPES:
            if VcfKeys.bnd_end_2 not in info:
                info[VcfKeys.bnd_end_2] = info[VcfKeys.end]
            info[VcfKeys.end] = f"{pos + 1}"
            modified_info = True
        else:
            # not a break-end, figure out what's going on here...
            if VcfKeys.source in info:
                chr2, end2 = _str_to_loc(info[VcfKeys.source])
                if info[VcfKeys.bnd_contig_2] != chr2:
                    record_id = columns[2]
                    message = (
                        f"Variant {record_id} with {VcfKeys.svtype}={info[VcfKeys.svtype]} has "
                        f"{VcfKeys.bnd_contig_2}!={VcfKeys.source} chrom "
                        f"({info[VcfKeys.bnd_contig_2]}!={chr2})"
                    )
                    if error_on_mismatched_chr2:
                        raise ValueError(message)
                    else:
                        warnings.warn(message)
                        info[VcfKeys.bnd_contig_2] = chr2

                info[VcfKeys.bnd_end_2] = end2
                modified_info = True
            else:
                # no alternate source and not a break-end
                chrom = columns[0]
                if info[VcfKeys.bnd_contig_2] != chrom:
                    # use this somehow?
                    # if record.info[VcfKeys.svtype] in "INS" and record.stop == record.pos:
                    #     record.stop = record.pos + 1
                    record_id = columns[2]
                    raise ValueError(
                        f"Variant {record_id} has {VcfKeys.svtype}={info[VcfKeys.svtype]} but "
                        f"{VcfKeys.bnd_contig_2}!=chrome ({info[VcfKeys.bnd_contig_2]}!={chrom})"
                    )
                info.pop(VcfKeys.bnd_contig_2)
                modified_info = True

    end = int(info[VcfKeys.end])
    if pos >= end:
        record_id = columns[2]
        if pos == end:
            # empty interval, valid for INS, otherwise an error
            if info[VcfKeys.svtype] not in VALID_EMPTY_INTERVAL_TYPES:
                raise ValueError(
                    f"Variant {record_id} with {VcfKeys.svtype}={info[VcfKeys.svtype]} has "
                    f"pos=stop={pos}"
                )
        else:
            # sometimes this error occurs for INS or similar SV types
            warnings.warn(f"{record_id} had END < POS")
            if flipped_interval_strategy == FlippedIntervalStrategy.end_to_pos:
                end = pos
            else:
                # just swap end and pos
                pos, end = end, pos
                modified_pos = True
            info[VcfKeys.end] = end
            modified_info = True

    if modified_info:
        if modified_pos:
            columns[1] = f"{pos}"
        columns[7] = ";".join(f"{key}={value}" for key, value in info.items())
        record_line = "\t".join(columns)
    return SortableRecordLine(
        record_line=record_line,
        chrom=columns[0],
        pos=pos,
        end=end,
    )


def _str_to_loc(loc: str) -> (str, int):
    f"""
    Extract contig and end position from {VcfKeys.source} INFO field
    Args:
        loc: str
            Text of {VcfKeys.source} INFO field
    Returns:
        contig: str
            contig of source
        end: int
            end position of source
    """
    contig, remainder = loc.split('_', 1)[1].split(':', 1)
    start, end = remainder.split('-', 1)
    return contig, int(end)


def _low_mem_cheat_sort(
    record_iterator: Iterator[SortableRecordLine]
) -> Iterator[SortableRecordLine]:
    """
    Sort variants without storing whole VCF in memory. Fixing the records can alter start and end coordinates, but it
    cannot alter the contig. This is much faster than sorting the whole VCF afterwards, because most records will be
    emitted as they are produced, with only a handful of records moving slightly in the order (and it is done in memory
    without touching disk).
    Args:
        record_iterator: Iterator[pysam.VariantRecord]
            Iterator of fixed records, almost in order
    Returns:
        sorted_record_iterator: Iterator[pysam.VariantRecord]
            Iterator of fixed records in order
    """
    chrom = "Dummy Start Contig"
    buffer: list[SortableRecordLine] = []
    for sortable_record_line in record_iterator:
        if sortable_record_line.chrom == chrom:
            # still on same contig, sort
            buffer.append(sortable_record_line)
        else:
            # finished the current contig
            yield from sorted(
                buffer, key=lambda r: (genomics_io.contig_sort_key(r.chrom), r.pos, r.end)
            )
            buffer = [sortable_record_line]
            chrom = sortable_record_line.chrom
    # finished last contig
    yield from sorted(buffer, key=lambda r: (genomics_io.contig_sort_key(r.chrom), r.pos, r.end))


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Put sv pipeline VCF into spec compliance. Maintains sorted VCFs, does not sort unsorted VCFs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("input_vcf", type=str, help=".vcf file to fix")
    parser.add_argument("output_vcf", type=str, help='in-spec output vcf')
    parser.add_argument("--num_threads", "-@", type=int, default=Default.num_threads,
                        help="number of threads for compressing output vcf")
    parser.add_argument("--index-output-vcf", type=common.argparse_bool, default=Default.index_output_vcf,
                        help="if true, create tabix index for output vcf")
    parser.add_argument("--error-on-mismatched-chr2", type=common.argparse_bool,
                        default=Default.error_on_mismatched_chr2,
                        help=f"If {VcfKeys.bnd_contig_2} does not match info in SOURCE, then throw an error if true, "
                             "otherwise warn")
    parser.add_argument("--flipped-interval-strategy", type=str, default=Default.flipped_interval_strategy.name,
                        help="If END < POS, fix the interval with the defined strategy.",
                        choices={name for name, __ in FlippedIntervalStrategy.__members__.items()})
    parser.add_argument("--gq-scale", type=float, default=None,
                        help="multiplicative scale to GQ")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[List[Text]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    fix_vcf(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        error_on_mismatched_chr2=args.error_on_mismatched_chr2,
        flipped_interval_strategy=FlippedIntervalStrategy[args.flipped_interval_strategy],
        index_output_vcf=args.index_output_vcf
    )


if __name__ == "__main__":
    main()
