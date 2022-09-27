#!/usr/bin/env python

import sys
import os
import warnings
import enum
import argparse
import tempfile
import pysam
from typing import List, Text, Optional, Iterator


from sv_utils import genomics_io


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


class VcfKeys:
    gt = genomics_io.VcfKeys.gt
    gq = genomics_io.VcfKeys.gq
    svtype = genomics_io.VcfKeys.svtype
    bnd_contig_2 = genomics_io.VcfKeys.bnd_contig_2
    bnd_end_2 = genomics_io.VcfKeys.bnd_end_2
    source = genomics_io.VcfKeys.source


BND_TYPES = frozenset({"BND", "CTX"})
GAIN_TYPES = frozenset({"INS", "CPX", "MEI", "DUP"})
VALID_EMPTY_INTERVAL_TYPES = frozenset({"INS", "BND", "CPX", "CTX"})


def str_to_loc(loc: str) -> (str, int):
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


def fix_bad_end2_header(input_vcf: str, encoding=Default.encoding) -> pysam.VariantHeader:
    f"""
    Fix {VcfKeys.bnd_end_2} header by writing to new temporary VCF
    Args:
        input_vcf: str
            Path to VCF with bad {VcfKeys.bnd_end_2} header
        encoding: str
            Encoding for python strings
    Returns:
        temporary_vcf: str
            Path to temporary VCF with fixed {VcfKeys.bnd_end_2} header
    """
    print(f"Fixing out of order {VcfKeys.bnd_end_2} header line in {input_vcf}", file=sys.stderr)
    bad_start = f"##INFO=<ID={VcfKeys.bnd_end_2}".encode(encoding=encoding)
    comment = '#'.encode(encoding=encoding)
    newline = '\n'.encode(encoding=encoding)
    fixed = (
        f"##INFO=<ID={VcfKeys.bnd_end_2},Number=1,Type=Integer,Description=\"Position of breakpoint on "
        f"{VcfKeys.bnd_contig_2}\">"
    ).encode(encoding=encoding)
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=True) as f_out:
        # load the old header as a plain tabix file, and fix the bad END2 header line
        with pysam.BGZFile(input_vcf, "rb", index=None) as f_in:
            for line in f_in:
                if not line.startswith(comment):
                    break  # only take the header
                write_line = fixed if line.startswith(bad_start) else line
                f_out.write(write_line + newline)
        f_out.flush()
        # now load the fixed header
        with pysam.VariantFile(f_out.name, "r") as f_in:
            return f_in.header


def fix_record(
        record: pysam.VariantRecord,
        error_on_mismatched_chr2: bool = Default.error_on_mismatched_chr2,
        flipped_interval_strategy: FlippedIntervalStrategy = Default.flipped_interval_strategy,
        gq_scale: Optional[float] = None
) -> pysam.VariantRecord:
    f"""
    For supplied VCF record, detect and fix errors with pos/end, and allele designation
    Args:
        record: pysam.VariantRecord
            Variant line in VCF
        error_on_mismatched_chr2: bool (Default={Default.error_on_mismatched_chr2})
            If {VcfKeys.bnd_contig_2} does not match info in SOURCE, then throw an error if true, otherwise warn
        flipped_interval_strategy: FlippedIntervalStrategy (Default={Default.flipped_interval_strategy})
            If END < POS, fix the interval with the defined strategy
        gq_scale: Optional[float] (Default=None)
            If not None, multiply GQ by this scale
    Returns:
        fixed_record: pysam.VariantRecord
            Corrected variant line
    """
    if VcfKeys.bnd_contig_2 in record.info:
        if record.info[VcfKeys.svtype] in BND_TYPES:
            if VcfKeys.bnd_end_2 not in record.info:
                record.info[VcfKeys.bnd_end_2] = record.stop
            record.stop = record.pos + 1
        else:
            # not a break-end, figure out what's going on here...
            if VcfKeys.source in record.info.keys():
                chr2, end2 = str_to_loc(record.info[VcfKeys.source])
                if record.info[VcfKeys.bnd_contig_2] != chr2:
                    message = f"Variant {record.id} with {VcfKeys.svtype}={record.info[VcfKeys.svtype]} has " \
                        f"{VcfKeys.bnd_contig_2}!={VcfKeys.source} chrom ({record.info[VcfKeys.bnd_contig_2]}!={chr2})"
                    if error_on_mismatched_chr2:
                        raise ValueError(message)
                    else:
                        warnings.warn(message)
                        record.info[VcfKeys.bnd_contig_2] = chr2

                record.info[VcfKeys.bnd_end_2] = end2
            else:
                # no alternate source and not a break-end
                if record.info[VcfKeys.bnd_contig_2] != record.chrom:
                    # use this somehow?
                    # if record.info[VcfKeys.svtype] in "INS" and record.stop == record.pos:
                    #     record.stop = record.pos + 1
                    raise ValueError(
                        f"Variant {record.id} has {VcfKeys.svtype}={record.info[VcfKeys.svtype]} but "
                        f"{VcfKeys.bnd_contig_2}!=chrome ({record.info[VcfKeys.bnd_contig_2]}!={record.chrom})"
                    )
                record.info.pop(VcfKeys.bnd_contig_2)

    if record.pos >= record.stop:
        if record.pos == record.stop:
            # empty interval, valid for INS, otherwise an error
            if record.info[VcfKeys.svtype] not in ["INS", "BND", "CPX", "CTX"]:
                raise ValueError(
                    f"Variant {record.id} with {VcfKeys.svtype}={record.info[VcfKeys.svtype]} has pos=stop={record.pos}"
                )
        else:
            # sometimes this error occurs for INS or similar SV types
            warnings.warn(f"{record.id} had END < POS")
            if flipped_interval_strategy == FlippedIntervalStrategy.end_to_pos:
                record.stop = record.pos
            else:
                # just swap end and pos
                record.pos, record.stop = record.stop, record.pos

    if record.alts == ('<DUP>',):
        # some older VCFs had DUP alleles marked as 2 instead of 1
        num_bad_gt = 0
        for sample in record.samples.values():
            gt = sample[VcfKeys.gt]
            if gt is not None and any(allele is not None and allele > 1 for allele in gt):
                # bad sample GT, DUP allele should be 0 or 1
                num_bad_gt += 1
                sample[VcfKeys.gt] = tuple(
                    a - 2 if (a is not None and a > 1) else a
                    for a in gt
                )
        if num_bad_gt > 0:
            warnings.warn(f"{record.id} had {num_bad_gt} bad samples")
    if gq_scale is not None:
        for sample in record.samples.values():
            gq = sample[VcfKeys.gq]
            if gq is not None:
                sample[VcfKeys.gq] = int(round(gq * gq_scale))

    return record


def low_mem_cheat_sort(record_iterator: Iterator[pysam.VariantRecord]) -> Iterator[pysam.VariantRecord]:
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
    buffer = []
    for record in record_iterator:
        if record.chrom == chrom:
            # still on same contig, sort
            buffer.append(record)
        else:
            # finished the current contig
            yield from sorted(buffer, key=lambda r: (genomics_io.contig_sort_key(r.chrom), r.pos, r.stop))
            buffer = [record]
            chrom = record.chrom
    # finished last contig
    yield from sorted(buffer, key=lambda r: (genomics_io.contig_sort_key(r.chrom), r.pos, r.stop))


def fix_vcf(
        input_vcf: str,
        output_vcf: str,
        num_threads: int = Default.num_threads,
        index_output_vcf: bool = Default.index_output_vcf,
        error_on_mismatched_chr2: bool = Default.error_on_mismatched_chr2,
        flipped_interval_strategy: FlippedIntervalStrategy = Default.flipped_interval_strategy,
        gq_scale: Optional[float] = None,
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
        num_threads: int (default={Default.num_threads})
            Number of threads for compressing/decompressing VCFs
        index_output_vcf: bool (default={Default.index_output_vcf})
            If True, create tabix index for output_vcf
        error_on_mismatched_chr2: bool (Default={Default.error_on_mismatched_chr2})
            If {VcfKeys.bnd_contig_2} does not match info in SOURCE, then throw an error if true, otherwise warn
        flipped_interval_strategy: FlippedIntervalStrategy (Default={Default.flipped_interval_strategy})
            If END < POS, fix the interval with the defined strategy
        gq_scale: Optional[float] (Default=None)
            If not None, multiply GQ by this scale
        encoding: str (default={Default.encoding})
            Format for encoding bytes for header processing
    """
    with pysam.VariantFile(input_vcf, 'r') as f_in:
        output_header = f_in.header
        if VcfKeys.bnd_end_2 in output_header.info:
            print(f"{VcfKeys.bnd_end_2} already in header")
            end2_rec = str(output_header.info[VcfKeys.bnd_end_2].record)
            if end2_rec.index("Type") < end2_rec.index("Number"):
                # bad order, need to fix this tag. Unfortunately GATK can't ignore this "problem" and
                # pysam can't fix it. So need to use text processing to fix this header into a temp VCF,
                # then process *THAT* for the main fix
                output_header = fix_bad_end2_header(input_vcf=input_vcf, encoding=encoding)
        else:
            print(f"Adding {VcfKeys.bnd_end_2} to header")
            # add END2 tag to header
            output_header.info.add(
                VcfKeys.bnd_end_2, '1', "Integer", "Distal position of breakend"
            )
        print("Fixing problems with records", file=sys.stderr)
        with pysam.VariantFile(output_vcf, 'w', header=output_header, threads=num_threads) as f_out:
            for record in low_mem_cheat_sort(
                fix_record(_record, error_on_mismatched_chr2=error_on_mismatched_chr2,
                           flipped_interval_strategy=flipped_interval_strategy,
                           gq_scale=gq_scale)
                for _record in f_in.fetch()
            ):
                f_out.write(record)

    if index_output_vcf:
        pysam.tabix_index(output_vcf, preset="vcf", force=True)


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
    parser.add_argument("--index-output-vcf", type=bool, default=Default.index_output_vcf,
                        help="if true, create tabix index for output vcf")
    parser.add_argument("--error-on-mismatched-chr2", type=bool, default=Default.error_on_mismatched_chr2,
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
    fix_vcf(input_vcf=args.input_vcf, output_vcf=args.output_vcf,
            num_threads=args.num_threads, error_on_mismatched_chr2=args.error_on_mismatched_chr2,
            flipped_interval_strategy=FlippedIntervalStrategy[args.flipped_interval_strategy],
            gq_scale=args.gq_scale, index_output_vcf=args.index_output_vcf)


if __name__ == "__main__":
    main()
