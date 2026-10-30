import gzip
import hashlib
import logging
import os
import pathlib
import pysam
import shutil
import re
import subprocess
import tarfile
import tempfile
import uuid

from typing import Callable, List
from dataclasses import dataclass
from pathlib import Path
from tqdm import tqdm


logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


@dataclass
class Region:
    chr: str
    start: int
    end: int

    @staticmethod
    def to_file(working_dir, regions):
        filename = str(uuid.uuid4())
        filename = os.path.join(working_dir, filename + ".bed")
        with open(filename, "w") as regions_file:
            for r in regions:
                regions_file.write("\t".join([str(r.chr), str(r.start), str(r.end)]) + "\n")
        return filename


class MiniReference:
    """
    A reference genome made of only the sequences of the target regions;
    the target regions on each contig are concatenated (in the order of their start positions)
    to make a contig with the same name. For instance, given the regions chr1:100-200 and chr1:500-600,
    the mini reference contains `chr1` of length 200 where chr1:100-200 is placed at chr1:0-100
    and chr1:500-600 is placed at chr1:100-200.

    The class also implements methods to lift over coordinates from the original reference to
    the mini reference. All the coordinates are 0-based and half-open (i.e., BED-style).
    """
    def __init__(self, reference_fasta: str, regions: List[Region], working_dir: str):
        with pysam.FastaFile(reference_fasta) as fasta:
            contigs_order = {name: i for i, name in enumerate(fasta.references)}
            self.regions = sorted(regions, key=lambda r: (contigs_order[r.chr], r.start))
            self.offsets = []
            self.contig_lengths = {}
            sequences = {}
            for r in self.regions:
                self.offsets.append(self.contig_lengths.get(r.chr, 0))
                self.contig_lengths[r.chr] = self.contig_lengths.get(r.chr, 0) + r.end - r.start
                sequences[r.chr] = sequences.get(r.chr, "") + fasta.fetch(r.chr, r.start, r.end)

        # The regions (and their offsets on the mini reference) per contig, for fast lookups.
        self.regions_by_chr = {}
        for region, offset in zip(self.regions, self.offsets):
            self.regions_by_chr.setdefault(region.chr, []).append((region, offset))

        self.contig_md5s = {k: hashlib.md5(v.upper().encode()).hexdigest() for k, v in sequences.items()}
        self.working_dir = tempfile.mkdtemp(prefix="tmp_mini_reference_", dir=working_dir)
        self.fasta_filename = os.path.join(self.working_dir, Path(reference_fasta).name)
        line_width = 100
        with open(self.fasta_filename, "w") as f:
            for name, sequence in sequences.items():
                f.write(f">{name}\n")
                for i in range(0, len(sequence), line_width):
                    f.write(sequence[i:i + line_width] + "\n")
        pysam.faidx(self.fasta_filename)

    def cleanup(self):
        shutil.rmtree(self.working_dir)

    def get_sam_header_sequences(self) -> List[dict]:
        return [{"SN": name, "LN": length, "M5": self.contig_md5s[name]}
                for name, length in self.contig_lengths.items()]

    def get_sam_header_sequence_lines(self) -> List[str]:
        return [f"@SQ\tSN:{sq['SN']}\tLN:{sq['LN']}\tM5:{sq['M5']}\n" for sq in self.get_sam_header_sequences()]

    def write_sequence_dict(self, filename: str):
        """
        Writes the sequence dictionary of the mini reference. The dictionary is written
        from the mini reference (rather than its FASTA file), so the dictionaries created
        for the inputs of different workflows are identical.
        """
        with open(filename, "w") as f:
            f.write("@HD\tVN:1.6\n")
            f.writelines(self.get_sam_header_sequence_lines())

    def lift(self, chr: str, start: int, end: int):
        """
        Returns the coordinates of the given interval on the mini reference,
        or None if the interval is not entirely within a target region.
        """
        for region, offset in self.regions_by_chr.get(chr, []):
            if region.start <= start and end <= region.end:
                return offset + start - region.start, offset + end - region.start
        return None

    def clip_and_lift(self, chr: str, start: int, end: int) -> List[tuple]:
        """
        Returns the coordinates on the mini reference of the parts of the given interval
        overlapping target regions; the parts adjacent on the mini reference are merged.
        """
        if start == end:
            # A zero-length interval (e.g., the site of an insertion or a breakend) is lifted over as a point.
            lifted_interval = self.lift(chr, start, start + 1)
            return [] if lifted_interval is None else [(lifted_interval[0], lifted_interval[0])]
        lifted = []
        for region, offset in self.regions_by_chr.get(chr, []):
            if max(start, region.start) < min(end, region.end):
                lifted_start = offset + max(start, region.start) - region.start
                lifted_end = offset + min(end, region.end) - region.start
                if lifted and lifted[-1][1] == lifted_start:
                    lifted[-1] = (lifted[-1][0], lifted_end)
                else:
                    lifted.append((lifted_start, lifted_end))
        return lifted


class BaseTransformer:
    def __init__(self, working_dir: str, callback: Callable[[str, ...], dict], input_root: str = None,
                 mini_reference: MiniReference = None, output_subdir: str = None,
                 ref_panel_samples: List[str] = None, **kwargs):
        # Convert the string to an ABS path and make sure it exists.
        self.working_dir = Path(working_dir).resolve(strict=True)
        self.callback = callback
        # When set, the output files are written to the working directory following
        # the same directory hierarchy as the input files have relative to this directory.
        self.input_root = Path(input_root).resolve() if input_root else None
        # When set, the coordinates in the output files are lifted over to this reference.
        self.mini_reference = mini_reference
        # When set, the output files are written to this subdirectory of the working directory
        # (e.g., `intermediates_HG00512/gather_sample_evidence` for an output of an upstream task),
        # regardless of the directory of the input file.
        self.output_subdir = output_subdir
        # When set, the reference panel is subset to these samples (see `is_excluded_sample_file`).
        self.ref_panel_samples = set(ref_panel_samples) if ref_panel_samples else None

    def is_excluded_sample(self, sample_id: str) -> bool:
        return self.ref_panel_samples is not None and sample_id not in self.ref_panel_samples

    def is_excluded_sample_file(self, filename: str) -> bool:
        """
        Returns True if the reference panel is subset, and the given per-sample file (e.g., `HG00096.pe.txt.gz`
        or `std_000.manta.HG00096.vcf.gz`) has none of the kept samples in its dot-separated name.
        """
        return self.ref_panel_samples is not None \
            and self.ref_panel_samples.isdisjoint(Path(filename).name.split("."))

    @staticmethod
    def get_supported_file_types() -> List[str]:
        """
        The file types should include the '.' prefix to match with Path().suffix output
        (e.g., it should return '.cram' instead of 'cram').
        """
        raise NotImplementedError()

    def get_output_filename(self, input_filename, output_prefix):
        output_dir = self.working_dir
        input_dir = Path(input_filename).resolve().parent
        if self.output_subdir is not None:
            output_dir = output_dir.joinpath(self.output_subdir)
            output_dir.mkdir(parents=True, exist_ok=True)
        elif self.input_root is not None and input_dir.is_relative_to(self.input_root):
            output_dir = output_dir.joinpath(input_dir.relative_to(self.input_root))
            output_dir.mkdir(parents=True, exist_ok=True)
        return str(output_dir.joinpath(f"{output_prefix}{Path(input_filename).name}"))

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs):
        raise NotImplementedError()


class BedToIntervalListConverter(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], sequence_dict_filename: str, picard_path: str, **kwargs):
        super().__init__(working_dir, callback, **kwargs)
        self.sequence_dict_filename = sequence_dict_filename
        self.picard_path = picard_path

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".interval_list"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        regions_filename = Region.to_file(self.working_dir, regions)

        subprocess.run(
            ["java", "-jar", self.picard_path, "BedToIntervalList",
             "-I", regions_filename, "-O", output_filename, "-SD", self.sequence_dict_filename],
            check=True)

        os.remove(regions_filename)
        return self.callback(output_filename)


class CramDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str, str], dict], reference_fasta: str, reference_index: str, **kwargs):
        super().__init__(working_dir, callback, **kwargs)
        self.reference_fasta = reference_fasta
        self.reference_index = reference_index
        # Decoding/encoding CRAM files requires the reference; when a local copy of the reference
        # is not available, htslib attempts to fetch the reference sequences by their MD5 from EBI.
        self.local_reference_fasta = reference_fasta if os.path.isfile(reference_fasta) else None

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".cram"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        include_discordant_reads = kwargs.get("include_discordant_reads", False)
        output_filename_unsorted = self.get_output_filename(input_filename, f"unsorted_{output_prefix}")
        with pysam.AlignmentFile(input_filename, "rc", reference_filename=self.local_reference_fasta) as cram:
            header = cram.header

        output_reference_fasta = self.local_reference_fasta
        if self.mini_reference is not None:
            header_dict = header.to_dict()
            header_dict["SQ"] = self.mini_reference.get_sam_header_sequences()
            header = pysam.AlignmentHeader.from_dict(header_dict)
            output_reference_fasta = self.mini_reference.fasta_filename

        # Shared across regions so that pairs with mates in different regions
        # (e.g., discordant pairs supporting an SV) are paired without a linear scan of the file.
        reads_for_second_pass = {}
        with pysam.AlignmentFile(output_filename_unsorted, "wc", header=header,
                                 reference_filename=output_reference_fasta) as output_cram_file:
            for region in regions:
                for r1, r2 in self.read_pair_generator(
                        input_filename, region=f"{region.chr}:{region.start + 1}-{region.end}",
                        reference_filename=self.local_reference_fasta, read_dict=reads_for_second_pass):
                    self.write_pair(output_cram_file, r1, r2)

            if include_discordant_reads:
                for r1, r2 in self.get_discordant_reads(
                        input_filename, reads_for_second_pass, reference_filename=self.local_reference_fasta):
                    self.write_pair(output_cram_file, r1, r2)

        output_filename = self.get_output_filename(input_filename, output_prefix)
        reference_args = ["--reference", output_reference_fasta] if output_reference_fasta else []
        # Newer htslib versions write CRAM v3.1 by default, which is not supported by
        # tools linked against older htslib versions (e.g., Manta), hence explicitly writing CRAM v3.0.
        pysam.sort("-O", "cram,version=3.0", "-o", output_filename, *reference_args, output_filename_unsorted)
        os.remove(output_filename_unsorted)
        index_filename = f"{output_filename}.crai"
        pysam.index(output_filename, index_filename)
        return self.callback(output_filename, index_filename)

    def write_pair(self, output_file: pysam.AlignmentFile, r1: pysam.AlignedSegment, r2: pysam.AlignedSegment):
        if self.mini_reference is None:
            output_file.write(r1)
            output_file.write(r2)
            return

        # Lift over the reads to the mini reference; skip the pair if either of the reads
        # is not entirely within a target region.
        lifted = []
        for read in [r1, r2]:
            # An unmapped read with a mapped mate is placed at the position of its mate.
            end = read.reference_start + 1 if read.is_unmapped else read.reference_end
            lifted_interval = self.mini_reference.lift(read.reference_name, read.reference_start, end)
            if lifted_interval is None:
                return
            lifted.append(lifted_interval[0])

        shifts = [lifted[0] - r1.reference_start, lifted[1] - r2.reference_start]
        for read, mate, read_start, mate_start in [(r1, r2, lifted[0], lifted[1]), (r2, r1, lifted[1], lifted[0])]:
            d = read.to_dict()
            d["ref_pos"] = str(read_start + 1)
            d["next_ref_pos"] = str(mate_start + 1)
            if shifts[0] != shifts[1] and read.reference_name == mate.reference_name \
                    and not read.is_unmapped and not mate.is_unmapped:
                # The mates are lifted from different target regions, hence the template length is updated.
                left = min(read_start, mate_start)
                right = max(read_start + read.reference_length, mate_start + mate.reference_length)
                is_leftmost = read_start < mate_start or (read_start == mate_start and read.is_read1)
                d["length"] = str(right - left if is_leftmost else left - right)
            d["tags"] = [self.lift_tag(tag) for tag in d.get("tags", [])]
            d["tags"] = [tag for tag in d["tags"] if tag is not None]
            output_file.write(pysam.AlignedSegment.from_dict(d, output_file.header))

    def lift_tag(self, tag: str):
        """
        Lifts over the coordinates in a tag of a read, formatted as in `pysam.AlignedSegment.to_dict`;
        returns None if the tag should be dropped.
        """
        if tag.startswith("XA:"):
            # Alternative alignments are mostly outside the target regions.
            return None
        if tag.startswith("SA:"):
            # Format: SA:Z:(rname,pos,strand,CIGAR,mapQ,NM;)+
            lifted_alignments = []
            for alignment in tag[len("SA:Z:"):].rstrip(";").split(";"):
                cols = alignment.split(",")
                lifted_interval = self.mini_reference.lift(cols[0], int(cols[1]) - 1, int(cols[1]))
                if lifted_interval is not None:
                    cols[1] = str(lifted_interval[0] + 1)
                    lifted_alignments.append(",".join(cols))
            return f"SA:Z:{';'.join(lifted_alignments)};" if lifted_alignments else None
        return tag

    @staticmethod
    def read_pair_generator(filename: str, region: str = None, reference_filename: str = None, read_dict: dict = None):
        """
        Generate read pairs in a CRAM file, including discordant (i.e., not properly paired) pairs.
        If the `region` string is provided, it only generates the reads overlapping this region.
        Reads whose mates are not seen yet are kept in `read_dict`, which can be shared between
        calls on different regions to pair reads whose mates are in different regions.
        """
        read_dict = {} if read_dict is None else read_dict
        with pysam.AlignmentFile(filename, "rc", reference_filename=reference_filename) as cram:
            for read in cram.fetch(region=region):
                q_name = read.query_name
                if read.is_secondary or read.is_supplementary:
                    continue
                if q_name not in read_dict:
                    read_dict[q_name] = [None, None]
                    read_dict[q_name][0 if read.is_read1 else 1] = read
                elif read_dict[q_name][0 if read.is_read1 else 1] is not None:
                    # It is the same read seen before (e.g., overlapping regions).
                    continue
                else:
                    if read.is_read1:
                        yield read, read_dict[q_name][1]
                    else:
                        yield read_dict[q_name][0], read
                    del read_dict[q_name]

        return read_dict

    @staticmethod
    def get_discordant_reads(filename: str, discordant_reads: dict, reference_filename: str = None):
        logging.info(f"Linearly scanning {filename} for discordant pairs of {len(discordant_reads)} reads.")
        with pysam.AlignmentFile(filename, "rc", reference_filename=reference_filename) as cram:
            for read in tqdm(cram.fetch(), desc="Iterating on reads", unit=" read", dynamic_ncols=True,
                             bar_format="{desc}: {n:,} [{elapsed}] {rate_fmt}"):
                if read.is_secondary or read.is_supplementary:
                    continue
                query_name = read.query_name
                pair = discordant_reads.get(query_name)
                if pair is not None:
                    if (read.is_read1 and pair[0] is not None) or (not read.is_read1 and pair[1] is not None):
                        # It is the same read as the one seen before.
                        continue
                    if read.is_read1:
                        yield read, pair[1]
                    else:
                        yield pair[0], read
                    del discordant_reads[query_name]
        if len(discordant_reads) > 0:
            logging.warning(f"Did not find discordant pairs for {len(discordant_reads)} reads.")
        logging.info("Finished linear scanning.")


class VcfDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".vcf"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        with pysam.VariantFile(input_filename) as input_file:
            header = input_file.header
            if self.mini_reference is not None:
                header = pysam.VariantHeader()
                for record in input_file.header.records:
                    if record.type != "CONTIG":
                        header.add_record(record)
                for name, length in self.mini_reference.contig_lengths.items():
                    header.contigs.add(name, length=length)
                for sample in input_file.header.samples:
                    header.add_sample(sample)

            with pysam.VariantFile(output_filename, "w", header=header) as output_file:
                for record in input_file:
                    if self.mini_reference is not None:
                        lifted_interval = self.mini_reference.lift(record.contig, record.start, record.stop)
                        if lifted_interval is None:
                            continue
                        has_end = "END" in record.info
                        record.translate(header)
                        record.start = lifted_interval[0]
                        if has_end:
                            record.stop = lifted_interval[1]
                        output_file.write(record)
                        continue
                    for region in regions:
                        if record.contig == region.chr and region.start <= record.start < region.end:
                            output_file.write(record)
                            break
        return self.callback(output_filename)


class IntervalListDownsampler(BaseTransformer):
    """
    Downsamples files with a SAM-style header (i.e., lines starting with `@`) and 1-based closed
    intervals in their first three columns; e.g., interval lists, and read counts files
    (`*.counts.tsv.gz`, created by GATK CollectReadCounts) which also have a columns header line.
    """
    # If True, the intervals partially overlapping the target regions are clipped to the regions;
    # otherwise, only the intervals entirely within the target regions are kept.
    clip_intervals = True

    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".interval_list", ".tsv.gz"]

    def is_kept(self, chr: str, start: int, end: int, regions: List[Region]) -> bool:
        if self.clip_intervals:
            return any(chr == r.chr and max(start, r.start) < min(end, r.end) for r in regions)
        return any(chr == r.chr and r.start <= start and end <= r.end for r in regions)

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        is_gz = Path(input_filename).suffix == ".gz"
        if is_gz:
            output_filename = output_filename.removesuffix(".gz")
        header_lines_count = 0
        # Note that this algorithm is not efficient.
        with gzip.open(input_filename, "rt") if is_gz else open(input_filename, "r") as input_file, \
                open(output_filename, "w") as output_file:
            for line in input_file:
                cols = line.rstrip("\n").split("\t")
                if line.startswith("@") or not cols[1].isdigit():
                    # A SAM-style header line, or a columns header line (e.g., `CONTIG START END COUNT`).
                    if self.mini_reference is None or not line.startswith("@SQ"):
                        output_file.write(line)
                        header_lines_count += 1
                    if self.mini_reference is not None and line.startswith("@HD"):
                        sequence_lines = self.mini_reference.get_sam_header_sequence_lines()
                        output_file.writelines(sequence_lines)
                        header_lines_count += len(sequence_lines)
                elif self.mini_reference is not None:
                    # Intervals are 1-based and closed.
                    start, end = int(cols[1]) - 1, int(cols[2])
                    if self.clip_intervals:
                        lifted_intervals = self.mini_reference.clip_and_lift(cols[0], start, end)
                    else:
                        lifted_interval = self.mini_reference.lift(cols[0], start, end)
                        lifted_intervals = [lifted_interval] if lifted_interval is not None else []
                    for start, end in lifted_intervals:
                        output_file.write("\t".join([cols[0], str(start + 1), str(end)] + cols[3:]) + "\n")
                elif self.is_kept(cols[0], int(cols[1]) - 1, int(cols[2]), regions):
                    output_file.write(line)

        if is_gz:
            output_gz_filename = output_filename + ".gz"
            pysam.tabix_compress(output_filename, output_gz_filename, force=True)
            os.remove(output_filename)
            pysam.tabix_index(output_gz_filename, seq_col=0, start_col=1, end_col=2,
                              line_skip=header_lines_count, zerobased=False, force=True)
            return self.callback(output_gz_filename, output_gz_filename + ".tbi")
        return self.callback(output_filename)


class BedDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".bed", ".bed.gz"]

    def read_lines(self, input_filename: str, regions: List[Region]):
        """
        Yields the lines of the given BED file. If the file is bgzipped and indexed,
        it only yields the header lines and the lines overlapping the given regions.
        """
        if not (Path(input_filename).suffix == ".gz" and os.path.isfile(f"{input_filename}.tbi")):
            # Note that this algorithm is not efficient.
            with gzip.open(input_filename, "rt") if Path(input_filename).suffix == ".gz" \
                    else open(input_filename, "r") as input_file:
                yield from input_file
            return

        with pysam.TabixFile(input_filename) as tabix_file:
            for line in tabix_file.header:
                yield line + "\n"
            for i, region in enumerate(regions):
                if region.chr not in tabix_file.contigs:
                    continue
                for line in tabix_file.fetch(region.chr, region.start, region.end):
                    cols = line.split("\t")
                    start, end = int(cols[1]), int(cols[2])
                    # Skip the lines already yielded for a previous region they overlap with.
                    if any(r.chr == region.chr and max(start, r.start) < min(end, r.end) for r in regions[:i]):
                        continue
                    yield line + "\n"

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        is_gz = Path(input_filename).suffix == ".gz"
        if is_gz:
            output_filename = output_filename.removesuffix(".gz")
        # Regions are processed in the order of the mini reference, so the output is sorted accordingly.
        ordered_regions = self.mini_reference.regions if self.mini_reference is not None else regions
        with open(output_filename, "w") as output_file:
            for line in self.read_lines(input_filename, ordered_regions):
                if line.startswith("#"):
                    output_file.write(line)
                    continue
                cols = line.rstrip("\n").split("\t")
                chr, start, end = cols[0], int(cols[1]), int(cols[2])
                if self.mini_reference is not None:
                    for lifted_start, lifted_end in self.mini_reference.clip_and_lift(chr, start, end):
                        output_file.write("\t".join([chr, str(lifted_start), str(lifted_end)] + cols[3:]) + "\n")
                    continue
                for region in regions:
                    if chr == region.chr and max(start, region.start) < min(end, region.end):
                        output_file.write(line)
                        break

        if is_gz:
            output_gz_filename = output_filename + ".gz"
            pysam.tabix_compress(output_filename, output_gz_filename, force=True)
            os.remove(output_filename)
            pysam.tabix_index(output_gz_filename, preset="bed", force=True)
            return self.callback(output_gz_filename, output_gz_filename + ".tbi")
        return self.callback(output_filename)


class PrimaryContigsDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str, str], dict], delimiter: str = "\t", **kwargs):
        super().__init__(working_dir, callback, **kwargs)
        self.delimiter = delimiter

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return []

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        include_chrs = set([r.chr for r in regions])

        is_gz = True if pathlib.Path(input_filename).suffix == ".gz" else False
        input_file = gzip.open(input_filename, "rt") if is_gz else open(input_filename, "r")

        if is_gz:
            output_filename = output_filename.removesuffix(".gz")

        with open(output_filename, "w") as out_file:
            for line in input_file:
                cols = line.strip().split(self.delimiter)
                if cols[0] in include_chrs:
                    out_file.write(self.delimiter.join(cols) + "\n")
        input_file.close()

        if is_gz:
            _, tmp_sorted_output_filename = tempfile.mkstemp(prefix="tmp_", dir=self.working_dir)
            sort_command = ["sort", "-k1,1", "-k2,2n", output_filename, "-o", tmp_sorted_output_filename]
            _ = subprocess.run(sort_command, capture_output=True, text=True, check=True)

            os.remove(output_filename)
            os.rename(tmp_sorted_output_filename, output_filename)

            output_gz_filename = output_filename + ".gz"
            bgzip_command = ["bgzip", "-f", output_filename, "-o", output_gz_filename]
            _ = subprocess.run(bgzip_command, check=True, capture_output=True)
            pysam.tabix_index(output_gz_filename, preset="bed", force=True)

            return self.callback(output_gz_filename, output_gz_filename + ".tbi")
        else:
            return self.callback(output_filename, None)


class FaiDownsampler(BaseTransformer):
    """
    Downsamples a FASTA index (or a file in the same format, e.g., a list of primary contigs and their lengths)
    to the contigs of the target regions. If a mini reference is given, the contig lengths are updated accordingly.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".fai"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        include_chrs = set([r.chr for r in regions])
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                cols = line.rstrip("\n").split("\t")
                if cols[0] not in include_chrs:
                    continue
                if self.mini_reference is not None:
                    cols[1] = str(self.mini_reference.contig_lengths[cols[0]])
                output_file.write("\t".join(cols) + "\n")
        return self.callback(output_filename)


class CnmopsContigsFaiDownsampler(FaiDownsampler):
    """
    Downsamples a FASTA index of the contigs that cn.mops runs on (e.g., the `chrom_file` of sv-shell's cnmops),
    as `FaiDownsampler` does, and if `cnmops_contigs` is given, further restricts it to those contigs.
    Since cn.mops runs on each contig separately, restricting its contigs shortens its runtime.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], cnmops_contigs: List[str] = None, **kwargs):
        super().__init__(working_dir, callback, **kwargs)
        self.cnmops_contigs = cnmops_contigs

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        if self.cnmops_contigs:
            regions = [r for r in regions if r.chr in self.cnmops_contigs]
        return super().transform(input_filename, output_prefix, regions, **kwargs)


class ReferenceDownsampler(BaseTransformer):
    """
    Writes the mini reference built from the target regions, and creates its index,
    sequence dictionary, and BWA index (requires `bwa` installed and added to PATH).
    """
    def __init__(self, working_dir, callback: Callable[[str, str, str, str, str, str, str, str, str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".fasta", ".fa"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        if self.mini_reference is None:
            raise ValueError("Downsampling a reference requires a mini reference.")
        output_filename = self.get_output_filename(input_filename, output_prefix)
        shutil.copyfile(self.mini_reference.fasta_filename, output_filename)
        pysam.faidx(output_filename)

        dict_filename = str(Path(output_filename).with_suffix(".dict"))
        self.mini_reference.write_sequence_dict(dict_filename)

        # The `-6` option names the index files as `<fasta>.64.*`, consistent with the GATK hg38 resource bundle.
        subprocess.run(["bwa", "index", "-6", output_filename], check=True, capture_output=True)
        # None of the contigs of the mini reference is an ALT contig, hence an empty ALT file.
        bwa_prefix = f"{output_filename}.64"
        open(f"{bwa_prefix}.alt", "w").close()

        return self.callback(
            output_filename, f"{output_filename}.fai", dict_filename,
            f"{bwa_prefix}.alt", f"{bwa_prefix}.amb", f"{bwa_prefix}.ann",
            f"{bwa_prefix}.bwt", f"{bwa_prefix}.pac", f"{bwa_prefix}.sa")


class SequenceDictDownsampler(BaseTransformer):
    """
    Writes the sequence dictionary of the mini reference built from the target regions;
    the output filename follows the same naming as that of `ReferenceDownsampler`.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".dict"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        if self.mini_reference is None:
            raise ValueError("Downsampling a sequence dictionary requires a mini reference.")
        output_filename = self.get_output_filename(input_filename, output_prefix)
        self.mini_reference.write_sequence_dict(output_filename)
        return self.callback(output_filename)


class GcnvReadCountsDownsampler(IntervalListDownsampler):
    """
    Downsamples read counts files that are inputs to gCNV in case mode. Since gCNV requires
    the intervals of the read counts to match the intervals of its models, only the intervals
    entirely within the target regions are kept (as in `GcnvModelDownsampler`), instead of clipping them.
    """
    clip_intervals = False


class GcnvModelDownsampler(BaseTransformer):
    """
    Downsamples a gCNV model tarball; i.e., a model shard created by GATK GermlineCNVCaller,
    or a contig-ploidy model created by GATK DetermineGermlineContigPloidy, in cohort mode.

    The model's intervals (`interval_list.tsv`) are subset to the intervals entirely within the target regions,
    and the model's arrays are subset accordingly; i.e., the arrays whose first dimension is the number of
    intervals (e.g., `mu_log_mean_bias_t.tsv`) are subset to the kept intervals, and the arrays whose first
    dimension is the number of contigs (e.g., `mu_psi_j_log__.tsv`) are subset to the contigs of the kept
    intervals, in the order they appear in `interval_list.tsv` (as gcnvkernel orders them).

    Returns None (i.e., no output) if none of the model's intervals is within the target regions.
    """
    interval_list_filename = "interval_list.tsv"
    contig_ploidy_prior_filename = "contig_ploidy_prior.tsv"

    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".tar.gz"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs):
        with tempfile.TemporaryDirectory(prefix="tmp_gcnv_model_", dir=self.working_dir) as tmp_dir:
            with tarfile.open(input_filename, "r:gz") as tar:
                tar.extractall(tmp_dir)

            kept_indices, kept_contigs, contigs = self.subset_interval_list(
                os.path.join(tmp_dir, self.interval_list_filename), regions)
            if len(kept_indices) == 0:
                logging.info(f"Skipping {input_filename}, since none of its intervals is within the target regions.")
                return None
            kept_contig_indices = [contigs.index(c) for c in kept_contigs]

            for filename in os.listdir(tmp_dir):
                path = os.path.join(tmp_dir, filename)
                if filename == self.contig_ploidy_prior_filename:
                    self.subset_contig_ploidy_prior(path, kept_contigs)
                elif filename.endswith(".tsv") and filename != self.interval_list_filename:
                    self.subset_array(path, {len(contigs): kept_contig_indices, self.intervals_count: kept_indices})

            output_filename = self.get_output_filename(input_filename, output_prefix)
            with tarfile.open(output_filename, "w:gz") as tar:
                tar.add(tmp_dir, arcname=".")
        return self.callback(output_filename)

    def subset_interval_list(self, filename: str, regions: List[Region]):
        """
        Subsets (and lifts over if a mini reference is given) the intervals in place, and returns the indices
        of the kept intervals, the contigs of the kept intervals, and all the contigs, in the order they appear.
        """
        with open(filename, "r") as f:
            lines = f.readlines()
        kept_indices, kept_contigs, contigs = [], [], []
        index = 0
        with open(filename, "w") as f:
            for line in lines:
                cols = line.rstrip("\n").split("\t")
                if line.startswith("@") or not cols[1].isdigit():
                    if self.mini_reference is None or not line.startswith("@SQ"):
                        f.write(line)
                    if self.mini_reference is not None and line.startswith("@HD"):
                        f.writelines(self.mini_reference.get_sam_header_sequence_lines())
                    continue
                if cols[0] not in contigs:
                    contigs.append(cols[0])
                # Intervals are 1-based and closed.
                start, end = int(cols[1]) - 1, int(cols[2])
                if any(cols[0] == r.chr and r.start <= start and end <= r.end for r in regions):
                    if self.mini_reference is not None:
                        lifted_interval = self.mini_reference.lift(cols[0], start, end)
                        cols[1], cols[2] = str(lifted_interval[0] + 1), str(lifted_interval[1])
                    f.write("\t".join(cols) + "\n")
                    kept_indices.append(index)
                if kept_indices and kept_indices[-1] == index and cols[0] not in kept_contigs:
                    kept_contigs.append(cols[0])
                index += 1
        self.intervals_count = index
        if self.intervals_count == len(contigs):
            raise ValueError(f"Cannot distinguish per-interval and per-contig arrays in {filename}.")
        return kept_indices, kept_contigs, contigs

    @staticmethod
    def subset_array(filename: str, indices_by_size: dict):
        """
        Subsets the rows of an array serialized by gcnvkernel (i.e., a header with `@CO shape:(...)`,
        `@CO dtype:...`, and a columns header line, followed by the rows) in place, given the indices
        of the rows to keep for each possible size of the array's first dimension.
        """
        with open(filename, "r") as f:
            lines = f.readlines()
        shape_line_index = next(i for i, line in enumerate(lines) if line.startswith("@CO\tshape:"))
        size = int(re.search(r"shape:\((\d+)", lines[shape_line_index]).group(1))
        if size not in indices_by_size:
            return
        header_lines_count = next(i for i, line in enumerate(lines) if not line.startswith("@")) + 1
        rows = lines[header_lines_count:]
        if len(rows) != size:
            raise ValueError(f"Expected {size} rows in {filename}, found {len(rows)}.")
        indices = indices_by_size[size]
        lines[shape_line_index] = re.sub(r"shape:\(\d+", f"shape:({len(indices)}", lines[shape_line_index])
        with open(filename, "w") as f:
            f.writelines(lines[:header_lines_count])
            f.writelines(rows[i] for i in indices)

    @staticmethod
    def subset_contig_ploidy_prior(filename: str, contigs: List[str]):
        with open(filename, "r") as f:
            lines = f.readlines()
        with open(filename, "w") as f:
            f.write(lines[0])
            f.writelines(line for line in lines[1:] if line.split("\t")[0] in contigs)


class SvVcfDownsampler(BaseTransformer):
    """
    Downsamples a bgzipped SV VCF file (e.g., the genotyped PESR or depth VCFs of a batch) to the records
    entirely within the target regions, and lifts over their coordinates if a mini reference is given.
    The VCFs are processed as text, since the END of a breakend (BND) record is a position on another
    contig (i.e., the INFO field CHR2), which may also be encoded in the ALT allele (e.g., `N[chr2:100[`).
    """
    bnd_alt_regex = re.compile(r"([\[\]])([^:\[\]]+):(\d+)([\[\]])")
    # The intervals of complex SVs in the INFO fields CPX_INTERVALS and SOURCE (e.g., `DUP_chr1:100-200`).
    cpx_interval_regex = re.compile(r"^([A-Za-z]+)_([^:]+):(\d+)-(\d+)$")
    cpx_interval_info_fields = ["CPX_INTERVALS", "SOURCE"]

    def __init__(self, working_dir, callback: Callable[[str, str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".vcf.gz"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        self.downsample_vcf(input_filename, output_filename, regions)
        pysam.tabix_index(output_filename, preset="vcf", force=True)
        return self.callback(output_filename, output_filename + ".tbi")

    def lift_interval(self, chr: str, start: int, end: int, regions: List[Region]):
        """
        Returns the given 0-based half-open interval on the mini reference (or unchanged if no mini reference
        is given) if it is entirely within a target region, otherwise None.
        """
        if self.mini_reference is not None:
            return self.mini_reference.lift(chr, start, end)
        if any(chr == r.chr and r.start <= start and end <= r.end for r in regions):
            return start, end
        return None

    def downsample_vcf(self, input_filename: str, output_filename: str, regions: List[Region]):
        """
        Downsamples the given bgzipped VCF file, and writes the output bgzipped;
        the input and output files can be the same (i.e., downsampling in place).
        """
        uncompressed_filename = output_filename.removesuffix(".gz") + ".tmp"
        with gzip.open(input_filename, "rt") as input_file, open(uncompressed_filename, "w") as output_file:
            for line in input_file:
                if line.startswith("##contig") and self.mini_reference is not None:
                    continue
                if line.startswith("#CHROM") and self.mini_reference is not None:
                    for name, length in self.mini_reference.contig_lengths.items():
                        output_file.write(f"##contig=<ID={name},length={length}>\n")
                if line.startswith("#"):
                    output_file.write(line)
                    continue
                lifted_line = self.lift_record(line.rstrip("\n").split("\t"), regions)
                if lifted_line is not None:
                    output_file.write(lifted_line + "\n")
        pysam.tabix_compress(uncompressed_filename, output_filename, force=True)
        os.remove(uncompressed_filename)

    def lift_record(self, cols: List[str], regions: List[Region]):
        """
        Returns the given VCF record with its coordinates lifted over, or None
        if any of its positions are not within the target regions.
        """
        chr, pos = cols[0], int(cols[1])
        info = [x.split("=", maxsplit=1) for x in cols[7].split(";")]
        info_dict = {x[0]: x[1] for x in info if len(x) == 2}
        end = int(info_dict.get("END", pos))
        updated_info = {}

        for field in self.cpx_interval_info_fields:
            if field in info_dict:
                lifted_intervals = []
                for interval in info_dict[field].split(","):
                    m = self.cpx_interval_regex.match(interval)
                    if m is None:
                        # e.g., the SOURCE of an insertion that is not a complex SV.
                        lifted_intervals.append(interval)
                        continue
                    lifted_interval = self.lift_interval(m.group(2), int(m.group(3)) - 1, int(m.group(4)), regions)
                    if lifted_interval is None:
                        return None
                    lifted_intervals.append(f"{m.group(1)}_{m.group(2)}:{lifted_interval[0] + 1}-{lifted_interval[1]}")
                updated_info[field] = ",".join(lifted_intervals)

        if "END2" in info_dict:
            # GATK-style records with a second breakpoint at END2 on CHR2 (e.g., BND and CTX);
            # the END of these records is on the same contig as POS.
            lifted_interval = self.lift_interval(chr, pos - 1, max(pos, end), regions)
            end2 = int(info_dict["END2"])
            lifted_end2 = self.lift_interval(info_dict.get("CHR2", chr), end2 - 1, end2, regions)
            if lifted_interval is None or lifted_end2 is None:
                return None
            new_pos = lifted_interval[0] + 1
            new_end = new_pos + (end - pos)
            updated_info["END2"] = str(lifted_end2[1])
        elif info_dict.get("SVTYPE") == "BND":
            # The END of a breakend is a position on CHR2 if set, otherwise, on the contig of POS.
            lifted_pos = self.lift_interval(chr, pos - 1, pos, regions)
            lifted_end = self.lift_interval(info_dict.get("CHR2", chr), end - 1, end, regions)
            if lifted_pos is None or lifted_end is None:
                return None
            new_pos, new_end = lifted_pos[1], lifted_end[1]
            # The mate in the ALT allele (e.g., `N[chr2:100[`) is lifted on its own, since a breakend
            # can encode its mate only in ALT (e.g., the output of `update_variant_representations.py`).
            alt_match = self.bnd_alt_regex.search(cols[4])
            if alt_match is not None:
                alt_pos = int(alt_match.group(3))
                lifted_alt = self.lift_interval(alt_match.group(2), alt_pos - 1, alt_pos, regions)
                if lifted_alt is None:
                    return None
                cols[4] = self.bnd_alt_regex.sub(
                    lambda m: f"{m.group(1)}{m.group(2)}:{lifted_alt[1]}{m.group(4)}", cols[4])
        else:
            lifted_interval = self.lift_interval(chr, pos - 1, max(pos, end), regions)
            if lifted_interval is None:
                return None
            new_pos = lifted_interval[0] + 1
            new_end = new_pos + (end - pos)

        cols[1] = str(new_pos)
        if "END" in info_dict:
            updated_info["END"] = str(new_end)
        cols[7] = ";".join(f"{x[0]}={updated_info[x[0]]}" if x[0] in updated_info else "=".join(x) for x in info)
        return "\t".join(cols)


class SvVcfTarDownsampler(SvVcfDownsampler):
    """
    Downsamples a tarball of SV VCF files (e.g., the standardized Manta, Wham, or Scramble VCFs of a batch);
    see `SvVcfDownsampler` for details. If the reference panel is subset, the VCFs of the other samples are dropped.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".tar.gz"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        with tempfile.TemporaryDirectory(prefix="tmp_sv_vcf_tar_", dir=self.working_dir) as tmp_dir:
            with tarfile.open(input_filename, "r:gz") as tar:
                tar.extractall(tmp_dir)
            for root, _, filenames in os.walk(tmp_dir):
                for filename in filenames:
                    if not filename.endswith(".vcf.gz"):
                        continue
                    vcf_filename = os.path.join(root, filename)
                    if self.is_excluded_sample_file(filename):
                        os.remove(vcf_filename)
                        continue
                    self.downsample_vcf(vcf_filename, vcf_filename, regions)
            output_filename = self.get_output_filename(input_filename, output_prefix)
            with tarfile.open(output_filename, "w:gz") as tar:
                tar.add(tmp_dir, arcname=".")
        return self.callback(output_filename)


class PeEvidenceDownsampler(BaseTransformer):
    """
    Downsamples a bgzipped and indexed discordant pair (PE) evidence file (i.e., created by GATK
    CollectSVEvidence, or merged for a batch), with the columns `contig, position, strand, mate contig,
    mate position, mate strand, sample` (positions are 0-based), to the pairs with both ends within the target regions, and lifts
    over both ends if a mini reference is given. Only the lines overlapping the target regions are read
    using the file's index, hence, it is efficient for large (e.g., multi-gigabyte) files.
    """
    def __init__(self, working_dir, callback: Callable[[str, str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".txt.gz"]

    def lift_position(self, chr: str, pos: int, regions: List[Region]):
        if self.mini_reference is not None:
            lifted_interval = self.mini_reference.lift(chr, pos, pos + 1)
            return None if lifted_interval is None else lifted_interval[0]
        return pos if any(chr == r.chr and r.start <= pos < r.end for r in regions) else None

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        uncompressed_filename = output_filename.removesuffix(".gz")
        # Regions are processed in the order of the mini reference, so the output is sorted accordingly.
        ordered_regions = self.mini_reference.regions if self.mini_reference is not None else regions
        with pysam.TabixFile(input_filename) as tabix_file, open(uncompressed_filename, "w") as output_file:
            for region in ordered_regions:
                if region.chr not in tabix_file.contigs:
                    continue
                for line in tabix_file.fetch(region.chr, region.start, region.end):
                    cols = line.split("\t")
                    pos = self.lift_position(cols[0], int(cols[1]), regions)
                    mate_pos = self.lift_position(cols[3], int(cols[4]), regions)
                    if pos is None or mate_pos is None:
                        continue
                    cols[1], cols[4] = str(pos), str(mate_pos)
                    output_file.write("\t".join(cols) + "\n")
        pysam.tabix_compress(uncompressed_filename, output_filename, force=True)
        os.remove(uncompressed_filename)
        pysam.tabix_index(output_filename, seq_col=0, start_col=1, end_col=1, zerobased=True, force=True)
        return self.callback(output_filename, output_filename + ".tbi")


class IntronTableDownsampler(BaseTransformer):
    """
    Downsamples a gzipped table of introns (e.g., `gencode.v39.CDS.intron.tsv.gz`, used for annotating
    SVs in CleanVcf) with the columns `contig, source, start, end, gene` (1-based closed intervals),
    to the introns entirely within the target regions, and lifts them over if a mini reference is given.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".tsv.gz"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        uncompressed_filename = output_filename.removesuffix(".gz")
        with gzip.open(input_filename, "rt") as input_file, open(uncompressed_filename, "w") as output_file:
            for line in input_file:
                cols = line.rstrip("\n").split("\t")
                # Some coordinates are written in scientific notation (e.g., `4.2e+07`).
                start, end = int(float(cols[2])) - 1, int(float(cols[3]))
                if self.mini_reference is not None:
                    lifted_interval = self.mini_reference.lift(cols[0], start, end)
                    if lifted_interval is None:
                        continue
                    cols[2], cols[3] = str(lifted_interval[0] + 1), str(lifted_interval[1])
                elif not any(cols[0] == r.chr and r.start <= start and end <= r.end for r in regions):
                    continue
                output_file.write("\t".join(cols) + "\n")
        pysam.tabix_compress(uncompressed_filename, output_filename, force=True)
        os.remove(uncompressed_filename)
        return self.callback(output_filename)


class PositionEvidenceDownsampler(BaseTransformer):
    """
    Downsamples a bgzipped and indexed evidence file whose first two columns are a contig and a 0-based
    position, e.g., split read (SR) or site depth (SD) evidence files created by GATK CollectSVEvidence,
    to the lines within the target regions, and lifts over the positions if a mini reference is given.
    Only the lines overlapping the target regions are read using the file's index.
    """
    def __init__(self, working_dir, callback: Callable[[str, str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".txt.gz"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        uncompressed_filename = output_filename.removesuffix(".gz")
        # Regions are processed in the order of the mini reference, so the output is sorted accordingly.
        ordered_regions = self.mini_reference.regions if self.mini_reference is not None else regions
        with pysam.TabixFile(input_filename) as tabix_file, open(uncompressed_filename, "w") as output_file:
            for region in ordered_regions:
                if region.chr not in tabix_file.contigs:
                    continue
                for line in tabix_file.fetch(region.chr, region.start, region.end):
                    cols = line.split("\t")
                    pos = int(cols[1])
                    if not region.start <= pos < region.end:
                        continue
                    if self.mini_reference is not None:
                        cols[1] = str(self.mini_reference.lift(cols[0], pos, pos + 1)[0])
                    output_file.write("\t".join(cols) + "\n")
        pysam.tabix_compress(uncompressed_filename, output_filename, force=True)
        os.remove(uncompressed_filename)
        pysam.tabix_index(output_filename, seq_col=0, start_col=1, end_col=1, zerobased=True, force=True)
        return self.callback(output_filename, output_filename + ".tbi")


class GtfDownsampler(BaseTransformer):
    """
    Downsamples a GTF file (e.g., `gencode.v47.basic.protein_coding.canonical.gtf`) to the genes entirely
    within the target regions, keeping all the records (e.g., transcripts, exons, and CDS) of those genes,
    and lifts them over if a mini reference is given. Genes are kept or dropped as a whole, so the
    hierarchy of the records (gene, transcript, exon) of the kept genes is intact.
    """
    gene_id_regex = re.compile(r'gene_id "([^"]+)"')

    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".gtf"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        kept_genes = set()
        with open(input_filename, "r") as input_file:
            for line in input_file:
                cols = line.split("\t")
                if line.startswith("#") or cols[2] != "gene":
                    continue
                start, end = int(cols[3]) - 1, int(cols[4])
                if any(cols[0] == r.chr and r.start <= start and end <= r.end for r in regions):
                    kept_genes.add(self.gene_id_regex.search(cols[8]).group(1))

        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                if line.startswith("#"):
                    output_file.write(line)
                    continue
                cols = line.split("\t")
                if self.gene_id_regex.search(cols[8]).group(1) not in kept_genes:
                    continue
                if self.mini_reference is not None:
                    # Intervals are 1-based and closed.
                    start, end = self.mini_reference.lift(cols[0], int(cols[3]) - 1, int(cols[4]))
                    cols[3], cols[4] = str(start + 1), str(end)
                output_file.write("\t".join(cols))
        return self.callback(output_filename)


class GcnvModelShardDownsampler(GcnvModelDownsampler):
    """
    Downsamples a gCNV model shard, as `GcnvModelDownsampler` does, and if `gcnv_contigs` is given, keeps only
    the shard's intervals on those contigs; hence, the shards with no intervals on those contigs are dropped.
    Since gCNV runs on each model shard separately, with a mostly fixed runtime per shard, fewer shards
    shorten its runtime. Not to be used for the contig-ploidy model, which needs all the contigs.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], gcnv_contigs: List[str] = None, **kwargs):
        super().__init__(working_dir, callback, **kwargs)
        self.gcnv_contigs = gcnv_contigs

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs):
        if self.gcnv_contigs:
            regions = [r for r in regions if r.chr in self.gcnv_contigs]
        return super().transform(input_filename, output_prefix, regions, **kwargs)


class FileListDownsampler(BaseTransformer):
    """
    Downsamples a text file listing files (one path per line), e.g., the evidence files or the gCNV model
    shards of a reference panel, by downsampling each of the listed files using `item_transformer`,
    and writes a list of the downsampled files. The listed files that are not on the local filesystem
    (e.g., paths in a container) are resolved by their names in the directory of the list file.
    The listed files with no data within the target regions (i.e., `item_transformer` returns None)
    are dropped from the list. If `per_sample` and the reference panel is subset, the files of the other
    samples are dropped too.
    """
    item_transformer: type = None
    per_sample: bool = False

    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)
        self.item_transformer_instance = self.item_transformer(
            working_dir=working_dir, callback=lambda *x: x[0], **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".list"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        with open(input_filename, "r") as f:
            filenames = [line.strip() for line in f if line.strip()]
        downsampled_filenames = []
        for filename in filenames:
            if self.per_sample and self.is_excluded_sample_file(filename):
                continue
            if not os.path.isfile(filename):
                filename = os.path.join(os.path.dirname(input_filename), os.path.basename(filename))
            downsampled_filename = self.item_transformer_instance.transform(filename, output_prefix, regions)
            if downsampled_filename is not None:
                downsampled_filenames.append(downsampled_filename)
        with open(output_filename, "w") as f:
            f.writelines(x + "\n" for x in downsampled_filenames)
        return self.callback(output_filename)


class PeEvidenceListDownsampler(FileListDownsampler):
    item_transformer = PeEvidenceDownsampler
    per_sample = True


class PositionEvidenceListDownsampler(FileListDownsampler):
    item_transformer = PositionEvidenceDownsampler
    per_sample = True


class GcnvModelListDownsampler(FileListDownsampler):
    item_transformer = GcnvModelShardDownsampler


class SampleListDownsampler(BaseTransformer):
    """
    Subsets a text file with a sample ID per line (e.g., the samples of a reference panel), or a tab-delimited
    file with a sample ID per row in column `sample_column` (e.g., a PED file), to the kept reference panel samples.
    """
    sample_column: int = 0

    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".list", ".txt"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                cols = line.rstrip("\n").split("\t")
                if line.startswith("#") or not self.is_excluded_sample(cols[self.sample_column]):
                    output_file.write(line)
        return self.callback(output_filename)


class PedDownsampler(SampleListDownsampler):
    sample_column = 1

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".ped"]


class MedianCovDownsampler(BaseTransformer):
    """
    Subsets a table of the median coverage of samples (e.g., `medianCov.transposed.bed`, with a row of
    sample IDs followed by a row of their median coverages) to the kept reference panel samples.
    """
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".bed"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        with open(input_filename, "r") as f:
            rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
        kept_indices = [i for i, sample_id in enumerate(rows[0]) if not self.is_excluded_sample(sample_id)]
        with open(output_filename, "w") as f:
            for row in rows:
                f.write("\t".join(row[i] for i in kept_indices) + "\n")
        return self.callback(output_filename)


class BincovMatrixDownsampler(BedDownsampler):
    """
    Downsamples a bin coverage matrix (i.e., a BED file with a column of read counts per sample),
    as `BedDownsampler` does, and subsets its columns to the kept reference panel samples.
    """
    def read_lines(self, input_filename: str, regions: List[Region]):
        kept_indices = None
        for line in super().read_lines(input_filename, regions):
            cols = line.rstrip("\n").split("\t")
            if kept_indices is None:
                # The header line, e.g., `#Chr Start End HG00096 HG00129 ...`.
                kept_indices = [i for i, x in enumerate(cols) if i < 3 or not self.is_excluded_sample(x)]
            yield "\t".join(cols[i] for i in kept_indices) + "\n"


class RefPanelDepthCallsDownsampler(BedDownsampler):
    """
    Downsamples the depth calls (DEL or DUP) of a reference panel, as `BedDownsampler` does,
    and keeps only the calls of the kept reference panel samples (in the `sample` column).
    """
    def read_lines(self, input_filename: str, regions: List[Region]):
        for line in super().read_lines(input_filename, regions):
            if line.startswith("#") or not self.is_excluded_sample(line.split("\t")[4]):
                yield line
