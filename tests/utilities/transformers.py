import gzip
import hashlib
import logging
import os
import pathlib
import pysam
import shutil
import subprocess
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

    def lift(self, chr: str, start: int, end: int):
        """
        Returns the coordinates of the given interval on the mini reference,
        or None if the interval is not entirely within a target region.
        """
        for region, offset in zip(self.regions, self.offsets):
            if chr == region.chr and region.start <= start and end <= region.end:
                return offset + start - region.start, offset + end - region.start
        return None

    def clip_and_lift(self, chr: str, start: int, end: int) -> List[tuple]:
        """
        Returns the coordinates on the mini reference of the parts of the given interval
        overlapping target regions; the parts adjacent on the mini reference are merged.
        """
        lifted = []
        for region, offset in zip(self.regions, self.offsets):
            if chr == region.chr and max(start, region.start) < min(end, region.end):
                lifted_start = offset + max(start, region.start) - region.start
                lifted_end = offset + min(end, region.end) - region.start
                if lifted and lifted[-1][1] == lifted_start:
                    lifted[-1] = (lifted[-1][0], lifted_end)
                else:
                    lifted.append((lifted_start, lifted_end))
        return lifted


class BaseTransformer:
    def __init__(self, working_dir: str, callback: Callable[[str, ...], dict], input_root: str = None,
                 mini_reference: MiniReference = None, **kwargs):
        # Convert the string to an ABS path and make sure it exists.
        self.working_dir = Path(working_dir).resolve(strict=True)
        self.callback = callback
        # When set, the output files are written to the working directory following
        # the same directory hierarchy as the input files have relative to this directory.
        self.input_root = Path(input_root).resolve() if input_root else None
        # When set, the coordinates in the output files are lifted over to this reference.
        self.mini_reference = mini_reference

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
        if self.input_root is not None and input_dir.is_relative_to(self.input_root):
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
                        input_filename, region=f"{region.chr}:{region.start}-{region.end}",
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
            logging.warning(f"Did not find discordant pairs for {len(read)} reads.")
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
                        if record.contig == region.chr and region.start <= record.pos <= region.end:
                            output_file.write(record)
                            break
        return self.callback(output_filename)


class IntervalListDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".interval_list"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        # Note that this algorithm is not efficient.
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                if line.startswith("@"):
                    if self.mini_reference is None or not line.startswith("@SQ"):
                        output_file.write(line)
                    if self.mini_reference is not None and line.startswith("@HD"):
                        for sq in self.mini_reference.get_sam_header_sequences():
                            output_file.write(f"@SQ\tSN:{sq['SN']}\tLN:{sq['LN']}\tM5:{sq['M5']}\n")
                elif self.mini_reference is not None:
                    # Intervals are 1-based and closed.
                    cols = line.rstrip("\n").split("\t")
                    for start, end in self.mini_reference.clip_and_lift(cols[0], int(cols[1]) - 1, int(cols[2])):
                        output_file.write("\t".join([cols[0], str(start + 1), str(end)] + cols[3:]) + "\n")
                else:
                    cols = line.rstrip().split()
                    chr, start, end = cols[0], int(cols[1]), int(cols[2])
                    for region in regions:
                        if chr == region.chr and max(start, region.start) < min(end, region.end):
                            output_file.write(line)
                            break
        return self.callback(output_filename)


class BedDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback, **kwargs)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".bed", ".bed.gz"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        is_gz = Path(input_filename).suffix == ".gz"
        if is_gz:
            output_filename = output_filename.removesuffix(".gz")
        # Note that this algorithm is not efficient.
        with gzip.open(input_filename, "rt") if is_gz else open(input_filename, "r") as input_file, \
                open(output_filename, "w") as output_file:
            for line in input_file:
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
        pysam.dict("-o", dict_filename, output_filename)

        # The `-6` option names the index files as `<fasta>.64.*`, consistent with the GATK hg38 resource bundle.
        subprocess.run(["bwa", "index", "-6", output_filename], check=True, capture_output=True)
        # None of the contigs of the mini reference is an ALT contig, hence an empty ALT file.
        bwa_prefix = f"{output_filename}.64"
        open(f"{bwa_prefix}.alt", "w").close()

        return self.callback(
            output_filename, f"{output_filename}.fai", dict_filename,
            f"{bwa_prefix}.alt", f"{bwa_prefix}.amb", f"{bwa_prefix}.ann",
            f"{bwa_prefix}.bwt", f"{bwa_prefix}.pac", f"{bwa_prefix}.sa")
