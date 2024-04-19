import os
import pysam
import subprocess
import uuid

from collections import defaultdict
from typing import Callable, List
from dataclasses import dataclass
from pathlib import Path
from tqdm import tqdm


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


class BaseTransformer:
    def __init__(self, working_dir: str, callback: Callable[[str, ...], dict]):
        # Convert the string to an ABS path and make sure it exists.
        self.working_dir = Path(working_dir).resolve(strict=True)
        self.callback = callback

    @staticmethod
    def get_supported_file_types() -> List[str]:
        """
        The file types should include the '.' prefix to match with Path().suffix output
        (e.g., it should return '.cram' instead of 'cram').
        """
        raise NotImplementedError()

    def get_output_filename(self, input_filename, output_prefix):
        return str(self.working_dir.joinpath(f"{output_prefix}{Path(input_filename).name}"))

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs):
        raise NotImplementedError()


class BedToIntervalListConverter(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], sequence_dict_filename: str, picard_path: str, **kwargs):
        super().__init__(working_dir, callback)
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
        super().__init__(working_dir, callback)
        self.reference_fasta = reference_fasta
        self.reference_index = reference_index

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".cram"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename_unsorted = self.get_output_filename(input_filename, f"unsorted_{output_prefix}")
        with pysam.AlignmentFile(input_filename, "rc") as cram:
            header = cram.header
            references = cram.references

        reads_for_second_pass = {}
        with pysam.AlignmentFile(output_filename_unsorted, "wc", header=header, reference_names=references) as output_cram_file:
            for region in regions:
                generator = self.read_pair_generator(input_filename, region=f"{region.chr}:{region.start}-{region.end}")
                try:
                    while True:
                        r1, r2 = next(generator)
                        output_cram_file.write(r1)
                        output_cram_file.write(r2)
                except StopIteration as e:
                    reads_for_second_pass = reads_for_second_pass | e.value

            for r1, r2 in self.get_distant_read_pairs(input_filename, reads_for_second_pass):
                output_cram_file.write(r1)
                output_cram_file.write(r2)


                # for r1, r2 in self.read_pair_generator(input_filename, region=f"{region.chr}:{region.start}-{region.end}"):
                #     output_cram_file.write(r1)
                #     output_cram_file.write(r2)

        output_filename = self.get_output_filename(input_filename, output_prefix)
        pysam.sort("-o", output_filename, output_filename_unsorted)
        os.remove(output_filename_unsorted)
        index_filename = f"{output_filename}.crai"
        pysam.index(output_filename, index_filename)
        return self.callback(output_filename, index_filename)

    @staticmethod
    def read_pair_generator(filename, region=None):
        """
        Generate read pairs in a CRAM file.
        If the `region` string is provided, it only generates the reads overlapping this region.


        """
        read_dict = {}
        secondary_dict = {}
        with pysam.AlignmentFile(filename, "rc") as cram:
            for read in cram.fetch(region=region):
                q_name = read.query_name
                if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                    # Maybe it would be better to exclude such reads
                    if q_name not in secondary_dict:
                        secondary_dict[q_name] = [None, None]
                    secondary_dict[q_name][0 if read.is_read1 else 1] = read
                    continue
                if q_name not in read_dict:
                    read_dict[q_name] = [None, None]
                    read_dict[q_name][0 if read.is_read1 else 1] = read
                else:
                    if read.is_read1:
                        yield read, read_dict[q_name][1]
                    else:
                        yield read_dict[q_name][0], read
                    del read_dict[q_name]

        return read_dict #| secondary_dict

    @staticmethod
    def get_distant_read_pairs(filename, reads):
        with pysam.AlignmentFile(filename, "rc") as cram:
            for read in tqdm(cram.fetch(), desc="Iterating on reads", unit=" read", dynamic_ncols=True,
                             bar_format="{desc}: {n:,} [{elapsed}] {rate_fmt}"):
                query_name = read.query_name
                pair = reads.get(query_name)
                if pair is not None:
                    if (read.is_read1 and pair[0] is not None) or (not read.is_read1 and pair[1] is not None):
                        # It is the same read as the one seen before.
                        continue
                    if read.is_read1:
                        yield read, pair[1]
                    else:
                        yield pair[0], read
                    del reads[query_name]
        print(f"\n\n\n size at the end: {len(reads)}")



class VcfDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".vcf"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        with \
                pysam.VariantFile(input_filename) as input_file, \
                pysam.VariantFile(output_filename, "w", header=input_file.header) as output_file:
            for record in input_file:
                for region in regions:
                    if record.contig == region.chr and region.start <= record.pos <= region.end:
                        output_file.write(record)
                        break
        return self.callback(output_filename)


class IntervalListDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".interval_list"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        # Note that this algorithm is not efficient.
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                if line.startswith("@"):
                    output_file.write(line)
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
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".bed"]

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        # Note that this algorithm is not efficient.
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                cols = line.rstrip().split()
                chr, start, end = cols[0], int(cols[1]), int(cols[2])
                for region in regions:
                    if chr == region.chr and max(start, region.start) < min(end, region.end):
                        output_file.write(line)
                        break
        return self.callback(output_filename)


class PrimaryContigsDownsampler(BaseTransformer):
    def __init__(self, working_dir, callback: Callable[[str], dict], delimiter: str = "\t", **kwargs):
        super().__init__(working_dir, callback)
        self.delimiter = delimiter

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return []

    def transform(self, input_filename: str, output_prefix: str, regions: List[Region], **kwargs) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        include_chrs = set([r.chr for r in regions])
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                cols = line.strip().split(self.delimiter)
                if cols[0] in include_chrs:
                    output_file.write(self.delimiter.join(cols) + "\n")
        return self.callback(output_filename)
