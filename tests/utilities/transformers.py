import os
import pysam
import subprocess

from typing import Callable, List
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Region:
    chr: str
    start: int
    end: int


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


class BaseConverter(BaseTransformer):
    def __init__(self, working_dir: str, callback: Callable[[str, ...], dict]):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        raise NotImplementedError()

    def convert(self, input_filename: str, output_prefix: str) -> dict:
        raise NotImplementedError()


class BedToIntervalListConverter(BaseConverter):
    def __init__(self, working_dir, callback: Callable[[str], dict], sequence_dict_filename: str, picard_path: str, **kwargs):
        super().__init__(working_dir, callback)
        self.sequence_dict_filename = sequence_dict_filename
        self.picard_path = picard_path

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".interval_list"]

    def convert(self, input_filename: str, output_prefix: str) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        subprocess.run(
            ["java", "-jar", self.picard_path, "BedToIntervalList",
             "-I", input_filename, "-O", output_filename, "-SD", self.sequence_dict_filename],
            check=True)
        return self.callback(output_filename)



class BaseDownsampler(BaseTransformer):
    def __init__(self, working_dir: str, callback: Callable[[str, ...], dict]):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        """
        The file types should include the '.' prefix to match with Path().suffix output
        (e.g., it should return '.cram' instead of 'cram').
        """
        raise NotImplementedError()

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
        raise NotImplementedError()


class CramDownsampler(BaseDownsampler):
    def __init__(self, working_dir, callback: Callable[[str, str], dict], reference_fasta: str, reference_index: str, **kwargs):
        super().__init__(working_dir, callback)
        self.reference_fasta = reference_fasta
        self.reference_index = reference_index

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".cram"]

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
        # Implementation notes:
        # 1. This method calls `samtools` instead of using `pysam` since it needs to include
        #    distant pair-end reads in the downsampled regions (i.e., the `--fetch-pairs` flag).
        #    Such reads are needed for MELT to function properly.
        #
        # 2. The method writes target regions to a BED file. While taking the BED file containing
        #    the regions as input seems a better option, the current approach is implemented as
        #    taking a BED file as input instead of a list of regions would convolut the method signature.

        regions_filename = os.path.join(self.working_dir, "_tmp_regions.bed")
        with open(regions_filename, "w") as f:
            for region in regions:
                f.write("\t".join([str(region.chr), str(region.start), str(region.end)]) + "\n")

        output_filename = self.get_output_filename(input_filename, output_prefix)
        subprocess.run(
            ["samtools", "view", input_filename, "--reference", self.reference_fasta,
             "--targets-file", regions_filename, "--output", output_filename, "--cram", "--fetch-pairs"]
        )

        subprocess.run(["samtools", "index", output_filename])

        os.remove(regions_filename)
        return self.callback(output_filename, output_filename + ".fai")


class VcfDownsampler(BaseDownsampler):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".vcf"]

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
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


class IntervalListDownsampler(BaseDownsampler):
    # Implementation note:
    # An alternative to the down sampling approach implemented here is to take
    # a BED file containing target regions as input, and convert the BED to .interval_list
    # as the following.
    #
    # > java -jar picard.jar BedToIntervalList I=regions.bed O=regions.interval_list SD=/Homo_sapiens_assembly38.dict
    #
    # There are two downsides to converting BED to .interval_list:
    # 1. It needs the picard tool installed;
    # 2. It needs an additional "SD" input, which would make the `downsample` method signature complicated.

    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".interval_list"]

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
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


class BedDownsampler(BaseDownsampler):
    def __init__(self, working_dir, callback: Callable[[str], dict], **kwargs):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".bed"]

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
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


class PrimaryContigsDownsampler(BaseDownsampler):
    def __init__(self, working_dir, callback: Callable[[str], dict], delimiter: str = "\t", **kwargs):
        super().__init__(working_dir, callback)
        self.delimiter = delimiter

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return []

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        include_chrs = set([r.chr for r in regions])
        with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
            for line in input_file:
                cols = line.strip().split(self.delimiter)
                if cols[0] in include_chrs:
                    output_file.write(self.delimiter.join(cols) + "\n")
        return self.callback(output_filename)
