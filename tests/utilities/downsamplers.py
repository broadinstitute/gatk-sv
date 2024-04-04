import pysam

from typing import Callable, List
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Region:
    chr: str
    start: int
    end: int


class BaseDownsampler:
    def __init__(self, working_dir: str, callback: Callable[[str, ...], dict]):
        # Convert the string to an ABS path and make sure it exists.
        self.working_dir = Path(working_dir).resolve(strict=True)
        self.callback = callback

    def get_output_filename(self, input_filename, output_prefix):
        return str(self.working_dir.joinpath(f"{output_prefix}{Path(input_filename).name}"))

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
    def __init__(self, working_dir, callback: Callable[[str, str], dict]):
        super().__init__(working_dir, callback)

    @staticmethod
    def get_supported_file_types() -> List[str]:
        return [".cram"]

    def downsample(self, input_filename: str, output_prefix: str, regions: List[Region]) -> dict:
        output_filename = self.get_output_filename(input_filename, output_prefix)
        with \
                pysam.AlignmentFile(input_filename, "rc") as input_cram_file, \
                pysam.AlignmentFile(output_filename, "wc",
                                    header=input_cram_file.header,
                                    reference_names=input_cram_file.references) as output_cram_file:
            for region in regions:
                for read in input_cram_file.fetch(region=f"{region.chr}:{region.start}-{region.end}"):
                    output_cram_file.write(read)
        index_filename = f"{output_filename}.crai"
        pysam.index(output_filename, index_filename)
        return self.callback(output_filename, index_filename)


class VcfDownsampler(BaseDownsampler):
    def __init__(self, working_dir, callback: Callable[[str], dict]):
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

    def __init__(self, working_dir, callback: Callable[[str], dict]):
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
        return self.callback(output_filename)


class PrimaryContigsDownsampler(BaseDownsampler):
    def __init__(self, working_dir, callback: Callable[[str], dict], delimiter: str = "\t"):
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
