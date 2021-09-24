"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

import enum
import os
import svtk.utils as svu
import sys
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path
from pysam.libcbcf import VariantFile, VariantRecord
from subprocess import check_call
from typing import Iterator


SVTYPE = "SVTYPE"
BLANK_SAMPLES = "Blank_Samples"


class SVType:
    DUP = "DUP"
    DEL = "DEL"


class VariantFormatTypes:
    # Predicted copy state
    RD_CN = "RD_CN"
    # Classes of evidence supporting final genotype
    EV = "EV"


@enum.unique
class WindowTypes(Enum):

    """
    Revise Nested CNVs.
    """
    RNCNV = 1


class BaseWindow(ABC):
    """
    A 'window' encapsulates a set of overlapping intervals.
    In the following illustration, the five intervals to left
    form Window_1, and the two interval to right form Window_2
    (Window_1 and Window_2 are different instances of BaseWindow).

    ----------- Window_1 -------------       ---- Window_2 ----
    ░░░░░░░░░░░░░░░░░░░░░░░   ░░░░░░░░       ░░░░░░░░░░░░
       ░░░░░░  ░░░░░░  ░░░░░░░░░                   ░░░░░░░░░░░░
    """
    def __init__(self):
        self.intervals = []
        self.overlapping_pairs = []
        self.stop = 0

    @staticmethod
    def _determine_wider(x: VariantRecord, y: VariantRecord) \
            -> (VariantRecord, VariantRecord):
        """
        Returns the input intervals ordered as: wider, narrower
        """
        if x.stop - x.start >= y.stop - y.start:
            return x, y
        else:
            return y, x

    def get_coverage(self, x: VariantRecord, y: VariantRecord) -> float:
        """
        Returns the Jaccard index of the overlapping intervals;
        returns 0.0 if the intervals do not overlap.
        """
        w, n = self._determine_wider(x, y)
        if w.start <= n.stop and n.start <= w.stop:
            intersection = min(n.stop, w.stop) - max(n.start, w.start)
            return intersection / (n.stop - n.start)
        else:
            return 0.0

    def add(self, interval: VariantRecord) -> None:
        self.intervals.append(interval)
        self.stop = max(self.stop, interval.stop)

    @abstractmethod
    def add_overlapping_pairs(self, x: VariantRecord, y: VariantRecord) \
            -> None:
        """
        Implement application-specific logic for pairs of overlapping
        intervals in the window.
        Should be implemented in a derived type.
        """
        raise NotImplementedError


class RNCNVWindow(BaseWindow):
    def __init__(self, min_var_width=5000, min_coverage=0.5,
                 include_var_types=None):
        super().__init__()
        self.include_variant_types = include_var_types or [SVType.DEL,
                                                           SVType.DUP]
        self.min_width = min_var_width
        self.min_coverage = min_coverage
        self.rd_cn = {}
        self.ev = {}
        # A set of OverlappingVariantsDifferBySample
        self.ovds = set()

    def _shall_process(self, variant: VariantRecord) -> bool:
        """
        Returns true if the variation type matches the specified types
        and is wider than a given width threshold.
        """
        return \
            variant.info[SVTYPE] in self.include_variant_types and \
            variant.stop - variant.start >= self.min_width

    def add(self, variant: VariantRecord) -> None:
        super().add(variant)
        rd_cn = []
        ev = []
        for k, v in variant.samples.items():
            rd_cn.append(v[VariantFormatTypes.RD_CN])
            ev.append(v[VariantFormatTypes.EV])
        self.rd_cn[variant.id] = rd_cn
        self.ev[variant.id] = ev

    def add_overlapping_pairs(self, x: VariantRecord, y: VariantRecord) \
            -> None:
        if not self._shall_process(x) or not self._shall_process(y):
            return

        if x.info[SVTYPE] == y.info[SVTYPE]:
            return

        wider, narrower = self._determine_wider(x, y)
        w_samples = set(svu.get_called_samples(wider))
        n_samples = set(svu.get_called_samples(narrower))
        if not w_samples:
            return

        non_common_samples = []
        for sample in w_samples:
            if sample not in n_samples:
                non_common_samples.append(sample)

        coverage_percentage = self.get_coverage(wider, narrower)
        if coverage_percentage >= self.min_coverage:
            for sample in non_common_samples:
                self.ovds.add(OverlappingVariantsDifferBySample(
                    wider, narrower, sample))


class WindowFactory:
    def __init__(self):
        self.windows = {WindowTypes.RNCNV: RNCNVWindow}

    def get_window(self, window_type: WindowTypes) -> BaseWindow:
        return self.windows[window_type]()


class StreamIntersect:
    """
    Finds and processes overlapping intervals on a stream of intervals
    and provides a stream of processed intervals. In other words, it takes
    an iterator on input intervals and returns a generator of processed
    intervals.

    It requires the iterator to yield sorted intervals.

    It yields a set of overlapping intervals encapsulated in a "Window"
    (BaseWindow or derived type). A window contains overlapping intervals
    and implements any business logic on the overlapping pairs.
    """

    def __init__(self, window_type: WindowTypes):
        self.factory = WindowFactory()
        self.window_type = window_type

    def _get_window(self):
        return self.factory.get_window(self.window_type)

    def get_windows(self, intervals: VariantFile) -> Iterator[RNCNVWindow]:
        """
        A generator of "window"s, each encapsulating a
        set of overlapping intervals.
        """
        window = self._get_window()

        # window.add(next(intervals))
        # while interval := next(intervals, None):
        # Since this ^^ syntax is not supported on python versions earlier
        # than 3.8 (at the time of writing this, gatk-sv-base is using
        # python 3.6.5); then the following four lines are used instead.
        for interval in intervals:
            if not window.intervals:
                window.add(interval)
                continue

            if interval.start <= window.stop:
                for x in window.intervals:
                    if x.start <= interval.stop and x.stop >= interval.start:
                        window.add_overlapping_pairs(x, interval)
                window.add(interval)
            else:
                yield window
                window = self._get_window()
                window.add(interval)

        yield window


# TODO: any better descriptive/shorter name?!
class OverlappingVariantsDifferBySample:
    """
    This is container of two overlapping variants that differ
    in their called sample. The wider variant is called on
    the given sample, but the not the narrower variant.
    See the :func:`RNCNVWindow.add_overlapping_pair` for details.
    """
    def __init__(self, wider: VariantRecord, narrower: VariantRecord,
                 sample: str):
        self.wider = wider
        self.narrower = narrower
        self.sample = sample

    def __hash__(self):
        return hash((self.wider.id, self.narrower.id, self.sample))

    def __eq__(self, other):
        if isinstance(other, OverlappingVariantsDifferBySample):
            return self.__hash__() == other.__hash__()
        else:
            return False


class VCFReviser:
    def __init__(self):
        self.samples = {}
        self.multi_cnvs = []

    def get_revised_variants(self, variants: VariantFile) \
            -> Iterator[VariantRecord]:
        """
        This is a generator; it iterates through the input variants,
        revises them as necessary, and yields them.
        """
        si = StreamIntersect(WindowTypes.RNCNV)
        for window in si.get_windows(variants):
            geno_normal_revise_dict = {}

            for x in window.ovds:
                sample_index = self.samples[x.sample]
                new_val = None
                if x.wider.info[SVTYPE] == SVType.DUP and \
                        window.rd_cn[x.narrower.id][sample_index] == 2 and \
                        window.rd_cn[x.wider.id][sample_index] == 3:
                    new_val = 1
                elif x.wider.info[SVTYPE] == SVType.DEL and \
                        window.rd_cn[x.narrower.id][sample_index] == 2 and \
                        window.rd_cn[x.wider.id][sample_index] == 1:
                    new_val = 3

                if new_val:
                    if x.narrower.id not in geno_normal_revise_dict:
                        geno_normal_revise_dict[x.narrower.id] = {}
                    geno_normal_revise_dict[x.narrower.id][x.sample] = new_val

            for variant in window.intervals:
                if variant.id in geno_normal_revise_dict:
                    for sample_id in geno_normal_revise_dict[variant.id]:
                        # original sample
                        o = variant.samples[sample_id]
                        o.update({"GT": (0, 1)})
                        o.update({"GQ": o["RD_GQ"]})

                if variant.stop - variant.start >= 1000:
                    if variant.info[SVTYPE] in [SVType.DEL, SVType.DUP]:
                        is_del = variant.info[SVTYPE] == SVType.DEL
                        for k, v in variant.samples.items():
                            rd_cn = v[VariantFormatTypes.RD_CN]
                            if not rd_cn:
                                continue
                            if (is_del and rd_cn > 3) or \
                                    (not is_del and (rd_cn < 1 or rd_cn > 4)):
                                self.multi_cnvs.append(variant.id)
                                break
                yield variant

    def revise(self, input_vcf: str, output_vcf: str,
               multi_cnvs_filename: str) -> None:
        with VariantFile(input_vcf, "r") as f_in:
            header = f_in.header
            i = -1
            for sample in header.samples:
                i += 1
                self.samples[sample] = i

            with VariantFile(output_vcf, "w", header=header) as f_out:
                for revised_variant in self.get_revised_variants(f_in.fetch()):
                    f_out.write(revised_variant)

        with open(multi_cnvs_filename, "w") as f:
            for x in self.multi_cnvs:
                f.write(x + "\n")


def ensure_file(filename: str) -> str:
    filename = Path(filename)
    if filename.exists():
        os.remove(filename)
    return filename.name


def main(int_vcf_gz):
    reviser = VCFReviser()
    revised_vcf = ensure_file("./normal.revise.vcf")
    revised_vcf_zipped = ensure_file(revised_vcf + ".gz")
    multi_cnvs = ensure_file("./multi.cnvs.txt")

    reviser.revise(int_vcf_gz, revised_vcf, multi_cnvs)

    # Since the input is sorted and the script preserves
    # their order, hence sorting output should not be necessary.
    # sorted_output = revised_vcf + "sorted"
    # cmd = f"cat {revised_vcf} | vcf-sort > {sorted_output}"
    # ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    # ps.communicate()[0]
    check_call(["bgzip", revised_vcf])
    check_call(["bcftools", "index", revised_vcf_zipped])


if __name__ == '__main__':
    main(sys.argv[1])
