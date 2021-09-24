"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

import os
import pybedtools
import pysam
import svtk.utils as svu
import sys
from pathlib import Path
from subprocess import check_call, Popen, PIPE, STDOUT


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


class VCFReviser:
    def __init__(self):
        self.rd_cn = {}
        self.ev = {}
        self.samples = {}

    def _update_rd_cn_ev(self, variant):
        rd_cn = []
        ev = []
        for k, v in variant.samples.items():
            rd_cn.append(v[VariantFormatTypes.RD_CN])
            ev.append(v[VariantFormatTypes.EV])
        self.rd_cn[variant.id] = rd_cn
        self.ev[variant.id] = ev

    def filter_variant_types(self, variants, filter_types=None,
                             min_width=5000):
        for variant in variants:
            sv_type = variant.info[SVTYPE]
            if filter_types and sv_type not in filter_types:
                continue
            if variant.stop - variant.start < min_width:
                continue

            self._update_rd_cn_ev(variant)

            samples = ','.join(svu.get_called_samples(variant))
            if len(samples) == 0:
                samples = BLANK_SAMPLES

            yield variant.contig, variant.start, variant.stop, \
                variant.id, sv_type, samples

    @staticmethod
    def get_wider(f):
        if int(f[2]) - int(f[1]) >= int(f[8]) - int(f[7]):
            return f[0:6], f[6:12]
        else:
            return f[6:12], f[0:6]

    @staticmethod
    def get_coverage(wider, narrower):
        n_start = int(narrower[1])
        n_stop = int(narrower[2])
        w_start = int(wider[1])
        w_stop = int(wider[2])

        coverage = 0
        if w_start <= n_stop and n_start <= w_stop:
            intersection_size = min(n_stop, w_stop) - max(n_start, w_start)
            coverage = intersection_size / (n_stop - n_start)
        return coverage

    def get_geno_normal_revise(self, int_vcf_gz):
        overlap_test_text = {}
        with pysam.VariantFile(int_vcf_gz, "r") as f:
            header = f.header
            i = -1
            for sample in header.samples:
                i += 1
                self.samples[sample] = i

            dels_and_dups = pybedtools.BedTool(
                self.filter_variant_types(
                    f.fetch(), [SVType.DEL, SVType.DUP], 5000)).saveas()

            overlapping_variants = dels_and_dups.intersect(
                dels_and_dups, wa=True, wb=True)

            for interval in overlapping_variants.intervals:
                # IDs are identical, hence it is the overlap
                # of an interval with itself.
                if interval.fields[3] == interval.fields[9]:
                    continue

                # Variant types are identical.
                if interval.fields[4] == interval.fields[10]:
                    continue

                wider, narrower = self.get_wider(interval.fields)
                if wider[5] == BLANK_SAMPLES:
                    continue
                non_common_samples = []
                for x in wider[5].split(","):
                    if x not in narrower[5].split(","):
                        non_common_samples.append(x)

                coverage = self.get_coverage(wider, narrower)
                if coverage >= 0.5:
                    for x in non_common_samples:
                        overlap_test_text[f"{narrower[3]}@{x}"] = \
                            [f"{narrower[3]}@{x}", wider[3], wider[4]]

        geno_normal_revise_dict = {}
        for k, v in overlap_test_text.items():
            var_id = k.split("@")[0]
            sample_index = self.samples[k.split("@")[1]]
            new_val = None
            if v[2] == SVType.DUP and \
                    self.rd_cn[var_id][sample_index] == 2 and \
                    self.rd_cn[v[1]][sample_index] == 3:
                new_val = 1
            elif v[2] == SVType.DEL and \
                    self.rd_cn[var_id][sample_index] == 2 \
                    and self.rd_cn[v[1]][sample_index] == 1:
                new_val = 3

            if new_val:
                if var_id not in geno_normal_revise_dict:
                    geno_normal_revise_dict[var_id] = {}
                geno_normal_revise_dict[var_id][v[0].split("@")[1]] = new_val

        return geno_normal_revise_dict

    def modify_variants(self, int_vcf_gz, revised_vcf, multi_cnvs):
        geno_normal_revise_dict = self.get_geno_normal_revise(int_vcf_gz)

        with pysam.VariantFile(int_vcf_gz, "r") as f_in:
            header = f_in.header

            with pysam.VariantFile(revised_vcf, "w", header=header) as f_out, \
                    open(multi_cnvs, "w") as multi_cnvs_f:
                variants = f_in.fetch()
                for variant in variants:
                    if variant.id in geno_normal_revise_dict:
                        for sample_id in geno_normal_revise_dict[variant.id]:
                            o = variant.samples[sample_id]
                            o.update({"GT": (0, 1)})
                            o.update({"GQ": o["RD_GQ"]})

                    if variant.stop - variant.start >= 1000:
                        if variant.info[SVTYPE] in [SVType.DEL, SVType.DUP]:
                            is_del = variant.info[SVTYPE] == SVType.DEL
                            for k, v in variant.samples.items():
                                rd_cn = v[VariantFormatTypes.RD_CN]
                                if rd_cn is None:
                                    continue
                                if (is_del and rd_cn > 3) or \
                                   (not is_del and (rd_cn < 1 or rd_cn > 4)):
                                    multi_cnvs_f.write(variant.id + "\n")
                                    break

                    f_out.write(variant)


def ensure_file(filename):
    filename = os.path.join(".", filename)
    filename = Path(filename)
    if filename.exists():
        os.remove(filename)
    return filename.name


def main(int_vcf_gz):
    revised_vcf_filename = ensure_file("normal.revise.unsorted.vcf")
    multi_cnvs_filename = ensure_file("multi.cnvs.txt")

    reviser = VCFReviser()
    reviser.modify_variants(
        int_vcf_gz, revised_vcf_filename, multi_cnvs_filename)

    sorted_output = ensure_file("normal.revise.vcf")
    cmd = f"cat {revised_vcf_filename} | vcf-sort > {sorted_output}"
    ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    ps.communicate()[0]
    ensure_file("normal.revise.vcf.gz")
    check_call(["bgzip", sorted_output])
    output_name = sorted_output + ".gz"
    check_call(["bcftools", "index", output_name])


if __name__ == '__main__':
    main(sys.argv[1])
