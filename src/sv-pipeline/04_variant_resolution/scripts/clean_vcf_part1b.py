"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

# TODO: most operations compress intermediate files and processed
#  the compressed files; if there is no particular reason for compressing
#  the intermediate files, replace them with their uncompressed versions.

# TODO: many intermediate files are used, some algorithm changes are
#  required to aggregate operations, and maybe using stdio/stdout for
#  streaming data would be better alternatives.

import copy
import gzip
import json
import os
import pybedtools
import pysam
import shutil
import sys
import tempfile
import time
import svtk.utils as svu
from subprocess import check_call, Popen, PIPE, STDOUT
from typing import Iterable, Iterator, Tuple


DELIMITER = "\t"
VCF_DELIMITER = "\t"
BED_DELIMITER = "\t"


class Keys:
    SVTYPE = "SVTYPE"
    DUP = "DUP"
    DEL = "DEL"


BLANK_SAMPLES = "Blank_Samples"


class Filter:
    def __init__(self):
        self.rd_cn = {}
        self.rd_cn_idx = None
        self.ev = {}
        self.ev_idx = None
        self.samples = {}

    def _update_rd_cn_ev(self, variant):
        rd_cn = []
        ev = []
        for k, v in variant.samples.items():
            rd_cn.append(v.values()[self.rd_cn_idx])
            ev.append(v.values()[self.ev_idx])
        self.rd_cn[variant.id] = rd_cn
        self.ev[variant.id] = ev

    def filter_variant_types(self, records: Iterable[pysam.VariantRecord], filter_types=None, min_interval_width=1) -> Iterator[Tuple]:
        # TODO: rename records to variants
        for record in records:
            sv_type = record.info[Keys.SVTYPE]
            if filter_types and sv_type not in filter_types:
                continue
            if record.stop - record.start < min_interval_width:
                continue

            self._update_rd_cn_ev(record)

            samples = ','.join(svu.get_called_samples(record))
            if len(samples) == 0:
                samples = BLANK_SAMPLES
            yield record.contig, record.start, record.stop, record.id, sv_type, samples

    def get_overlapping_variants(self, int_vcf_gz):
        overlapping_variants_ids = set()
        overlap_test_text = {}
        with pysam.VariantFile(int_vcf_gz, "r") as f:
            header = f.header
            # print(header)
            # print("\n\n")
            # print(dir(header))
            # print("\n\n")
            # print(len(header.samples))
            # print(dir(header.samples))
            i = -1
            for sample in header.samples:
                i += 1
                self.samples[sample] = i
                # SAMPLES.append(sample)
            # print("\n\n----------")
            # print(dir(header.formats))
            formats = header.formats.keys()

            # TODO: this should be re-checked with every line
            #  since the lines can differ on this. This part
            #  should be moved to where the variants are traversed.
            self.rd_cn_idx = formats.index("RD_CN")
            self.ev_idx = formats.index("EV")
            # print(rd_cn_idx)
            # print(header.formats.keys())
            # print(header.formats.values())
            # exit()
            dels_and_dups = pybedtools.BedTool(self.filter_variant_types(f.fetch(), [Keys.DEL, Keys.DUP], 5000)).saveas()
            overlapping_variants = dels_and_dups.intersect(dels_and_dups, wa=True, wb=True)
            for interval in overlapping_variants.intervals:
                # IDs are identical, hence it is the overlap
                # of an interval with itself.
                var_a_id = interval.fields[3]
                var_b_id = interval.fields[9]
                if var_a_id == var_b_id:
                    continue

                # Variant types are identical.
                if interval.fields[4] == interval.fields[10]:
                    continue

                f = interval.fields
                non_common_samples = []
                if int(f[2]) - int(f[1]) >= int(f[8]) - int(f[7]):
                    if f[5] == BLANK_SAMPLES:
                        continue
                    overlapping_variants_ids.add(var_a_id)
                    overlapping_variants_ids.add(var_b_id)
                    wider_var_id = f[3]
                    wider_var_interval = [int(f[1]), int(f[2])]
                    wider_var_type = f[4]
                    narrower_var_id = f[9]
                    narrower_var_interval = [int(f[7]), int(f[8])]
                    for x in interval.fields[5].split(","):
                        if x not in interval.fields[11].split(","):
                            non_common_samples.append(x)
                else:
                    if f[11] == BLANK_SAMPLES:
                        continue
                    overlapping_variants_ids.add(var_a_id)
                    overlapping_variants_ids.add(var_b_id)
                    narrower_var_id = f[3]
                    narrower_var_interval = [int(f[1]), int(f[2])]
                    wider_var_id = f[9]
                    wider_var_interval = [int(f[7]), int(f[8])]
                    wider_var_type = f[10]
                    for x in interval.fields[11].split(","):
                        if x not in interval.fields[5].split(","):
                            non_common_samples.append(x)

                ss = narrower_var_interval[0]
                se = narrower_var_interval[1]
                ls = wider_var_interval[0]
                le = wider_var_interval[1]

                coverage_percentage = 0
                if ls <= se and ss <= le:
                    intersection_size = min(se, le) - max(ss, ls)
                    coverage_percentage = intersection_size / (se-ss)

                if coverage_percentage >= 0.5:
                    for x in non_common_samples:
                        overlap_test_text[f"{narrower_var_id}@{x}"] = [f"{narrower_var_id}@{x}", wider_var_id, wider_var_type]

        geno_normal_revise_txt = []
        geno_normal_revise_dict = {}
        for k, v in overlap_test_text.items():
            var_id = k.split("@")[0]
            sample_index = self.samples[k.split("@")[1]]
            if v[2] == "DUP" and self.rd_cn[var_id][sample_index] == 2 and self.rd_cn[v[1]][sample_index] == 3:
                geno_normal_revise_txt.append(v[0].replace("@", "\t") + "\t" + "1" + "\n")
                if var_id not in geno_normal_revise_dict:
                    geno_normal_revise_dict[var_id] = {}
                geno_normal_revise_dict[var_id][v[0].split("@")[1]] = 1
            elif v[2] == "DEL" and self.rd_cn[var_id][sample_index] == 2 and self.rd_cn[v[1]][sample_index] == 1:
                geno_normal_revise_txt.append(v[0].replace("@", "\t") + "\t" + "3" + "\n")
                if var_id not in geno_normal_revise_dict:
                    geno_normal_revise_dict[var_id] = {}
                geno_normal_revise_dict[var_id][v[0].split("@")[1]] = 3

        return overlapping_variants_ids, overlap_test_text, geno_normal_revise_txt, geno_normal_revise_dict, header

    def modify_variants(self, int_vcf_gz, geno_normal_revise_dict, header):
        output_vcf = get_tmp_file(get_working_dir())
        multi_cnvs_txt = get_tmp_file(get_working_dir())
        with pysam.VariantFile(int_vcf_gz, "r") as f_in, pysam.VariantFile(output_vcf.name, "w", header=header) as f_out, multi_cnvs_txt as multi_cnvs_f:
            variants = f_in.fetch()
            for variant in variants:
                # TODO: minimize dict access here.
                if variant.id in geno_normal_revise_dict:
                    for sample_id in geno_normal_revise_dict[variant.id]:
                        original_sample = variant.samples[sample_id]
                        original_sample.update({"GT": (0, 1)})
                        original_sample.update({"GQ": original_sample["RD_GQ"]})
                        # TODO: double-check above with shell script & Harison.

                if variant.stop - variant.start >= 1000:
                    if variant.info[Keys.SVTYPE] == Keys.DEL or variant.info[Keys.SVTYPE] == Keys.DUP:
                        is_del = variant.info[Keys.SVTYPE] == Keys.DEL
                        for k, v in variant.samples.items():
                            rd_cn_val = v.values()[self.rd_cn_idx]
                            if (is_del and rd_cn_val > 3) or (not is_del and (rd_cn_val < 1 or rd_cn_val > 4)):
                                multi_cnvs_f.write(variant.id + "\n")
                                break

                f_out.write(variant)

        # TODO: redo the following in a better way.
        sorted_output = output_vcf.name + "sorted"
        cmd = f"cat {output_vcf.name} | vcf-sort > {sorted_output}"
        ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        ps.communicate()[0]
        check_call(["bgzip", sorted_output])
        output_name = sorted_output + ".gz"
        check_call(["bcftools", "index", output_name])
        return output_name, multi_cnvs_txt.name


def get_columns_headers(filename):
    with gzip.open(filename, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            else:
                return line


def get_working_dir(base_dir="."):
    while True:
        working_dir = os.path.join(base_dir, "tmp_" + str(int(time.time())))
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
            return working_dir


def get_tmp_file(working_dir):
    return tempfile.NamedTemporaryFile(dir=working_dir, mode="w", delete=False)


def filter_variant_types(records: Iterable[pysam.VariantRecord], filter_types=None, min_interval_width=1) -> Iterator[Tuple]:
    for record in records:
        sv_type = record.info[Keys.SVTYPE]
        if filter_types and sv_type not in filter_types:
            continue
        if record.stop - record.start < min_interval_width:
            continue

        rd_cn = []
        for k, v in record.samples.items():
            print("\n-----------")
            print(f"{k}\t{v}")
            print(v.values())
            print(rd_cn_idx)
            print(v.values()[rd_cn_idx])
            print(dir(v))
            exit()
        exit()

        samples = ','.join(svu.get_called_samples(record))
        if len(samples) == 0:
            samples = BLANK_SAMPLES


        yield record.contig, record.start, record.stop, record.id, sv_type, samples

RD_CN = {}
rd_cn_idx = None
SAMPLES = []

def get_overlapping_variants(int_vcf_gz):
    overlapping_variants_ids = set()
    overlap_test_text = {}
    with pysam.VariantFile(int_vcf_gz, "r") as f:
        header = f.header
        # print(header)
        # print("\n\n")
        # print(dir(header))
        # print("\n\n")
        # print(len(header.samples))
        # print(dir(header.samples))
        for sample in header.samples:
            SAMPLES.append(sample)
        # print("\n\n----------")
        # print(dir(header.formats))
        formats = header.formats.keys()

        rd_cn_idx = formats.index("RD_CN")
        # print(rd_cn_idx)
        # print(header.formats.keys())
        # print(header.formats.values())
        # exit()
        dels_and_dups = pybedtools.BedTool(filter_variant_types(f.fetch(), [Keys.DEL, Keys.DUP], 5000)).saveas()
        overlapping_variants = dels_and_dups.intersect(dels_and_dups, wa=True, wb=True)
        for interval in overlapping_variants.intervals:
            # IDs are identical, hence it is the overlap
            # of an interval with itself.
            var_a_id = interval.fields[3]
            var_b_id = interval.fields[9]
            if var_a_id == var_b_id:
                continue

            # Variant types are identical.
            if interval.fields[4] == interval.fields[10]:
                continue

            f = interval.fields
            non_common_samples = []
            if int(f[2]) - int(f[1]) >= int(f[8]) - int(f[7]):
                if f[5] == BLANK_SAMPLES:
                    continue
                overlapping_variants_ids.add(var_a_id)
                overlapping_variants_ids.add(var_b_id)
                wider_var_id = f[3]
                wider_var_interval = [int(f[1]), int(f[2])]
                wider_var_type = f[4]
                narrower_var_id = f[9]
                narrower_var_interval = [int(f[7]), int(f[8])]
                for x in interval.fields[5].split(","):
                    if x not in interval.fields[11].split(","):
                        non_common_samples.append(x)
            else:
                if f[11] == BLANK_SAMPLES:
                    continue
                overlapping_variants_ids.add(var_a_id)
                overlapping_variants_ids.add(var_b_id)
                narrower_var_id = f[3]
                narrower_var_interval = [int(f[1]), int(f[2])]
                wider_var_id = f[9]
                wider_var_interval = [int(f[7]), int(f[8])]
                wider_var_type = f[10]
                for x in interval.fields[11].split(","):
                    if x not in interval.fields[5].split(","):
                        non_common_samples.append(x)

            ss = narrower_var_interval[0]
            se = narrower_var_interval[1]
            ls = wider_var_interval[0]
            le = wider_var_interval[1]

            coverage_percentage = 0
            if ls <= se and ss <= le:
                intersection_size = min(se, le) - max(ss, ls)
                coverage_percentage = intersection_size / (se-ss)

            if coverage_percentage >= 0.5:
                for x in non_common_samples:
                    overlap_test_text[f"{narrower_var_id}@{x}"] = [f"{narrower_var_id}@{x}", wider_var_id, wider_var_type]

    return overlapping_variants_ids, overlap_test_text


def get_depth_based_copy_number_variant_2(vcf_gz, overlapping_variants_id):
    pass


def get_filtered_vcf_in_bed(working_dir, int_vcf_gz):
    # TODO: This script uses intermediate files for multiple manipulations,
    #  can all manipulations implemented in a single pass?
    int_bed = get_tmp_file(working_dir)
    filtered_vcf = get_tmp_file(working_dir)
    with gzip.open(int_vcf_gz, "rt") as f, filtered_vcf as tf:
        for line in f:
            columns = line.split(VCF_DELIMITER)
            if "#" in columns[0] or "DEL" in columns[4] or "DUP" in columns[4]:
                tf.write(line)

    tmp_bed = get_tmp_file(working_dir)
    check_call(["svtk", "vcf2bed", filtered_vcf.name, tmp_bed.name])

    with open(tmp_bed.name, "r") as f, open(int_bed.name, "w") as m:
        for line in f:
            columns = line.strip().split(DELIMITER)
            if len(columns) < 6:
                # TODO: seems like we can ignore writing this line
                #  as it is skipped downstream.
                m.write(f"{line.strip()}{DELIMITER}blanksample\n")
            else:
                m.write(line)

    check_call(["gzip", "--force", int_bed.name])
    os.remove(filtered_vcf.name)
    os.remove(tmp_bed.name)
    return int_bed.name + ".gz"


def get_normal_overlap(int_bed_gz):
    """
    list of potential overlaps with a normal copy state variant
    (>5kb variants require depth but nested events could be missed;
    i.e a duplication with a nest deletion will have a normal copy
    state for the deletion).

    Flip bed intersect so largest is CNV is always first.
    """
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp2 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp3 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    normaloverlap = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(int_bed_gz, "rt") as f, tmp as t:
        for line in f:
            # TODO: maybe skip writing header line in the first place.
            if line.startswith("#"):
                continue
            sline = line.split(DELIMITER)
            if int(sline[2]) - int(sline[1]) >= 5000:
                t.write(line)

    check_call(["bedtools", "intersect", "-wa", "-wb", "-a", tmp.name, "-b", tmp.name], stdout=tmp2)

    with open(tmp2.name, "r") as f, tmp3 as t:
        for line in f:
            line = line.strip()
            s = line.strip().split(BED_DELIMITER)
            if s[3] != s[9] and int(s[2]) - int(s[1]) >= int(s[8]) - int(s[7]) and s[4] != s[10]:
                if s[5] == "blanksample":
                    continue
                t.write(line + "\n")
            elif s[3] != s[9] and s[4] != s[10]:
                x = [s[6], s[7], s[8], s[9], s[10], s[11], s[0], s[1], s[2], s[3], s[4], s[5]]
                if x[5] == "blanksample":
                    continue
                t.write(BED_DELIMITER.join(x) + "\n")

    check_call(["sort", "-u", tmp3.name], stdout=normaloverlap)
    os.remove(tmp.name)
    os.remove(tmp2.name)
    os.remove(tmp3.name)
    return normaloverlap.name


def get_depth_based_copy_number_variant(vcf_gz, normaloverlap):
    # A set, implemented as a hashtable, so lookup
    # should be O(log n) or O(1) depending on how it is implemented.
    ids = set()
    with open(normaloverlap) as f:
        for line in f:
            sline = line.strip().split(BED_DELIMITER)
            ids.add(sline[3])
            ids.add(sline[9])

    # keep only those lines in vcf whose ID is in the ids set.
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(vcf_gz, "rt") as f, tmp as t:
        for line in f:
            if line.startswith("#"):
                t.write(line)
            else:
                sline = line.split(VCF_DELIMITER)
                if sline[2] in ids:
                    t.write(line)

    tmp2 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with open(tmp.name) as t1, tmp2 as t2:
        for line in t1:
            sline = line.strip().split(VCF_DELIMITER)
            if "#" not in sline[0]:
                sline[0] = sline[2]
            if "#" in sline[0] or sline[4] == "<DEL>" or sline[4] == "<DUP>":
                t2.write(VCF_DELIMITER.join(sline) + "\n")

    tmp3 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp3.close()
    format_id = "RD_CN"
    check_call(["vcftools", "--vcf", tmp2.name, "--extract-FORMAT-info", format_id, "--out", tmp3.name])
    vcftools_output = f"{tmp3.name}.{format_id}.FORMAT"

    c = 0
    header = []
    tmp4 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with open(vcftools_output, "r") as f, tmp4 as t:
        for line in f:
            sline = line.strip().split(VCF_DELIMITER)
            c += 1
            if c == 1:
                for i in range(2, len(sline)):
                    header.append(sline[i])
            else:
                for i in range(2, len(sline)):
                    t.write(f"{sline[0]}@{header[i-2]}\t{sline[i]}\n")

    tmp5 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    check_call(["sort", "-k1,1", tmp4.name], stdout=tmp5)
    check_call(["gzip", tmp5.name])
    return tmp5.name + ".gz"


def get_evidence_supporting_each_normal_overlapping_variant(vcf_gz, normaloverlap):
    # Probably the only difference between this method and the previous one
    # is the VCF filter they use. another diff is at creating tmp2,
    # double-check the above method, it seems that method is doing an
    # unnecessary additional step.
    ids = set()
    with open(normaloverlap) as f:
        for line in f:
            sline = line.strip().split(BED_DELIMITER)
            ids.add(sline[3])
            ids.add(sline[9])

    # keep only those lines in vcf whose ID is in the ids set.
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(vcf_gz, "rt") as f, tmp as t:
        for line in f:
            if line.startswith("#"):
                t.write(line)
            else:
                sline = line.split(VCF_DELIMITER)
                if sline[2] in ids:
                    t.write(line)

    tmp2 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with open(tmp.name) as t1, tmp2 as t2:
        for line in t1:
            sline = line.strip().split(VCF_DELIMITER)
            if "#" not in sline[0]:
                sline[0] = sline[2]
            t2.write(VCF_DELIMITER.join(sline) + "\n")

    tmp3 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tmp3.close()
    format_id = "EV"
    check_call(["vcftools", "--vcf", tmp2.name, "--extract-FORMAT-info", format_id, "--out", tmp3.name])
    vcftools_output = f"{tmp3.name}.{format_id}.FORMAT"

    c = 0
    header = []
    tmp4 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with open(vcftools_output, "r") as f, tmp4 as t:
        for line in f:
            sline = line.strip().split(VCF_DELIMITER)
            c += 1
            if c == 1:
                for i in range(2, len(sline)):
                    header.append(sline[i])
            else:
                for i in range(2, len(sline)):
                    t.write(f"{sline[0]}@{header[i-2]}\t{sline[i]}\n")

    tmp5 = tempfile.NamedTemporaryFile(mode="w", delete=False)
    check_call(["sort", "-k1,1", tmp4.name], stdout=tmp5)
    check_call(["gzip", tmp5.name])
    return tmp5.name + ".gz"


def check_if_nested_is_incorrectly_classified_as_normal(normaloverlap):
    overlap_test_txt = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with open(normaloverlap, "r") as f, overlap_test_txt as o:
        for line in f:
            line = line.strip()
            line = line.replace(" ", "\t")
            sline = line.split("\t")
            # lfile.write("\t".join(sline[0:6]) + "\n") # maybe 5
            # sfile.write("\t".join(sline[6:12]) + "\n") # maybe 11

            ss = int(sline[7])
            se = int(sline[8])
            ls = int(sline[1])
            le = int(sline[2])

            coverage_percentage = 0
            if ls <= se and ss <= le:
                intersection_size = min(se, le) - max(ss, ls)
                coverage_percentage = intersection_size / (se-ss)

            alist = []
            if coverage_percentage >= 0.5:
                fgrep_part_a = sline[11].split(",")
                smallid = sline[9]
                # wrong variable name
                large_ids = sline[5].split(",")
                for x in large_ids:
                    alist.append(f"{smallid}@{x}\t{sline[3]}@{x}\t{sline[4]}")

                not_match_count = 0
                # not sure if this check is what it was intended.
                for x in alist:
                    for y in fgrep_part_a:
                        if y in x:
                            break
                    else:
                        o.write(x + "\n")
    return overlap_test_txt.name


def get_variants_to_be_revised_from_normal_copy_state_into_cnv(
        overlap_test_txt, rd_cn_normalcheck_format_gz,
        ev_normalcheck_format_gz):
    data = {}
    with open(overlap_test_txt) as f:
        for line in f:
            sline = line.strip().split("\t")
            data[sline[0]] = [sline[0], sline[1], sline[2]]

    rd_cn_dict = {}
    with gzip.open(rd_cn_normalcheck_format_gz, "rt") as f:
        for line in f:
            sline = line.strip().split("\t")
            if sline[0] in data:
                data[sline[0]].append(int(sline[1]))
            rd_cn_dict[sline[0]] = int(sline[1])

    with gzip.open(ev_normalcheck_format_gz, "rt") as f:
        for line in f:
            sline = line.strip().split("\t")
            if sline[0] in data:
                data[sline[0]].append(sline[1])

    joined = []
    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with output as t:
        for k, v in data.items():
            if v[0] in rd_cn_dict:
                joined.append([v[0], v[1], v[2], v[3], v[4], rd_cn_dict[v[1]]])
                if v[2] == "DUP" and v[3] == 2 and rd_cn_dict[v[1]] == 3:
                    aaa = v[0].replace("@", "\t") + "\t" + "1" + "\n"
                    t.write(aaa)
                elif v[2] == "DEL" and v[3] == 2 and rd_cn_dict[v[1]] == 1:
                    t.write(v[0].replace("@", "\t") + "\t" + "3" + "\n")

    return output.name


def get_subset_vcf(geno_normal_revise_txt, int_vcf_gz):
    ids = set()
    with open(geno_normal_revise_txt, "r") as f:
        for line in f:
            ids.add(line.strip().split("\t")[0])

    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(int_vcf_gz, "rt") as f, output as o:
        for line in f:
            if line.startswith("#"):
                continue
            id = line.split("\t")[2]
            if id in ids:
                o.write(line)

    return output.name


def pull_out_and_revise_vcf_line_that_needs_to_be_edited(geno_normal_revise_txt, subset_vcf, col_txt, cols):
    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    variants_dict = {}
    col_headers = {}
    for x in range(len(cols)):
        col_headers[cols[x][0]] = x

    with open(subset_vcf) as f:
        for line in f:
            sline = line.strip().split("\t")
            variants_dict[sline[2]] = sline

    ids = []
    with open(geno_normal_revise_txt) as f:
        for line in f:
            sline = line.strip().split("\t")
            ids.append(sline)

    for _id in ids:  # variant, v in ids.items():
        variant = _id[0]
        # copy the dict to its state is reset at every iteration.
        tmp_cols = copy.deepcopy(cols)

        sline_txt = variants_dict[variant]
        for x in range(len(sline_txt)):
            tmp_cols[x].extend(sline_txt[x].split(":"))

        i = col_headers[_id[1]]
        tmp_cols[i][1] = "0/1"
        tmp_cols[i][2] = tmp_cols[i][4]  # ??? should not be the value in the third column of _id? also are you sure to modify the 2nd item and not the one related to RD_CN or EV?

        for i in range(len(sline_txt)):
            if i < 9:
                continue
            sline_txt[i] = ":".join(tmp_cols[i][1:])

    with output as o:
        sorted_keys = list(variants_dict.keys())
        sorted_keys.sort()
        for key in sorted_keys:
            o.write("\t".join(variants_dict[key]) + "\n")
    return output.name


def modify_vcf(int_vcf_gz, normal_revise_vcf_lines_txt):
    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(int_vcf_gz, "rt") as in_file, open(output.name, "wt") as out_file:
        for l in in_file:
            if not l.startswith("#"):
                sl = l.strip().split("\t")
                variant_id = sl[2]
                with open(normal_revise_vcf_lines_txt) as m_file:
                    for l2 in m_file:
                        sl2 = l2.strip().split("\t")
                        if sl2[2] == variant_id:
                            sl = sl2
                            break
                l = "\t".join(sl) + "\n"
            out_file.write(l)

    sorted_output = output.name + "sorted"
    cmd = f"cat {output.name} | vcf-sort > {sorted_output}"
    ps = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
    ps.communicate()[0]
    check_call(["bgzip", sorted_output])
    output_name = sorted_output + ".gz"
    check_call(["bcftools", "index", output_name])
    return output_name


def get_copystate_per_variant(normal_revise_vcf_gz):
    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(normal_revise_vcf_gz, "rt") as f, tmp as t:
        for l in f:
            if l.startswith("#"):
                t.write(l)
                continue
            sline = l.strip().split("\t")
            sline[0] = sline[2]
            t.write("\t".join(sline) + "\n")

    format_id = "RD_CN"
    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    check_call(["vcftools", "--vcf", tmp.name, "--extract-FORMAT-info", format_id, "--out", output.name])
    output_name = output.name + "." + format_id + ".FORMAT"
    check_call(["gzip", output_name])
    return output_name + ".gz"


def get_copystate_per_variant_2(copystate_rd_cn_format_gz):
    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    output_data = set()
    with gzip.open(copystate_rd_cn_format_gz, "rt") as f, output as o:
        c = 0
        for l in f:
            c += 1
            if c <= 1:
                continue
            sl = l.strip().split("\t")
            for i in range(2, len(sl)):
                x = f"{sl[0]}\t{sl[i]}"
                if x not in output_data:
                    output_data.add(x)

        for x in output_data:
            o.write(x + "\n")
    check_call(["gzip", output.name])
    return output.name + ".gz"


def find_multi_allelic_for_del(copystate_per_variant_txt_gz, int_bed_gz):
    variants = set()
    with gzip.open(copystate_per_variant_txt_gz, "rt") as f:
        for l in f:
            sl = l.strip().split("\t")
            if sl[1] != "." and int(sl[1]) > 3:
                variants.add(sl[0])

    output = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(int_bed_gz, "rt") as f, output as o:
        for l in f:
            sl = l.strip().split("\t")
            if sl[3] in variants:
                if sl[4] == "DEL" and int(sl[2]) - int(sl[1]) >= 1000:
                    o.write(sl[3] + "\n")

    return output.name


def find_multi_allelic_for_dup(copystate_per_variant_txt_gz, int_bed_gz, multi_cnvs_txt):
    variants = set()
    with gzip.open(copystate_per_variant_txt_gz, "rt") as f:
        for l in f:
            sl = l.strip().split("\t")
            if sl[1] != "." and (int(sl[1]) < 1 or int(sl[1]) > 4):
                variants.add(sl[0])

    with gzip.open(int_bed_gz, "rt") as f, open(multi_cnvs_txt, "a") as o:
        for l in f:
            sl = l.strip().split("\t")
            if sl[3] in variants:
                if sl[4] == "DUP" and int(sl[2]) - int(sl[1]) >= 1000:
                    o.write(sl[3] + "\n")

    return multi_cnvs_txt


def main(int_vcf_gz):
    working_dir = get_working_dir()

    f = Filter()
    overlapping_variants_ids, overlap_test_text, geno_normal_revise_txt, geno_normal_revise_dict, header = f.get_overlapping_variants(int_vcf_gz)
    new_normal_revise_vcf_gz, multi_cnvs_txt = f.modify_variants(int_vcf_gz, geno_normal_revise_dict, header)
    shutil.move(new_normal_revise_vcf_gz, "./normal.revise.vcf.gz")
    shutil.move(new_normal_revise_vcf_gz + ".csi", "./normal.revise.vcf.gz.csi")
    shutil.move(multi_cnvs_txt, "./multi.cnvs.txt")
    exit()



    headers = get_columns_headers(int_vcf_gz)
    headers = headers.strip().split("\t")
    col_txt = "col.txt"
    cols = {}
    with open(col_txt, "w") as f:
        for x in range(len(headers)):
            f.write(f"{x}\t{headers[x]}\n")
            cols[x] = [headers[x]]
    int_bed_gz = get_filtered_vcf_in_bed(working_dir, int_vcf_gz)
    normaloverlap_txt = get_normal_overlap(int_bed_gz)
    rd_cn_normalcheck_FORMAT_gz = get_depth_based_copy_number_variant(int_vcf_gz, normaloverlap_txt)
    ev_normalcheck_format_gz = get_evidence_supporting_each_normal_overlapping_variant(int_vcf_gz, normaloverlap_txt)
    overlap_test_txt = check_if_nested_is_incorrectly_classified_as_normal(normaloverlap_txt)
    geno_normal_revise_txt = get_variants_to_be_revised_from_normal_copy_state_into_cnv(overlap_test_txt, rd_cn_normalcheck_FORMAT_gz, ev_normalcheck_format_gz)
    subset_vcf = get_subset_vcf(geno_normal_revise_txt, int_vcf_gz)
    normal_revise_vcf_lines_txt = pull_out_and_revise_vcf_line_that_needs_to_be_edited(geno_normal_revise_txt, subset_vcf, col_txt, cols)
    normal_revise_vcf_gz = modify_vcf(int_vcf_gz, normal_revise_vcf_lines_txt)
    copystate_rd_cn_format_gz = get_copystate_per_variant(normal_revise_vcf_gz)
    copystate_per_variant_txt_gz = get_copystate_per_variant_2(copystate_rd_cn_format_gz)
    multi_cnvs_txt = find_multi_allelic_for_del(copystate_per_variant_txt_gz, int_bed_gz)
    multi_cnvs_txt = find_multi_allelic_for_dup(copystate_per_variant_txt_gz, int_bed_gz, multi_cnvs_txt)

    shutil.move(int_bed_gz, "/Users/vahid/Desktop/mine/int_bed.gz")
    shutil.move(normaloverlap_txt, "/Users/vahid/Desktop/mine/normaloverlap")
    shutil.move(rd_cn_normalcheck_FORMAT_gz, "/Users/vahid/Desktop/mine/rd_cn_normalcheck_FORMAT.gz")
    shutil.move(ev_normalcheck_format_gz, "/Users/vahid/Desktop/mine/ev_normalcheck_format.gz")
    shutil.move(overlap_test_txt, "/Users/vahid/Desktop/mine/overlap_test_txt")
    shutil.move(geno_normal_revise_txt, "/Users/vahid/Desktop/mine/geno_normal_revise_txt")
    shutil.move(subset_vcf, "/Users/vahid/Desktop/mine/subset_vcf")
    shutil.move(normal_revise_vcf_lines_txt, "/Users/vahid/Desktop/mine/normal_revise_vcf_lines_txt")
    shutil.copy(normal_revise_vcf_gz, "/Users/vahid/Desktop/mine/normal_revise_vcf.gz")
    shutil.move(copystate_rd_cn_format_gz, "/Users/vahid/Desktop/mine/copystate_rd_cn_format.gz")
    shutil.move(copystate_per_variant_txt_gz, "/Users/vahid/Desktop/mine/copystate_per_variant_txt.gz")
    shutil.copy(multi_cnvs_txt, "/Users/vahid/Desktop/mine/multi_cnvs_txt")

    assert os.path.getsize(multi_cnvs_txt) == 5854
    assert os.path.getsize(normal_revise_vcf_gz) == 5814655

    return multi_cnvs_txt, normal_revise_vcf_gz, normal_revise_vcf_gz + "csi"


if __name__ == '__main__':
    main(sys.argv[1])
