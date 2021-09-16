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
import os
import sys
import tempfile
from subprocess import check_call, Popen, PIPE, STDOUT


VCF_DELIMITER = "\t"
BED_DELIMITER = "\t"


def get_columns_headers(filename):
    with gzip.open(filename, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            else:
                return line


def get_filtered_vcf_in_bed(filename):
    # TODO: This script uses intermediate files for multiple manipulations,
    #  can all manipulations implemented in a single pass?
    int_bed = "./int.bed"
    tmp_filtered_vcf = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with gzip.open(filename, "rt") as f, tmp_filtered_vcf as tf:
        for line in f:
            sline = line.split(VCF_DELIMITER)
            if "#" in sline[0] or "DEL" in sline[4] or "DUP" in sline[4]:
                tf.write(line)

    tmp_bed = tempfile.NamedTemporaryFile(mode="w", delete=False)
    check_call(["svtk", "vcf2bed", tmp_filtered_vcf.name, tmp_bed.name])

    with open(tmp_bed.name, "r") as f, open(int_bed, "w") as m:
        for line in f:
            sline = line.strip().split(BED_DELIMITER)
            if len(sline) < 6:
                # TODO: seems like we can ignore writing this line
                #  as it is skipped downstream.
                m.write(f"{line.strip()}{BED_DELIMITER}blanksample\n")
            else:
                m.write(line)

    check_call(["gzip", "--force", int_bed])
    os.remove(tmp_filtered_vcf.name)
    os.remove(tmp_bed.name)
    return int_bed + ".gz"


def get_normal_overlap(filename):
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
    with gzip.open(filename, "rt") as f, tmp as t:
        for line in f:
            # TODO: maybe skip writing header line in the first place.
            if line.startswith("#"):
                continue
            sline = line.split(BED_DELIMITER)
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
        tmp_cols[i][2] = tmp_cols[i][4]

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
            for i in range(3, len(sl)):
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


def main(int_vcf_gz):
    headers = get_columns_headers(int_vcf_gz)
    headers = headers.strip().split("\t")
    col_txt = "col.txt"
    cols = {}
    with open(col_txt, "w") as f:
        for x in range(len(headers)):
            f.write(f"{x}\t{headers[x]}\n")
            cols[x] = [headers[x]]
    int_bed_gz = get_filtered_vcf_in_bed(int_vcf_gz)
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
    print(multi_cnvs_txt)


if __name__ == '__main__':
    main(sys.argv[1])
