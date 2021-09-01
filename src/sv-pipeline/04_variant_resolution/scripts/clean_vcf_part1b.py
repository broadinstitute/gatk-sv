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

import gzip
import os
import sys
import tempfile
from subprocess import check_call


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
    os.remove(tmp)
    os.remove(tmp2)
    os.remove(tmp3)
    return normaloverlap.name


def main(input_vcf):
    headers = get_columns_headers(input_vcf)
    headers = headers.replace("\t", "\n")
    with open("col.txt", "w") as f:
        f.writelines(headers)
    int_bed_gz = get_filtered_vcf_in_bed(input_vcf)
    normaloverlap_txt = get_normal_overlap(int_bed_gz)


if __name__ == '__main__':
    main(sys.argv[1])
