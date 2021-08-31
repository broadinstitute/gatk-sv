"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

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
                m.write(f"{line.strip()}{BED_DELIMITER}blanksample\n")
            else:
                m.write(line)

    check_call(["gzip", "--force", int_bed])
    os.remove(tmp_filtered_vcf.name)
    os.remove(tmp_bed.name)
    return int_bed + ".gz"


def main(input_vcf):
    headers = get_columns_headers(input_vcf)
    headers = headers.replace("\t", "\n")
    with open("col.txt", "w") as f:
        f.writelines(headers)
    int_bed_gz = get_filtered_vcf_in_bed(input_vcf)


if __name__ == '__main__':
    main(sys.argv[1])
