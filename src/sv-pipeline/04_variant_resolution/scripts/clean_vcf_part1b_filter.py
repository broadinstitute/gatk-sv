"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

import os
import logging
import pysam
import sys
from pathlib import Path
import json
import gzip

SVTYPE = "SVTYPE"
BLANK_SAMPLES = "B"


class SVType:
    DUP = "DUP"
    DEL = "DEL"


class VariantFormatTypes:
    # Predicted copy state
    RD_CN = "RD_CN"
    # Classes of evidence supporting final genotype
    EV = "EV"


def modify_variants(dict_file_gz, vcf, multi_cnvs):
    logging.info('Loading dictionary')
    with gzip.open(dict_file_gz, 'rt') as f:
        geno_normal_revise_dict = json.load(f)

    logging.info('Filtering variants')
    with pysam.VariantFile(vcf, "r") as f_in:
        header = f_in.header
        sys.stdout.write(str(header))
        with open(multi_cnvs, "w") as multi_cnvs_f:
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

                sys.stdout.write(str(variant))


def ensure_file(filename):
    filename = os.path.join(".", filename)
    filename = Path(filename)
    if filename.exists():
        os.remove(filename)
    return filename.name


def main(args):
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    logging.info('Starting script')
    multi_cnvs_filename = ensure_file("multi.cnvs.txt")
    dict_file_gz = args[1]
    vcf_file = args[2]
    modify_variants(dict_file_gz, vcf_file, multi_cnvs_filename)
    logging.info('Done')


if __name__ == '__main__':
    main(sys.argv)
