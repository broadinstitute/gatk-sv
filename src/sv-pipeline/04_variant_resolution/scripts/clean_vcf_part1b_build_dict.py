"""
Remove CNVs that are improperly genotyped by depth because they are nested
within a real CNV
"""

import logging
import pybedtools
import pysam
import sys
import json

from collections import defaultdict

SVTYPE = "SVTYPE"
BLANK_SAMPLES = "blanksample"


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
        self.sample_indices_dict = {}
        self.sample_list = []

    def _update_rd_cn(self, variant, sample_indices):
        self.rd_cn[variant.id] = {s: variant.samples[s][VariantFormatTypes.RD_CN] for s in sample_indices}

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

    def get_geno_normal_revise(self, vcf_file, bed_file):
        overlap_test_text = defaultdict(dict)
        with pysam.VariantFile(vcf_file, "r") as f:
            header = f.header
            i = -1
            for sample in header.samples:
                i += 1
                self.sample_indices_dict[sample] = i
                self.sample_list.append(sample)

            logging.info("Filtering intersect results")
            bed = pybedtools.BedTool(bed_file)
            for interval in bed.intervals:
                wider, narrower = self.get_wider(interval.fields)
                if wider[5] == BLANK_SAMPLES:
                    continue

                coverage = self.get_coverage(wider, narrower)
                if coverage >= 0.5:
                    wider_samples = set(wider[5].split(","))
                    narrower_samples = set(narrower[5].split(","))
                    non_common_samples = [self.sample_indices_dict[s] for s in wider_samples - narrower_samples]
                    for x in non_common_samples:
                        vid = narrower[3]
                        overlap_test_text[vid][x] = (wider[3], wider[4])

            # Determine for which vid/sample pairs we need RD_CN
            # Substantially reduces memory
            logging.info('Getting revised variant IDs')
            revise_vids = defaultdict(set)
            for var_id, samples_dict in overlap_test_text.items():
                for sample_index, v in samples_dict.items():
                    if v[1] == SVType.DUP or v[1] == SVType.DEL:
                        revise_vids[var_id].add(sample_index)
                        revise_vids[v[0]].add(sample_index)

            logging.info('Getting RD_CN/EV')
            for variant in f:
                if variant.id in revise_vids:
                    sample_indices = revise_vids[variant.id]
                    self._update_rd_cn(variant, sample_indices)

        logging.info('Generating geno_normal_revise_dict')
        geno_normal_revise_dict = {}
        for var_id, samples_dict in overlap_test_text.items():
            for sample_index, v in samples_dict.items():
                new_val = None
                if sample_index not in revise_vids[v[0]]:
                    sys.stderr.write("{} {}\n".format(sample_index, v[0]))
                if v[1] == SVType.DUP and \
                        self.rd_cn[var_id][sample_index] == 2 and \
                        self.rd_cn[v[0]][sample_index] == 3:
                    new_val = 1
                elif v[1] == SVType.DEL and \
                        self.rd_cn[var_id][sample_index] == 2 \
                        and self.rd_cn[v[0]][sample_index] == 1:
                    new_val = 3

                if new_val:
                    if var_id not in geno_normal_revise_dict:
                        geno_normal_revise_dict[var_id] = {}
                    sample_id = self.sample_list[sample_index]
                    geno_normal_revise_dict[var_id][sample_id] = new_val

        return geno_normal_revise_dict


def main(args):
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    logging.info('Starting script')
    reviser = VCFReviser()
    filtered_vcf = args[1]
    intersected_bed = args[2]
    geno_normal_revise_dict = reviser.get_geno_normal_revise(filtered_vcf, intersected_bed)
    logging.info('Dumping dictionary')
    sys.stdout.write(json.dumps(geno_normal_revise_dict))
    logging.info('Done')


if __name__ == '__main__':
    main(sys.argv)
