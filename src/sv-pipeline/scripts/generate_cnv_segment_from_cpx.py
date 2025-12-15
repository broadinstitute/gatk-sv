import os
import sys
import argparse

# Implementation note:
# This code implements the same python code in the GenerateCnvSegmentFromCpx task without alterations.


def split_cpx_interval(info):
    # eg of info: DUP_chrY:3125606-3125667
    out = [info.split('_')[0], info.split('_')[1].split(':')[0]] + [int(i) for i in
                                                                    info.split('_')[1].split(':')[1].split('-')]
    return out[1:4] + [out[0]]


def get_sample_batch_map(sample_batch_pe_map):
    sample_to_batch = {}
    with open(sample_batch_pe_map, 'r') as inp:
        for line in inp:
            sample, batch, pe_file = line.strip().split("\t")
            sample_to_batch[sample] = batch
    return sample_to_batch


def readin_cpx_cnv(bed_input, sample_to_batch, min_depth_size):
    CPX_CNV = []
    f_bed = os.popen(r'''zcat %s''' % (bed_input))
    for line in f_bed:
        pin = line.strip().split('\t')
        if pin[0][0] == '#':
            pos_CPX_TYPE = pin.index('CPX_TYPE')
            pos_CPX_INTERVALS = pin.index('CPX_INTERVALS')
            pos_SOURCE = pin.index('SOURCE')
            pos_SAMPLES = pin.index('samples')
        else:
            interval = []
            if 'DEL_' in pin[pos_CPX_INTERVALS] or 'DUP_' in pin[pos_CPX_INTERVALS]:
                interval += [split_cpx_interval(i) + [pin[3]] for i in pin[pos_CPX_INTERVALS].split(",") if
                             "DEL_" in i or "DUP_" in i]
            if 'DEL_' in pin[pos_SOURCE] or 'DUP_' in pin[pos_SOURCE]:
                interval += [split_cpx_interval(i) + [pin[3]] for i in pin[pos_SOURCE].split(",") if
                             "DEL_" in i or "DUP_" in i]
            if len(interval) > 0:
                for i in interval:
                    if i[2] - i[1] >= min_depth_size:
                        sample_names = pin[pos_SAMPLES].split(',')
                        if i[3] == "DEL":
                            for j in sample_names:
                                CPX_CNV.append(i + [j, sample_to_batch[j]])
                        if i[3] == "DUP":
                            for j in sample_names:
                                CPX_CNV.append(i + [j, sample_to_batch[j]])
    f_bed.close()
    return CPX_CNV


def write_cpx_cnv(CPX_CNV, fileout):
    fo = open(fileout, 'w')
    for info in CPX_CNV:
        print('\t'.join([str(i) for i in info]), file=fo)
    fo.close()


def main(bed, sample_batch_pe_map, min_depth_size, prefix):
    fileout = f"{prefix}.lg_CNV.bed"
    sample_to_batch = get_sample_batch_map(sample_batch_pe_map)
    cpx_cnv = readin_cpx_cnv(bed, sample_to_batch, min_depth_size)
    write_cpx_cnv(cpx_cnv, fileout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True)
    parser.add_argument("--sample-batch-pe-map", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--min-depth-size", default=5000, type=int)
    args = parser.parse_args()

    main(args.bed, args.sample_batch_pe_map, args.min_depth_size, args.prefix)
