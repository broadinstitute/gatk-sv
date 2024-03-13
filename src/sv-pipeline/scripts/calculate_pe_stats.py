#!/bin/env python
import argparse
import gzip


def count_discontinuous(pe, sort_on = 0):
    DIR1 = 2
    DIR2 = 5
    # sort by position then chromosome (first in pair if sort_on = 0, second in pair if sort_on = 3)
    # column indexes: sort_on = contig, sort_on + 1 = position, sort_on + 2 = +/-
    pe_sort = sorted(pe, key=lambda x: (x[sort_on + 1], x[sort_on]))
    # count changes in direction (expect 2 for each sorting method, one per column)
    dir1_discont = len([i for i in range(1, len(pe_sort)) if pe_sort[i][DIR1]!=pe_sort[i-1][DIR1]])
    dir2_discont = len([i for i in range(1, len(pe_sort)) if pe_sort[i][DIR2]!=pe_sort[i-1][DIR2]])
    return dir1_discont + dir2_discont


def count_discontinous_all(pe, pe_header):
    if len(pe) == 0:
        return "NA", "NA"
    return count_discontinuous(pe, 0), count_discontinuous(pe, 3)


def count_patterns(pe, pe_header):
    counts = {"++": 0, "+-": 0, "-+": 0, "--": 0}
    for pair in pe:
        dir1 = pair[pe_header["dir1"]]
        dir2 = pair[pe_header["dir2"]]
        counts[dir1 + dir2] += 1
    return counts.values()


def determine_expected_patterns(descriptor, descriptor_header):
    if "INV" in descriptor[descriptor_header["SVTYPE"]] or \
        "INV" in descriptor[descriptor_header["CPX_TYPE"]] or \
        descriptor[descriptor_header["CPX_TYPE"]] == "CTX_PQ/QP":
        # "INV" in descriptor[descriptor_header["SOURCE"]] or \
        return ["++", "--"]
    else:
        return ["+-", "-+"]


def calc_dist_all(pe, pe_header, expected_patterns):
    pattern1_chrom1 = []
    pattern1_chrom2 = []
    pattern2_chrom1 = []
    pattern2_chrom2 = []

    # get positions of all pairs matching expected patterns
    for pair in pe:
        dir1 = pair[pe_header["dir1"]]
        dir2 = pair[pe_header["dir2"]]
        pattern = dir1 + dir2
        if pattern == expected_patterns[0]:
            pattern1_chrom1.append(int(pair[pe_header["pos1"]]))
            pattern1_chrom2.append(int(pair[pe_header["pos2"]]))
        elif pattern == expected_patterns[1]:
            pattern2_chrom1.append(int(pair[pe_header["pos1"]]))
            pattern2_chrom2.append(int(pair[pe_header["pos2"]]))

    # sort and calculate distances
    dists = []
    for pos_list in [pattern1_chrom1, pattern1_chrom2, pattern2_chrom1, pattern2_chrom2]:
        if len(pos_list) < 1:
            dists.append("NA")
        else:
            pos_list.sort()
            dists.append(pos_list[-1] - pos_list[0])
    return dists


def evaluate(descriptor, pe, descriptor_header, pe_header, background):
    svid = descriptor[descriptor_header["name"]]
    sample = descriptor[descriptor_header["sample"]]
    cpx_type = descriptor[descriptor_header["CPX_TYPE"]]
    # if len(pe) == 0:
    #     return [svid, sample, cpx_type] + ["0"]*4 + ["NA"]*6
        # return [svid, sample, cpx_type] + ["0"]*4 + ["NA"]*2
    discont_sortonfirst, discont_sortonsecond = count_discontinous_all(pe, pe_header)
    pp, pm, mp, mm = count_patterns(pe, pe_header)
    dist1, dist2, dist3, dist4 = calc_dist_all(pe, pe_header, determine_expected_patterns(descriptor, descriptor_header))
    background_fields = []
    if background is not None:
        background_fields = background[svid]
    return [svid, sample, cpx_type] + [str(x) for x in [pp, pm, mp, mm, discont_sortonfirst, discont_sortonsecond, dist1, dist2, dist3, dist4]] + background_fields


def process(pe_evidence, out_file, background):
    with gzip.open(pe_evidence, 'rt') as pe, open(out_file, 'w') as out:
        if background is not None:
            out.write("\t".join("#SVID sample CPX_TYPE ++ +- -+ -- discont_sortonfirst discont_sortonsecond dist1 dist2 dist3 dist4 background_min1 background_min4 background_min10".split()) + "\n")
        else:
            out.write("\t".join("#SVID sample CPX_TYPE ++ +- -+ -- discont_sortonfirst discont_sortonsecond dist1 dist2 dist3 dist4".split()) + "\n")
        first = True
        descriptor_header = None
        pe_header = {x:i for i,x in enumerate("chrom1 pos1 dir1 chrom2 pos2 dir2 sample".split())}
        curr_descriptor = None
        curr_pe = None
        for line in pe:
            fields = line.strip().lstrip("#").split("\t")
            if first:
                descriptor_header = {x:i for i, x in enumerate(fields)}
                first = False
            elif line.startswith("#"):
                if curr_pe is not None:
                    curr_out = evaluate(curr_descriptor, curr_pe, descriptor_header, pe_header, background)
                    out.write("\t".join(curr_out) + "\n")
                curr_descriptor = fields
                curr_pe = []
            else:
                curr_pe.append(fields)
        # handle last variant
        if len(curr_pe) > 0:
            curr_out = evaluate(curr_descriptor, curr_pe, descriptor_header, pe_header, background)
            out.write("\t".join(curr_out) + "\n")


def load_background(background_file):
    background = dict()  # {svid: [background_min1, background_min4, background_min10]}
    with gzip.open(background_file, 'rt') as bg:
        for line in bg:
            fields = line.strip().lstrip("#").split("\t")
            if not line.startswith("#"):
                svid = fields[0]
                background[svid] = fields[1:]
    return background


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pe-evidence", required=True, help="PE evidence file for manual review")
    parser.add_argument("-b", "--background-pe", required=False, help="Background PE stats")
    parser.add_argument("-o", "--out-file", required=True, help="Name for output table")
    args = parser.parse_args()

    background = None
    if args.background_pe is not None:
        background = load_background(args.background_pe)
    process(args.pe_evidence, args.out_file, background)


if __name__ == '__main__':
    main()
