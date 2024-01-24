#!/bin/env python

def count_discontinuous(seq, reverse = False):
    # sort by position
    seq.sort(key=lambda a: a[0], reverse=reverse)
    # count changes in direction (expect 1)
    return len([i for i in range(1, len(seq)) if seq[i][0]!=seq[i-1][0]])


def count_discontinous_all(pe, pe_header):
    dir1 = [(x[pe_header["pos1"]], x[pe_header["dir1"]]) for x in pe]
    fwd = count_discontinuous(dir1)
    rev = count_discontinuous(dir1, reverse=True)
    dir2 = [(x[pe_header["pos2"]], x[pe_header["dir2"]]) for x in pe]
    fwd += count_discontinuous(dir2)
    rev += count_discontinuous(dir2, reverse=True)
    return fwd, rev


def count_patterns(pe, pe_header):
    counts = {"++": 0, "+-": 0, "-+": 0, "--": 0}
    for pair in pe:
        dir1 = pair[pe_header["dir1"]]
        dir2 = pair[pe_header["dir2"]]
        counts[dir1 + dir2] += 1
    return counts.values()


def calc_dist(chrom_seq, pos_seq):
    max_jump = 0
    max_jump_idx = 0
    prev_chrom = None
    prev_pos = None
    curr_idx = 0
    for chrom, pos in zip(chrom_seq, pos_seq):
        if curr_idx = 0:
            prev_chrom = chrom
            prev_pos = pos
        else:
            if chrom != prev_chrom:
                max_jump = 248956422
                max_jump_idx = curr_idx
            elif pos - prev_pos > max_jump:
                max_jump = pos - prev_pos
                max_jump_idx = curr_idx
        curr_idx += 1
    return pos_seq[max_jump_idx - 1] - pos_seq[0], pos_seq[-1] - pos_seq[max_jump_idx]


def calc_dist_all(pe, pe_header):
    dist1, dist3 = calc_dist([x[pe_header["chrom1"]] for x in pe], [int(x[pe_header["pos1"]]) for x in pe])
    dist2, dist4 = calc_dist([x[pe_header["chrom2"]] for x in pe], [int(x[pe_header["pos2"]]) for x in pe])
    return dist1, dist2, dist3, dist4


def evaluate(descriptor, pe, descriptor_header, pe_header):
    svid = descriptor[descriptor_header["name"]]
    sample = descriptor[descriptor_header["sample"]]
    if len(pe) == 0:
        return [svid, sample] + [0]*10
    fwd, rev = count_discontinous_all(pe, pe_header)
    pp, pm, mp, mm = count_patterns(pe, pe_header)
    dist1, dist2, dist3, dist4 = calc_dist_all(pe, pe_header)
    return [svid, sample, pp, pm, mp, mm, fwd, rev, dist1, dist2, dist3, dist4]


def process(pe_evidence, out_file)
    with open(pe_evidence, 'r') as pe, open(out_file, 'w') as out:
        out.write("\t".join("#SVID sample ++ +- -+ -- discont_fwd discont_rev dist1 dist1 dist2 dist3 dist4".split()) + "\n")
        first = True
        descriptor_header = None
        pe_header = {x:i for i,x in enumerate("chrom1 pos1 dir1 chrom2 pos2 dir2 sample".split())}
        curr_descriptor = None
        curr_pe = None
        for line in pe:
            fields = line.strip().lstrip("#").split("\t")
            if first:
                descriptor_header = {x:i for i, x in enumerate(fields)}
            elif line.startswith("#"):
                if curr_pe is not None:
                    curr_out = evaluate(curr_descriptor, curr_pe, descriptor_header, pe_header)
                    out.write("\t".join(curr_out) + "\n")
                curr_descriptor = fields
                curr_pe = []
            else:
                curr_pe.append(fields)
        # handle last variant
        if len(curr_pe) > 0:
            curr_out = evaluate(curr_descriptor, curr_pe, descriptor_header, pe_header)
            out.write("\t".join(curr_out) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "pe-evidence", required=True, help="PE evidence file for manual review")
    parser.add_argument("-o", "--out-file", required=True, help="Name for output table")
    args = parser.parse_args()

    process(args.pe_evidence, args.out_file)


if __name__ == '__main__':
    main()
