import os
import argparse

# Implementation note:
# This code implements the same python code in the CalculateCpxEvidences task without alterations.


def sample_svid_readin(filein):
    out = {}
    fin = open(filein)
    for line in fin:
        pin = line.strip().split()
        if not pin[1] in out.keys():
            out[pin[1]] = []
        out[pin[1]].append(pin[0])
    fin.close()
    return out


def PE_supp_readin(filein):
    fin = os.popen(r'''zcat %s''' % (filein))
    out = {}
    for line in fin:
        pin = line.strip().split()
        if not pin[4] in out.keys():
            out[pin[4]] = {}
        if not pin[3] in out[pin[4]].keys():
            out[pin[4]][pin[3]] = [[], []]
        if pin[1] == "+":
            out[pin[4]][pin[3]][0] = [int(pin[0])] + pin[1:3]
        elif pin[1] == "-":
            out[pin[4]][pin[3]][1] = [int(pin[0])] + pin[1:3]
    fin.close()
    return out


def Depth_supp_readin(filein):
    fin = os.popen(r'''zcat %s''' % (filein))
    out = {}
    for line in fin:
        pin = line.strip().split()
        if not pin[0][0] == '#':
            if not pin[4] in out.keys():
                out[pin[4]] = {}
            if not pin[5] in out[pin[4]].keys():
                out[pin[4]][pin[5]] = []
            out[pin[4]][pin[5]].append(float(pin[7]))
    fin.close()
    return out


def add_PE_supp(sample_svid, PE_supp, min_pe_cpx):
    out = {}
    for svid in sample_svid.keys():
        out[svid] = {}
        if svid in PE_supp.keys():
            for sample in sample_svid[svid]:
                if sample in PE_supp[svid].keys():
                    pe_supp = 'no_PE'
                    if len(PE_supp[svid][sample][0]) == 0:
                        PE_supp[svid][sample][0].append(0)
                    if len(PE_supp[svid][sample][1]) == 0:
                        PE_supp[svid][sample][1].append(0)
                    if PE_supp[svid][sample][0][0] > 0 or PE_supp[svid][sample][1][0] > 0:
                        pe_supp = 'partial_PE'
                    if PE_supp[svid][sample][0][0] > 0 and PE_supp[svid][sample][1][0] > 0:
                        pe_supp = 'low_PE'
                    if PE_supp[svid][sample][0][0] >= min_pe_cpx and PE_supp[svid][sample][1][0] >= min_pe_cpx:
                        pe_supp = 'high_PE'
                    out[svid][sample] = [pe_supp]
                else:
                    out[svid][sample] = ['no_PE']
        else:
            for sample in sample_svid[svid]:
                out[svid][sample] = ['no_PE']
    return out


def add_depth_supp(sample_svid_pe, Depth_supp):
    for svid in sample_svid_pe.keys():
        if svid in Depth_supp.keys():
            for sample in sample_svid_pe[svid].keys():
                if sample in Depth_supp[svid].keys():
                    depth_supp = 'lack_depth'
                    if Depth_supp[svid][sample][0] > .5:
                        depth_supp = 'depth'
                        if len(Depth_supp[svid][sample]) > 1:
                            if not Depth_supp[svid][sample][1] > .5:
                                depth_supp += ',lack_depth'
                    else:
                        if len(Depth_supp[svid][sample]) > 1:
                            if Depth_supp[svid][sample][1] > .5:
                                depth_supp += ',depth'
                    sample_svid_pe[svid][sample].append(depth_supp)
                else:
                    sample_svid_pe[svid][sample].append('NA')
        else:
            for sample in sample_svid_pe[svid].keys():
                sample_svid_pe[svid][sample].append('NA')
    return sample_svid_pe


def write_pe_depth_supp(sample_svid_pe_depth, fileout):
    fo = open(fileout, 'w')
    print('\t'.join(['VID', 'sample', 'supportive_PE_counts', 'depth_supp']), file=fo)
    for svid in sample_svid_pe_depth.keys():
        for sample in sample_svid_pe_depth[svid].keys():
            print('\t'.join([svid, sample] + sample_svid_pe_depth[svid][sample]), file=fo)
    fo.close()


def main(sample_svid, pe_supp, depth_supp, prefix, min_pe_cpx):
    sample_svid = sample_svid_readin(sample_svid)
    PE_supp = PE_supp_readin(pe_supp)
    Depth_supp = Depth_supp_readin(depth_supp)
    sample_svid_pe = add_PE_supp(sample_svid, PE_supp, min_pe_cpx)
    sample_svid_pe_depth = add_depth_supp(sample_svid_pe, Depth_supp)
    write_pe_depth_supp(sample_svid_pe_depth, f"{prefix}.manual_revise.CPX_results")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-svid", required=True)
    parser.add_argument("--pe-supp", required=True)
    parser.add_argument("--depth-supp", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--min-pe-cpx", type=int, required=True)
    args = parser.parse_args()

    main(args.sample_svid, args.pe_supp, args.depth_supp, args.prefix, args.min_pe_cpx)
