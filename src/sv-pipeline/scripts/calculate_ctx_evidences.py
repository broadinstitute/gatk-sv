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
    fin = open(filein)
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


def add_PE_supp(sample_svid, PE_supp, min_pe_ctx):
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
                    if PE_supp[svid][sample][0][0] >= min_pe_ctx and PE_supp[svid][sample][1][0] >= min_pe_ctx:
                        pe_supp = 'high_PE'
                    out[svid][sample] = [pe_supp]
                else:
                    out[svid][sample] = ['no_PE']
        else:
            for sample in sample_svid[svid]:
                out[svid][sample] = ['no_PE']
    return out


def write_pe_depth_supp(sample_svid_pe_depth, fileout):
    fo = open(fileout, 'w')
    print('\t'.join(['VID', 'sample', 'supportive_PE_counts']), file=fo)
    for svid in sample_svid_pe_depth.keys():
        for sample in sample_svid_pe_depth[svid].keys():
            print('\t'.join([svid, sample] + sample_svid_pe_depth[svid][sample]), file=fo)
    fo.close()


def main(sample_svid, pe_supp, prefix, min_pe_ctx):
    sample_svid = sample_svid_readin(sample_svid)
    PE_supp = PE_supp_readin(pe_supp)
    sample_svid_pe = add_PE_supp(sample_svid, PE_supp, min_pe_ctx)
    write_pe_depth_supp(sample_svid_pe, f"{prefix}.manual_revise.CTX_results")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-svid", required=True)
    parser.add_argument("--pe-supp", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--min-pe-ctx", type=int, required=True)
    args = parser.parse_args()

    main(args.sample_svid, args.pe_supp, args.prefix, args.min_pe_ctx)
