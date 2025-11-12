import os
import argparse


def main(args):
    fin = os.popen(r'''zcat %s''' % (args.vcf))
    header = []
    body = {}
    SVID_key = []
    for line in fin:
        pin = line.strip().split()
        if pin[0][:2] == '##':
            header.append(pin)
        else:
            body[pin[2]] = pin
            SVID_key.append(pin[2])
    header.append(['##INFO=<ID=' + args.ref_prefix + '_SVID' + ',Number=1,Type=String,Description="SVID of an overlapping event in gnomad used for external allele frequency annotation.">'])

    fin.close()
    fin = open(args.labeled_bed)
    colname = fin.readline().strip().split()

    for j in range(len(colname) - 1):
        if j > 1:
            header.append(['##INFO=<ID=' + args.ref_prefix + '_' + colname[j] + ',Number=1,Type=Float,Description="Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad.">'])

    for line in fin:
        pin = line.strip().split()
        if pin[0] == 'query_svid': continue
        if not pin[0] in body.keys(): continue
        info_add = [args.ref_prefix + '_SVID' + '=' + pin[1]]
        for j in range(len(colname) - 1):
            if j > 1:
                info_add.append(args.ref_prefix + '_' + colname[j] + '=' + pin[j])
        body[pin[0]][7] += ';' + ';'.join(info_add)
    fin.close()

    output_filename = args.output_filename
    fo = open(output_filename, 'w')
    for i in header:
        print(' '.join(i), file=fo)
    for i in SVID_key:
        print('\t'.join(body[i]), file=fo)
    fo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--ref-prefix", required=True)
    parser.add_argument("--labeled-bed", required=True)
    parser.add_argument("--output-filename", required=True)
    args = parser.parse_args()
    main(args)
