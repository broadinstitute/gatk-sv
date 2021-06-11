import pysam
import argparse


def calcu_inheri_hash(vcf_file, fam_file):
    fam_info = trio_info_readin(fam_file)
    fvcf = pysam.VariantFile(vcf_file)
    inheri_hash = {}
    for child in fam_info.keys():
        inheri_hash[child] = []
    for record in fvcf:
        print(record.id)
        for child in fam_info.keys():
            trio = fam_info[child] + [child]
            trio_len = sum(
                [1 if i in record.samples.keys() else 0 for i in trio])
            if trio_len == 3:
                gt = [record.samples[i]['GT'] for i in trio]
                if (None, None) in gt:
                    continue
                if gt == [(0, 0), (0, 0), (0, 0)]:
                    continue
                else:
                    if gt[1] == (0, 0) and not gt[0] == (0, 0) and not gt[2] == (0, 0):
                        inheri_hash[child].append(
                            ['fa_pb', record.info['SVTYPE']])
                    if gt[0] == (0, 0) and not gt[1] == (0, 0) and not gt[2] == (0, 0):
                        inheri_hash[child].append(
                            ['mo_pb', record.info['SVTYPE']])
                    if not gt[0] == (0, 0) and not gt[1] == (0, 0) and not gt[2] == (0, 0):
                        inheri_hash[child].append(
                            ['fa_mo_pb', record.info['SVTYPE']])
                    if gt[0] == (0, 0) and gt[1] == (0, 0) and not gt[2] == (0, 0):
                        inheri_hash[child].append(
                            ['denovo', record.info['SVTYPE']])
    fvcf.close()
    return inheri_hash


def inheri_hash_to_stat(inheri_hash):
    inheri_stat = {}
    for child in inheri_hash.keys():
        inheri_stat[child] = {}
        for rec in inheri_hash[child]:
            if not rec[1] in inheri_stat[child].keys():
                inheri_stat[child][rec[1]] = {}
            if not rec[0] in inheri_stat[child][rec[1]].keys():
                inheri_stat[child][rec[1]][rec[0]] = 0
            inheri_stat[child][rec[1]][rec[0]] += 1
    return inheri_stat


def trio_info_readin(fam_file):
    fam_info = {}
    fin = open(fam_file)
    for line in fin:
        pin = line.strip().split()
        if pin[2] == '0' and pin[3] == '0':
            continue
        if not pin[1] in fam_info.keys():
            fam_info[pin[1]] = pin[2:4]
    fin.close()
    return fam_info


def unique_list(list):
    out = []
    for i in list:
        if i not in out:
            out.append(i)
    return out


def write_output_stat(fileout, inheri_stat):
    fo = open(fileout, 'w')
    print('\t'.join(['sample', 'svtype', 'fa_mo_pb',
                     'fa_pb', 'mo_pb', 'denovo']), file=fo)
    for samp in inheri_stat.keys():
        for svt in inheri_stat[samp].keys():
            tmp = []
            for inh in ['fa_mo_pb', 'fa_pb', 'mo_pb', 'denovo']:
                if inh in inheri_stat[samp][svt].keys():
                    tmp.append(inheri_stat[samp][svt][inh])
                else:
                    tmp.append(0)
            print('\t'.join([str(i) for i in [samp, svt] + tmp]), file=fo)
    fo.close()


def main():
    parser = argparse.ArgumentParser("GATK-SV.S1.vcf2bed.py")
    parser.add_argument('fam_file', type=str, help='fam / ped file')
    parser.add_argument('vcf_file', type=str, help='vcf file')
    parser.add_argument('inheri_stat', type=str, help='name of output stat')
    args = parser.parse_args()
    # read_write_basic_vcf(args.vcfname,args.bedname)
    fam_file = args.fam_file
    vcf_file = args.vcf_file
    fileout = args.inheri_stat
    # readin fam information
    # only complete trios would be read in here
    inheri_hash = calcu_inheri_hash(vcf_file, fam_file)
    inheri_stat = inheri_hash_to_stat(inheri_hash)
    write_output_stat(fileout, inheri_stat)


if __name__ == '__main__':
    main()
