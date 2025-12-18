import pysam
import argparse

# Implementation note:
# This code implements the same python code in the CalculateCpxEvidences task without alterations.


def CPX_manual_readin(CPX_manual):
    out = {}
    fin = open(CPX_manual)
    for line in fin:
        pin = line.strip().split()
        if not pin[0] in out.keys():
            out[pin[0]] = {}
        out[pin[0]][pin[1]] = pin[2:]
    fin.close()
    return out


def CTX_manual_readin(CTX_manual):
    fin = open(CTX_manual)
    out = {}
    for line in fin:
        pin = line.strip().split('\t')
        if not pin[0] in out.keys():
            out[pin[0]] = {}
        if not pin[1] in out[pin[0]].keys():
            out[pin[0]][pin[1]] = pin[2:]
    fin.close()
    return out


def unresolved_readin(unresolved_svids):
    svids = set()
    with open(unresolved_svids, 'r') as inp:
        for line in inp:
            svids.add(line.strip())
    return svids


def revise_vcf(vcf_input, vcf_output, hash_CPX_manual, unresolved_svids, hash_CTX_manual):
    fin = pysam.VariantFile(vcf_input)
    # revise vcf header
    header = fin.header
    fo = pysam.VariantFile(vcf_output, 'w', header=header)
    for record in fin:
        if record.id in unresolved_svids:
            if 'PASS' in record.filter:
                record.filter.clear()
            record.filter.add('UNRESOLVED')
        # label CPX with manual review results:
        elif record.id in hash_CPX_manual.keys():
            unresolve_rec = 0
            for sample in hash_CPX_manual[record.id].keys():
                if sample in record.samples.keys():
                    if hash_CPX_manual[record.id][sample][0] == 'no_PE':
                        if hash_CPX_manual[record.id][sample][1] in ["NA", "lack_depth", "lack_depth,lack_depth", "lack_depth,depth", "depth,lack_depth"]:
                            record.samples[sample]['GT'] = [None, None]
                    elif hash_CPX_manual[record.id][sample][0] == 'low_PE':
                        if hash_CPX_manual[record.id][sample][1] in ["NA", "lack_depth", "lack_depth,lack_depth", "lack_depth,depth", "depth,lack_depth"]:
                            record.samples[sample]['GT'] = [None, None]
                    elif hash_CPX_manual[record.id][sample][0] == 'partial_PE':
                        if hash_CPX_manual[record.id][sample][1] in ["NA", "lack_depth", "lack_depth,lack_depth", "lack_depth,depth", "depth,lack_depth"]:
                            record.samples[sample]['GT'] = [None, None]
                            unresolve_rec += 1
                    elif hash_CPX_manual[record.id][sample][0] == 'high_PE':
                        if hash_CPX_manual[record.id][sample][1] in ["lack_depth", "lack_depth,lack_depth", "lack_depth,depth", "depth,lack_depth"]:
                            record.samples[sample]['GT'] = [None, None]
            if not unresolve_rec / len(hash_CPX_manual[record.id].keys())<.5:
                if 'PASS' in record.filter:
                    record.filter.clear()
                record.filter.add('UNRESOLVED')
        # label CTX with manual review results:
        elif record.id in hash_CTX_manual.keys():
            if 'NA' in hash_CTX_manual[record.id].keys():
                record.filter.add('UNRESOLVED')
            else:
                for sample in hash_CTX_manual[record.id].keys():
                    if sample in record.samples.keys():
                        if hash_CTX_manual[record.id][sample][0] in ['no_PE', 'low_PE', 'partial_PE']:
                            record.samples[sample]['GT'] = [None, None]
        # revisions to insertion-type complex events
        if record.info['SVTYPE'] == "CPX" and record.info['CPX_TYPE'] in ['dDUP', 'dDUP_iDEL', 'INS_iDEL']:
            if record.info['CPX_TYPE'] == 'INS_iDEL':
                record.info['CPX_INTERVALS'] = ','.join([x for x in record.info['CPX_INTERVALS']] + [record.info['SOURCE']])
            del record.info['SOURCE']
        elif record.info['SVTYPE'] == "INS" and 'SOURCE' in record.info.keys() and "INV" in record.info['SOURCE']:
            record.info['SVTYPE'] = "CPX"
            record.alts = ('<CPX>',)
            del_section = record.stop - record.pos
            if del_section < 50:
                record.info['CPX_TYPE'] = "dDUP"
                record.info['CPX_INTERVALS'] = record.info['SOURCE'].replace('INV', 'DUP') + ',' + record.info['SOURCE']
            else:
                record.info['CPX_TYPE'] = "dDUP_iDEL"
                record.info['CPX_INTERVALS'] = record.info['SOURCE'].replace('INV', 'DUP') + ',' + record.info['SOURCE'] + ',' + "DEL_" + record.chrom + ":" + str(record.pos) + '-' + str(record.stop)
            del record.info['SOURCE']
        fo.write(record)  # write out every record that was in the input - NCR will remove ones with no carriers left
    fin.close()
    fo.close()


def main(unresolved_svids, CPX_manual, CTX_manual, vcf_file, prefix):
    unresolved_svids = unresolved_readin(unresolved_svids)
    hash_CPX_manual = CPX_manual_readin(CPX_manual)
    hash_CTX_manual = CTX_manual_readin(CTX_manual)
    print(len(hash_CPX_manual.keys()))
    revise_vcf(vcf_file, f"{prefix}.Manual_Revised.vcf.gz", hash_CPX_manual, unresolved_svids, hash_CTX_manual)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--unresolved-svids", required=True)
    parser.add_argument("--CPX-manual", required=True)
    parser.add_argument("--CTX-manual", required=True)
    parser.add_argument("--vcf-file", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    main(args.unresolved_svids, args.CPX_manual, args.CTX_manual, args.vcf_file, args.prefix)
