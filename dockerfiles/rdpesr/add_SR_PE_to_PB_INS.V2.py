import os


def INS_readin(filein):
    fin = open(filein)
    out = []
    for line in fin:
        pin = line.strip().split()
        if pin[0][0] == '#':
            continue
        # if pin[4]=='INS':
        out.append(pin)
    fin.close()
    return out


def add_Num_SR_le(sr_index, info, flank_length=100):
    # eg of info:    ['chr1', '137221', '137339', 'HOM', 'INS']
    fin = os.popen(r'''tabix %s %s:%d-%d''' % (sr_index,
                                               info[0], int(info[1]) - flank_length, int(info[1]) + flank_length))
    tmp = []
    for line in fin:
        pin = line.strip().split()
        tmp.append(pin)
    fin.close()
    if len(tmp) == 0:
        return 0
    else:
        return max([int(i[3]) for i in tmp])


def add_Num_SR_ri(sr_index, info, flank_length=100):
    # eg of info:    ['chr1', '137221', '137339', 'HOM', 'INS']
    fin = os.popen(r'''tabix %s %s:%d-%d''' % (sr_index,
                                               info[0], int(info[2]) - flank_length, int(info[2]) + flank_length))
    tmp = []
    for line in fin:
        pin = line.strip().split()
        tmp.append(pin)
    fin.close()
    if len(tmp) == 0:
        return 0
    else:
        return max([int(i[3]) for i in tmp])


def add_Num_PE_le(pe_index, info, flank_length=300):
    fin = os.popen(r'''tabix %s %s:%d-%d''' % (pe_index,
                                               info[0], int(info[1]) - 2 * flank_length, int(info[1]) + flank_length))
    tmp = []
    for line in fin:
        pin = line.strip().split()
        if 'INS' in pin[4] or pin[4] in ['INS', 'ALU', 'LINE1', 'SVA']:
            tmp.append(pin)
        else:
            if pin[0] == pin[3]:
                if abs(int(pin[4]) - int(pin[1])) > 100 * (int(info[2]) - int(info[1])):
                    continue
                else:
                    tmp.append(pin)
    fin.close()
    # if len(tmp)==0:
    #    return 0
    # else:
    #    cluster_hash= cluster_pe_mate(tmp)
    #    return cluster_hash[0]
    return len(tmp)


def add_Num_PE_ri(pe_index, info, flank_length=300):
    fin = os.popen(r'''tabix %s %s:%d-%d''' % (pe_index,
                                               info[0], int(info[2]) - flank_length, int(info[2]) + 2 * flank_length))
    tmp = []
    for line in fin:
        pin = line.strip().split()
        if 'INS' in pin[4] or pin[4] in ['INS', 'ALU', 'LINE1', 'SVA']:
            tmp.append(pin)
        else:
            if pin[0] == pin[3]:
                if abs(int(pin[4]) - int(pin[1])) > 100 * (int(info[2]) - int(info[1])):
                    continue
                else:
                    tmp.append(pin)
    fin.close()
    # if len(tmp)==0:
    #    return 0
    # else:
    #    cluster_hash= cluster_pe_mate(tmp)
    #    return cluster_hash[0]
    return len(tmp)


def cluster_pe_mate(tmp):
    out = {}
    for i in tmp:
        if not i[3] in out.keys():
            out[i[3]] = []
        out[i[3]].append(int(i[4]))
    key_name = [i for i in out.keys()]
    key_lengh = [len(out[i]) for i in key_name]
    most_abundant = key_name[key_lengh.index(max(key_lengh))]
    return [most_abundant, sorted(out[most_abundant])]


def write_Num_SR(info_list, fileout):
    fo = open(fileout, 'w')
    for i in info_list:
        print('\t'.join([str(j) for j in i]), file=fo)
    fo.close()


def main():
    import argparse
    parser = argparse.ArgumentParser("add_SR_PE_to_PB_INS.py")
    parser.add_argument('PB_bed', type=str,
                        help='name of input PacBio bed file')
    parser.add_argument('pe_file', type=str,
                        help='name of pe files with index')
    parser.add_argument('sr_file', type=str,
                        help='name of sr files with index')
    parser.add_argument('output', type=str,
                        help='name of output files with index')
    args = parser.parse_args()
    filein = args.PB_bed
    pe_index = args.pe_file
    sr_index = args.sr_file
    info_list = INS_readin(filein)
    for i in info_list:
        i += [add_Num_PE_le(pe_index, i)]
        i += [add_Num_PE_ri(pe_index, i)]
        if i[4] == 'INS' or i[4] == 'MEI':
            i += [add_Num_SR_le(sr_index, i, 50)]
            i += [add_Num_SR_ri(sr_index, i, 50)]
        else:
            if int(i[5]) < 300:
                i += [add_Num_SR_le(sr_index, i, int(i[5]) / 2)]
                i += [add_Num_SR_ri(sr_index, i, int(i[5]) / 2)]
            else:
                i += [add_Num_SR_le(sr_index, i, 150)]
                i += [add_Num_SR_ri(sr_index, i, 150)]
        i += [add_Num_SR_le(sr_index, i, 0)]
        i += [add_Num_SR_ri(sr_index, i, 0)]
    write_Num_SR(info_list, args.output)


if __name__ == '__main__':
    main()
