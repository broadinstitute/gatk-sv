import os


def INS_readin(filein):
    fin = open(filein)
    out = []
    for line in fin:
        pin = line.strip().split()
        # if pin[4]=='INS':
        out.append(pin)
    fin.close()
    return out


def add_Num_SR(sr_index, info, flank_length=50):
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


def add_Num_PE(pe_index, info, flank_length=100):
    fin = os.popen(r'''tabix %s %s:%d-%d''' % (pe_index,
                                               info[0], int(info[1]) - flank_length, int(info[1]) + flank_length))
    tmp = []
    for line in fin:
        pin = line.strip().split()
        # if int(pin[2])-int(pin[1])>100*(int(info[2])-int(info[1])): continue
        tmp.append(pin)
    fin.close()
    # if len(tmp)==0:
    #    return 0
    # else:
    #    cluster_hash= cluster_pe_mate(tmp)
    #    return cluster_hash[0]
    return len(tmp)


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
    # parser.add_argument('sr_file', type=str, help='name of sr files with index')
    args = parser.parse_args()
    filein = args.PB_bed
    pe_index = args.pe_file
    # sr_index = args.sr_file
    info_list = INS_readin(filein)
    for i in info_list:
        # i+=[add_Num_SR(sr_index,i)]
        i += [add_Num_PE(pe_index, i)]
        print(i)
    write_Num_SR(info_list, filein + '.with_INS_PE')


if __name__ == '__main__':
    main()
