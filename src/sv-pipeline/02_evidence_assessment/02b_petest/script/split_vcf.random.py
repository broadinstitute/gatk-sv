import os
import argparse
import numpy
import random


def filein_readin(filein, file_size):
    fin = os.popen(r'''zcat %s''' % (filein))
    header = []
    header2 = []
    info = []
    for line in fin:
        pin = line.strip().split()
        if pin[0][:2] == '##':
            header.append(pin)
        elif pin[0][0] == '#':
            header2.append(pin)
        else:
            info.append(pin)
    info_index = [i for i in range(len(info))]
    random.shuffle(info_index)
    info_new = []
    for i in info_index:
        if info_new == []:
            info_new.append([info[i]])
        elif len(info_new[-1]) < file_size:
            info_new[-1].append(info[i])
        else:
            info_new.append([info[i]])
    return [header, header2, info_new]


def fileout_write(header, header2, info_sub, fileout):
    fo = open(fileout, 'w')
    for i in header:
        print(' '.join(i), file=fo)
    for i in header2:
        print('\t'.join(i), file=fo)
    for i in info_sub:
        print('\t'.join(i), file=fo)
    fo.close()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser = argparse.ArgumentParser("split_pesrtest_random.py")
    parser.add_argument("input", type=str,
                        help="namd of input vcf.gz to be splited")
    parser.add_argument("output", type=str, help="prefix of output")
    parser.add_argument("-s", '--size', type=str, help="size of outputs")
    args = parser.parse_args()
    filein = args.input
    out_prefix = args.output
    file_size = int(args.size)
    [header, header2, info] = filein_readin(filein, file_size)
    out_digits = int(numpy.log10(len(info))) + 1
    rec = -1
    for i in info:
        rec += 1
        fileout = out_prefix + format(rec, '0' + str(out_digits))
        fileout_write(header, header2, i, fileout)


if __name__ == '__main__':
    main()
