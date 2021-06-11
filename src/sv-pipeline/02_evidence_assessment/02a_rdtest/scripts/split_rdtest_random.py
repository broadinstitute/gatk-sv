#!python
# !python script to random shuffle and split rdtest big files:
import random
import numpy
import argparse
parser = argparse.ArgumentParser("split_rdtest_random.py")
parser.add_argument("input", type=str, help="namd of input to be splited")
parser.add_argument("output", type=str, help="prefix of output")
parser.add_argument("-s", '--size', type=str, help="size of outputs")

# define workding and file locations:
args = parser.parse_args()
file_name = args.input
out_prefix = args.output
out_size = 100
if args.size:
    out_size = int(args.size)
print(out_size)
info = []
fin = open(file_name)
for line in fin:
    pin = line.strip().split()
    if not pin[0][0] == '#':
        info.append(pin)
fin.close()
out = [i for i in range(len(info))]
random.shuffle(out)
out2 = []
for i in out:
    if out2 == []:
        out2.append([i])
    elif len(out2[-1]) < out_size:
        out2[-1].append(i)
    else:
        out2.append([i])

out_digits = int(numpy.log10(len(out2))) + 1
rec = -1
for i in out2:
    rec += 1
    fileout = out_prefix + format(rec, '0' + str(out_digits))
    print(rec, out_digits, format(rec, '0' + str(out_digits)))
    print(fileout)
    tmp = [info[x] for x in i]
    fo = open(fileout, 'w')
    for j in tmp:
        if int(j[2]) - int(j[1]) < 1000000 and int(j[2]) - int(j[1]) > 0:
            print('\t'.join(j), file=fo)
    fo.close()
