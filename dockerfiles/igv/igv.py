import sys
[_, varfile] = sys.argv
plotdir = "plots"
igvfile = "igv.txt"
igvsh = "igv.sh"
with open(varfile, 'r') as f:
    for line in f:
        dat = line.split('\t')
        chr = dat[0]
        start = dat[1]
        end = dat[2]
        data = dat[3].split(',')
