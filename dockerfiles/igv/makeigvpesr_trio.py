import os
import argparse
# [_,varfile,buff,fasta]=sys.argv #assume the varfile has *.bed in the end
# Usage
# python makeigvpesr_cram.py varfile fasta sample ped cram_list buffer chromosome
# bash IL.DUP.HG00514.V2.sh
# bash igv.sh -b IL.DUP.HG00514.V2.txt


parser = argparse.ArgumentParser("makeigvsplit_cram.py")
parser.add_argument('varfile', type=str,
                    help='name of variant file in bed format, with cram and SVID in last two columns')
parser.add_argument('fasta', type=str, help='reference sequences')
parser.add_argument('sample', type=str, help='name of sample to make igv on')
parser.add_argument('cram_list', type=str,
                    help='a file including sample and cram path')
parser.add_argument('outdir', type=str, help='output folder')
parser.add_argument('-b', '--buff', type=str,
                    help='length of buffer to add around variants', default=500)
parser.add_argument('-c', '--chromosome', type=str,
                    help='name of chromosome to make igv on', default='all')

args = parser.parse_args()


buff = int(args.buff)
fasta = args.fasta
varfile = args.varfile

outstring = os.path.basename(varfile)[0:-4]
bamdir = "pe_bam"
outdir = args.outdir
igvfile = "pe.txt"
bamfiscript = "pe.sh"
###################################
sample = args.sample
chromosome = args.chromosome


def ped_info_readin(ped_file):
    out = {}
    fin = open(ped_file)
    for line in fin:
        pin = line.strip().split()
        if not pin[1] in out.keys():
            out[pin[1]] = [pin[1]]
        if not(pin[2]) == 0:
            out[pin[1]].append(pin[2])
        if not(pin[3]) == 0:
            out[pin[1]].append(pin[3])
    fin.close()
    return out


def cram_info_readin(cram_file):
    out = {}
    fin = open(cram_file)
    for line in fin:
        pin = line.strip().split()
        if not pin[0] in out.keys():
            out[pin[0]] = pin[1:]
    fin.close()
    return(out)


# ped_info = ped_info_readin(args.ped)
# cram_info = cram_info_readin(args.cram_list)
cram_list = args.cram_list.split(',')

with open(bamfiscript, 'w') as h:
    h.write("#!/bin/bash\n")
    h.write("set -e\n")
    h.write("mkdir -p {}\n".format(bamdir))
    h.write("mkdir -p {}\n".format(outdir))
    with open(igvfile, 'w') as g:
        g.write('new\n')
        g.write('genome {}\n'.format(fasta))
        with open(varfile, 'r') as f:
            for line in f:
                dat = line.rstrip().split("\t")
                Chr = dat[0]
                if not chromosome == 'all':
                    if not Chr == chromosome:
                        continue
                Start = str(int(dat[1]) - buff)
                End = str(int(dat[2]) + buff)
                ID = dat[4]
                for cram in cram_list:
                    g.write('load ' + cram + '\n')
                if int(End) - int(Start) < 10000:
                    g.write('goto ' + Chr + ":" + Start + '-' + End + '\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('squish\n')
                    g.write('snapshotDirectory ' + outdir + '\n')
                    g.write('snapshot ' + sample + '_' + ID + '.png\n')
                else:
                    # Extra 1kb buffer if variant large
                    g.write('goto ' + Chr + ":" + Start +
                            '-' + str(int(Start) + 1000) + '\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('squish\n')
                    g.write('snapshotDirectory ' + outdir + '\n')
                    g.write('snapshot ' + sample + '_' + ID + '.left.png\n')
                    g.write('goto ' + Chr + ":" +
                            str(int(End) - 1000) + '-' + End + '\n')
                    g.write('sort base\n')
                    g.write('squish\n')
                    g.write('snapshotDirectory ' + outdir + '\n')
                    g.write('snapshot ' + sample + '_' + ID + '.right.png\n')
                # g.write('goto '+Chr+":"+Start+'-'+End+'\n')
                # g.write('sort base\n')
                # g.write('viewaspairs\n')
                # g.write('squish\n')
                # g.write('snapshotDirectory '+outdir+'\n')
                # g.write('snapshot '+ID+'.png\n' )
                g.write('new\n')
        g.write('exit\n')
# with open(bamfiscript,'w') as g:
