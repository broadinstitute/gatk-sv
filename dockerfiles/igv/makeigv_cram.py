import os
import argparse
# [_,varfile,buff,fasta]=sys.argv #assume the varfile has *.bed in the end
# Example
# python makeigv.py /data/talkowski/xuefang/local/src/IGV_2.4.14/IL_DUP/IL.DUP.HG00514.V2.bed /data/talkowski/Samples/1000Genomes/HGSV_Illumina_Alignment_GRCh38 400
# bash IL.DUP.HG00514.V2.sh
# bash igv.sh -b IL.DUP.HG00514.V2.txt


parser = argparse.ArgumentParser("makeigvsplit_cram.py")
parser.add_argument('varfile', type=str,
                    help='name of variant file in bed format, with cram and SVID in last two columns')
parser.add_argument(
    'buff', type=str, help='length of buffer to add around variants')
parser.add_argument('fasta', type=str, help='reference sequences')

parser.add_argument('sample', type=str, help='name of sample to make igv on')
parser.add_argument('chromosome', type=str,
                    help='name of chromosome to plot igv on')
args = parser.parse_args()


buff = int(args.buff)
fasta = args.fasta
varfile = args.varfile

outstring = os.path.basename(varfile)[0:-4]
bamdir = "bam"
outdir = "screenshot"
igvfile = "igv.txt"
bamfiscript = "igv.sh"
###################################
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
                if not Chr == args.chromosome:
                    continue
                Start = str(int(dat[1]) - buff)
                End = str(int(dat[2]) + buff)
                Dat = dat[3].split(',')
                ID = dat[4]
                for cram in Dat:
                    # sample=cram.split("/")[-1].split('.')[0]
                    g.write('load ' + bamdir + '/' + args.sample +
                            '_' + args.chromosome + '.bam\n')
                if int(End) - int(Start) < 10000:
                    g.write('goto ' + Chr + ":" + Start + '-' + End + '\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('collapse\n')
                    g.write('snapshotDirectory ' + outdir + '\n')
                    g.write('snapshot ' + args.sample + '_' + ID + '.png\n')
                else:
                    # Extra 1kb buffer if variant large
                    g.write('goto ' + Chr + ":" + Start +
                            '-' + str(int(Start) + 1000) + '\n')
                    g.write('sort base\n')
                    g.write('viewaspairs\n')
                    g.write('collapse\n')
                    g.write('snapshotDirectory ' + outdir + '\n')
                    g.write('snapshot ' + args.sample +
                            '_' + ID + '.left.png\n')
                    g.write('goto ' + Chr + ":" +
                            str(int(End) - 1000) + '-' + End + '\n')
                    g.write('sort base\n')
                    g.write('collapse\n')
                    g.write('snapshotDirectory ' + outdir + '\n')
                    g.write('snapshot ' + args.sample +
                            '_' + ID + '.right.png\n')
                # g.write('goto '+Chr+":"+Start+'-'+End+'\n')
                # g.write('sort base\n')
                # g.write('viewaspairs\n')
                # g.write('collapse\n')
                # g.write('snapshotDirectory '+outdir+'\n')
                # g.write('snapshot '+ID+'.png\n' )
                g.write('new\n')
        g.write('exit\n')
# with open(bamfiscript,'w') as g:
