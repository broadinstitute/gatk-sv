#script to add cov to SVs

def add_ILL_cov(pb_uni_svs,bincov):
        for i in pb_uni_svs.keys():
                for j in pb_uni_svs[i]:
                        cov_list=cov_SV_readin(j, bincov)
                        j+=[len(cov_list),np.median(cov_list), np.mean(cov_list),np.std(cov_list)]
                        print(j)
        return pb_uni_svs

def bed_info_readin(input):
        fin=open(input)
        out={}
        for line in fin:
                pin=line.strip().split()
                if pin[0][0]=='#': continue
                if not pin[0] in out.keys():
                        out[pin[0]]=[]
                out[pin[0]].append([pin[0],int(pin[1]),int(pin[2])]+pin[3:])
        fin.close()
        return out

def cov_SV_readin(svpos, bincov):
        fin=os.popen(r'''tabix %s %s:%d-%d'''%(bincov, svpos[0],svpos[1],svpos[2]))
        normCov_list=[]
        for line in fin:
                pin=line.strip().split()
                normCov_list.append(float(pin[-1]))
        fin.close() 
        return normCov_list
          
def path_modify(path):
        if not path[-1]=='/':
                path+='/'
        return path 

def write_output(input,pb_uni_svs):
        fo=open(input+'.Seq_Cov','w') 
        for k1 in pb_uni_svs.keys(): 
                for k2 in pb_uni_svs[k1]:
                        print('\t'.join([str(i) for i in k2]),file=fo)
        fo.close() 

def main():
        parser = argparse.ArgumentParser(description='S2a.calcu.Seq_Cov.of.PB_Uni.py')
        parser.add_argument('input', help='name of input file containing PacBio unique SVs in bed format')
        parser.add_argument('bincov',help='name of bincov metrics of the sample to be processed')
        args = parser.parse_args()
        pb_uni_svs=bed_info_readin(args.input)
        pb_uni_svs=add_ILL_cov(pb_uni_svs,args.bincov)
        write_output(args.input,pb_uni_svs)

import os
import numpy as np
import argparse
main()



