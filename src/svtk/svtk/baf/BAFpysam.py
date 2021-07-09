#!/usr/bin/env python
from scipy import stats
import numpy as np
from sklearn import mixture


# calculate the Del statistic given a FME combo in het files
def Deltest(F, M, E, length, crit=0.01, thres1=0.0005):
    # if True:
    thres1 = min(50 / length, thres1)
    if F / length < thres1 and M / length < thres1  \
            or E / length < thres1 and M / length < thres1:
        return "ROH"
    else:
        flank = min(F, E)
        ratio = np.log10((M + thres1 * length) / (flank + thres1 * length))
#         if m<flank:
#             print()
#         if ratio<0:
#             print(ratio)
        return ratio


def ROH(F, M, E, length, thres=0.0001):
    if min(F, M, E) < length * thres:
        return True
    else:
        return False


class DeletionTest:

    def __init__(self, obj, probands, length):
        self.length = length  # length of SV
        self.obj = obj  # het file
        # self.homobkgrd=0
        # self.hetbkgrd=0
        self.probands = probands  # python list of proband IDs
        self.count = {}  # record of FME count and Deltest statistic for everyone
        self.nullratio = []  # list of del statistic for non-ROH controls
        # if not os.path.isfile(self.regionfile):
        # raise ValueError("file not found")

        if self.obj.shape[0] == 0:
            self.nullavg = 'nan'
            self.ns = 0
        else:
            nsROH = 0  # total number of SNP in nonROH controls in SV region
            ns = 0  # total number of SNPs in SV region
            for index, row in self.obj.iterrows():
                F = row['before']
                M = row['inside']
                E = row['after']
                # print(dat[-1])
                self.count[row['sample']] = {
                    'F': F, 'M': M, 'E': E, 'Ratio': Deltest(F, M, E, self.length)}
                if row['sample'] not in self.probands and Deltest(F, M, E, self.length) != 'ROH':
                    self.nullratio.append(Deltest(F, M, E, self.length))
                    # print(F,M,E)
                    nsROH += M
                ns += M

            self.nullavg = nsROH / (len(self.nullratio) + 1)
            self.ns = ns
            self.nullratio = np.array(self.nullratio).reshape(-1, 1)
            if len(self.nullratio) > 10:
                self.gmm = mixture.BayesianGaussianMixture(
                    n_components=3, covariance_type='spherical').fit(self.nullratio)
    # def Ttest(self,sample):

        # testlist=[self.count[x]['Ratio'] for x in sample if self.count[x]['Ratio']!='ROH']
        # if len(self.nullratio)<10 or np.std(self.nullratio)==0:
            # return 'nan',"ROHregion"
        # elif len(testlist)==0:
            # return 'nan',"ROH"
        # elif len(testlist)==1:
            # stat=testlist[0]
            # if stat=="ROH":
                # return 'nan',"ROH"
            # else:
                # stat1=(stat-np.mean(self.nullratio))/(np.std(self.nullratio))
                # ans=stats.norm.cdf(stat1)
                # return 10**-stat,ans
        # else:
            # tstat,pvalue=stats.ttest_ind(testlist,self.nullratio)
            # mean=np.mean([10**-x for x in testlist])
            # if tstat<0:
                # return mean,pvalue
            # else:
                # return mean,1-pvalue
    # def Ttest(self,sample):

        # testlist=[self.count[x]['Ratio'] for x in sample if self.count[x]['Ratio']!='ROH']
        # if len(self.nullratio)<10 or max(self.nullratio)-min(self.nullratio)<0.0001:
            # return 'nan',"ROHregion"
        # elif len(testlist)==0:
            # return 'nan',"ROH"
        # elif len(testlist)>len(self.nullratio) or self.ns<10:
            # return 'nan',"Potential ROHregion or reference error"
        # elif len(testlist)==1:
            # stat=testlist[0]
            # if stat=="ROH":
                # return 'nan',"ROH"
            # else:

                # _,ans=stats.mannwhitneyu(testlist, self.nullratio, use_continuity=False,alternative='less')

                # return 10**-stat,ans
        # else:

            # _,pvalue=stats.mannwhitneyu(testlist, self.nullratio, use_continuity=False, alternative='less')
            # mean=np.mean([10**-x for x in testlist])
            # return mean,pvalue

    def Ttest(self, sample):

        testlist = [self.count[x]['Ratio']
                    for x in sample if self.count[x]['Ratio'] != 'ROH']
        if len(self.nullratio) <= 10 or max(self.nullratio) - min(self.nullratio) < 0.0001:
            return 'nan', "ROHregion"
        elif len(testlist) == 0:
            return 'nan', "ROH"
        elif len(testlist) > len(self.nullratio) or self.ns < 10:
            return 'nan', "Potential ROHregion or reference error"
        elif len(testlist) == 1:
            stat = testlist[0]
            if stat == "ROH":
                return 'nan', "ROH"
            else:
                # stat1=(stat-np.mean(self.nullratio))/(np.std(self.nullratio))
                # gmm=mixture.BayesianGaussianMixture(n_components=3, covariance_type='spherical').fit(a.reshape(-1,1))
                # ans=stats.norm.cdf(stat1)
                ans = self.gmm.score(np.array([stat]).reshape(-1, 1))
                return 10**-stat, ans
        else:
            # tstat,pvalue=stats.ttest_ind(testlist,self.nullratio)
            # _,pvalue=stats.mannwhitneyu(testlist, self.nullratio, use_continuity=False, alternative='less')
            ans = self.gmm.score(np.array(testlist).reshape(-1, 1))
            mean = np.mean([10**-x for x in testlist])
            return mean, ans
            # if tstat<0:
            #   return mean,pvalue
            # else:
            #   return mean,1-pvalue

    def stats(self, sample):
        nsnp = 0
        for x in sample:
            nsnp += self.count[x]['M']
        testlist = [self.count[x]['Ratio']
                    for x in sample if self.count[x]['Ratio'] != 'ROH']
        nsamplenullratio = len(self.nullratio)
        nonrohsample = len(testlist)
        nsample = len(sample)
        nnorm = len(self.count.keys()) - nsample
        return str(nsnp) + ',' + str(self.ns) + '\t' + str(nonrohsample) + ',' + str(nsample) + '\t' + \
            str(self.nullavg) + ',' + str(nsamplenullratio) + ',' + str(nnorm)
#       return Deltest(self.count[sample]['F'], self.count[sample]['M'], self.count[sample]['E'], self.length)
#       with open(self.regionfile,'r') as f:
#       for line in f:
#           dat=line.rstrip().split("\t")
#           hom=int(dat[3])
#           he=int(dat[4])
#           if dat[-1]==sample:
#               homo=hom
#               het=he
#       oddsratio, pvalue =stats.fisher_exact([[self.count[sample][0],s


class KS2sample:

    def __init__(self, obj, probands):
        self.obj = obj
        self.probands = probands
        self.controlst = []
        self.dct = {}
        # if not os.path.isfile(self.regionfile):
        # raise ValueError("file not found")
        if obj.shape[0] == 0:
            self.mean = ''
        for index, row in self.obj.iterrows():
            if row['sample'] not in probands:
                self.controlst.append(row['baf'])
            else:
                if row['sample'] not in self.dct.keys():
                    self.dct[row['sample']] = [row['baf']]
                else:
                    self.dct[row['sample']].append(row['baf'])

        # self.mean=np.mean(self.controlst)##
        # self.sd=np.std(self.controlst)##
        # print(self.mean,self.sd)

    def test(self, samples):
        testset = []
        for sample in samples:
            if sample in self.dct.keys():
                testset += self.dct[sample]
        if len(testset) < 1:
            return 'nan', "lowSNPs"
        elif len(self.controlst) < 1:
            return 'nan', "noBG"
        else:
            # testset=[(x-self.mean)/self.sd for x in testset]
            ks = stats.ks_2samp(testset, self.controlst)
            # ks=stats.kstest(testset,'norm')
            return ks
#############
# import sys
# [_,txt,het,chr,start,end,cnvid,sample,type]=sys.argv
# samplelst=sample.split(",")
# Del=DeletionTest(het,samplelst,int(end)-int(start))
# delp=Del.Ttest(samplelst)
# KS=KS2sample(txt,samplelst)
# ksp=KS.test(samplelst)
# stats=Del.stats(samplelst)
# print(chr+'\t'+start+'\t'+end+'\t'+cnvid+'\t'+sample+'\t'+type+'\t'+str(delp)+"\t"+str(ksp)+'\t'+stats)
