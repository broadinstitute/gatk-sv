# This assumes that RDTest splits are already merged
# This assume for each batch+source, there is one X file for each sex, and one Y file for each sex
# The usage is ```python rd-merge.py {batch} {source} {chrom} {inputfolder} {outputfolder}```
import argparse
import pandas as pd
parser = argparse.ArgumentParser("merge_allosomes.py")
parser.add_argument("batch", type=str, help="")
parser.add_argument("source", type=str, help="")
parser.add_argument("chrom", type=str, help="")
parser.add_argument("inputfolder", type=str, help="")
parser.add_argument("outputfolder", type=str, help="")
# inputfolder="split_rdtest/"
# outputfolder="rdtest/"

args = parser.parse_args()

inputmales = args.inputfolder + '/' + \
    '.'.join([args.batch, args.source, args.chrom, 'males.metrics'])
inputfemales = args.inputfolder + '/' + \
    '.'.join([args.batch, args.source, args.chrom, 'females.metrics'])
output = args.outputfolder + '/' + \
    '.'.join([args.batch, args.source, args.chrom, 'metrics'])

males = pd.read_table(inputmales)
females = pd.read_table(inputfemales)
if args.chrom == "X":
    if males.shape[0] != females.shape[0]:
        raise Exception('mismatched table sizes')
    male_only = females.P == 'No_samples_for_analysis'
    females.loc[male_only] = males
    females.to_csv(output, sep='\t', index=False)
elif args.chrom == "Y":
    if males.shape[0] != females.shape[0]:
        raise Exception('mismatched table sizes')
    males.to_csv(output, sep='\t', index=False)
else:
    raise Exception("Non-sex chromosome")
