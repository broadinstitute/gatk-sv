import glob
import pandas as pd
import sys

suffix = 'AllMetrices'
file_name = sys.argv[1]

files = []
with open(file_name, "r") as fin:
    files = fin.readlines()
output_file = file_name + "_REViewer_AllMetrices.txt"

# files = glob.glob('*' + suffix)
df = pd.concat((pd.read_csv(f, sep='\t') for f in files))
na_free = df.dropna()
only_na = df[df.isna().any(axis=1)]
na_free.to_csv(output_file,sep='\t',index=None)
only_na.to_csv("missing_metrices.txt",sep='\t',index=None)
