import pandas as pd
import sys

file_name = sys.argv[1]

files = []
with open(file_name, "r") as fin:
    files = [line.rstrip() for line in fin]
output_file = "reviewer_all_metrics.txt"
missing_metrics_filename = "missing_metrics.txt"

dfs = []
for filename in files:
    try:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        if not df.empty:
            dfs.append(df)
    except pd.errors.EmptyDataError:
        print(f"Empty file: {filename}")

if len(dfs) > 0:
    df = pd.concat(dfs, axis=0, ignore_index=True)
    na_free = df.dropna()
    only_na = df[df.isna().any(axis=1)]
    na_free.to_csv(output_file,sep='\t',index=None)
    only_na.to_csv(missing_metrics_filename,sep='\t',index=None)
else:
    open(output_file, "w").close()
    open(missing_metrics_filename, "w").close()
