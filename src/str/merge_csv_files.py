import argparse
import pandas as pd


def merge_files(input_filename: str, metrics_filename: str, missing_metrics_filename: str):
    with open(input_filename, "r") as fin:
        files = [line.rstrip() for line in fin]

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
        na_free.to_csv(metrics_filename, sep='\t', index=None)
        only_na.to_csv(missing_metrics_filename, sep='\t', index=None)
    else:
        open(metrics_filename, "w").close()
        open(missing_metrics_filename, "w").close()


def main():
    parser = argparse.ArgumentParser(
        description="Merges the CSV files listed in a text file into two files: only non-NaN, only NaN.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-i", "--input-filename",
        help="A text file containing the list of CSV files to merge, with one file per line."
    )

    parser.add_argument(
        "-o", "--metrics",
        help="The output filename containing all the metrics (excluding the NaN values)."
    )

    parser.add_argument(
        "-m", "--missing-metrics",
        help="The output filename containing only the missing metrics."
    )

    args = parser.parse_args()
    merge_files(args.input_filename, args.metrics, args.missing_metrics)


if __name__ == '__main__':
    main()
