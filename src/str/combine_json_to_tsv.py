"""Combines json files into a single .tsv table for analysis."""

import argparse
import json
import logging
import os
import pathlib
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def parse_args(args_list=None):
    """Parse command-line args and return the ArgumentParser args object.

    Args:
        arg_list (list): optional artificial list of command-line args to use for testing.

    Returns:
        args
    """

    p = argparse.ArgumentParser()
    p.add_argument(
        "-d",
        "--add-dirname-column",
        action="store_true",
        help="Add Dirname column containing the relative path of directory containing the json file."
    )
    p.add_argument(
        "-f",
        "--add-filename-column",
        action="store_true",
        help="Add Filename column containing the json filename."
    )
    p.add_argument(
        "-m",
        "--sample-metadata",
        help="Table of sample annotations. If specified, all columns from this table will be added to the output table."
    )
    p.add_argument(
        "--json-sample-id-key",
        help="The json field that contains a sample id to use for joining with the --sample-metadata table. "
             "If not specified, 'sample_id' and other variations like 'SampleId', 'sample', and 'ParticipantId' "
             "will be tried."
    )
    p.add_argument(
        "--sample-metadata-key",
        help="The column name in the --sample-metdata table that contains a sample id. "
             "If not specified, 'sample_id' and other variations like 'SampleId', 'sample', and 'ParticipantId' "
             "will be tried."
    )
    p.add_argument(
        "-o",
        "--output-prefix",
        help="Combined table output filename prefix",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print additional logging messages",
    )
    p.add_argument(
        "json_paths",
        help="json path(s). If not specified, this script will retrieve all json files in the current directory and subdirectories",
        type=pathlib.Path,
        nargs="*"
    )
    args = p.parse_args(args=args_list)

    if args.json_paths:
        # check if files exist
        for json_path in args.json_paths:
            if not os.path.isfile(json_path):
                p.error(f"File not found: {json_path}")
    else:
        # find all .json files underneath the current directory
        args.json_paths = [p for p in pathlib.Path(".").glob("**/*.json")]
        if args.verbose and len(args.json_paths) > 0:
            print("Found json files: ")
            for p in args.json_paths:
                print(f"  {p}")
        print(f"Found {len(args.json_paths)} .json files under {os.getcwd()}")
        if len(args.json_paths) == 0:
            sys.exit(0)

    if args.sample_metadata_key and not args.sample_metadata:
        p.error("--sample-metadata-key should only be specified along with --sample-metadata")

    if args.json_sample_id_key and not args.sample_metadata:
        p.error("--json-sample-id-key should only be specified along with --sample-metadata")

    return args


SAMPLE_ID_COLUMN_ALAISES = {
    "sampleid",
    "sample_id",
    "sample",
    "participantid",
}


def get_sample_id_column_index(df, column_name=None):
    """Try to find a column in df that contains sample ids

    Args:
         df (DataFrame): pandas DataFrame
         column_name (str): if specified this function will get the index of this specific column name

    Return:
        Index of sample id column, or -1 if not found
    """
    for i, column in enumerate(df.columns):
        if column_name is not None:
            if column == column_name:
                return i
        elif column.lower() in SAMPLE_ID_COLUMN_ALAISES:
            return i

    return -1


class ParseError(Exception):
    pass


def parse_json_files(json_paths, add_dirname_column=False, add_filename_column=False):
    """Takes json file paths and yields the contents of each one as a dictionary or list"""
    for json_path in json_paths:
        if not os.path.isfile(json_path):
            raise ValueError(f"{json_path} not found")

        with open(json_path, "rt", encoding="UTF-8") as f:
            try:
                json_contents = json.load(f)
            except Exception as e:
                raise ParseError(f"Unable to parse {json_path}: {e}")

            if isinstance(json_contents, dict):
                if add_dirname_column:
                    json_contents["Dirname"] = os.path.dirname(json_path)
                if add_filename_column:
                    json_contents["Filename"] = os.path.basename(json_path)

                json_contents_excluding_complex_values = {
                    key: value for key, value in sorted(json_contents.items()) if isinstance(value, (int, str, bool, float, tuple))
                }

                yield json_contents_excluding_complex_values


def join_with_sample_metadata(df, df_sample_id_columns, sample_metadata_df, sample_metadata_df_sample_id_column, verbose=False):
    """Performs a LEFT join between df and sample_metadata_df.

    Args:
        df (pandas.DataFrame): DataFrame representing the contents of the .json files.
        df_sample_id_columns (str): Name of the column in df that should be used as the join key.
        sample_metadata_df (pandas.DataFrame): DataFrame representing additional sample-level metadata.
        sample_metadata_df_sample_id_column (str): Name of the column in sample_metadata_df that should be used as the join key.
        verbose (bool): Whether to print additional log statements.

    Return:
        pandas.DataFrame: The DataFrame resulting from a LEFT join between df and sample_metadata_df
    """
    sample_id_column_idx1 = get_sample_id_column_index(df, column_name=df_sample_id_columns)
    if sample_id_column_idx1 == -1:
        raise ValueError(f"'sample_id' field not found in json files. The fields found were: {df.columns}")

    sample_id_column_idx2 = get_sample_id_column_index(sample_metadata_df, column_name=sample_metadata_df_sample_id_column)
    if sample_id_column_idx2 == -1:
        raise ValueError(f"'sample_id' column not found in sample metadata table. The columns found were: {sample_metadata_df.columns}")

    sample_id_column1 = df.columns[sample_id_column_idx1]
    sample_id_column2 = sample_metadata_df.columns[sample_id_column_idx2]
    print(f"Doing a LEFT JOIN with sample metadata table using keys {sample_id_column1} and {sample_id_column2}")

    set1 = set(df[sample_id_column1])
    set2 = set(sample_metadata_df[sample_id_column2])
    shared_sample_ids = set1 & set2
    df_unique_sample_ids = sorted(set1 - set2)
    sample_metadata_df_unique_sample_ids = sorted(set2 - set1)

    print(f"{len(shared_sample_ids)} out of {len(df)} ({100*len(shared_sample_ids)/len(df):0.1f}%) sample ids "
          f"in the json files have a matching sample id in the sample metadata table")

    if len(df_unique_sample_ids) > 0 and len(df_unique_sample_ids) < 100:
        print("Sample ids found only in the json files: ", ", ".join(df_unique_sample_ids))

    if verbose:
        print(f"{len(shared_sample_ids)} out of {len(sample_metadata_df)} "
              f"({100*len(shared_sample_ids)/len(sample_metadata_df):0.1f}%) sample ids in the sample metadata table "
              f"have a matching sample id in the json files")

        if len(sample_metadata_df_unique_sample_ids) > 0 and len(sample_metadata_df_unique_sample_ids) < 100:
            print("Sample ids found only in the sample metadata table: ", ", ".join(sample_metadata_df_unique_sample_ids))

    if len(set2) < len(sample_metadata_df):
        print("WARNING: the sample metadata table has", len(sample_metadata_df) - len(set2), "duplicate sample id(s)")
        sample_ids = sorted(list(sample_metadata_df[sample_id_column2]))
        print(", ".join([str(i) for i in sample_ids if sample_ids.count(i) > 1]))

    sample_metadata_df = sample_metadata_df.rename(columns={
        c: f"Sample_{c}" for c in sample_metadata_df.columns if c != sample_id_column2
    })
    df = pd.merge(df, sample_metadata_df, how="left", left_on=sample_id_column1, right_on=sample_id_column2)
    return df


def main():
    """Combine json to tsv"""

    args = parse_args()

    output_prefix = args.output_prefix or "combined"
    output_prefix += f".{len(args.json_paths)}_json_files"

    df = pd.DataFrame(parse_json_files(
        args.json_paths,
        add_dirname_column=args.add_dirname_column,
        add_filename_column=args.add_filename_column))

    if args.sample_metadata:
        sample_metadata_df = pd.read_table(args.sample_metadata)
        df = join_with_sample_metadata(df, args.json_sample_id_key, sample_metadata_df, args.sample_metadata_key, args.verbose)

    output_filename = f"{output_prefix}.tsv"
    df.to_csv(output_filename, index=False, header=True, sep="\t")
    print(f"Wrote {len(df)} rows to {output_filename}")


if __name__ == "__main__":
    main()
