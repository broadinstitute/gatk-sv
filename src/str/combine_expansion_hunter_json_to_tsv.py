import argparse
import collections
import json
import logging
import os
import pathlib
import re

import numpy as np
import pandas as pd

from combine_json_to_tsv import get_sample_id_column_index

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
ALREADY_WARNED_ABOUT = set()  # used for logging

MISSING_LOCUS_RESULTS_WARNING_MESSAGE = \
    "Skipping {json_path} as it does not contain `LocusResults`. " \
    "One possible cause is running ExpansionHunter on male samples with default " \
    "value of the `--sex` argument, which defaults to female. " \
    "Rerunning EH with a correct value of `--sex` argument may fix this issue. " \
    "Another issue could be running ExpansionHunter on female samples with a variant " \
    "catalog that only contains sites on ChrY; or in more general sense, running " \
    "ExpansionHunter on a sample that contains no reads on any of the chromosomes " \
    "specified in the variant catalog (e.g., the sample only contains reads on Chr1 " \
    "while variant catalog contains sites on Chr2 only)."


def parse_args(args_list=None):
    """Parse command line args and return the argparse args object"""
    p = argparse.ArgumentParser()
    p.add_argument(
        "-c",
        "--variant-catalog",
        help="path of json variant catalog file. If specified, annotations for each locus in this variant catalog will "
             "be added to the output table as additional columns. This option can be specified more than once to "
             "use multiple variant catalogs as an annotation source.",
        action="append",
        default=[],
    )
    p.add_argument(
        "-m",
        "--sample-metadata",
        help="Table of sample annotations. If specified, all columns from this table will be added to the output table."
    )
    p.add_argument(
        "--sample-metadata-key",
        help="The column name in the --sample-metdata table that contains a sample id. "
             "If not specified, 'sample_id' and other variations like 'SampleId', 'sample', and 'ParticipantId' "
             "will be tried."
    )
    p.add_argument(
        "--include-all-fields",
        action="store_true",
        help="If specified, all fields from the ExpansionHunter will be added."
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
        help="ExpansionHunter output json path(s). If not specified, this script will retrieve all json files in the "
             "current directory and subdirectories",
        type=pathlib.Path,
        nargs="*"
    )

    args = p.parse_args(args=args_list)

    for variant_catalog in args.variant_catalog:
        if not os.path.isfile(variant_catalog):
            p.error(f"{variant_catalog} not found")

    if not args.json_paths:
        args.json_paths = [p for p in pathlib.Path(".").glob("**/*.json")]
        logging.info(f"Found {len(args.json_paths)} json files in .")

    return args


def main():
    """Main"""

    args = parse_args()

    output_prefix = args.output_prefix or "combined_expansion_hunter"
    output_prefix += f".{len(args.json_paths)}_json_files"

    combined_variant_catalog_contents = []
    if args.variant_catalog:
        for variant_catalog in args.variant_catalog:
            variant_catalog_contents = parse_json_file(variant_catalog)
            logging.info(f"Parsed {len(variant_catalog_contents)} records from {variant_catalog}")
            combined_variant_catalog_contents.extend(variant_catalog_contents)

    sample_metadata_lookup = collections.defaultdict(list)
    if args.sample_metadata:
        sample_metadata_df = pd.read_table(args.sample_metadata)
        sample_id_column_idx = get_sample_id_column_index(sample_metadata_df, column_name=args.sample_metadata_key)
        if sample_id_column_idx == -1:
            raise ValueError(f"'sample_id' column not found in sample metadata table. The columns found were: {sample_metadata_df.columns}")
        sample_id_column = sample_metadata_df.columns[sample_id_column_idx]
        for _, row in sample_metadata_df.iterrows():
            if args.verbose and len(sample_metadata_lookup[row[sample_id_column]]) > 0:
                logging.info(f"  {row[sample_id_column]} is a duplicate sample id in {args.sample_metadata}")
            sample_metadata_lookup[row[sample_id_column]].append(row)
        logging.info(f"Parsed {len(sample_metadata_df)} rows from {args.sample_metadata}")

    variant_table_columns = []
    allele_table_columns = []
    wrote_variant_table_header = wrote_allele_table_header = False
    variant_records_counter = allele_records_counter = 0
    variant_output_file = open(f"{output_prefix}_variants.tsv", "wt")
    allele_output_file = open(f"{output_prefix}_alleles.tsv", "wt")
    sample_metadata_lookup_counters = {}

    for just_get_header in True, False:
        if just_get_header:
            json_paths = args.json_paths[:20]
        else:
            json_paths = args.json_paths

        for json_path in json_paths:
            try:
                json_contents = parse_json_file(json_path)
            except Exception as e:
                logging.warning(f"Skipping {json_path}... Unable to parse json: {e}")
                continue

            if "LocusResults" not in json_contents:
                logging.warning(MISSING_LOCUS_RESULTS_WARNING_MESSAGE.format(json_path=json_path))
                continue

            if not isinstance(json_contents, dict) or "SampleParameters" not in json_contents:
                logging.info(f"Skipping {json_path}... Expected key 'SampleParameters' not found.")
                continue

            for record in convert_expansion_hunter_json_to_tsv_columns(
                json_contents,
                variant_catalog_contents=combined_variant_catalog_contents,
                sample_metadata_lookup=sample_metadata_lookup,
                sample_metadata_lookup_counters=sample_metadata_lookup_counters,
                json_file_path=json_path,
                return_allele_records=False,
                include_all_fields=args.include_all_fields,
            ):
                if just_get_header:
                    variant_table_columns.extend([k for k in record.keys() if k not in variant_table_columns])
                else:
                    if not wrote_variant_table_header:
                        variant_output_file.write("\t".join(variant_table_columns) + "\n")
                        wrote_variant_table_header = True
                    variant_records_counter += 1
                    variant_output_file.write("\t".join([str(record[c] if c in record else "") for c in variant_table_columns]) + "\n")

            for record in convert_expansion_hunter_json_to_tsv_columns(
                json_contents,
                variant_catalog_contents=combined_variant_catalog_contents,
                sample_metadata_lookup=sample_metadata_lookup,
                json_file_path=json_path,
                return_allele_records=True,
                include_all_fields=args.include_all_fields,
            ):
                if just_get_header:
                    allele_table_columns.extend([k for k in record.keys() if k not in allele_table_columns])
                else:
                    if not wrote_allele_table_header:
                        allele_output_file.write("\t".join(allele_table_columns) + "\n")
                        wrote_allele_table_header = True
                    allele_records_counter += 1
                    allele_output_file.write("\t".join([str(record[c] if c in record else "") for c in allele_table_columns]) + "\n")

    if sample_metadata_lookup_counters:
        logging.info(
            f"Found matches for {sample_metadata_lookup_counters['sample_ids_found']} out of "
            f"{sample_metadata_lookup_counters['total_sample_ids']} "
            f"({100*sample_metadata_lookup_counters['sample_ids_found']/sample_metadata_lookup_counters['total_sample_ids']:0.1f}%) "
            f"ExpansionHunter json sample ids in {args.sample_metadata}")

        if len(sample_metadata_lookup) != len(sample_metadata_df):
            logging.info(f"{len(sample_metadata_df) - len(sample_metadata_lookup)} duplicate sample ids in {args.sample_metadata}")
    logging.info(f"Wrote {variant_records_counter} records to {output_prefix}_variants.tsv")
    logging.info(f"Wrote {allele_records_counter} records to {output_prefix}_alleles.tsv")


class ParseError(Exception):
    """Represents an error that occurs while parsing"""
    pass


def parse_json_file(json_path):
    """Takes a json file path and loads its contents into a dictionary or list"""

    if not os.path.isfile(json_path):
        raise ValueError(f"{json_path} not found")

    with open(json_path, "rt", encoding="UTF-8") as f:
        try:
            return json.load(f)
        except Exception as e:
            raise ParseError(f"Unable to parse {json_path}: {e}")


def parse_read_count_tuples(tuple_string):
    """Parse the string value of one of the CountsOf* ExpansionHunter output fields
    (CountsOfSpanningReads, CountsOfFlankingReads, or CountsOfInrepeatReads). For example:

        "CountsOfSpanningReads": "(6, 1), (17, 10), (19, 16)"

    Return: a list of 2-tuples where t[0] = number of repeats, t[1] = number of supporting reads for this repeat size
    """
    json_string = tuple_string.replace("(", "[").replace(")", "]")

    return [tuple(t) for t in json.loads(f"[{json_string}]") if t]


def weighted_avg_and_std(values, weights):
    """Return the weighted average and standard deviation for a given array of values, and weights.
    Copied from: https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values - average) ** 2, weights=weights)
    return average, np.sqrt(variance)


def convert_expansion_hunter_json_to_tsv_columns(
    json_contents,
    variant_catalog_contents=None,
    sample_metadata_lookup=None,
    sample_metadata_lookup_counters=None,
    variant_info=None,
    json_file_path="",
    return_allele_records=True,
    include_all_fields=False,
):
    """Converts a dictionary that represents the contents of an ExpansionHunter v3 or v4 json output file to
    a dictionary of tsv column values.

    Args:
        json_contents (dict): a dict with the contents of an ExpansionHunter v3 or v4 json output file
        variant_catalog_contents (dict): optional dict with the contents of the variant catalog used when running
            expansion hunter. If provided, the fields will be added to the output table.
        sample_metadata_lookup (dict): maps sample id to row of sample metadata
        sample_metadata_lookup_counters (dict): dictionary for accumulating counters and stats about the
            ability to find samples in sample_metadata_lookup
        variant_info (dict): if provided, results will be added to this dict. Otherwise, a new dict will be created.
        json_file_path (str): if provided, it will be added as a field to the output table, and also used for logging
        return_allele_records (bool): if True, the returned list will have one record per allele rather than per variant
        include_all_fields (bool): if True, include all fields in the returned records, not just the basic ones.

    Returns:
        list: a list of tsv rows
    """

    if variant_info is None:
        variant_info = collections.OrderedDict()

    if variant_catalog_contents:
        variant_catalog_contents = {c['LocusId']: c for c in variant_catalog_contents}

    if json_file_path:
        variant_info["Dirname"] = os.path.dirname(json_file_path)
        variant_info["Filename"] = os.path.basename(json_file_path)

    try:
        variant_info["Sex"] = json_contents["SampleParameters"]["Sex"]
        variant_info["SampleId"] = json_contents["SampleParameters"]["SampleId"]
        locus_results_list = json_contents["LocusResults"].values()
    except KeyError as e:
        if "LocusResults" not in json_contents:
            logging.warning(MISSING_LOCUS_RESULTS_WARNING_MESSAGE.format(json_path=json_file_path))
            return []
        raise ValueError(f"Key {e} not found in json file {json_file_path}")

    if sample_metadata_lookup:
        sample_metadata_rows = sample_metadata_lookup.get(variant_info["SampleId"])
        if sample_metadata_rows is None:
            sample_id_from_filename = re.sub("([._]expansion_hunter[0-9]?)?.json.*$", "", variant_info["Filename"])
            sample_metadata_rows = sample_metadata_lookup.get(sample_id_from_filename)

        if sample_metadata_rows is None:
            if sample_metadata_lookup_counters is not None:
                sample_metadata_lookup_counters["total_sample_ids"] = sample_metadata_lookup_counters.get(
                    "total_sample_ids", 0) + 1
                sample_metadata_lookup_counters["sample_ids_found"] = sample_metadata_lookup_counters.get(
                    "sample_ids_found", 0)
        else:
            if sample_metadata_lookup_counters is not None:
                sample_metadata_lookup_counters["total_sample_ids"] = sample_metadata_lookup_counters.get(
                    "total_sample_ids", 0) + len(sample_metadata_rows)
                sample_metadata_lookup_counters["sample_ids_found"] = sample_metadata_lookup_counters.get(
                    "sample_ids_found", 0) + len(sample_metadata_rows)

            for sample_metadata_row in sample_metadata_rows:
                row_dict = sample_metadata_row.to_dict()
                for key, value in row_dict.items():
                    key = f"Sample_{key}"

                    if value is not None and not pd.isna(value):
                        output_value = str(value)
                    else:
                        output_value = ""

                    if key not in variant_info:
                        variant_info[key] = output_value
                    else:
                        variant_info[key] += f"; {output_value}"

    records_to_return = []
    for locus_json in locus_results_list:
        locus_record = collections.OrderedDict(variant_info)
        for key in "LocusId", "ReadLength", "Coverage", "AlleleCount":
            locus_record[key] = locus_json[key]

        # if variant catalog specified, transfer info from it
        locus_id = locus_record["LocusId"]
        if variant_catalog_contents:
            locus_record["InVariantCatalog"] = locus_id in variant_catalog_contents
            if locus_record["InVariantCatalog"]:
                excluded_keys = "LocusId", "ReferenceRegion", "VariantId", "VariantType", "OfftargetRegions"
                locus_record.update({f"VariantCatalog_{k}": v for k, v in variant_catalog_contents[locus_id].items() if k not in excluded_keys})
                locus_record["VariantCatalog_NumOfftargetRegions"] = len(variant_catalog_contents[locus_id].get("OfftargetRegions", []))
            else:
                if locus_id not in ALREADY_WARNED_ABOUT:
                    ALREADY_WARNED_ABOUT.add(locus_id)
                    logging.warning(f"{json_file_path} locus id {locus_id} not found in variant catalog.")

        for variant_json in locus_json.get("Variants", {}).values():
            if "Genotype" not in variant_json:
                continue

            # docs @ https://github.com/Illumina/ExpansionHunter/blob/master/docs/05_OutputJsonFiles.md
            variant_record = collections.OrderedDict(locus_record)
            variant_record["VariantId"] = variant_json.get("VariantId", "")
            variant_record["VariantSubtype"] = variant_json.get("VariantSubtype", "")
            variant_record["RepeatUnit"] = variant_json["RepeatUnit"]
            variant_record["RepeatUnitLength"] = len(variant_json["RepeatUnit"])
            variant_record["ReferenceRegion"] = variant_json["ReferenceRegion"]

            if include_all_fields:
                variant_record["CountsOfSpanningReads"] = variant_json["CountsOfSpanningReads"]
                variant_record["CountsOfFlankingReads"] = variant_json["CountsOfFlankingReads"]
                variant_record["CountsOfInrepeatReads"] = variant_json["CountsOfInrepeatReads"]

                spanning_read_tuples = parse_read_count_tuples(variant_json["CountsOfSpanningReads"])
                flanking_read_tuples = parse_read_count_tuples(variant_json["CountsOfFlankingReads"])
                inrepeat_read_tuples = parse_read_count_tuples(variant_json["CountsOfInrepeatReads"])

                variant_record["NumSpanningReads"] = sum(t[1] for t in spanning_read_tuples)
                variant_record["NumFlankingReads"] = sum(t[1] for t in flanking_read_tuples)
                variant_record["NumInrepeatReads"] = sum(t[1] for t in inrepeat_read_tuples)
                variant_record["NumReadsTotal"] = sum([variant_record[k] for k in ("NumSpanningReads", "NumFlankingReads", "NumInrepeatReads")])

                variant_record["NumAllelesSupportedBySpanningReads"] = len(spanning_read_tuples)
                variant_record["NumAllelesSupportedByFlankingReads"] = len(flanking_read_tuples)
                variant_record["NumAllelesSupportedByInrepeatReads"] = len(inrepeat_read_tuples)
                variant_record["NumAllelesSupportedTotal"] = sum([variant_record[k] for k in ("NumAllelesSupportedBySpanningReads", "NumAllelesSupportedByFlankingReads", "NumAllelesSupportedByInrepeatReads")])

            variant_record["Genotype"] = variant_json["Genotype"]
            variant_record["GenotypeConfidenceInterval"] = variant_json["GenotypeConfidenceInterval"]
            genotype_tuples = list(zip(
                variant_json["Genotype"].split("/"),
                variant_json["GenotypeConfidenceInterval"].split("/"),
            ))

            for i, (genotype, genotypeCI) in enumerate(genotype_tuples):
                if return_allele_records:
                    suffix = ""
                    allele_record = collections.OrderedDict(variant_record)
                else:
                    suffix = f": Allele {i+1}"
                    allele_record = variant_record

                confidence_interval_start, confidence_interval_end = genotypeCI.split("-")
                allele_record[f"Allele Number{suffix}"] = i + 1
                allele_record[f"Num Repeats{suffix}"] = int(genotype)
                allele_record[f"Repeat Size (bp){suffix}"] = int(genotype) * len(variant_json.get("RepeatUnit", ""))
                allele_record[f"CI start{suffix}"] = int(confidence_interval_start)
                allele_record[f"CI end{suffix}"] = int(confidence_interval_end)
                allele_record[f"CI size{suffix}"] = int(confidence_interval_end) - int(confidence_interval_start)

                if include_all_fields:
                    is_homozygous = len(genotype_tuples) > 1 and genotype_tuples[0][0] == genotype_tuples[1][0]
                    divisor = 2 if is_homozygous else 1
                    allele_record[f"NumSpanningReadsThatSupportGenotype{suffix}"] = sum(t[1] for t in spanning_read_tuples if t[0] == int(genotype)) / divisor
                    allele_record[f"NumFlankingReadsThatSupportGenotype{suffix}"] = sum(t[1] for t in flanking_read_tuples if t[0] == int(genotype)) / divisor
                    allele_record[f"NumInrepeatReadsThatSupportGenotype{suffix}"] = sum(t[1] for t in inrepeat_read_tuples if t[0] == int(genotype)) / divisor
                    allele_record[f"NumReadsTotalThatSupportGenotype{suffix}"] = (
                        allele_record[f"NumSpanningReadsThatSupportGenotype{suffix}"] +
                        allele_record[f"NumFlankingReadsThatSupportGenotype{suffix}"] +
                        allele_record[f"NumInrepeatReadsThatSupportGenotype{suffix}"]
                    )
                    allele_record[f"FractionOfReadsThatSupportsGenotype{suffix}"] = (
                        allele_record[f"NumReadsTotalThatSupportGenotype{suffix}"] / float(allele_record["NumReadsTotal"]) if int(allele_record["NumReadsTotal"]) > 0 else 0
                    )

                if return_allele_records:
                    records_to_return.append(allele_record)

            if not return_allele_records:
                records_to_return.append(allele_record)

    return records_to_return


if __name__ == "__main__":
    main()
