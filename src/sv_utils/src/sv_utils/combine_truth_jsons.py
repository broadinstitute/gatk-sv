#!/usr/bin/env python

import sys
from pathlib import Path
import argparse
import json
import dask.array
import dask.dataframe

from sv_utils.get_truth_overlap import (
    ConfidentVariants, ConfidentVariantsCombineStrategy
)
from typing import Optional
from collections.abc import Iterator

DaskDataFrame = dask.dataframe.DataFrame
DaskSeries = dask.dataframe.Series
DaskArray = dask.array.core.Array


class Default:
    combine_strategy = ConfidentVariantsCombineStrategy.OmitConflicting


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Form truth JSON by finding genotypes with clear adherence / violation of "
                    "mendelian inheritance",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--input-json", '-i', type=Path, action="extend", nargs="+", required=True,
                        help="Input truth JSON to be combined")
    parser.add_argument("--output-json", "-o", type=Path, required=True,
                        help="path to output combined JSON file")
    # noinspection PyTypeChecker
    parser.add_argument("--combine-strategy", "-c", type=ConfidentVariantsCombineStrategy,
                        choices=ConfidentVariantsCombineStrategy.choices,
                        default=Default.combine_strategy,
                        help="Rule for combining potentially-conflicting JSONs")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    return parsed_arguments


def main(argv: Optional[list[str]] = None) -> None:
    arguments = __parse_arguments(sys.argv if argv is None else argv)
    combine_truth_jsons(
        input_truth_jsons=arguments.input_json,
        output_truth_json=arguments.output_json,
        combine_strategy=arguments.combine_strategy
    )


def combine_truth_jsons(
        input_truth_jsons: list[Path],
        output_truth_json: Path,
        combine_strategy: ConfidentVariantsCombineStrategy = Default.combine_strategy
):
    output_confident_variants = combine_confident_variants(
        (load_confident_variants(input_truth_json) for input_truth_json in input_truth_jsons),
        combine_strategy=combine_strategy
    )
    with open(output_truth_json, 'w') as f_out:
        json.dump(output_confident_variants, f_out, indent="  ", default=lambda x: x.__dict__)


def load_confident_variants(truth_json: Path) -> ConfidentVariants:
    with open(truth_json, 'r') as f_in:
        return json.load(f_in, object_hook=ConfidentVariants.from_json)


def combine_confident_variants(
        confident_variants_iter: Iterator[ConfidentVariants],
        combine_strategy: ConfidentVariantsCombineStrategy
) -> ConfidentVariants:
    confident_variants = next(confident_variants_iter)
    for other_confident_variants in confident_variants_iter:
        confident_variants = confident_variants.combine(
            other_confident_variants, combine_strategy=combine_strategy
        )
    return confident_variants
