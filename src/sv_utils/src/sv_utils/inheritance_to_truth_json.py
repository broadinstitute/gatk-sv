#!/usr/bin/env python

import sys
from pathlib import Path
import time
import argparse
import json

import pysam

from sv_utils import common, genomics_io, pedigree_tools, combine_truth_jsons
from sv_utils.get_truth_overlap import (
    ConfidentVariants, SampleConfidentVariants, ConfidentVariantsCombineStrategy
)
from typing import Optional
from collections.abc import Iterator


class VcfKeys:
    gt = genomics_io.VcfKeys.gt
    rd_cn = genomics_io.VcfKeys.rd_cn
    cn = genomics_io.VcfKeys.cn


class Default:
    inheritance_af_rareness = 0.05
    num_threads = common.num_physical_cpus
    use_copy_number = False
    use_cn = False


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Form truth JSON by finding genotypes with clear adherence / violation of "
                    "mendelian inheritance",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--vcf", '-v', type=Path, action="extend", nargs="+", required=True,
                        help="VCF with genotypes of unknown truth")
    parser.add_argument("--ped-file", "-p", type=Path, action="extend", nargs='*', required=True,
                        help="Pedigree file for samples.")
    parser.add_argument("--output", "-o", type=Path, required=True,
                        help="path to output JSON file")
    parser.add_argument(
        "--inheritance-af-rareness", type=float, default=Default.inheritance_af_rareness,
        help="Maximum allele frequency for a variant to use trio inheritance as a truth signal."
    )
    parser.add_argument(
        "--use-copy-number", type=common.argparse_bool, default=Default.use_copy_number,
        help="Where genotype is insufficient, use copy number for estimating allele frequency and "
             "carrier status"
    )
    parser.add_argument("--num_threads", "-@", type=int, default=Default.num_threads,
                        help="number of threads for compressing output vcf")
    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    if parsed_arguments.vcf is None:
        raise ValueError("Must supply one or more --vcf")
    if parsed_arguments.ped_file is None:
        raise ValueError("Must supply one or more --ped-file")
    return parsed_arguments


def main(argv: Optional[list[str]] = None) -> None:
    arguments = __parse_arguments(sys.argv if argv is None else argv)
    inheritance_to_truth_json(
        vcfs=arguments.vcf,
        pedigree_files=arguments.ped_file,
        output_truth_json=arguments.output,
        inheritance_af_rareness=arguments.inheritance_af_rareness,
        use_copy_number=arguments.use_copy_number,
        num_threads=arguments.num_threads
    )


def inheritance_to_truth_json(
        vcfs: list[Path],
        pedigree_files: list[Path],
        output_truth_json: Path,
        inheritance_af_rareness: float = Default.inheritance_af_rareness,
        use_copy_number: bool = Default.use_copy_number,
        num_threads: int = Default.num_threads
):
    confident_variants = inheritance_to_confident_variants(
        vcfs=vcfs, pedigree_files=pedigree_files, inheritance_af_rareness=inheritance_af_rareness,
        use_copy_number=use_copy_number, num_threads=num_threads
    )
    with open(output_truth_json, "w") as f_out:
        json.dump(confident_variants, f_out, default=lambda x: x.__dict__)


def inheritance_to_confident_variants(
        vcfs: list[Path],
        pedigree_files: list[Path],
        inheritance_af_rareness: float = Default.inheritance_af_rareness,
        use_copy_number: bool = Default.use_copy_number,
        num_threads: int = Default.num_threads
) -> ConfidentVariants:
    pedigree_file_info = pedigree_tools.PedigreeFileInfo.load(pedigree_files=pedigree_files)
    confident_variants: Iterator[ConfidentVariants] = (
        get_vcf_inheritance_confident_variants(
            vcf=vcf,
            pedigree_file_info=pedigree_file_info,
            inheritance_af_rareness=inheritance_af_rareness,
            use_copy_number=use_copy_number,
            num_threads=num_threads
        ) for vcf in vcfs
    )
    return combine_truth_jsons.combine_confident_variants(
        confident_variants,
        combine_strategy=ConfidentVariantsCombineStrategy.OmitConflicting
    )


def _compute_if_unknown(dask_object) -> int:
    return dask_object.compute() if hasattr(dask_object, "compute") else dask_object


def _dask_shape(dask_object) -> tuple[int, ...]:
    return tuple(
        _compute_if_unknown(s) for s in dask_object.shape
    )


def get_vcf_inheritance_confident_variants(
        vcf: Path,
        pedigree_file_info: pedigree_tools.PedigreeFileInfo,
        inheritance_af_rareness: float = Default.inheritance_af_rareness,
        use_copy_number: bool = Default.use_copy_number,
        use_cn: bool = Default.use_cn,
        num_threads: int = Default.num_threads
) -> ConfidentVariants:
    t0 = time.time()
    num_variants = 0
    with pysam.VariantFile(f"{vcf}", "r", threads=num_threads) as vcf_in:
        # subset trios to samples present in the VCF
        # noinspection PyTypeChecker
        pedigree_file_info = pedigree_file_info.subset_participants(
            vcf_in.header.samples, allow_unknown=False
        )
        # construct dict to hold allele counts for all samples in VCF. This may include non-trio
        # samples, which we want because we need to estimate cohort allele_frequency
        sample_acs = {sample_id: -1 for sample_id in vcf_in.header.samples}
        # construct ConfidentVariants object with no initial confident variants
        confident_variants = ConfidentVariants({
            sample_id: SampleConfidentVariants(good_variant_ids=tuple(), bad_variant_ids=tuple())
            for sample_id in pedigree_file_info.participant_ids
        })
        # iterate over VCF adding confident variants as they are discovered
        for vcf_line in vcf_in:
            num_variants += 1
            for sample_id, sample_format in vcf_line.samples.items():
                ac = genomics_io.genotype_to_allele_count(
                    sample_format.get(VcfKeys.gt), no_call_ac=-1
                )
                if use_copy_number and ac < 0:
                    # Genotype is no-call, set ac to 1 (HET) if copy-number is called and != 2
                    cn = sample_format.get(VcfKeys.rd_cn, -1)
                    if cn < 0 and use_cn:
                        cn = sample_format.get(VcfKeys.cn, -1)
                    ac = -1 if cn < 0 else 0 if cn == 2 else 1
                sample_acs[sample_id] = ac

            num_called_alleles = sum((2 for ac in sample_acs.values() if ac >= 0), start=0)
            if num_called_alleles == 0:
                continue  # this variant is not useful for inheritance
            called_allele_count = sum((ac for ac in sample_acs.values() if ac > 0), start=0)
            allele_frequency = called_allele_count / num_called_alleles
            if called_allele_count == 0 or allele_frequency > inheritance_af_rareness:
                continue  # this variant is not useful for inheritance

            variant_id = vcf_line.id
            for pedigree_line in pedigree_file_info.pedigree_lines:
                proband_ac = sample_acs[pedigree_line.proband_id]
                if proband_ac == 0:
                    continue  # this variant x trio is not useful for inheritance
                mother_ac = sample_acs[pedigree_line.mother_id]
                father_ac = sample_acs[pedigree_line.father_id]
                max_inherit_ac = (mother_ac > 0) + (father_ac > 0)
                if max_inherit_ac == 0:
                    # trio is denovo, consider all calls unreliable
                    confident_variants[pedigree_line.proband_id].append_bad(variant_id)
                    confident_variants[pedigree_line.mother_id].append_bad(variant_id)
                    confident_variants[pedigree_line.father_id].append_bad(variant_id)
                min_inherit_ac = mother_ac // 2 + father_ac // 2
                if min_inherit_ac <= proband_ac <= max_inherit_ac:
                    # trio is compatible with mendelian inheritance, consider all calls reliable
                    confident_variants[pedigree_line.proband_id].append_good(variant_id)
                    confident_variants[pedigree_line.mother_id].append_good(variant_id)
                    confident_variants[pedigree_line.father_id].append_good(variant_id)

    t1 = time.time()
    print(
        f"{vcf}: processed {num_variants} variants in "
        f"{common.elapsed_time(t1-t0, seconds_precision=0)}"
    )
    return confident_variants
