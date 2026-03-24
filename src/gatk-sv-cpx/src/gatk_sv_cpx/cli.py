"""
``gatk-sv-cpx`` — unified CLI for complex SV resolution, evidence evaluation,
and genotype refinement.

Usage::

    gatk-sv-cpx resolve           -i in.vcf.gz -o resolved.vcf.gz [OPTIONS]
    gatk-sv-cpx evaluate-evidence -i resolved.vcf.gz -o evidence.json [OPTIONS]
    gatk-sv-cpx genotype-and-refine -i resolved.vcf.gz --evidence evidence.json -o final.vcf.gz [OPTIONS]
"""

from __future__ import annotations

import argparse
import logging
import sys


def _add_resolve_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "resolve",
        help="Resolve complex SV linkages, filter breakpoint overlaps, and rename variants.",
        description=(
            "Single-pass complex SV resolution replacing the scattered "
            "ResolveComplexVariants (RCV) workflow.  Performs breakpoint "
            "overlap filtering, post-CPX cleanup, and variant renaming."
        ),
    )
    p.add_argument(
        "-i", "--input-vcf", required=True,
        help="Clustered input VCF (bgzipped + tabix-indexed).",
    )
    p.add_argument(
        "-o", "--output-vcf", required=True,
        help="Path for resolved output VCF (.vcf.gz).",
    )
    p.add_argument(
        "--bothside-pass",
        help="SR bothside-pass list (one variant ID per line).",
    )
    p.add_argument(
        "--background-fail",
        help="SR background-fail list (one variant ID per line).",
    )
    p.add_argument(
        "--variant-prefix", default="CPX",
        help="Prefix for renamed variant IDs (default: CPX).",
    )
    p.set_defaults(func=_run_resolve)


def _run_resolve(args: argparse.Namespace) -> None:
    from gatk_sv_cpx.resolve import resolve
    resolve(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        bothside_pass=args.bothside_pass,
        background_fail=args.background_fail,
        variant_prefix=args.variant_prefix,
    )


def _add_evaluate_evidence_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "evaluate-evidence",
        help="Gather PE and RD evidence for CPX/CTX variants.",
        description=(
            "Replaces the fragmented PE and RD evidence gathering in "
            "GenotypeComplexVariants and RefineComplexVariants.  Reads a "
            "resolved VCF and evidence files, emits a JSON evidence bundle."
        ),
    )
    p.add_argument(
        "-i", "--input-vcf", required=True,
        help="Resolved VCF (output of 'resolve').",
    )
    p.add_argument(
        "-o", "--output-json", required=True,
        help="Path for the evidence JSON output.",
    )
    p.add_argument(
        "--pe-metrics",
        help="Concatenated/merged PE metrics file.",
    )
    p.add_argument(
        "--depth-support",
        help="Depth support BED file.",
    )
    p.add_argument(
        "--min-pe-cpx", type=int, default=3,
        help="Min PE count for CPX 'high_PE' (default: 3).",
    )
    p.add_argument(
        "--min-pe-ctx", type=int, default=3,
        help="Min PE count for CTX 'high_PE' (default: 3).",
    )
    p.add_argument(
        "--depth-threshold", type=float, default=0.5,
        help="Fraction threshold for depth support (default: 0.5).",
    )
    p.set_defaults(func=_run_evaluate_evidence)


def _run_evaluate_evidence(args: argparse.Namespace) -> None:
    from gatk_sv_cpx.evaluate_evidence import evaluate_evidence
    evaluate_evidence(
        input_vcf=args.input_vcf,
        output_json=args.output_json,
        pe_metrics=args.pe_metrics,
        depth_support=args.depth_support,
        min_pe_cpx=args.min_pe_cpx,
        min_pe_ctx=args.min_pe_ctx,
        depth_threshold=args.depth_threshold,
    )


def _add_genotype_and_refine_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "genotype-and-refine",
        help="Apply depth genotyping and PE/RD refinement.",
        description=(
            "Replaces GCV genotyping and RefCV VCF revision.  Reads the "
            "evidence JSON from 'evaluate-evidence', applies no-call rules, "
            "marks UNRESOLVED variants, and converts INS+INV to CPX."
        ),
    )
    p.add_argument(
        "-i", "--input-vcf", required=True,
        help="Resolved VCF (output of 'resolve').",
    )
    p.add_argument(
        "--evidence", required=True,
        help="Evidence JSON (output of 'evaluate-evidence').",
    )
    p.add_argument(
        "-o", "--output-vcf", required=True,
        help="Path for the finalized output VCF (.vcf.gz).",
    )
    p.add_argument(
        "--unresolved-ids",
        help="Optional file listing variant IDs to force-mark UNRESOLVED.",
    )
    p.set_defaults(func=_run_genotype_and_refine)


def _run_genotype_and_refine(args: argparse.Namespace) -> None:
    from gatk_sv_cpx.genotype_and_refine import genotype_and_refine
    genotype_and_refine(
        input_vcf=args.input_vcf,
        evidence_json=args.evidence,
        output_vcf=args.output_vcf,
        unresolved_ids_file=args.unresolved_ids,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Main entry point
# ═══════════════════════════════════════════════════════════════════════════

def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        prog="gatk-sv-cpx",
        description=(
            "Consolidated complex SV resolution, evidence evaluation, and "
            "genotype refinement for the GATK-SV pipeline."
        ),
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_resolve_parser(subparsers)
    _add_evaluate_evidence_parser(subparsers)
    _add_genotype_and_refine_parser(subparsers)

    args = parser.parse_args(argv)

    # Configure logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        stream=sys.stderr,
    )

    args.func(args)


if __name__ == "__main__":
    main()
