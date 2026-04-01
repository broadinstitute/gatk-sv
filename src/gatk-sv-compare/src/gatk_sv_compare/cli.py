"""CLI entrypoint for gatk-sv-compare."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from .config import AnalysisConfig
from .preprocess import parse_reference_dict, read_contig_list, run_preprocess
from .validate import validate_and_render


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser."""
    parser = argparse.ArgumentParser(prog="gatk-sv-compare")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser("validate", help="Check VCF format")
    validate_parser.add_argument("--vcf", required=True, type=Path, help="Path to the input VCF")
    validate_parser.set_defaults(handler=_handle_validate)

    preprocess_parser = subparsers.add_parser("preprocess", help="Run SVConcordance + SVRegionOverlap")
    preprocess_parser.add_argument("--vcf-a", required=True, type=Path)
    preprocess_parser.add_argument("--vcf-b", required=True, type=Path)
    preprocess_parser.add_argument("--reference-dict", required=True, type=Path)
    preprocess_parser.add_argument("--contig-list", required=True, type=Path)
    preprocess_parser.add_argument("--output-dir", required=True, type=Path)
    preprocess_parser.add_argument("--seg-dup-track", type=Path)
    preprocess_parser.add_argument("--simple-repeat-track", type=Path)
    preprocess_parser.add_argument("--repeatmasker-track", type=Path)
    preprocess_parser.add_argument("--gatk-path", default="gatk")
    preprocess_parser.add_argument("--java-options", default="-Xmx4g")
    preprocess_parser.add_argument("--num-workers", type=int, default=1)
    preprocess_parser.set_defaults(handler=_handle_preprocess)

    analyze_parser = subparsers.add_parser("analyze", help="Run analysis modules on preprocessed VCFs")
    analyze_parser.add_argument("--vcf-a", required=True, type=Path)
    analyze_parser.add_argument("--vcf-b", required=True, type=Path)
    analyze_parser.add_argument("--output-dir", required=True, type=Path)
    analyze_parser.set_defaults(handler=_handle_not_implemented)

    run_parser = subparsers.add_parser("run", help="Run preprocess + analyze end-to-end")
    run_parser.add_argument("--vcf-a", required=True, type=Path)
    run_parser.add_argument("--vcf-b", required=True, type=Path)
    run_parser.add_argument("--reference-dict", required=True, type=Path)
    run_parser.add_argument("--contig-list", required=True, type=Path)
    run_parser.add_argument("--output-dir", required=True, type=Path)
    run_parser.set_defaults(handler=_handle_not_implemented)

    return parser


def _handle_validate(args: argparse.Namespace) -> int:
    summary, output = validate_and_render(args.vcf)
    print(output)
    return 1 if summary.has_errors else 0


def _handle_preprocess(args: argparse.Namespace) -> int:
    contigs = read_contig_list(args.contig_list)
    contig_lengths = parse_reference_dict(args.reference_dict)
    config = AnalysisConfig(
        vcf_a_path=args.vcf_a,
        vcf_b_path=args.vcf_b,
        output_dir=args.output_dir,
        reference_dict=args.reference_dict,
        contigs=contigs,
        contig_lengths=contig_lengths,
        n_workers=args.num_workers,
        seg_dup_track=args.seg_dup_track,
        simple_repeat_track=args.simple_repeat_track,
        repeatmasker_track=args.repeatmasker_track,
        gatk_path=args.gatk_path,
        java_options=args.java_options,
    )
    annotated_a, annotated_b = run_preprocess(config)
    print(f"annotated_a={annotated_a}")
    print(f"annotated_b={annotated_b}")
    return 0


def _handle_not_implemented(_: argparse.Namespace) -> int:
    raise NotImplementedError("This Phase 1 CLI command is scaffolded but not implemented yet.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
