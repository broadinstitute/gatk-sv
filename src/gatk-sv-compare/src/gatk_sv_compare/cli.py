"""CLI entrypoint for gatk-sv-compare."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

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
    preprocess_parser.set_defaults(handler=_handle_not_implemented)

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


def _handle_not_implemented(_: argparse.Namespace) -> int:
    raise NotImplementedError("This Phase 1 CLI command is scaffolded but not implemented yet.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
