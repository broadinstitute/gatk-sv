"""CLI entrypoint for gatk-sv-compare."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Type

import pysam

from .aggregate import aggregate
from .config import AnalysisConfig
from .modules import ALL_MODULES, AnalysisModule
from .preprocess import parse_reference_dict, read_contig_list, run_preprocess
from .validate import fix_and_render, validate_and_render


def _configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s", force=True)


def _resolve_num_workers(requested_workers: Optional[int], contig_count: int) -> int:
    cpu_count = os.cpu_count() or 1
    if requested_workers is None:
        return max(1, min(contig_count, cpu_count, 4))
    return max(1, min(requested_workers, contig_count))


def _context_overlap_value(raw_value: str) -> float:
    value = float(raw_value)
    if not 0.0 <= value <= 1.0:
        raise argparse.ArgumentTypeError("--context-overlap must be between 0 and 1")
    return value


def _default_fix_output_path(vcf_path: Path) -> Path:
    if vcf_path.name.endswith(".vcf.gz"):
        return vcf_path.with_name(vcf_path.name[:-7] + ".fixed.vcf.gz")
    if vcf_path.suffix == ".vcf":
        return vcf_path.with_name(vcf_path.stem + ".fixed.vcf")
    return vcf_path.with_name(vcf_path.name + ".fixed.vcf")


def _parse_module_names(raw_value: Optional[str]) -> Optional[List[str]]:
    if raw_value is None:
        return None
    module_names = [item.strip() for item in raw_value.split(",") if item.strip()]
    return list(dict.fromkeys(module_names))


def _module_registry() -> Dict[str, Type[AnalysisModule]]:
    return {module_type().name: module_type for module_type in ALL_MODULES}


def _resolve_module_types(requested_names: Optional[Iterable[str]]) -> List[Type[AnalysisModule]]:
    registry = _module_registry()
    if requested_names is None:
        return list(ALL_MODULES)
    requested_names_list = list(requested_names)
    unknown = [name for name in requested_names_list if name not in registry]
    if unknown:
        available = ", ".join(sorted(registry))
        raise ValueError(f"Unknown module(s): {', '.join(unknown)}. Available modules: {available}")
    return [registry[name] for name in requested_names_list]


def _read_vcf_contigs(vcf_path: Path) -> Tuple[List[str], Dict[str, int]]:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        contigs = list(vcf.header.contigs)
        lengths = {name: int(vcf.header.contigs[name].length or 0) for name in contigs}
    return contigs, lengths


def _infer_common_contigs(vcf_a_path: Path, vcf_b_path: Path) -> Tuple[List[str], Dict[str, int]]:
    contigs_a, lengths_a = _read_vcf_contigs(vcf_a_path)
    contigs_b, lengths_b = _read_vcf_contigs(vcf_b_path)
    contigs_b_set = set(contigs_b)
    common_contigs = [contig for contig in contigs_a if contig in contigs_b_set]
    if not common_contigs:
        raise ValueError(f"No shared contigs found between {vcf_a_path} and {vcf_b_path}")
    contig_lengths = {
        contig: lengths_a.get(contig) or lengths_b.get(contig) or 0 for contig in common_contigs
    }
    return common_contigs, contig_lengths


def _resolve_requested_preprocess_contigs(
    contig_list_path: Path,
    reference_dict_path: Path,
    requested_contig: Optional[str] = None,
) -> Tuple[List[str], Dict[str, int]]:
    contigs = read_contig_list(contig_list_path)
    contig_lengths = parse_reference_dict(reference_dict_path)
    if requested_contig is None:
        return contigs, contig_lengths
    if requested_contig not in contigs:
        raise ValueError(f"Requested contig {requested_contig} is not present in {contig_list_path}")
    if requested_contig not in contig_lengths:
        raise ValueError(f"Requested contig {requested_contig} is not present in {reference_dict_path}")
    return [requested_contig], {requested_contig: contig_lengths[requested_contig]}


def _build_analysis_config(
    *,
    vcf_a_path: Path,
    vcf_b_path: Path,
    output_dir: Path,
    label_a: str,
    label_b: str,
    module_names: Optional[List[str]],
    pass_only: bool,
    context_overlap: float,
    per_chrom: bool,
    enable_site_match_table: bool,
    per_sample_counts_table: bool,
    ped_file: Optional[Path],
    requested_workers: Optional[int],
    contigs: Optional[List[str]] = None,
    contig_lengths: Optional[Dict[str, int]] = None,
) -> AnalysisConfig:
    resolved_contigs = list(contigs) if contigs is not None else None
    resolved_lengths = dict(contig_lengths) if contig_lengths is not None else None
    if resolved_contigs is None or resolved_lengths is None:
        resolved_contigs, resolved_lengths = _infer_common_contigs(vcf_a_path, vcf_b_path)
    resolved_workers = _resolve_num_workers(requested_workers, len(resolved_contigs))
    return AnalysisConfig(
        vcf_a_path=vcf_a_path,
        vcf_b_path=vcf_b_path,
        vcf_a_label=label_a,
        vcf_b_label=label_b,
        output_dir=output_dir,
        contigs=resolved_contigs,
        contig_lengths=resolved_lengths,
        n_workers=resolved_workers,
        modules=module_names,
        pass_only=pass_only,
        context_overlap=context_overlap,
        per_chrom=per_chrom,
        enable_site_match_table=enable_site_match_table,
        per_sample_counts_table=per_sample_counts_table,
        ped_file=ped_file,
    )


def _skip_reason(module: AnalysisModule, data, config: AnalysisConfig) -> Optional[str]:
    if module.requires_shared_samples and not data.shared_samples:
        return "no shared samples"
    if module.requires_ped_file and config.ped_file is None:
        return "no --ped provided"
    return None


def _run_analysis(config: AnalysisConfig) -> int:
    logger = logging.getLogger(__name__)
    module_types = _resolve_module_types(config.modules)
    logger.info(
        "Starting analysis: %s vs %s across %s contig(s) with %s worker(s)",
        config.vcf_a_label,
        config.vcf_b_label,
        len(config.contigs),
        config.n_workers,
    )
    logger.info("Selected modules: %s", ", ".join(module_type().name for module_type in module_types))
    data = aggregate(config)
    logger.info(
        "Aggregation complete: %s sites in %s, %s sites in %s, %s matched pair(s), %s shared sample(s)",
        len(data.sites_a),
        data.label_a,
        len(data.sites_b),
        data.label_b,
        len(data.matched_pairs),
        len(data.shared_samples),
    )

    ran_modules: List[str] = []
    skipped_modules: List[str] = []
    for module_type in module_types:
        module = module_type()
        reason = _skip_reason(module, data, config)
        if reason is not None:
            logger.info("Skipping module %s: %s", module.name, reason)
            skipped_modules.append(module.name)
            continue
        logger.info("Running module %s", module.name)
        module.run(data, config)
        logger.info("Completed module %s", module.name)
        ran_modules.append(module.name)

    print(f"analysis_output_dir={config.output_dir}")
    print(f"modules_ran={','.join(ran_modules)}")
    if skipped_modules:
        print(f"modules_skipped={','.join(skipped_modules)}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser."""
    parser = argparse.ArgumentParser(prog="gatk-sv-compare")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser("validate", help="Check VCF format")
    validate_parser.add_argument("--vcf", required=True, type=Path, help="Path to the input VCF")
    validate_parser.add_argument("--fix", action="store_true", help="Write an automatically corrected VCF when only fixable issues are present")
    validate_parser.add_argument("--out", type=Path, help="Output VCF path for --fix mode")
    validate_parser.add_argument("--ploidy-table", type=Path, help="Tab-delimited sample ploidy table used to repair missing ECN fields in --fix mode")
    validate_parser.add_argument("--drop-bnd", action="store_true", help="In --fix mode, drop all BND records from the VCF")
    validate_parser.add_argument("--drop-ctx", action="store_true", help="In --fix mode, drop all CTX records from the VCF")
    validate_parser.set_defaults(handler=_handle_validate)

    preprocess_parser = subparsers.add_parser("preprocess", help="Run SVConcordance + SVRegionOverlap")
    preprocess_parser.add_argument("--vcf-a", required=True, type=Path)
    preprocess_parser.add_argument("--vcf-b", required=True, type=Path)
    preprocess_parser.add_argument("--reference-dict", required=True, type=Path)
    preprocess_parser.add_argument("--contig-list", required=True, type=Path)
    preprocess_parser.add_argument("--contig", help="Restrict preprocessing to a single contig from --contig-list")
    preprocess_parser.add_argument("--output-dir", required=True, type=Path)
    preprocess_parser.add_argument("--seg-dup-track", type=Path)
    preprocess_parser.add_argument("--simple-repeat-track", type=Path)
    preprocess_parser.add_argument("--repeatmasker-track", type=Path)
    preprocess_parser.add_argument("--gatk-path", default="gatk")
    preprocess_parser.add_argument("--java-options", default="-Xmx4g")
    preprocess_parser.add_argument(
        "--num-workers",
        type=int,
        default=None,
        help="Number of contig shards to run in parallel (default: auto, capped at 4 and the contig count)",
    )
    preprocess_parser.set_defaults(handler=_handle_preprocess)

    analyze_parser = subparsers.add_parser("analyze", help="Run analysis modules on preprocessed VCFs")
    analyze_parser.add_argument("--vcf-a", required=True, type=Path)
    analyze_parser.add_argument("--vcf-b", required=True, type=Path)
    analyze_parser.add_argument("--label-a", default="VCF_A")
    analyze_parser.add_argument("--label-b", default="VCF_B")
    analyze_parser.add_argument("--output-dir", required=True, type=Path)
    analyze_parser.add_argument("--modules", help="Comma-separated module list to run")
    analyze_parser.add_argument("--pass-only", action="store_true")
    analyze_parser.add_argument(
        "--context-overlap",
        type=_context_overlap_value,
        default=0.5,
        help="Minimum overlap fraction required to assign span-based genomic contexts (default: 0.5)",
    )
    analyze_parser.add_argument("--per-chrom", action="store_true")
    analyze_parser.add_argument("--enable-site-match-table", action="store_true")
    analyze_parser.add_argument("--per-sample-counts-table", action="store_true")
    analyze_parser.add_argument("--ped", dest="ped_file", type=Path)
    analyze_parser.add_argument(
        "--num-workers",
        type=int,
        default=None,
        help="Number of contigs to aggregate in parallel (default: auto, capped at 4 and the contig count)",
    )
    analyze_parser.set_defaults(handler=_handle_analyze)

    run_parser = subparsers.add_parser("run", help="Run preprocess + analyze end-to-end")
    run_parser.add_argument("--vcf-a", required=True, type=Path)
    run_parser.add_argument("--vcf-b", required=True, type=Path)
    run_parser.add_argument("--label-a", default="VCF_A")
    run_parser.add_argument("--label-b", default="VCF_B")
    run_parser.add_argument("--reference-dict", required=True, type=Path)
    run_parser.add_argument("--contig-list", required=True, type=Path)
    run_parser.add_argument("--output-dir", required=True, type=Path)
    run_parser.add_argument("--seg-dup-track", type=Path)
    run_parser.add_argument("--simple-repeat-track", type=Path)
    run_parser.add_argument("--repeatmasker-track", type=Path)
    run_parser.add_argument("--gatk-path", default="gatk")
    run_parser.add_argument("--java-options", default="-Xmx4g")
    run_parser.add_argument("--modules", help="Comma-separated module list to run")
    run_parser.add_argument("--pass-only", action="store_true")
    run_parser.add_argument(
        "--context-overlap",
        type=_context_overlap_value,
        default=0.5,
        help="Minimum overlap fraction required to assign span-based genomic contexts (default: 0.5)",
    )
    run_parser.add_argument("--per-chrom", action="store_true")
    run_parser.add_argument("--enable-site-match-table", action="store_true")
    run_parser.add_argument("--per-sample-counts-table", action="store_true")
    run_parser.add_argument("--ped", dest="ped_file", type=Path)
    run_parser.add_argument(
        "--num-workers",
        type=int,
        default=None,
        help="Number of contigs to process in parallel (default: auto, capped at 4 and the contig count)",
    )
    run_parser.set_defaults(handler=_handle_run)

    return parser


def _handle_validate(args: argparse.Namespace) -> int:
    if args.fix:
        out_path = args.out or _default_fix_output_path(args.vcf)
        result, output = fix_and_render(
            args.vcf,
            out_path,
            ploidy_table_path=args.ploidy_table,
            drop_bnd=bool(args.drop_bnd),
            drop_ctx=bool(args.drop_ctx),
        )
        print(output)
        return 1 if result.has_errors else 0
    summary, output = validate_and_render(args.vcf)
    print(output)
    return 1 if summary.has_errors else 0


def _handle_preprocess(args: argparse.Namespace) -> int:
    _configure_logging()
    contigs, contig_lengths = _resolve_requested_preprocess_contigs(
        args.contig_list,
        args.reference_dict,
        requested_contig=args.contig,
    )
    resolved_workers = _resolve_num_workers(args.num_workers, len(contigs))
    logging.getLogger(__name__).info(
        "Loaded preprocess inputs: %s contigs from %s; using %s worker(s)",
        len(contigs),
        args.contig_list,
        resolved_workers,
    )
    config = AnalysisConfig(
        vcf_a_path=args.vcf_a,
        vcf_b_path=args.vcf_b,
        output_dir=args.output_dir,
        reference_dict=args.reference_dict,
        contigs=contigs,
        contig_lengths=contig_lengths,
        n_workers=resolved_workers,
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


def _handle_analyze(args: argparse.Namespace) -> int:
    _configure_logging()
    config = _build_analysis_config(
        vcf_a_path=args.vcf_a,
        vcf_b_path=args.vcf_b,
        output_dir=args.output_dir,
        label_a=args.label_a,
        label_b=args.label_b,
        module_names=_parse_module_names(args.modules),
        pass_only=bool(args.pass_only),
        context_overlap=float(args.context_overlap),
        per_chrom=bool(args.per_chrom),
        enable_site_match_table=bool(args.enable_site_match_table),
        per_sample_counts_table=bool(args.per_sample_counts_table),
        ped_file=args.ped_file,
        requested_workers=args.num_workers,
    )
    return _run_analysis(config)


def _handle_run(args: argparse.Namespace) -> int:
    _configure_logging()
    contigs = read_contig_list(args.contig_list)
    contig_lengths = parse_reference_dict(args.reference_dict)
    resolved_workers = _resolve_num_workers(args.num_workers, len(contigs))
    logging.getLogger(__name__).info(
        "Starting end-to-end run with %s contig(s) and %s worker(s)",
        len(contigs),
        resolved_workers,
    )
    preprocess_config = AnalysisConfig(
        vcf_a_path=args.vcf_a,
        vcf_b_path=args.vcf_b,
        output_dir=args.output_dir,
        reference_dict=args.reference_dict,
        contigs=contigs,
        contig_lengths=contig_lengths,
        n_workers=resolved_workers,
        seg_dup_track=args.seg_dup_track,
        simple_repeat_track=args.simple_repeat_track,
        repeatmasker_track=args.repeatmasker_track,
        gatk_path=args.gatk_path,
        java_options=args.java_options,
    )
    annotated_a, annotated_b = run_preprocess(preprocess_config)
    print(f"annotated_a={annotated_a}")
    print(f"annotated_b={annotated_b}")

    analyze_config = _build_analysis_config(
        vcf_a_path=annotated_a,
        vcf_b_path=annotated_b,
        output_dir=args.output_dir,
        label_a=args.label_a,
        label_b=args.label_b,
        module_names=_parse_module_names(args.modules),
        pass_only=bool(args.pass_only),
        context_overlap=float(args.context_overlap),
        per_chrom=bool(args.per_chrom),
        enable_site_match_table=bool(args.enable_site_match_table),
        per_sample_counts_table=bool(args.per_sample_counts_table),
        ped_file=args.ped_file,
        requested_workers=resolved_workers,
        contigs=contigs,
        contig_lengths=contig_lengths,
    )
    return _run_analysis(analyze_config)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    raise SystemExit(main())
