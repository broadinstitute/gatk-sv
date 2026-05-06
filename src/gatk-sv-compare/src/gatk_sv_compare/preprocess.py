"""Preprocess orchestration for gatk-sv-compare."""

from __future__ import annotations

from dataclasses import dataclass
import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pysam

from .config import AnalysisConfig
from .dimensions import STATUS_MATCHED, STATUS_UNMATCHED
from .vcf_format import resolve_record_svtype

logger = logging.getLogger(__name__)


class PreprocessError(RuntimeError):
    """Raised when preprocessing fails."""


StrSiteKey = Tuple[str, int, Tuple[Tuple[str, str], ...]]
_STR_MATCH_EXCLUDED_INFO = frozenset({"STATUS", "TRUTH_VID"})


@dataclass(frozen=True)
class SplitVcfPaths:
    sv_vcf: Path
    str_vcf: Path
    sv_count: int
    str_count: int


def read_contig_list(contig_list_path: Path) -> List[str]:
    """Read contigs from a plain-text list, ignoring blanks and comments."""
    contigs: List[str] = []
    for raw_line in contig_list_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        contigs.append(line)
    return contigs


def parse_reference_dict(reference_dict: Path) -> Dict[str, int]:
    """Parse a Picard/GATK sequence dictionary into contig lengths."""
    contig_lengths: Dict[str, int] = {}
    for line in reference_dict.read_text().splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(item.split(":", 1) for item in line.split("\t")[1:] if ":" in item)
        if "SN" in fields and "LN" in fields:
            contig_lengths[fields["SN"]] = int(fields["LN"])
    return contig_lengths


def get_tracks(config: AnalysisConfig) -> Dict[str, Path]:
    """Return configured genome tracks in stable CLI order."""
    tracks: Dict[str, Path] = {}
    if config.seg_dup_track is not None:
        tracks["segdup"] = config.seg_dup_track
    if config.simple_repeat_track is not None:
        tracks["simple_repeat"] = config.simple_repeat_track
    if config.repeatmasker_track is not None:
        tracks["repeatmasker"] = config.repeatmasker_track
    return tracks


def _record_identifier(record: pysam.VariantRecord) -> str:
    return record.id or f"{record.contig}:{record.pos}"


def _ensure_info_header(header: pysam.VariantHeader, key: str, number: str | int, value_type: str, description: str) -> None:
    if key not in header.info:
        header.add_meta("INFO", items=[("ID", key), ("Number", number), ("Type", value_type), ("Description", description)])


def _ensure_str_passthrough_headers(header: pysam.VariantHeader) -> None:
    _ensure_info_header(header, "SVTYPE", 1, "String", "SV type")
    _ensure_info_header(header, "STATUS", 1, "String", "Match status")
    _ensure_info_header(header, "TRUTH_VID", 1, "String", "Matched truth variant ID")


def _translate_record(record: pysam.VariantRecord, header: pysam.VariantHeader) -> pysam.VariantRecord:
    translated = record.copy()
    translated.translate(header)
    return translated


def _set_missing_svtype(record: pysam.VariantRecord, svtype: Optional[str]) -> None:
    if svtype is None:
        return
    current = record.info.get("SVTYPE") if "SVTYPE" in record.header.info else None
    if current in (None, "."):
        record.info["SVTYPE"] = svtype


def _split_str_sv_vcfs(input_vcf: Path, sv_output: Path, str_output: Path) -> SplitVcfPaths:
    """Split one input VCF into GATK-compatible SV records and STR pass-through records."""
    sv_output.parent.mkdir(parents=True, exist_ok=True)
    sv_count = 0
    str_count = 0
    with pysam.VariantFile(str(input_vcf)) as in_vcf:
        header = in_vcf.header.copy()
        _ensure_str_passthrough_headers(header)
        with pysam.VariantFile(str(sv_output), "wz", header=header) as sv_vcf, pysam.VariantFile(str(str_output), "wz", header=header) as str_vcf:
            for record in in_vcf:
                svtype = resolve_record_svtype(record)
                if svtype == "STR":
                    output_record = _translate_record(record, str_vcf.header)
                    _set_missing_svtype(output_record, svtype)
                    str_vcf.write(output_record)
                    str_count += 1
                else:
                    output_record = _translate_record(record, sv_vcf.header)
                    _set_missing_svtype(output_record, svtype)
                    sv_vcf.write(output_record)
                    sv_count += 1
    pysam.tabix_index(str(sv_output), preset="vcf", force=True)
    pysam.tabix_index(str(str_output), preset="vcf", force=True)
    return SplitVcfPaths(sv_vcf=sv_output, str_vcf=str_output, sv_count=sv_count, str_count=str_count)


def _canonical_info_value(value: object) -> str:
    if value is True:
        return ""
    if isinstance(value, (tuple, list)):
        return ",".join(str(item) for item in value)
    return str(value)


def _str_site_key(record: pysam.VariantRecord) -> StrSiteKey:
    info_items: List[Tuple[str, str]] = []
    observed_svtype = False
    for key in sorted(str(key) for key in record.info.keys() if str(key) not in _STR_MATCH_EXCLUDED_INFO):
        value = record.info.get(key)
        if key == "SVTYPE":
            observed_svtype = True
            value = "STR"
        info_items.append((key, _canonical_info_value(value)))
    if not observed_svtype:
        info_items.append(("SVTYPE", "STR"))
        info_items.sort()
    return (record.contig, int(record.pos), tuple(info_items))


def _load_str_site_ids(str_vcf: Path) -> Dict[StrSiteKey, str]:
    site_ids: Dict[StrSiteKey, str] = {}
    with pysam.VariantFile(str(str_vcf)) as vcf:
        for record in vcf:
            key = _str_site_key(record)
            if key in site_ids:
                raise PreprocessError(f"Duplicate STR fixed site in {str_vcf}: {_record_identifier(record)} and {site_ids[key]}")
            site_ids[key] = _record_identifier(record)
    return site_ids


def _delete_info_if_present(record: pysam.VariantRecord, key: str) -> None:
    try:
        del record.info[key]
    except (KeyError, ValueError):
        pass


def _annotate_str_vcf(str_vcf: Path, output_path: Path, truth_ids_by_key: Dict[StrSiteKey, str]) -> Path:
    with pysam.VariantFile(str(str_vcf)) as in_vcf:
        header = in_vcf.header.copy()
        _ensure_str_passthrough_headers(header)
        with pysam.VariantFile(str(output_path), "wz", header=header) as out_vcf:
            for record in in_vcf:
                output_record = _translate_record(record, out_vcf.header)
                _set_missing_svtype(output_record, "STR")
                key = _str_site_key(output_record)
                truth_id = truth_ids_by_key.get(key)
                if truth_id is None:
                    output_record.info["STATUS"] = STATUS_UNMATCHED
                    _delete_info_if_present(output_record, "TRUTH_VID")
                else:
                    output_record.info["STATUS"] = STATUS_MATCHED
                    output_record.info["TRUTH_VID"] = truth_id
                out_vcf.write(output_record)
    pysam.tabix_index(str(output_path), preset="vcf", force=True)
    return output_path


def _annotate_str_pair(str_a: Path, str_b: Path, output_a: Path, output_b: Path) -> Tuple[Path, Path]:
    ids_a = _load_str_site_ids(str_a)
    ids_b = _load_str_site_ids(str_b)
    return (
        _annotate_str_vcf(str_a, output_a, ids_b),
        _annotate_str_vcf(str_b, output_b, ids_a),
    )


def _write_empty_vcf_like(source_path: Path, output_path: Path) -> Path:
    with pysam.VariantFile(str(source_path)) as source_vcf:
        header = source_vcf.header.copy()
        with pysam.VariantFile(str(output_path), "wz", header=header):
            pass
    pysam.tabix_index(str(output_path), preset="vcf", force=True)
    return output_path


def _add_missing_header_records(header: pysam.VariantHeader, extra_header: pysam.VariantHeader) -> None:
    for contig_name, contig in extra_header.contigs.items():
        if contig_name not in header.contigs:
            header.contigs.add(contig_name, length=contig.length)
    for key, value in extra_header.info.items():
        if key not in header.info:
            header.add_line(str(value.record).strip())
    for key, value in extra_header.formats.items():
        if key not in header.formats:
            header.add_line(str(value.record).strip())
    for key, value in extra_header.alts.items():
        if key not in header.alts:
            header.add_line(str(value).strip())
    for key, value in extra_header.filters.items():
        if key != "PASS" and key not in header.filters:
            header.add_line(str(value.record).strip())
    for sample in extra_header.samples:
        if sample not in header.samples:
            header.add_sample(sample)


def _iter_records_for_contig(vcf: pysam.VariantFile, contig: str) -> Iterable[pysam.VariantRecord]:
    try:
        yield from vcf.fetch(contig)
    except ValueError:
        for record in vcf:
            if record.contig == contig:
                yield record


def _merge_vcfs(input_paths: Sequence[Path], output_path: Path, contigs: Sequence[str]) -> Path:
    if not input_paths:
        raise PreprocessError("No input VCFs provided for merge.")
    with pysam.VariantFile(str(input_paths[0])) as first_vcf:
        header = first_vcf.header.copy()
    for input_path in input_paths[1:]:
        with pysam.VariantFile(str(input_path)) as vcf:
            _add_missing_header_records(header, vcf.header)

    opened_vcfs = [pysam.VariantFile(str(input_path)) for input_path in input_paths]
    try:
        with pysam.VariantFile(str(output_path), "wz", header=header) as out_vcf:
            for contig in contigs:
                records: List[Tuple[int, int, int, str, pysam.VariantRecord]] = []
                for source_index, vcf in enumerate(opened_vcfs):
                    for record in _iter_records_for_contig(vcf, contig):
                        records.append((int(record.pos), int(record.stop), source_index, _record_identifier(record), record.copy()))
                for _, _, _, _, record in sorted(records, key=lambda item: item[:4]):
                    output_record = _translate_record(record, out_vcf.header)
                    out_vcf.write(output_record)
    finally:
        for vcf in opened_vcfs:
            vcf.close()
    pysam.tabix_index(str(output_path), preset="vcf", force=True)
    return output_path


def build_svconcordance_command(
    *,
    eval_vcf: Path,
    truth_vcf: Path,
    contig: str,
    reference_dict: Path,
    output_path: Path,
    gatk_path: str,
    java_options: str,
    clustering_config: Optional[Path] = None,
    stratification_config: Optional[Path] = None,
) -> List[str]:
    """Build a GATK SVConcordance command."""
    command = [
        gatk_path,
        "--java-options",
        java_options,
        "SVConcordance",
        "-L",
        contig,
        "--sequence-dictionary",
        str(reference_dict),
        "--eval",
        str(eval_vcf),
        "--truth",
        str(truth_vcf),
        "-O",
        str(output_path),
    ]
    if clustering_config is not None:
        command.extend(["--clustering-config", str(clustering_config)])
    if stratification_config is not None:
        command.extend(["--stratify-config", str(stratification_config)])
    return command


def build_svregionoverlap_command(
    *,
    vcf: Path,
    reference_dict: Path,
    output_path: Path,
    gatk_path: str,
    java_options: str,
    tracks: Dict[str, Path],
) -> List[str]:
    """Build a GATK SVRegionOverlap command."""
    command = [
        gatk_path,
        "--java-options",
        java_options,
        "SVRegionOverlap",
        "-V",
        str(vcf),
        "-O",
        str(output_path),
        "--sequence-dictionary",
        str(reference_dict),
    ]
    for track_path in tracks.values():
        command.extend(["--track-intervals", str(track_path)])
    for track_name in tracks.keys():
        command.extend(["--track-name", track_name])
    return command


def run_command(command: Sequence[str], timeout_seconds: int = 3600) -> subprocess.CompletedProcess:
    """Run an external command and surface actionable errors."""
    logger.info("Running command: %s", " ".join(command))
    result = subprocess.run(
        list(command),
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
        check=False,
    )
    if result.returncode == 0:
        return result

    stderr = result.stderr or ""
    hints: List[str] = []
    if "OutOfMemoryError" in stderr:
        hints.append("Increase --java-options, e.g. '-Xmx8g'.")
    if "A USER ERROR has occurred" in stderr or "Unknown argument" in stderr:
        hints.append("Check the installed GATK version and command-line arguments.")
    hint_text = f" Hints: {' '.join(hints)}" if hints else ""
    raise PreprocessError(
        f"Command failed with exit code {result.returncode}: {' '.join(command)}\n{stderr}{hint_text}"
    )


def get_gatk_version(gatk_path: str) -> str:
    """Return the installed GATK version string if available."""
    result = run_command([gatk_path, "--version"], timeout_seconds=60)
    return (result.stdout or result.stderr).strip()


def _run_sv_concordance_contig(
    *,
    eval_vcf: Path,
    truth_vcf: Path,
    contig: str,
    reference_dict: Path,
    output_path: Path,
    gatk_path: str,
    java_options: str,
    clustering_config: Optional[Path],
    stratification_config: Optional[Path],
) -> Path:
    """Run SVConcordance for one contig."""
    command = build_svconcordance_command(
        eval_vcf=eval_vcf,
        truth_vcf=truth_vcf,
        contig=contig,
        reference_dict=reference_dict,
        output_path=output_path,
        gatk_path=gatk_path,
        java_options=java_options,
        clustering_config=clustering_config,
        stratification_config=stratification_config,
    )
    run_command(command)
    return output_path


def _run_sv_region_overlap_contig(
    *,
    vcf: Path,
    reference_dict: Path,
    output_path: Path,
    gatk_path: str,
    java_options: str,
    tracks: Dict[str, Path],
) -> Path:
    """Run SVRegionOverlap for one contig shard."""
    command = build_svregionoverlap_command(
        vcf=vcf,
        reference_dict=reference_dict,
        output_path=output_path,
        gatk_path=gatk_path,
        java_options=java_options,
        tracks=tracks,
    )
    run_command(command)
    return output_path


def concatenate_vcfs(input_paths: Sequence[Path], output_path: Path) -> Path:
    """Concatenate VCF shards into one bgzipped, tabix-indexed VCF."""
    if not input_paths:
        raise PreprocessError("No input VCF shards provided for concatenation.")

    with pysam.VariantFile(str(input_paths[0])) as first_vcf:
        header = first_vcf.header.copy()
        with pysam.VariantFile(str(output_path), "wz", header=header) as out_vcf:
            for shard_path in input_paths:
                with pysam.VariantFile(str(shard_path)) as shard_vcf:
                    for record in shard_vcf:
                        out_vcf.write(record)
    pysam.tabix_index(str(output_path), preset="vcf", force=True)
    return output_path


def _copy_with_index(source: Path, destination: Path) -> Path:
    """Copy a VCF and its tabix index."""
    shutil.copy2(source, destination)
    source_index = Path(str(source) + ".tbi")
    destination_index = Path(str(destination) + ".tbi")
    if source_index.exists():
        shutil.copy2(source_index, destination_index)
    else:
        pysam.tabix_index(str(destination), preset="vcf", force=True)
    return destination


def _scatter_concordance(
    *,
    eval_vcf: Path,
    truth_vcf: Path,
    contigs: Sequence[str],
    reference_dict: Path,
    output_dir: Path,
    gatk_path: str,
    java_options: str,
    n_workers: int,
    clustering_config: Optional[Path],
    stratification_config: Optional[Path],
) -> List[Path]:
    """Run SVConcordance across contigs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    max_workers = max(1, n_workers)
    total = len(contigs)
    logger.info(
        "Starting SVConcordance scatter for %s contigs: %s vs %s",
        total,
        eval_vcf.name,
        truth_vcf.name,
    )
    results: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_sv_concordance_contig,
                eval_vcf=eval_vcf,
                truth_vcf=truth_vcf,
                contig=contig,
                reference_dict=reference_dict,
                output_path=output_dir / f"{contig}.vcf.gz",
                gatk_path=gatk_path,
                java_options=java_options,
                clustering_config=clustering_config,
                stratification_config=stratification_config,
            ): contig
            for contig in contigs
        }
        for completed, future in enumerate(as_completed(futures), start=1):
            contig = futures[future]
            results.append(future.result())
            logger.info("Completed SVConcordance shard %s/%s: %s", completed, total, contig)
    return results


def _scatter_region_overlap(
    *,
    vcfs: Sequence[Path],
    reference_dict: Path,
    output_dir: Path,
    gatk_path: str,
    java_options: str,
    n_workers: int,
    tracks: Dict[str, Path],
) -> List[Path]:
    """Run SVRegionOverlap across per-contig concordance shards."""
    output_dir.mkdir(parents=True, exist_ok=True)
    max_workers = max(1, n_workers)
    total = len(vcfs)
    logger.info("Starting SVRegionOverlap scatter for %s contig shards", total)
    results: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_sv_region_overlap_contig,
                vcf=vcf,
                reference_dict=reference_dict,
                output_path=output_dir / vcf.name,
                gatk_path=gatk_path,
                java_options=java_options,
                tracks=tracks,
            ): vcf.name
            for vcf in vcfs
        }
        for completed, future in enumerate(as_completed(futures), start=1):
            shard_name = futures[future]
            results.append(future.result())
            logger.info("Completed SVRegionOverlap shard %s/%s: %s", completed, total, shard_name)
    return results


def run_preprocess(config: AnalysisConfig) -> Tuple[Path, Path]:
    """Run the Phase 2 preprocess workflow and return annotated VCF paths."""
    if config.vcf_a_path is None or config.vcf_b_path is None:
        raise PreprocessError("Both vcf_a_path and vcf_b_path are required for preprocessing.")
    if config.reference_dict is None:
        raise PreprocessError("reference_dict is required for preprocessing.")
    if not config.contigs:
        raise PreprocessError("At least one contig is required for preprocessing.")

    preprocess_dir = config.output_dir / "preprocess"
    preprocess_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting preprocess into %s", preprocess_dir)
    logger.info(
        "Inputs: vcf_a=%s, vcf_b=%s, contigs=%s, workers=%s",
        config.vcf_a_path,
        config.vcf_b_path,
        len(config.contigs),
        max(1, config.n_workers),
    )

    version = get_gatk_version(config.gatk_path)
    logger.info("Using GATK: %s", version)

    tracks = get_tracks(config)
    if tracks:
        logger.info("Region-overlap tracks enabled: %s", ", ".join(tracks.keys()))
    else:
        logger.info("No region-overlap tracks configured; annotated outputs will mirror concordance outputs")

    split_a = _split_str_sv_vcfs(
        config.vcf_a_path,
        preprocess_dir / "input_a.sv.vcf.gz",
        preprocess_dir / "input_a.str.vcf.gz",
    )
    split_b = _split_str_sv_vcfs(
        config.vcf_b_path,
        preprocess_dir / "input_b.sv.vcf.gz",
        preprocess_dir / "input_b.str.vcf.gz",
    )
    logger.info(
        "Split inputs: A has %s SV records and %s STR records; B has %s SV records and %s STR records",
        split_a.sv_count,
        split_a.str_count,
        split_b.sv_count,
        split_b.str_count,
    )

    logger.info("Running concordance with VCF A as eval and VCF B as truth")
    if split_a.sv_count:
        concordance_a_shards = _scatter_concordance(
            eval_vcf=split_a.sv_vcf,
            truth_vcf=split_b.sv_vcf,
            contigs=config.contigs,
            reference_dict=config.reference_dict,
            output_dir=preprocess_dir / "concordance_a_eval.per_contig",
            gatk_path=config.gatk_path,
            java_options=config.java_options,
            n_workers=config.n_workers,
            clustering_config=config.clustering_config,
            stratification_config=config.stratification_config,
        )
    else:
        logger.info("Skipping VCF A concordance because A has no non-STR SV records")
        concordance_a_shards = []
    logger.info("Running concordance with VCF B as eval and VCF A as truth")
    if split_b.sv_count:
        concordance_b_shards = _scatter_concordance(
            eval_vcf=split_b.sv_vcf,
            truth_vcf=split_a.sv_vcf,
            contigs=config.contigs,
            reference_dict=config.reference_dict,
            output_dir=preprocess_dir / "concordance_b_eval.per_contig",
            gatk_path=config.gatk_path,
            java_options=config.java_options,
            n_workers=config.n_workers,
            clustering_config=config.clustering_config,
            stratification_config=config.stratification_config,
        )
    else:
        logger.info("Skipping VCF B concordance because B has no non-STR SV records")
        concordance_b_shards = []

    logger.info("Concatenating %s concordance shards into %s", len(concordance_a_shards), preprocess_dir / "concordance_a.vcf.gz")
    concordance_a = concatenate_vcfs(sorted(concordance_a_shards), preprocess_dir / "concordance_a.vcf.gz") if concordance_a_shards else _write_empty_vcf_like(split_a.sv_vcf, preprocess_dir / "concordance_a.vcf.gz")
    logger.info("Concatenating %s concordance shards into %s", len(concordance_b_shards), preprocess_dir / "concordance_b.vcf.gz")
    concordance_b = concatenate_vcfs(sorted(concordance_b_shards), preprocess_dir / "concordance_b.vcf.gz") if concordance_b_shards else _write_empty_vcf_like(split_b.sv_vcf, preprocess_dir / "concordance_b.vcf.gz")

    annotated_a = preprocess_dir / "annotated_a.vcf.gz"
    annotated_b = preprocess_dir / "annotated_b.vcf.gz"
    annotated_sv_a = preprocess_dir / "annotated_a.sv.vcf.gz"
    annotated_sv_b = preprocess_dir / "annotated_b.sv.vcf.gz"
    if tracks:
        logger.info("Running region-overlap annotation for annotated_a")
        annotated_a_shards = _scatter_region_overlap(
            vcfs=sorted(concordance_a_shards),
            reference_dict=config.reference_dict,
            output_dir=preprocess_dir / "annotated_a.per_contig",
            gatk_path=config.gatk_path,
            java_options=config.java_options,
            n_workers=config.n_workers,
            tracks=tracks,
        ) if concordance_a_shards else []
        logger.info("Running region-overlap annotation for annotated_b")
        annotated_b_shards = _scatter_region_overlap(
            vcfs=sorted(concordance_b_shards),
            reference_dict=config.reference_dict,
            output_dir=preprocess_dir / "annotated_b.per_contig",
            gatk_path=config.gatk_path,
            java_options=config.java_options,
            n_workers=config.n_workers,
            tracks=tracks,
        ) if concordance_b_shards else []
        logger.info("Concatenating %s annotated shards into %s", len(annotated_a_shards), annotated_sv_a)
        concatenate_vcfs(sorted(annotated_a_shards), annotated_sv_a) if annotated_a_shards else _write_empty_vcf_like(concordance_a, annotated_sv_a)
        logger.info("Concatenating %s annotated shards into %s", len(annotated_b_shards), annotated_sv_b)
        concatenate_vcfs(sorted(annotated_b_shards), annotated_sv_b) if annotated_b_shards else _write_empty_vcf_like(concordance_b, annotated_sv_b)
    else:
        logger.info("Copying concordance outputs to annotated outputs (no region-overlap stage)")
        _copy_with_index(concordance_a, annotated_sv_a)
        _copy_with_index(concordance_b, annotated_sv_b)

    annotated_str_a, annotated_str_b = _annotate_str_pair(
        split_a.str_vcf,
        split_b.str_vcf,
        preprocess_dir / "annotated_a.str.vcf.gz",
        preprocess_dir / "annotated_b.str.vcf.gz",
    )
    _merge_vcfs([annotated_sv_a, annotated_str_a], annotated_a, config.contigs)
    _merge_vcfs([annotated_sv_b, annotated_str_b], annotated_b, config.contigs)

    logger.info("Preprocess complete: annotated_a=%s annotated_b=%s", annotated_a, annotated_b)
    return annotated_a, annotated_b
