"""Preprocess orchestration for gatk-sv-compare."""

from __future__ import annotations

import logging
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pysam

from .config import AnalysisConfig

logger = logging.getLogger(__name__)


class PreprocessError(RuntimeError):
    """Raised when preprocessing fails."""


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
    track_names: Optional[Sequence[str]] = None,
    track_intervals: Optional[Sequence[Path]] = None,
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
    for track_name in track_names or ():
        command.extend(["--track-name", track_name])
    for track_interval in track_intervals or ():
        command.extend(["--track-intervals", str(track_interval)])
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
    track_names: Optional[Sequence[str]],
    track_intervals: Optional[Sequence[Path]],
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
        track_names=track_names,
        track_intervals=track_intervals,
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
    track_names: Optional[Sequence[str]],
    track_intervals: Optional[Sequence[Path]],
) -> List[Path]:
    """Run SVConcordance across contigs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    max_workers = max(1, n_workers)
    results: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
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
                track_names=track_names,
                track_intervals=track_intervals,
            )
            for contig in contigs
        ]
        for future in futures:
            results.append(future.result())
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
    results: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _run_sv_region_overlap_contig,
                vcf=vcf,
                reference_dict=reference_dict,
                output_path=output_dir / vcf.name,
                gatk_path=gatk_path,
                java_options=java_options,
                tracks=tracks,
            )
            for vcf in vcfs
        ]
        for future in futures:
            results.append(future.result())
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

    version = get_gatk_version(config.gatk_path)
    logger.info("Using GATK: %s", version)

    tracks = get_tracks(config)
    track_names = list(tracks.keys()) if tracks else None
    track_intervals = list(tracks.values()) if tracks else None

    concordance_a_shards = _scatter_concordance(
        eval_vcf=config.vcf_a_path,
        truth_vcf=config.vcf_b_path,
        contigs=config.contigs,
        reference_dict=config.reference_dict,
        output_dir=preprocess_dir / "concordance_a_eval.per_contig",
        gatk_path=config.gatk_path,
        java_options=config.java_options,
        n_workers=config.n_workers,
        clustering_config=config.clustering_config,
        stratification_config=config.stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
    )
    concordance_b_shards = _scatter_concordance(
        eval_vcf=config.vcf_b_path,
        truth_vcf=config.vcf_a_path,
        contigs=config.contigs,
        reference_dict=config.reference_dict,
        output_dir=preprocess_dir / "concordance_b_eval.per_contig",
        gatk_path=config.gatk_path,
        java_options=config.java_options,
        n_workers=config.n_workers,
        clustering_config=config.clustering_config,
        stratification_config=config.stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
    )

    concordance_a = concatenate_vcfs(
        sorted(concordance_a_shards), preprocess_dir / "concordance_a.vcf.gz"
    )
    concordance_b = concatenate_vcfs(
        sorted(concordance_b_shards), preprocess_dir / "concordance_b.vcf.gz"
    )

    annotated_a = preprocess_dir / "annotated_a.vcf.gz"
    annotated_b = preprocess_dir / "annotated_b.vcf.gz"
    if tracks:
        annotated_a_shards = _scatter_region_overlap(
            vcfs=sorted(concordance_a_shards),
            reference_dict=config.reference_dict,
            output_dir=preprocess_dir / "annotated_a.per_contig",
            gatk_path=config.gatk_path,
            java_options=config.java_options,
            n_workers=config.n_workers,
            tracks=tracks,
        )
        annotated_b_shards = _scatter_region_overlap(
            vcfs=sorted(concordance_b_shards),
            reference_dict=config.reference_dict,
            output_dir=preprocess_dir / "annotated_b.per_contig",
            gatk_path=config.gatk_path,
            java_options=config.java_options,
            n_workers=config.n_workers,
            tracks=tracks,
        )
        concatenate_vcfs(sorted(annotated_a_shards), annotated_a)
        concatenate_vcfs(sorted(annotated_b_shards), annotated_b)
    else:
        _copy_with_index(concordance_a, annotated_a)
        _copy_with_index(concordance_b, annotated_b)

    return annotated_a, annotated_b
