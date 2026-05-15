"""Privacy-safe logging helpers for gatk-sv-ploidy CLI tools."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys
import time
from importlib import metadata
from typing import Any, Iterable, Iterator, Sequence


_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "PYTORCH_NUM_THREADS",
    "CUDA_VISIBLE_DEVICES",
    "PYTHONHASHSEED",
)
_DEPENDENCY_DISTS = (
    "gatk-sv-ploidy",
    "numpy",
    "pandas",
    "torch",
    "pyro-ppl",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "pysam",
)
_PATH_ARG_PARTS = (
    "artifacts",
    "data",
    "depth",
    "dir",
    "file",
    "ht",
    "input",
    "json",
    "list",
    "output",
    "path",
    "predictions",
    "regions",
    "stats",
    "tsv",
)
_SAFE_STRING_ARGS = {
    "autosome_prior_mode",
    "binq_field",
    "cn_inference_method",
    "device",
    "site_af_estimator",
}


class _UtcIsoFormatter(logging.Formatter):
    """Formatter that emits millisecond-resolution UTC ISO timestamps."""

    def formatTime(
        self,
        record: logging.LogRecord,
        datefmt: str | None = None,
    ) -> str:
        del datefmt
        timestamp = datetime.fromtimestamp(record.created, tz=timezone.utc)
        return timestamp.isoformat(timespec="milliseconds").replace(
            "+00:00",
            "Z",
        )


def _safe_path_label(value: str | os.PathLike[str]) -> str:
    """Return a path basename suitable for privacy-aware logs."""
    path_str = str(value)
    stripped = path_str.rstrip(os.sep)
    if not stripped:
        return path_str
    return os.path.basename(stripped) or stripped


def _is_path_arg(name: str, value: Any) -> bool:
    if name in _SAFE_STRING_ARGS:
        return False
    if isinstance(value, os.PathLike):
        return True
    if not isinstance(value, str):
        return False
    normalized = name.replace("-", "_").lower()
    tokens = set(normalized.split("_"))
    if tokens.intersection(_PATH_ARG_PARTS):
        return True
    return os.sep in value or value.endswith(
        (
            ".bed",
            ".gz",
            ".ht",
            ".json",
            ".npz",
            ".tsv",
            ".txt",
            ".vcf",
        )
    )


def _sanitize_value(name: str, value: Any) -> Any:
    """Sanitize an argparse value while preserving useful run metadata."""
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, os.PathLike):
        return {"path_label": _safe_path_label(value)}
    if isinstance(value, str):
        if _is_path_arg(name, value):
            return {"path_label": _safe_path_label(value)}
        if "sample" in name.lower() and value:
            return "<provided>"
        return value
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        values = list(value)
        if all(isinstance(item, (int, float, bool)) or item is None for item in values):
            return values
        if any(_is_path_arg(name, item) for item in values):
            return {
                "n_values": len(values),
                "path_labels": [_safe_path_label(item) for item in values[:5]],
                "truncated": len(values) > 5,
            }
        if "sample" in name.lower():
            return {"n_values": len(values)}
        return values
    return str(value)


def sanitize_namespace(args: argparse.Namespace) -> dict[str, Any]:
    """Return a privacy-safe dictionary of parsed command-line arguments."""
    return {
        key: _sanitize_value(key, value)
        for key, value in sorted(vars(args).items())
    }


def _json_dumps(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _dependency_versions() -> dict[str, str]:
    versions = {}
    for dist_name in _DEPENDENCY_DISTS:
        try:
            versions[dist_name] = metadata.version(dist_name)
        except metadata.PackageNotFoundError:
            versions[dist_name] = "not-installed"
        except Exception:
            versions[dist_name] = "unknown"
    return versions


def _git_commit() -> str:
    try:
        result = subprocess.run(
            [
                "git",
                "-C",
                str(Path(__file__).resolve().parent),
                "rev-parse",
                "--short",
                "HEAD",
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return "unknown"
    if result.returncode != 0:
        return "unknown"
    return result.stdout.strip() or "unknown"


def _environment_summary() -> dict[str, Any]:
    hostname = socket.gethostname()
    hostname_hash = hashlib.sha256(hostname.encode("utf-8")).hexdigest()[:12]
    return {
        "cpu_count": os.cpu_count(),
        "cwd_label": _safe_path_label(Path.cwd()),
        "hostname_sha256_12": hostname_hash,
        "pid": os.getpid(),
        "platform": platform.platform(),
        "python_executable_label": _safe_path_label(sys.executable),
        "python_version": platform.python_version(),
        "thread_env": {
            name: os.environ[name]
            for name in _THREAD_ENV_VARS
            if name in os.environ
        },
    }


def _reset_ploidy_handlers(logger: logging.Logger) -> None:
    for handler in list(logger.handlers):
        if getattr(handler, "_gatk_sv_ploidy_handler", False):
            logger.removeHandler(handler)
            handler.close()


def configure_tool_logging(
    *,
    tool_name: str,
    output_dir: str | os.PathLike[str],
    log_filename: str | None = None,
) -> Path:
    """Configure package logging to stderr and a per-tool text log file."""
    log_dir = Path(output_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / (log_filename or f"{tool_name}.log")

    package_logger = logging.getLogger("gatk_sv_ploidy")
    _reset_ploidy_handlers(package_logger)
    package_logger.setLevel(logging.INFO)
    package_logger.propagate = False

    formatter = _UtcIsoFormatter(
        "%(asctime)s %(levelname)s %(name)s: %(message)s"
    )
    stream_handler = logging.StreamHandler(stream=sys.stderr)
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    stream_handler._gatk_sv_ploidy_handler = True  # type: ignore[attr-defined]

    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    file_handler._gatk_sv_ploidy_handler = True  # type: ignore[attr-defined]

    package_logger.addHandler(stream_handler)
    package_logger.addHandler(file_handler)
    return log_path


def log_startup_metadata(
    logger: logging.Logger,
    *,
    tool_name: str,
    args: argparse.Namespace,
    log_path: str | os.PathLike[str],
    random_seeds: dict[str, Any] | None = None,
) -> None:
    """Log privacy-safe startup metadata for reproducibility."""
    logger.info("Starting gatk-sv-ploidy %s", tool_name)
    logger.info("Log file: %s", _safe_path_label(log_path))
    logger.info("Git commit: %s", _git_commit())
    logger.info("Dependency versions: %s", _json_dumps(_dependency_versions()))
    logger.info("Runtime environment: %s", _json_dumps(_environment_summary()))
    logger.info("Command arguments: %s", _json_dumps(sanitize_namespace(args)))
    if random_seeds:
        logger.info("Random seeds: %s", _json_dumps(random_seeds))


@contextmanager
def tool_logging_context(
    *,
    tool_name: str,
    output_dir: str | os.PathLike[str],
    args: argparse.Namespace,
    log_filename: str | None = None,
    random_seeds: dict[str, Any] | None = None,
) -> Iterator[logging.Logger]:
    """Configure logging and record start, failure, and completion events."""
    log_path = configure_tool_logging(
        tool_name=tool_name,
        output_dir=output_dir,
        log_filename=log_filename,
    )
    logger = logging.getLogger(f"gatk_sv_ploidy.{tool_name}")
    start_time = time.perf_counter()
    log_startup_metadata(
        logger,
        tool_name=tool_name,
        args=args,
        log_path=log_path,
        random_seeds=random_seeds,
    )
    try:
        yield logger
    except BaseException:
        elapsed = time.perf_counter() - start_time
        logger.exception("gatk-sv-ploidy %s failed after %.2f seconds", tool_name, elapsed)
        raise
    else:
        elapsed = time.perf_counter() - start_time
        logger.info(
            "gatk-sv-ploidy %s completed successfully in %.2f seconds",
            tool_name,
            elapsed,
        )


def log_output_artifacts(
    logger: logging.Logger,
    artifacts: Iterable[str | os.PathLike[str]],
) -> None:
    """Log output artifact basenames without exposing full paths."""
    labels = [_safe_path_label(artifact) for artifact in artifacts]
    logger.info(
        "Output artifacts: %s",
        _json_dumps({"n_artifacts": len(labels), "path_labels": labels}),
    )


def progress_iter(
    iterable: Iterable[Any],
    *,
    logger: logging.Logger,
    description: str,
    unit: str,
    enabled: bool = True,
    total: int | None = None,
    min_interval_seconds: float = 30.0,
    log_every: int | None = None,
) -> Iterator[Any]:
    """Yield items while logging bounded progress updates."""
    if not enabled:
        yield from iterable
        return

    start_time = time.perf_counter()
    last_log_time = start_time
    last_log_count = 0
    count = 0
    total_text = str(total) if total is not None else "unknown"
    logger.info("%s started: total_%s=%s", description, unit, total_text)

    for item in iterable:
        yield item
        count += 1
        now = time.perf_counter()
        should_log_count = log_every is not None and count % log_every == 0
        should_log_time = now - last_log_time >= min_interval_seconds
        should_log_final = total is not None and count >= total
        if should_log_count or should_log_time or should_log_final:
            elapsed = now - start_time
            if total:
                logger.info(
                    "%s progress: %d/%d %s (%.1f%%, %.1f seconds elapsed)",
                    description,
                    count,
                    total,
                    unit,
                    100.0 * count / total,
                    elapsed,
                )
            else:
                logger.info(
                    "%s progress: %d %s (%.1f seconds elapsed)",
                    description,
                    count,
                    unit,
                    elapsed,
                )
            last_log_time = now
            last_log_count = count

    if count != last_log_count:
        elapsed = time.perf_counter() - start_time
        logger.info(
            "%s complete: %d/%s %s (%.1f seconds elapsed)",
            description,
            count,
            total_text,
            unit,
            elapsed,
        )
