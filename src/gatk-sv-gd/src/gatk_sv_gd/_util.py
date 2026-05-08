"""
Shared configuration and utility functions for gatk_sv_gd.

Module-level state (e.g. VERBOSE flag) lives here so that every sub-module
can import it without creating circular dependencies.
"""

import atexit
import json
import logging
import os
import platform
import re
import subprocess
import sys
import threading
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd


# Module-level verbosity flag, set from --verbose in subcommand CLIs.
VERBOSE = False

_LOGGER_NAME = "gatk_sv_gd"
_MANAGED_HANDLER_ATTR = "_gatk_sv_gd_managed_handler"
_ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
_PATH_OR_URI_RE = re.compile(
    r"(?:[A-Za-z][A-Za-z0-9+.-]*://\S+|(?:~?/|/)[^\s'\"<>]+|[A-Za-z]:\\[^\s'\"<>]+)"
)
_FILENAME_RE = re.compile(
    r"\b[^\s'\"<>/\\]+\."
    r"(?:bed|bam|bai|cram|crai|vcf|bcf|gvcf|tsv|csv|txt|json|jsonl|gz|bgz|zip|"
    r"pdf|png|svg|npz|npy|hdf5|h5|pkl|pickle|wdl|py)"
    r"(?:\.[A-Za-z0-9]+)?\b"
)
_SENSITIVE_ARG_KEYWORDS = (
    "baf",
    "bed",
    "calls",
    "counts",
    "dir",
    "file",
    "gd_table",
    "gtf",
    "input",
    "locus",
    "loci",
    "manifest",
    "matrix",
    "output",
    "path",
    "ploidy",
    "posteriors",
    "prefix",
    "raw",
    "region",
    "report",
    "sample",
    "table",
    "truth",
    "vcf",
)
_logging_session = None


class _UTCTextFormatter(logging.Formatter):
    """Human-readable formatter with high-resolution UTC timestamps."""

    def formatTime(self, record, datefmt=None):  # noqa: N802 - logging API name
        timestamp = datetime.fromtimestamp(record.created, timezone.utc)
        return timestamp.isoformat(timespec="milliseconds").replace("+00:00", "Z")

    def format(self, record):
        formatted = super().format(record)
        fields = getattr(record, "fields", None)
        if fields:
            safe_fields = _sanitize_log_fields(fields)
            formatted = f"{formatted} {json.dumps(safe_fields, sort_keys=True, default=str)}"
        return _sanitize_log_text(formatted)


class _JSONFormatter(logging.Formatter):
    """JSON-lines formatter for workflow and sweep-scale log consumers."""

    def format(self, record):
        payload: Dict[str, Any] = {
            "timestamp": datetime.fromtimestamp(
                record.created, timezone.utc,
            ).isoformat(timespec="milliseconds").replace("+00:00", "Z"),
            "level": record.levelname,
            "logger": record.name,
            "message": _sanitize_log_text(record.getMessage()),
            "thread": record.threadName,
        }
        fields = getattr(record, "fields", None)
        if fields:
            payload.update(_sanitize_log_fields(fields))
        if record.exc_info:
            payload["exception"] = _sanitize_log_text(
                self.formatException(record.exc_info)
            )
        return json.dumps(payload, sort_keys=True, default=str)


class LoggingStream:
    """File-like stream that routes writes through ``logging``.

    ``print`` calls in existing modules are diagnostic output.  Replacing
    stdout/stderr with this stream keeps stdout clean for data while still
    preserving the diagnostic text in stderr and the per-run log file.
    Carriage-return progress updates are collapsed to their final line so
    log files are not filled with terminal redraws.
    """

    def __init__(self, logger: logging.Logger, level: int, passthrough_stream=None):
        self.logger = logger
        self.level = level
        self._passthrough_stream = passthrough_stream
        self.encoding = getattr(sys.__stderr__, "encoding", "utf-8")
        self._buffer = ""
        self._closed = False
        self._lock = threading.RLock()

    def writable(self):
        return True

    def isatty(self):
        return False

    def write(self, message: str):
        if self._closed:
            return 0
        text = str(message)
        if self._passthrough_stream is not None and "\r" in text:
            passthrough_text = _sanitize_log_text(_ANSI_ESCAPE_RE.sub("", text))
            self._passthrough_stream.write(passthrough_text)
            self._passthrough_stream.flush()
        with self._lock:
            for character in text:
                if character == "\r":
                    self._buffer = ""
                elif character == "\n":
                    self._emit_buffer()
                    self._buffer = ""
                else:
                    self._buffer += character
        return len(text)

    def flush(self):
        return None

    def close(self):
        with self._lock:
            self._emit_buffer()
            self._buffer = ""
            self._closed = True

    def _emit_buffer(self):
        cleaned = _sanitize_log_text(_ANSI_ESCAPE_RE.sub("", self._buffer).strip())
        if not cleaned:
            return
        self.logger.log(self._level_for_message(cleaned), cleaned)

    def _level_for_message(self, message: str) -> int:
        upper_message = message.lstrip().upper()
        if upper_message.startswith(("ERROR", "FATAL", "CRITICAL")):
            return logging.ERROR
        if upper_message.startswith(("WARNING", "WARN")) or " WARNING:" in upper_message:
            return logging.WARNING
        return self.level


class LoggingSession:
    """Handle returned by ``setup_logging`` for tests and explicit cleanup."""

    def __init__(
        self,
        log_path: str,
        original_stdout,
        original_stderr,
        handlers,
    ):
        self.log_path = log_path
        self._original_stdout = original_stdout
        self._original_stderr = original_stderr
        self._handlers = list(handlers)
        self._closed = False

    def flush(self):
        if self._closed:
            return
        for handler in self._handlers:
            handler.flush()

    def close(self):
        if self._closed:
            return
        if isinstance(sys.stdout, LoggingStream):
            sys.stdout.close()
            sys.stdout = self._original_stdout
        if isinstance(sys.stderr, LoggingStream):
            sys.stderr.close()
            sys.stderr = self._original_stderr
        root_logger = logging.getLogger()
        for handler in self._handlers:
            handler.flush()
            root_logger.removeHandler(handler)
            handler.close()
        logging.captureWarnings(False)
        self._closed = True


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a package logger."""
    if not name:
        return logging.getLogger(_LOGGER_NAME)
    if name.startswith(_LOGGER_NAME):
        return logging.getLogger(name)
    return logging.getLogger(f"{_LOGGER_NAME}.{name}")


def vlog(msg: str) -> None:
    """Emit a verbose developer diagnostic at DEBUG level."""
    get_logger().debug(msg)


def posterior_probability_to_qual(
    probability: Any,
    max_qual: float = 99.0,
) -> Any:
    """Convert posterior event probabilities into capped phred-scale QUAL.

    The returned QUAL is based on posterior event-vs-non-event log odds,
    computed as ``10 * log10(p / (1 - p))`` and clipped to ``[0, max_qual]``.
    Scalar inputs return a float; array-like inputs return a NumPy array.
    """
    probabilities = np.clip(np.asarray(probability, dtype=float), 0.0, 1.0)
    min_probability = 10.0 ** (-float(max_qual) / 10.0)
    event_probabilities = np.clip(probabilities, min_probability, 1.0)
    non_event_probabilities = np.clip(1.0 - probabilities, min_probability, 1.0)
    qual = np.clip(
        10.0 * np.log10(event_probabilities / non_event_probabilities),
        0.0,
        float(max_qual),
    )
    if qual.ndim == 0:
        return float(qual)
    return qual


def posterior_called_state_to_qual(
    event_probability: Any,
    called_event: Any,
    max_qual: float = 99.0,
) -> Any:
    """Convert posterior support for an expected state into capped QUAL.

    When ``called_event`` is truthy, the expected state is DEL/DUP on at least
    one haplotype and QUAL is based on ``p / (1 - p)``. When ``called_event`` is
    falsy, the expected state is the complement and QUAL is based on
    ``(1 - p) / p``. Scalar inputs return a float; array-like inputs return a
    NumPy array.
    """
    probabilities = np.clip(np.asarray(event_probability, dtype=float), 0.0, 1.0)
    called_mask = np.asarray(called_event, dtype=bool)
    min_probability = 10.0 ** (-float(max_qual) / 10.0)
    expected_probabilities = np.where(called_mask, probabilities, 1.0 - probabilities)
    alternative_probabilities = np.where(called_mask, 1.0 - probabilities, probabilities)
    expected_probabilities = np.clip(expected_probabilities, min_probability, 1.0)
    alternative_probabilities = np.clip(
        alternative_probabilities,
        min_probability,
        1.0,
    )
    qual = np.clip(
        10.0 * np.log10(expected_probabilities / alternative_probabilities),
        0.0,
        float(max_qual),
    )
    if qual.ndim == 0:
        return float(qual)
    return qual


def setup_logging(
    output_dir: str,
    filename: str = "log.txt",
    *,
    verbose: bool = False,
    command: Optional[str] = None,
    args: Any = None,
    seed_info: Optional[Dict[str, Any]] = None,
    log_format: Optional[str] = None,
) -> LoggingSession:
    """Configure timestamped diagnostics for a GD subcommand.

    Diagnostics are written to stderr and to ``output_dir/filename``.
    Existing diagnostic ``print`` calls are routed through logging so stdout
    remains clean for machine-readable data.  Set ``GATK_SV_GD_LOG_FORMAT``
    to ``json`` for JSON-lines logs.
    """
    global _logging_session

    os.makedirs(output_dir, exist_ok=True)
    if _logging_session is not None:
        _logging_session.close()

    log_path = os.path.join(output_dir, filename)
    selected_format = (log_format or os.getenv("GATK_SV_GD_LOG_FORMAT", "text")).lower()
    log_level = logging.DEBUG if verbose else logging.INFO

    root_logger = logging.getLogger()
    _remove_managed_handlers(root_logger)
    root_logger.setLevel(log_level)

    formatter: logging.Formatter
    if selected_format == "json":
        formatter = _JSONFormatter()
    else:
        formatter = _UTCTextFormatter(
            "%(asctime)s %(levelname)s %(name)s: %(message)s"
        )

    stream_handler = logging.StreamHandler(sys.__stderr__)
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    handlers = [stream_handler, file_handler]
    for handler in handlers:
        setattr(handler, _MANAGED_HANDLER_ATTR, True)
        handler.setLevel(log_level)
        handler.setFormatter(formatter)
        root_logger.addHandler(handler)

    logging.captureWarnings(True)
    logging.getLogger("py.warnings").setLevel(logging.WARNING)
    get_logger().setLevel(log_level)

    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = LoggingStream(get_logger("stdout"), logging.INFO)
    sys.stderr = LoggingStream(
        get_logger("stderr"),
        logging.ERROR,
        passthrough_stream=sys.__stderr__,
    )

    _logging_session = LoggingSession(
        log_path=log_path,
        original_stdout=original_stdout,
        original_stderr=original_stderr,
        handlers=handlers,
    )
    atexit.register(_logging_session.flush)

    log_startup_metadata(
        command=command,
        args=args,
        seed_info=seed_info,
        log_path=log_path,
    )
    return _logging_session


def log_startup_metadata(
    *,
    command: Optional[str],
    args: Any,
    seed_info: Optional[Dict[str, Any]],
    log_path: str,
) -> None:
    """Log enough run metadata to reproduce the invocation."""
    metadata_payload = {
        "event": "run_start",
        "command": command or _infer_command_name(),
        "software": {
            "package": "gatk-sv-gd",
            "version": _dependency_version("gatk-sv-gd"),
            "git": _git_metadata(),
        },
        "python": {
            "version": platform.python_version(),
        },
        "dependencies": _dependency_versions(),
        "invocation": _privacy_safe_invocation(args),
        "environment": _environment_metadata(),
        "random_seeds": seed_info or {},
    }
    get_logger().info("run_start", extra={"fields": metadata_payload})


def _sanitize_log_text(message: str) -> str:
    if not message:
        return message
    message = _PATH_OR_URI_RE.sub("<redacted-path>", str(message))
    return _FILENAME_RE.sub("<redacted-file>", message)


def _sanitize_log_fields(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _sanitize_log_fields(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize_log_fields(item) for item in value]
    if isinstance(value, tuple):
        return [_sanitize_log_fields(item) for item in value]
    if isinstance(value, str):
        return _sanitize_log_text(value)
    return value


def _privacy_safe_invocation(args: Any) -> Dict[str, Any]:
    return {
        "argument_count": max(len(sys.argv) - 1, 0),
        "flags": _argv_flags(sys.argv),
        "parsed_args": _privacy_safe_args(vars(args)) if args is not None else None,
    }


def _argv_flags(argv) -> list:
    return sorted({arg.split("=", 1)[0] for arg in argv[1:] if str(arg).startswith("--")})


def _privacy_safe_args(args_dict: Dict[str, Any]) -> Dict[str, Any]:
    return {
        str(key): _privacy_safe_arg_value(str(key), value)
        for key, value in sorted(args_dict.items())
    }


def _privacy_safe_arg_value(key: str, value: Any) -> Any:
    if value is None or isinstance(value, bool):
        return value
    if _is_sensitive_arg_key(key):
        return _redacted_summary(value)
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, (list, tuple, set)):
        return {"count": len(value)}
    if isinstance(value, str):
        if _PATH_OR_URI_RE.search(value) or _FILENAME_RE.search(value):
            return "<redacted>"
        return value
    return str(type(value).__name__)


def _is_sensitive_arg_key(key: str) -> bool:
    lowered = key.lower()
    return any(token in lowered for token in _SENSITIVE_ARG_KEYWORDS)


def _redacted_summary(value: Any) -> Any:
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, (list, tuple, set)):
        return {"count": len(value), "value": "<redacted>"}
    if isinstance(value, dict):
        return {"count": len(value), "value": "<redacted>"}
    return "<redacted>"


def _remove_managed_handlers(logger: logging.Logger) -> None:
    for handler in list(logger.handlers):
        if getattr(handler, _MANAGED_HANDLER_ATTR, False):
            logger.removeHandler(handler)
            handler.close()


def _infer_command_name() -> str:
    return "gatk-sv-gd"


def _dependency_version(distribution_name: str) -> str:
    try:
        return importlib_metadata.version(distribution_name)
    except importlib_metadata.PackageNotFoundError:
        return "unknown"


def _dependency_versions() -> Dict[str, str]:
    return {
        distribution_name: _dependency_version(distribution_name)
        for distribution_name in (
            "numpy",
            "pandas",
            "pysam",
            "torch",
            "pyro-ppl",
            "matplotlib",
            "intervaltree",
        )
    }


def _git_metadata() -> Dict[str, str]:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
    commit = _run_git_command(repo_root, ["rev-parse", "HEAD"])
    dirty_status = _run_git_command(repo_root, ["status", "--porcelain"])
    return {
        "commit": commit or "unknown",
        "dirty": "unknown" if dirty_status is None else str(bool(dirty_status)),
    }


def _run_git_command(repo_root: str, command_args) -> Optional[str]:
    try:
        result = subprocess.run(
            ["git", *command_args],
            cwd=repo_root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def _environment_metadata() -> Dict[str, Any]:
    env_keys = [
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "PYTHONHASHSEED",
    ]
    return {
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "env": {key: os.environ.get(key) for key in env_keys if key in os.environ},
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def get_sample_columns(df: pd.DataFrame) -> list:
    """
    Extract sample column names from a DataFrame by excluding metadata columns.

    Args:
        df: DataFrame with bins as rows and samples as columns

    Returns:
        List of sample column names
    """
    metadata_cols = ["Chr", "Start", "End", "source_file"]
    sample_cols = [col for col in df.columns if col not in metadata_cols]
    return sample_cols
