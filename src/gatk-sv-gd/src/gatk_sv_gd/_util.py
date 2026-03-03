"""
Shared configuration and utility functions for gatk_sv_gd.

Module-level state (e.g. VERBOSE flag) lives here so that every sub-module
can import it without creating circular dependencies.
"""

import os
import sys

import numpy as np
import pandas as pd

# Module-level verbosity flag, set from --verbose in the infer CLI.
VERBOSE = False

# Module-level handle for log-only (verbose) output.
# Set by setup_call_logging(); written to directly by vlog() so that
# verbose diagnostics appear in the log file but NOT on stdout.
_log_fh = None


def vlog(msg: str) -> None:
    """Write *msg* to the log file only (never to stdout).

    Used for verbose diagnostic output so that it doesn't flood the
    terminal while still being fully captured in the run log.
    """
    if _log_fh is not None:
        _log_fh.write(msg + "\n")
        _log_fh.flush()


class TeeStream:
    """Write to both the original stream and a log file."""

    def __init__(self, original_stream, log_file):
        self.original_stream = original_stream
        self.log_file = log_file

    def write(self, message: str):
        self.original_stream.write(message)
        self.log_file.write(message)

    def flush(self):
        self.original_stream.flush()
        self.log_file.flush()

    # Forward any other attribute lookups to the original stream so that
    # callers relying on e.g. ``fileno()`` or ``isatty()`` still work.
    def __getattr__(self, name):
        return getattr(self.original_stream, name)


def setup_logging(output_dir: str, filename: str = "log.txt"):
    """Redirect stdout and stderr so output goes to both the console and a log file.

    Call once at the beginning of ``main()``.  The log file is opened in
    write mode so each run produces a fresh log.

    Args:
        output_dir: Directory where the log file will be created.
        filename: Name of the log file (default ``log.txt``).

    Returns:
        The open file handle (kept so the caller can close it if desired).
    """
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    sys.stdout = TeeStream(sys.__stdout__, log_fh)
    sys.stderr = TeeStream(sys.__stderr__, log_fh)
    return log_fh


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
