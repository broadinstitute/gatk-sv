import os
from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files

__all__ = []

try:
    __version__ = version("svtk")
except PackageNotFoundError:  # running from a source tree with no installed metadata
    __version__ = "0.1"


def resource_path(*parts):
    """Absolute filesystem path to a file shipped inside the svtk package.

    Replaces ``pkg_resources.resource_filename('svtk', <relpath>)``. pkg_resources
    is vendored by setuptools and was removed in setuptools 82; the sv-pipeline
    image installs torch, which requires setuptools>=77.0.3, so the image cannot
    hold setuptools back at a version that still ships pkg_resources. svtk is
    imported by every ``svtk`` subcommand, so importing pkg_resources here (or in
    any module reachable from ``svtk.cli`` / ``svtk.utils``) breaks the tool
    outright rather than one code path.
    """
    return os.fspath(files("svtk").joinpath(*parts))
