"""osmotool — screen osmoadaptation genes in metagenomic datasets using DIAMOND and HMMER."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("osmotool")
except PackageNotFoundError:
    __version__ = "unknown"
