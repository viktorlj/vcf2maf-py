"""Shared utility functions."""

from __future__ import annotations

import gzip
import logging
import re
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Iterator

from vcf2maf_py.constants import AA3_TO_AA1

logger = logging.getLogger("vcf2maf_py")


def setup_logging(verbosity: int = 0) -> None:
    """Configure logging for the package."""
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )


@contextmanager
def open_file(path: str | Path, mode: str = "r") -> Iterator[IO]:
    """Open a file, transparently handling gzip compression."""
    path = Path(path)
    if path.suffix == ".gz" or path.suffixes[-2:] == [".vcf", ".gz"]:
        fh = gzip.open(path, mode + "t", encoding="utf-8")
    else:
        fh = open(path, mode, encoding="utf-8")
    try:
        yield fh
    finally:
        fh.close()


def hgvsp_short(hgvsp: str) -> str:
    """Convert 3-letter amino acid HGVS notation to 1-letter.

    Example: p.Val600Glu → p.V600E
    """
    if not hgvsp:
        return ""
    result = hgvsp
    for aa3, aa1 in AA3_TO_AA1.items():
        result = result.replace(aa3, aa1)
    return result


_CHROM_ORDER_MAP: dict[str, int] = {}


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a sort key for chromosome names (1-22, X, Y, MT, then others)."""
    c = chrom.replace("chr", "").replace("Chr", "")
    if c.isdigit():
        return (0, c.zfill(2))
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if c.upper() in order:
        return (0, str(order[c.upper()]).zfill(2))
    return (1, c)


def parse_info_field(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO field string into a dictionary.

    Flags (no '=') are stored with value ''.
    """
    if info_str == "." or not info_str:
        return {}
    result: dict[str, str] = {}
    for entry in info_str.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            result[key] = value
        else:
            result[entry] = ""
    return result


def parse_format_sample(format_str: str, sample_str: str) -> dict[str, str]:
    """Parse VCF FORMAT and sample column into a dict.

    Example: parse_format_sample("GT:AD:DP", "0/1:30,15:45")
             → {"GT": "0/1", "AD": "30,15", "DP": "45"}
    """
    if format_str == "." or sample_str == ".":
        return {}
    keys = format_str.split(":")
    values = sample_str.split(":")
    return dict(zip(keys, values))
