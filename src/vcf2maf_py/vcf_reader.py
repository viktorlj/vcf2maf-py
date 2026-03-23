"""Lightweight VCF parser — pure Python, no C dependencies.

Handles VCFv4.x files (plain text or gzipped) with full header parsing
and per-record INFO/FORMAT access.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator

from vcf2maf_py.utils import open_file, parse_format_sample, parse_info_field


# ---------------------------------------------------------------------------
# VCF header
# ---------------------------------------------------------------------------

@dataclass
class VCFInfoDef:
    """Parsed ##INFO header line."""
    id: str
    number: str
    type: str
    description: str


@dataclass
class VCFFormatDef:
    """Parsed ##FORMAT header line."""
    id: str
    number: str
    type: str
    description: str


@dataclass
class VCFHeader:
    """Parsed VCF header."""
    meta_lines: list[str] = field(default_factory=list)
    info_fields: dict[str, VCFInfoDef] = field(default_factory=dict)
    format_fields: dict[str, VCFFormatDef] = field(default_factory=dict)
    samples: list[str] = field(default_factory=list)
    contigs: list[str] = field(default_factory=list)
    column_header: str = ""

    @property
    def csq_fields(self) -> list[str] | None:
        """Extract VEP CSQ field names from the INFO header, or None."""
        for key in ("CSQ", "vep"):
            if key in self.info_fields:
                desc = self.info_fields[key].description
                match = re.search(r"Format:\s*(.+?)\"?\s*$", desc)
                if match:
                    return match.group(1).split("|")
        return None

    @property
    def ann_fields(self) -> list[str] | None:
        """Extract SnpEff ANN field names from the INFO header, or None."""
        if "ANN" in self.info_fields:
            desc = self.info_fields["ANN"].description
            # SnpEff description: "Functional annotations: 'Allele | Annotation | ...'"
            match = re.search(r"'([^']+)'", desc)
            if match:
                return [f.strip() for f in match.group(1).split("|")]
        return None

    @property
    def eff_fields(self) -> list[str] | None:
        """Extract legacy SnpEff EFF field names, or None."""
        if "EFF" in self.info_fields:
            return None  # EFF uses a different format; handled in annotation.py
        return None

    @property
    def annotation_source(self) -> str | None:
        """Detect which annotation source is present: 'VEP', 'SnpEff', or None."""
        if self.csq_fields is not None:
            return "VEP"
        if self.ann_fields is not None:
            return "SnpEff"
        if "ANN" in self.info_fields or "EFF" in self.info_fields:
            return "SnpEff"
        return None


# ---------------------------------------------------------------------------
# VCF record
# ---------------------------------------------------------------------------

@dataclass
class VCFRecord:
    """A single VCF data record."""
    chrom: str
    pos: int
    id: str
    ref: str
    alt: list[str]
    qual: str
    filter: str
    info_str: str
    format_str: str
    sample_strs: list[str]
    sample_names: list[str]

    _info: dict[str, str] | None = field(default=None, repr=False)

    @property
    def info(self) -> dict[str, str]:
        if self._info is None:
            self._info = parse_info_field(self.info_str)
        return self._info

    def sample_data(self, sample_index: int) -> dict[str, str]:
        """Parse FORMAT fields for a given sample index."""
        if sample_index >= len(self.sample_strs):
            return {}
        return parse_format_sample(self.format_str, self.sample_strs[sample_index])

    def sample_data_by_name(self, sample_name: str) -> dict[str, str]:
        """Parse FORMAT fields for a given sample name."""
        if sample_name not in self.sample_names:
            return {}
        idx = self.sample_names.index(sample_name)
        return self.sample_data(idx)

    @property
    def is_multiallelic(self) -> bool:
        return len(self.alt) > 1


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

_META_RE = re.compile(
    r"##(?P<key>\w+)=<(?P<content>.+)>"
)

_KV_RE = re.compile(
    r'(\w+)=(?:"([^"]*?)"|([^,>]*))'
)


def _parse_meta_line(line: str) -> tuple[str, dict[str, str]] | None:
    """Parse a structured ##KEY=<...> header line."""
    m = _META_RE.match(line)
    if not m:
        return None
    key = m.group("key")
    content = m.group("content")
    attrs: dict[str, str] = {}
    for km in _KV_RE.finditer(content):
        attrs[km.group(1)] = km.group(2) if km.group(2) is not None else km.group(3)
    return key, attrs


def _parse_record(line: str, sample_names: list[str]) -> VCFRecord:
    """Parse a VCF data line into a VCFRecord."""
    fields = line.rstrip("\n").split("\t")
    n_fixed = 8  # CHROM..INFO
    alt_str = fields[4]
    alts = alt_str.split(",") if alt_str != "." else []
    return VCFRecord(
        chrom=fields[0],
        pos=int(fields[1]),
        id=fields[2],
        ref=fields[3],
        alt=alts,
        qual=fields[5],
        filter=fields[6],
        info_str=fields[7] if len(fields) > 7 else ".",
        format_str=fields[8] if len(fields) > 9 else ".",
        sample_strs=fields[9:] if len(fields) > 9 else [],
        sample_names=sample_names,
    )


# ---------------------------------------------------------------------------
# Main reader
# ---------------------------------------------------------------------------

class VCFReader:
    """Iterate over records in a VCF file."""

    def __init__(self, path: str):
        self.path = path
        self.header = VCFHeader()
        self._parse_header()

    def _parse_header(self) -> None:
        with open_file(self.path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    self.header.meta_lines.append(line)
                    parsed = _parse_meta_line(line)
                    if parsed:
                        key, attrs = parsed
                        if key == "INFO" and "ID" in attrs:
                            self.header.info_fields[attrs["ID"]] = VCFInfoDef(
                                id=attrs["ID"],
                                number=attrs.get("Number", "."),
                                type=attrs.get("Type", "String"),
                                description=attrs.get("Description", ""),
                            )
                        elif key == "FORMAT" and "ID" in attrs:
                            self.header.format_fields[attrs["ID"]] = VCFFormatDef(
                                id=attrs["ID"],
                                number=attrs.get("Number", "."),
                                type=attrs.get("Type", "String"),
                                description=attrs.get("Description", ""),
                            )
                        elif key == "contig" and "ID" in attrs:
                            self.header.contigs.append(attrs["ID"])
                elif line.startswith("#CHROM"):
                    self.header.column_header = line
                    parts = line.split("\t")
                    if len(parts) > 9:
                        self.header.samples = parts[9:]
                    break

    def __iter__(self) -> Iterator[VCFRecord]:
        with open_file(self.path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                line = line.rstrip("\n")
                if not line:
                    continue
                yield _parse_record(line, self.header.samples)
