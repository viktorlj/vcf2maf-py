"""Parse VEP CSQ and SnpEff ANN annotations from VCF INFO fields.

Provides transcript-effect objects and logic for selecting the single
best effect per variant (biotype priority → effect severity → transcript
length → canonical flag).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Sequence

from vcf2maf_py.constants import (
    BIOTYPE_PRIORITY,
    DEFAULT_BIOTYPE_PRIORITY,
    DEFAULT_EFFECT_PRIORITY,
    EFFECT_PRIORITY,
    SNPEFF_TO_VEP_FIELDS,
    VEP_AF_COLUMN_REMAP,
)
from vcf2maf_py.vcf_reader import VCFHeader, VCFRecord

logger = logging.getLogger("vcf2maf_py")


# ---------------------------------------------------------------------------
# TranscriptEffect — one annotation entry (one transcript for one allele)
# ---------------------------------------------------------------------------

@dataclass
class TranscriptEffect:
    """A single transcript-level annotation for a variant allele."""

    allele: str = ""
    consequence: str = ""  # SO terms, "&"-separated
    impact: str = ""
    symbol: str = ""
    gene: str = ""
    feature_type: str = ""
    feature: str = ""  # transcript ID
    biotype: str = ""
    exon: str = ""
    intron: str = ""
    hgvsc: str = ""
    hgvsp: str = ""
    cdna_position: str = ""
    cds_position: str = ""
    protein_position: str = ""
    amino_acids: str = ""
    codons: str = ""
    existing_variation: str = ""
    distance: str = ""
    strand: str = ""
    canonical: str = ""
    sift: str = ""
    polyphen: str = ""
    # All raw fields preserved
    raw: dict[str, str] = field(default_factory=dict)

    @property
    def consequences(self) -> list[str]:
        """Split the consequence string into individual SO terms."""
        return self.consequence.replace(",", "&").split("&") if self.consequence else []

    @property
    def worst_consequence(self) -> str:
        """Return the most severe SO term from the consequence field."""
        terms = self.consequences
        if not terms:
            return ""
        return min(terms, key=lambda t: EFFECT_PRIORITY.get(t, DEFAULT_EFFECT_PRIORITY))

    @property
    def effect_priority(self) -> int:
        """Numeric severity of the worst consequence (lower = more severe)."""
        return EFFECT_PRIORITY.get(self.worst_consequence, DEFAULT_EFFECT_PRIORITY)

    @property
    def biotype_priority(self) -> int:
        """Numeric biotype rank (lower = preferred)."""
        return BIOTYPE_PRIORITY.get(self.biotype, DEFAULT_BIOTYPE_PRIORITY)

    @property
    def is_canonical(self) -> bool:
        return self.canonical.upper() == "YES"


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _unescape_vep(value: str) -> str:
    """VEP replaces commas with '&' inside CSQ values; undo that."""
    return value.replace("&", ",") if value else value


def parse_vep_csq(
    info: dict[str, str],
    field_names: list[str],
    alt_alleles: list[str] | None = None,
    target_allele: str | None = None,
) -> list[TranscriptEffect]:
    """Parse VEP CSQ entries from an INFO dict.

    Parameters
    ----------
    info : parsed INFO dict
    field_names : list of CSQ sub-field names (from the VCF header)
    alt_alleles : ALT alleles of the VCF record (for matching ALLELE_NUM)
    target_allele : if set, only return effects for this ALT allele
    """
    csq_str = info.get("CSQ") or info.get("vep") or ""
    if not csq_str:
        return []

    effects: list[TranscriptEffect] = []
    allele_num_idx = None
    if "ALLELE_NUM" in field_names:
        allele_num_idx = field_names.index("ALLELE_NUM")

    for entry in csq_str.split(","):
        values = entry.split("|")
        raw: dict[str, str] = {}
        for i, name in enumerate(field_names):
            raw[name] = values[i] if i < len(values) else ""

        # Apply old VEP AF column remapping
        for old_name, new_name in VEP_AF_COLUMN_REMAP.items():
            if old_name in raw and new_name not in raw:
                raw[new_name] = raw[old_name]

        # Filter by target allele if requested
        entry_allele = raw.get("Allele", "")
        if allele_num_idx is not None and alt_alleles and target_allele:
            allele_num = raw.get("ALLELE_NUM", "")
            if allele_num and allele_num.isdigit():
                idx = int(allele_num) - 1
                if 0 <= idx < len(alt_alleles) and alt_alleles[idx] != target_allele:
                    continue
        elif target_allele and entry_allele and entry_allele != target_allele:
            # Fall back to matching by allele string
            # VEP sometimes uses a shortened allele for indels (just the inserted bases)
            if target_allele not in (entry_allele, "-"):
                continue

        effect = TranscriptEffect(
            allele=entry_allele,
            consequence=raw.get("Consequence", ""),
            impact=raw.get("IMPACT", ""),
            symbol=raw.get("SYMBOL", "") or raw.get("SYMBOL_SOURCE", ""),
            gene=raw.get("Gene", ""),
            feature_type=raw.get("Feature_type", ""),
            feature=raw.get("Feature", ""),
            biotype=raw.get("BIOTYPE", ""),
            exon=raw.get("EXON", ""),
            intron=raw.get("INTRON", ""),
            hgvsc=raw.get("HGVSc", ""),
            hgvsp=raw.get("HGVSp", ""),
            cdna_position=raw.get("cDNA_position", ""),
            cds_position=raw.get("CDS_position", ""),
            protein_position=raw.get("Protein_position", ""),
            amino_acids=raw.get("Amino_acids", ""),
            codons=raw.get("Codons", ""),
            existing_variation=raw.get("Existing_variation", ""),
            distance=raw.get("DISTANCE", ""),
            strand=raw.get("STRAND", ""),
            canonical=raw.get("CANONICAL", ""),
            sift=raw.get("SIFT", ""),
            polyphen=raw.get("PolyPhen", ""),
            raw=raw,
        )
        effects.append(effect)

    return effects


def parse_snpeff_ann(
    info: dict[str, str],
    field_names: list[str] | None = None,
    target_allele: str | None = None,
) -> list[TranscriptEffect]:
    """Parse SnpEff ANN entries from an INFO dict."""
    ann_str = info.get("ANN", "")
    if not ann_str:
        return []

    # Default SnpEff ANN field order
    if field_names is None:
        field_names = [
            "Allele", "Annotation", "Annotation_Impact", "Gene_Name",
            "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
            "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length",
            "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance",
            "ERRORS / WARNINGS / INFO",
        ]

    effects: list[TranscriptEffect] = []
    for entry in ann_str.split(","):
        values = entry.split("|")
        raw: dict[str, str] = {}
        for i, name in enumerate(field_names):
            raw[name] = values[i] if i < len(values) else ""

        # Map SnpEff fields to VEP-like names
        mapped: dict[str, str] = {}
        for snpeff_name, vep_name in SNPEFF_TO_VEP_FIELDS.items():
            if snpeff_name in raw:
                mapped[vep_name] = raw[snpeff_name]

        # Filter by allele
        entry_allele = raw.get("Allele", "")
        if target_allele and entry_allele and entry_allele != target_allele:
            continue

        # Extract position part from "pos/length" fields
        cdna_pos = raw.get("cDNA.pos / cDNA.length", "")
        cds_pos = raw.get("CDS.pos / CDS.length", "")
        aa_pos = raw.get("AA.pos / AA.length", "")

        effect = TranscriptEffect(
            allele=entry_allele,
            consequence=mapped.get("Consequence", ""),
            impact=mapped.get("IMPACT", ""),
            symbol=mapped.get("SYMBOL", ""),
            gene=mapped.get("Gene", ""),
            feature_type=mapped.get("Feature_type", ""),
            feature=mapped.get("Feature", ""),
            biotype=mapped.get("BIOTYPE", ""),
            exon=mapped.get("EXON", ""),
            hgvsc=mapped.get("HGVSc", ""),
            hgvsp=mapped.get("HGVSp", ""),
            cdna_position=cdna_pos.split("/")[0] if "/" in cdna_pos else cdna_pos,
            cds_position=cds_pos.split("/")[0] if "/" in cds_pos else cds_pos,
            protein_position=aa_pos.split("/")[0] if "/" in aa_pos else aa_pos,
            distance=mapped.get("DISTANCE", ""),
            raw=raw,
        )
        effects.append(effect)

    return effects


def get_effects_for_record(
    record: VCFRecord,
    header: VCFHeader,
    target_allele: str | None = None,
) -> list[TranscriptEffect]:
    """Extract transcript effects from a VCF record, auto-detecting annotation source."""
    source = header.annotation_source

    if source == "VEP":
        field_names = header.csq_fields
        if field_names:
            return parse_vep_csq(
                record.info, field_names,
                alt_alleles=record.alt,
                target_allele=target_allele,
            )
    elif source == "SnpEff":
        field_names = header.ann_fields
        return parse_snpeff_ann(
            record.info,
            field_names=field_names,
            target_allele=target_allele,
        )

    return []


# ---------------------------------------------------------------------------
# Effect selection — pick the single "best" annotation per variant
# ---------------------------------------------------------------------------

def select_best_effect(
    effects: list[TranscriptEffect],
    isoform_overrides: dict[str, str] | None = None,
) -> TranscriptEffect | None:
    """Select the single best transcript effect using the priority cascade.

    Priority order (from mskcc/vcf2maf):
    1. Prefer custom isoform override transcript for the worst-affected gene
    2. Prefer VEP canonical transcript for the worst-affected gene
    3. Fall back to worst-affected custom ENST, then worst canonical, then worst overall
    Within each category: biotype priority → effect severity → canonical flag → feature length
    """
    if not effects:
        return None
    if len(effects) == 1:
        return effects[0]

    override_set: set[str] = set()
    if isoform_overrides:
        override_set = set(isoform_overrides.values())

    def sort_key(e: TranscriptEffect) -> tuple:
        return (
            e.biotype_priority,
            e.effect_priority,
            0 if e.is_canonical else 1,
            0 if e.feature in override_set else 1,
            # Prefer longer transcript IDs as a proxy for longer transcripts
            -len(e.feature),
        )

    # Group by gene symbol to find worst-affected gene
    sorted_effects = sorted(effects, key=sort_key)

    # First choice: custom override transcript for the worst-affected gene
    if override_set:
        worst_gene = sorted_effects[0].symbol
        for e in sorted_effects:
            if e.symbol == worst_gene and e.feature in override_set:
                return e

    # Second choice: canonical transcript for the worst-affected gene
    worst_gene = sorted_effects[0].symbol
    for e in sorted_effects:
        if e.symbol == worst_gene and e.is_canonical:
            return e

    # Third choice: any custom override, then any canonical, then overall worst
    for e in sorted_effects:
        if e.feature in override_set:
            return e
    for e in sorted_effects:
        if e.is_canonical:
            return e

    return sorted_effects[0]


def format_all_effects(effects: list[TranscriptEffect]) -> str:
    """Format all transcript effects as a semicolon-delimited summary string.

    Each entry: SYMBOL,Consequence,HGVSp_Short,Transcript,RefSeq
    """
    parts = []
    for e in effects:
        from vcf2maf_py.utils import hgvsp_short
        p_short = hgvsp_short(e.hgvsp)
        refseq = e.raw.get("RefSeq", "")
        parts.append(
            f"{e.symbol},{e.worst_consequence},{p_short},{e.feature},{refseq}"
        )
    return ";".join(parts)
