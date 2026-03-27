"""Core VCF → MAF conversion engine."""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import IO, TextIO

from vcf2maf_py.annotation import (
    TranscriptEffect,
    format_all_effects,
    get_effects_for_record,
    select_best_effect,
)
from vcf2maf_py.constants import (
    ALL_MAF_COLUMNS,
    CORE_MAF_COLUMNS,
    EFFECT_PRIORITY,
    EXTENDED_MAF_COLUMNS,
    SO_TO_MAF,
    VEP_MAF_COLUMNS,
)
from vcf2maf_py.utils import hgvsp_short, open_file, strip_hgvs_prefix
from vcf2maf_py.vcf_reader import VCFHeader, VCFReader, VCFRecord

logger = logging.getLogger("vcf2maf_py")


# ---------------------------------------------------------------------------
# Allele manipulation
# ---------------------------------------------------------------------------

def trim_common_bases(ref: str, alt: str, pos: int) -> tuple[str, str, int]:
    """Remove shared prefix and suffix bases from REF/ALT, adjusting POS.

    VCF pads indels with a reference base; MAF does not.
    Returns (trimmed_ref, trimmed_alt, adjusted_pos).
    """
    # Strip common prefix
    prefix_len = 0
    for r, a in zip(ref, alt):
        if r == a:
            prefix_len += 1
        else:
            break
    if prefix_len > 0:
        ref = ref[prefix_len:]
        alt = alt[prefix_len:]
        pos += prefix_len

    # Strip common suffix
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    return ref or "", alt or "", pos


def determine_variant_type(ref: str, alt: str) -> str:
    """Determine MAF Variant_Type from trimmed REF/ALT alleles."""
    ref_len = len(ref) if ref else 0
    alt_len = len(alt) if alt else 0

    if ref_len == 0 and alt_len > 0:
        return "INS"
    if ref_len > 0 and alt_len == 0:
        return "DEL"
    if ref_len == alt_len:
        if ref_len == 1:
            return "SNP"
        if ref_len == 2:
            return "DNP"
        if ref_len == 3:
            return "TNP"
        return "ONP"
    # Complex substitution — classify by net length change (matches Perl vcf2maf)
    if ref_len > alt_len:
        return "DEL"
    return "INS"


def get_maf_coordinates(
    vcf_pos: int, vcf_ref: str, vcf_alt: str
) -> tuple[str, str, int, int, str]:
    """Convert VCF alleles to MAF representation.

    Returns (maf_ref, maf_alt, start_pos, end_pos, variant_type).
    """
    ref, alt, pos = trim_common_bases(vcf_ref, vcf_alt, vcf_pos)
    variant_type = determine_variant_type(ref, alt)

    if variant_type == "INS":
        # Insertion: start = base before, end = base after
        maf_ref = "-"
        maf_alt = alt
        start = pos - 1
        end = pos
    elif variant_type == "DEL":
        # Pure deletion
        maf_ref = ref
        maf_alt = "-"
        start = pos
        end = pos + len(ref) - 1
    else:
        # SNP, DNP, TNP, ONP (including complex substitutions)
        maf_ref = ref if ref else "-"
        maf_alt = alt if alt else "-"
        start = pos
        end = pos + max(len(ref), 1) - 1

    return maf_ref, maf_alt, start, end, variant_type


def get_variant_classification(
    effect: TranscriptEffect | None, variant_type: str
) -> str:
    """Map a TranscriptEffect to a MAF Variant_Classification."""
    if effect is None:
        return ""

    for so_term in sorted(
        effect.consequences,
        key=lambda t: EFFECT_PRIORITY.get(t, 99),
    ):
        classification = SO_TO_MAF.get(so_term)
        if classification is not None:
            return classification
        # Handle terms that depend on variant type
        if so_term == "frameshift_variant":
            return "Frame_Shift_Del" if variant_type == "DEL" else "Frame_Shift_Ins"
        if so_term == "protein_altering_variant":
            if variant_type == "DEL":
                return "In_Frame_Del"
            if variant_type == "INS":
                return "In_Frame_Ins"
            return "Missense_Mutation"

    return "Targeted_Region"


# ---------------------------------------------------------------------------
# Genotype / depth extraction
# ---------------------------------------------------------------------------

def extract_depths(
    sample_data: dict[str, str],
    ref_allele: str,
    alt_allele: str,
    alt_index: int = 1,
) -> tuple[int, int, int]:
    """Extract (depth, ref_count, alt_count) from sample FORMAT fields.

    Handles multiple caller conventions (GATK, FreeBayes, Strelka, VarScan, etc.).
    Returns -1 for unknown values.
    """
    depth = ref_count = alt_count = -1

    # --- Method 1: AD (allelic depths) — most common (GATK, FreeBayes, etc.) ---
    ad = sample_data.get("AD", "")
    if ad and ad != ".":
        parts = ad.split(",")
        if len(parts) > alt_index:
            try:
                ref_count = int(parts[0]) if parts[0] != "." else -1
                alt_count = int(parts[alt_index]) if parts[alt_index] != "." else -1
            except ValueError:
                pass

    # --- Method 2: RO/AO (FreeBayes alternative) ---
    if ref_count == -1:
        ro = sample_data.get("RO", "")
        if ro and ro != ".":
            try:
                ref_count = int(ro)
            except ValueError:
                pass
    if alt_count == -1:
        ao = sample_data.get("AO", "")
        if ao and ao != ".":
            try:
                parts = ao.split(",")
                alt_count = int(parts[alt_index - 1]) if len(parts) >= alt_index else int(parts[0])
            except (ValueError, IndexError):
                pass

    # --- Method 3: Strelka SNV tier counts ({A,C,G,T}U) ---
    if alt_count == -1 and len(alt_allele) == 1 and alt_allele + "U" in sample_data:
        au = sample_data.get(alt_allele + "U", "")
        if au and au != ".":
            try:
                alt_count = int(au.split(",")[0])  # First tier
            except (ValueError, IndexError):
                pass
    if ref_count == -1 and len(ref_allele) == 1 and ref_allele + "U" in sample_data:
        ru = sample_data.get(ref_allele + "U", "")
        if ru and ru != ".":
            try:
                ref_count = int(ru.split(",")[0])
            except (ValueError, IndexError):
                pass

    # --- Method 4: Strelka indel (TIR/TAR) ---
    if alt_count == -1:
        tir = sample_data.get("TIR", "")
        if tir and tir != ".":
            try:
                alt_count = int(tir.split(",")[0])
            except (ValueError, IndexError):
                pass
    if ref_count == -1:
        tar = sample_data.get("TAR", "")
        if tar and tar != ".":
            try:
                ref_count = int(tar.split(",")[0])
            except (ValueError, IndexError):
                pass

    # --- Method 5: NR/NV ---
    if ref_count == -1 and alt_count == -1:
        nr = sample_data.get("NR", "")
        nv = sample_data.get("NV", "")
        if nr and nr != "." and nv and nv != ".":
            try:
                depth = int(nr)
                alt_count = int(nv)
                ref_count = depth - alt_count
            except ValueError:
                pass

    # --- Method 6: DP4 (ref-fwd, ref-rev, alt-fwd, alt-rev) ---
    if ref_count == -1 or alt_count == -1:
        dp4 = sample_data.get("DP4", "")
        if dp4 and dp4 != ".":
            parts = dp4.split(",")
            if len(parts) == 4:
                try:
                    ref_count = int(parts[0]) + int(parts[1])
                    alt_count = int(parts[2]) + int(parts[3])
                except ValueError:
                    pass

    # --- Method 7: DV (variant depth, e.g. bcftools) ---
    if alt_count == -1:
        dv = sample_data.get("DV", "")
        if dv and dv != ".":
            try:
                alt_count = int(dv)
            except ValueError:
                pass

    # --- DP (total depth) ---
    if depth == -1:
        dp = sample_data.get("DP", "")
        if dp and dp != ".":
            try:
                depth = int(dp)
            except ValueError:
                pass

    # Infer missing values
    if depth == -1 and ref_count >= 0 and alt_count >= 0:
        depth = ref_count + alt_count
    if ref_count == -1 and depth >= 0 and alt_count >= 0:
        ref_count = depth - alt_count
    if alt_count == -1 and depth >= 0 and ref_count >= 0:
        alt_count = depth - ref_count

    return depth, ref_count, alt_count


def get_genotype_alleles(
    sample_data: dict[str, str],
    ref: str,
    alts: list[str],
) -> tuple[str, str]:
    """Extract allele1, allele2 from a sample's GT field."""
    gt = sample_data.get("GT", "")
    if not gt or gt == ".":
        return ref, ref

    sep = "/" if "/" in gt else "|"
    parts = gt.split(sep)
    alleles = [ref] + alts

    def resolve(idx_str: str) -> str:
        if idx_str == ".":
            return ref
        try:
            idx = int(idx_str)
            return alleles[idx] if idx < len(alleles) else ref
        except ValueError:
            return ref

    a1 = resolve(parts[0])
    a2 = resolve(parts[1]) if len(parts) > 1 else a1
    return a1, a2


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def _build_maf_row(
    record: VCFRecord,
    alt_allele: str,
    alt_index: int,
    effect: TranscriptEffect | None,
    all_effects: list[TranscriptEffect],
    tumor_id: str,
    normal_id: str,
    tumor_idx: int | None,
    normal_idx: int | None,
    center: str,
    ncbi_build: str,
    output_columns: list[str],
) -> dict[str, str]:
    """Build a single MAF row dict for one variant allele."""
    maf_ref, maf_alt, start, end, variant_type = get_maf_coordinates(
        record.pos, record.ref, alt_allele
    )
    classification = get_variant_classification(effect, variant_type)

    # Extract depths
    t_depth = t_ref = t_alt = -1
    n_depth = n_ref = n_alt = -1

    trimmed_ref, trimmed_alt, _ = trim_common_bases(record.ref, alt_allele, record.pos)

    if tumor_idx is not None:
        t_data = record.sample_data(tumor_idx)
        t_depth, t_ref, t_alt = extract_depths(t_data, record.ref, alt_allele, alt_index)
        t_allele1_vcf, t_allele2_vcf = get_genotype_alleles(t_data, record.ref, record.alt)
    else:
        t_allele1_vcf, t_allele2_vcf = record.ref, alt_allele

    if normal_idx is not None:
        n_data = record.sample_data(normal_idx)
        n_depth, n_ref, n_alt = extract_depths(n_data, record.ref, alt_allele, alt_index)
        n_allele1_vcf, n_allele2_vcf = get_genotype_alleles(n_data, record.ref, record.alt)
    else:
        n_allele1_vcf, n_allele2_vcf = record.ref, record.ref

    # Tumor_Seq_Allele1: reference allele in MAF coords, unless homozygous alt
    is_hom_alt = (t_allele1_vcf != record.ref and t_allele2_vcf != record.ref)
    t_allele1_maf = maf_alt if is_hom_alt else maf_ref

    # Build dbSNP RS from existing variation
    dbsnp_rs = ""
    if effect:
        ev = effect.existing_variation or effect.raw.get("Existing_variation", "")
        rs_ids = [v for v in ev.replace("&", ",").split(",") if v.startswith("rs")]
        if rs_ids:
            dbsnp_rs = ";".join(rs_ids)
        elif ev:
            # Non-rs entries exist (e.g. COSMIC) — leave dbSNP_RS empty
            dbsnp_rs = ""
        else:
            dbsnp_rs = "novel"
    else:
        dbsnp_rs = record.id if record.id != "." else "novel"

    # Gene symbol and Entrez ID
    hugo = effect.symbol if effect else ""
    entrez = effect.raw.get("ENTREZ", "0") if effect else "0"
    if not entrez or entrez == "":
        entrez = "0"

    row: dict[str, str] = {
        "Hugo_Symbol": hugo,
        "Entrez_Gene_Id": entrez,
        "Center": center,
        "NCBI_Build": ncbi_build,
        "Chromosome": record.chrom,
        "Start_Position": str(start),
        "End_Position": str(end),
        "Strand": "+",
        "Variant_Classification": classification,
        "Variant_Type": variant_type,
        "Reference_Allele": maf_ref,
        "Tumor_Seq_Allele1": t_allele1_maf,
        "Tumor_Seq_Allele2": maf_alt,
        "dbSNP_RS": dbsnp_rs,
        "dbSNP_Val_Status": "",
        "Tumor_Sample_Barcode": tumor_id,
        "Matched_Norm_Sample_Barcode": normal_id,
        "Match_Norm_Seq_Allele1": maf_ref,
        "Match_Norm_Seq_Allele2": maf_ref,
        "Tumor_Validation_Allele1": "",
        "Tumor_Validation_Allele2": "",
        "Match_Norm_Validation_Allele1": "",
        "Match_Norm_Validation_Allele2": "",
        "Verification_Status": "",
        "Validation_Status": "",
        "Mutation_Status": "",
        "Sequencing_Phase": "",
        "Sequence_Source": "",
        "Validation_Method": "",
        "Score": "",
        "BAM_File": "",
        "Sequencer": "",
        "Tumor_Sample_UUID": "",
        "Matched_Norm_Sample_UUID": "",
    }

    # Extended columns — strip transcript/protein ID prefix from HGVS
    row["HGVSc"] = strip_hgvs_prefix(effect.hgvsc) if effect else ""
    row["HGVSp"] = strip_hgvs_prefix(effect.hgvsp) if effect else ""
    row["HGVSp_Short"] = hgvsp_short(strip_hgvs_prefix(effect.hgvsp)) if effect else ""
    row["Transcript_ID"] = effect.feature if effect else ""
    row["Exon_Number"] = effect.exon if effect else ""
    row["t_depth"] = str(t_depth) if t_depth >= 0 else ""
    row["t_ref_count"] = str(t_ref) if t_ref >= 0 else ""
    row["t_alt_count"] = str(t_alt) if t_alt >= 0 else ""
    row["n_depth"] = str(n_depth) if n_depth >= 0 else ""
    row["n_ref_count"] = str(n_ref) if n_ref >= 0 else ""
    row["n_alt_count"] = str(n_alt) if n_alt >= 0 else ""
    row["all_effects"] = format_all_effects(all_effects) if all_effects else ""

    # VEP fields where "&" separator should be converted to ","
    _VEP_MULTI_VALUE_FIELDS = {
        "Consequence", "CLIN_SIG", "SOMATIC", "PUBMED", "PHENO",
        "Existing_variation", "DOMAINS",
    }

    # VEP columns (pass through from annotation)
    if effect:
        raw = effect.raw
        for col in VEP_MAF_COLUMNS:
            if col not in row:
                # Map column names to VEP field names
                vep_key = col
                if col == "STRAND_VEP":
                    vep_key = "STRAND"
                val = raw.get(vep_key, "")
                if col in _VEP_MULTI_VALUE_FIELDS and "&" in val:
                    val = val.replace("&", ",")
                row[col] = val

    # Preserve original VCF fields
    row["vcf_id"] = record.id
    row["vcf_qual"] = record.qual
    row["vcf_pos"] = str(record.pos)
    row["vcf_ref"] = record.ref
    row["vcf_alt"] = alt_allele
    row["FILTER"] = record.filter

    return row


def convert_vcf_to_maf(
    input_vcf: str,
    output_maf: str | None = None,
    tumor_id: str | None = None,
    normal_id: str | None = None,
    ncbi_build: str = "GRCh38",
    center: str = ".",
    min_hom_vaf: float = 0.7,
    isoform_overrides: dict[str, str] | None = None,
    output_columns: list[str] | None = None,
) -> list[dict[str, str]]:
    """Convert a VCF file to MAF format.

    Parameters
    ----------
    input_vcf : Path to input VCF file (plain or gzipped).
    output_maf : Path to output MAF file. If None, results are only returned.
    tumor_id : Tumor sample ID. Auto-detected from VCF if not provided.
    normal_id : Normal sample ID. Auto-detected from VCF if not provided.
    ncbi_build : Reference genome build (default: GRCh38).
    center : Sequencing center name.
    min_hom_vaf : Minimum VAF to call homozygous when GT is missing.
    isoform_overrides : Dict mapping gene symbols to preferred transcript IDs.
    output_columns : List of MAF columns to include. Defaults to ALL_MAF_COLUMNS.

    Returns
    -------
    List of MAF row dicts.
    """
    reader = VCFReader(input_vcf)
    header = reader.header

    # Auto-detect annotation source
    ann_source = header.annotation_source
    if ann_source:
        logger.info("Detected %s annotations", ann_source)
    else:
        logger.warning(
            "No VEP (CSQ) or SnpEff (ANN) annotations found in VCF header. "
            "Variant_Classification and gene annotations will be empty. "
            "Consider annotating with VEP: "
            "https://www.ensembl.org/info/docs/tools/vep/"
        )

    # Determine sample indices
    samples = header.samples
    tumor_idx: int | None = None
    normal_idx: int | None = None

    if samples:
        if tumor_id:
            if tumor_id in samples:
                tumor_idx = samples.index(tumor_id)
            else:
                logger.warning("Tumor ID '%s' not found in VCF samples: %s", tumor_id, samples)
                tumor_idx = 0
                tumor_id = samples[0]
        else:
            tumor_idx = 0
            tumor_id = samples[0]

        if normal_id:
            if normal_id in samples:
                normal_idx = samples.index(normal_id)
            else:
                logger.warning("Normal ID '%s' not found in VCF samples: %s", normal_id, samples)
        elif len(samples) > 1:
            # Assume second sample is normal
            normal_idx = 1
            normal_id = samples[1]

    if not tumor_id:
        tumor_id = "TUMOR"
    if not normal_id:
        normal_id = "NORMAL"

    if output_columns is None:
        output_columns = ALL_MAF_COLUMNS

    rows: list[dict[str, str]] = []
    n_variants = 0
    n_skipped = 0

    for record in reader:
        if not record.alt:
            n_skipped += 1
            continue

        for alt_idx, alt_allele in enumerate(record.alt, start=1):
            if alt_allele == "*" or alt_allele == ".":
                continue

            n_variants += 1

            # Get all effects for this allele
            all_effects = get_effects_for_record(
                record, header, target_allele=alt_allele
            )
            best = select_best_effect(all_effects, isoform_overrides)

            row = _build_maf_row(
                record=record,
                alt_allele=alt_allele,
                alt_index=alt_idx,
                effect=best,
                all_effects=all_effects,
                tumor_id=tumor_id,
                normal_id=normal_id,
                tumor_idx=tumor_idx,
                normal_idx=normal_idx,
                center=center,
                ncbi_build=ncbi_build,
                output_columns=output_columns,
            )
            rows.append(row)

    logger.info(
        "Processed %d variant alleles (%d records skipped)",
        n_variants, n_skipped,
    )

    # Write output
    if output_maf:
        _write_maf(rows, output_maf, output_columns)
        logger.info("Wrote %d MAF records to %s", len(rows), output_maf)

    return rows


def _write_maf(
    rows: list[dict[str, str]],
    output_path: str,
    columns: list[str],
) -> None:
    """Write MAF rows to a file."""
    with open(output_path, "w", newline="") as fh:
        fh.write("#version 2.4\n")
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Isoform override loading
# ---------------------------------------------------------------------------

def load_isoform_overrides(path: str) -> dict[str, str]:
    """Load isoform override file (TSV: transcript_id<TAB>gene_symbol).

    Returns dict mapping gene_symbol → transcript_id.
    """
    overrides: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                transcript_id = parts[0].strip()
                gene_symbol = parts[1].strip()
                overrides[gene_symbol] = transcript_id
    return overrides
