"""MAF → VCF reverse conversion."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from vcf2maf_py.utils import open_file

logger = logging.getLogger("vcf2maf_py")


def _read_maf(path: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a MAF file, skipping comment/version lines.

    Returns (column_names, list_of_row_dicts).
    """
    rows: list[dict[str, str]] = []
    columns: list[str] = []

    with open_file(path) as fh:
        # Skip comment/version lines
        header_line = ""
        for line in fh:
            if line.startswith("#"):
                continue
            header_line = line
            break

        if not header_line:
            return columns, rows

        columns = header_line.strip().split("\t")
        reader = csv.DictReader(fh, fieldnames=columns, delimiter="\t")
        for row in reader:
            rows.append(row)

    return columns, rows


def _maf_to_vcf_alleles(
    ref_allele: str,
    alt_allele: str,
    start_pos: int,
    variant_type: str,
    ref_fasta_func=None,
    chrom: str = "",
) -> tuple[int, str, str]:
    """Convert MAF alleles to VCF representation.

    For indels, adds a padding base. If ref_fasta_func is provided,
    it fetches the actual reference base; otherwise uses 'N'.

    Returns (vcf_pos, vcf_ref, vcf_alt).
    """
    if ref_allele == "-" or not ref_allele:
        # Insertion: MAF has ref="-", alt=inserted bases
        # VCF needs: POS = start_pos, REF = padding_base, ALT = padding_base + alt
        pad_pos = start_pos
        if ref_fasta_func:
            pad_base = ref_fasta_func(chrom, pad_pos, pad_pos)
        else:
            pad_base = "N"
        return pad_pos, pad_base, pad_base + alt_allele

    if alt_allele == "-" or not alt_allele:
        # Deletion: MAF has ref=deleted bases, alt="-"
        # VCF needs: POS = start_pos - 1, REF = padding_base + deleted, ALT = padding_base
        pad_pos = start_pos - 1
        if ref_fasta_func:
            pad_base = ref_fasta_func(chrom, pad_pos, pad_pos)
        else:
            pad_base = "N"
        return pad_pos, pad_base + ref_allele, pad_base

    # SNP/DNP/ONP: no padding needed
    return start_pos, ref_allele, alt_allele


def convert_maf_to_vcf(
    input_maf: str,
    output_vcf: str,
    ref_fasta: str | None = None,
    per_sample_vcfs: bool = False,
) -> None:
    """Convert a MAF file to VCF format.

    Parameters
    ----------
    input_maf : Path to input MAF file.
    output_vcf : Path to output VCF file.
    ref_fasta : Path to reference FASTA (optional, for correct padding bases).
    per_sample_vcfs : If True, also write per-sample VCF files.
    """
    columns, rows = _read_maf(input_maf)
    if not rows:
        logger.warning("No records found in MAF file")
        return

    # Check required columns
    required = {"Chromosome", "Start_Position", "Reference_Allele", "Tumor_Sample_Barcode"}
    alt_col = None
    for candidate in ("Tumor_Seq_Allele2", "Tumor_Seq_Allele1"):
        if candidate in columns:
            alt_col = candidate
            break
    if alt_col is None:
        raise ValueError("MAF must contain Tumor_Seq_Allele2 or Tumor_Seq_Allele1")

    missing = required - set(columns)
    if missing:
        raise ValueError(f"MAF missing required columns: {missing}")

    # Optional reference FASTA via pysam
    ref_fasta_func = None
    if ref_fasta:
        try:
            import pysam
            fasta = pysam.FastaFile(ref_fasta)

            def _fetch_base(chrom: str, start: int, end: int) -> str:
                try:
                    return fasta.fetch(chrom, start - 1, end).upper()
                except (ValueError, KeyError):
                    # Try with/without "chr" prefix
                    alt_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else "chr" + chrom
                    try:
                        return fasta.fetch(alt_chrom, start - 1, end).upper()
                    except (ValueError, KeyError):
                        return "N"

            ref_fasta_func = _fetch_base
            logger.info("Using reference FASTA: %s", ref_fasta)
        except ImportError:
            logger.warning(
                "pysam not installed — using 'N' as padding base for indels. "
                "Install with: pip install 'vcf2maf-py[pysam]'"
            )

    # Collect all unique samples and tumor-normal pairs
    all_samples: dict[str, int] = {}
    pairs_file_rows: list[tuple[str, str]] = []

    for row in rows:
        tumor = row.get("Tumor_Sample_Barcode", "TUMOR")
        normal = row.get("Matched_Norm_Sample_Barcode", "NORMAL")
        if tumor not in all_samples:
            all_samples[tumor] = len(all_samples)
        if normal and normal not in all_samples:
            all_samples[normal] = len(all_samples)
        pairs_file_rows.append((tumor, normal))

    sample_names = sorted(all_samples.keys(), key=lambda s: all_samples[s])

    # Build VCF records
    vcf_records: list[dict] = []
    skipped = 0

    for row in rows:
        chrom = row["Chromosome"]
        start = int(row["Start_Position"])
        ref_allele = row["Reference_Allele"]
        alt_allele = row.get(alt_col, "")
        variant_type = row.get("Variant_Type", "")
        tumor = row.get("Tumor_Sample_Barcode", "TUMOR")
        normal = row.get("Matched_Norm_Sample_Barcode", "NORMAL")

        if not alt_allele or alt_allele == ref_allele:
            skipped += 1
            continue

        vcf_pos, vcf_ref, vcf_alt = _maf_to_vcf_alleles(
            ref_allele, alt_allele, start, variant_type,
            ref_fasta_func=ref_fasta_func, chrom=chrom,
        )

        # Build genotype info
        t_depth = row.get("t_depth", ".")
        t_ref_count = row.get("t_ref_count", ".")
        t_alt_count = row.get("t_alt_count", ".")
        n_depth = row.get("n_depth", ".")
        n_ref_count = row.get("n_ref_count", ".")
        n_alt_count = row.get("n_alt_count", ".")

        def _make_gt_ad_dp(depth: str, ref_ct: str, alt_ct: str, is_tumor: bool) -> str:
            gt = "0/1" if is_tumor else "0/0"
            ad_parts = []
            for v in (ref_ct, alt_ct):
                ad_parts.append(v if v and v != "" else ".")
            ad = ",".join(ad_parts)
            dp = depth if depth and depth != "" else "."
            return f"{gt}:{ad}:{dp}"

        vcf_records.append({
            "chrom": chrom,
            "pos": vcf_pos,
            "id": row.get("dbSNP_RS", ".") if row.get("dbSNP_RS", "") not in ("", "novel") else ".",
            "ref": vcf_ref,
            "alt": vcf_alt,
            "qual": ".",
            "filter": row.get("FILTER", ".") or ".",
            "info": ".",
            "tumor": tumor,
            "normal": normal,
            "tumor_gt": _make_gt_ad_dp(t_depth, t_ref_count, t_alt_count, True),
            "normal_gt": _make_gt_ad_dp(n_depth, n_ref_count, n_alt_count, False),
        })

    # Write VCF
    _write_vcf(vcf_records, output_vcf, sample_names, rows[0].get("NCBI_Build", ""))

    logger.info(
        "Wrote %d VCF records to %s (%d skipped)",
        len(vcf_records), output_vcf, skipped,
    )

    # Write pairs file
    pairs_path = str(output_vcf).replace(".vcf", ".pairs.tsv")
    unique_pairs = list(dict.fromkeys(pairs_file_rows))
    with open(pairs_path, "w") as fh:
        for tumor, normal in unique_pairs:
            fh.write(f"{tumor}\t{normal}\n")


def _write_vcf(
    records: list[dict],
    output_path: str,
    sample_names: list[str],
    ncbi_build: str,
) -> None:
    """Write VCF records to a file."""
    with open(output_path, "w") as fh:
        # Header
        fh.write("##fileformat=VCFv4.2\n")
        if ncbi_build:
            fh.write(f"##reference={ncbi_build}\n")
        fh.write("##source=vcf2maf-py\n")
        fh.write('##INFO=<ID=.,Number=.,Type=String,Description="No INFO fields">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,'
            'Description="Allelic depths for the ref and alt alleles">\n'
        )
        fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')

        # Determine unique samples from records
        sample_set: dict[str, None] = {}
        for rec in records:
            sample_set[rec["tumor"]] = None
            if rec["normal"]:
                sample_set[rec["normal"]] = None
        samples = list(sample_set.keys())

        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in samples:
            fh.write(f"\t{s}")
        fh.write("\n")

        # Records
        for rec in records:
            fh.write(
                f"{rec['chrom']}\t{rec['pos']}\t{rec['id']}\t{rec['ref']}\t"
                f"{rec['alt']}\t{rec['qual']}\t{rec['filter']}\t.\tGT:AD:DP"
            )
            for s in samples:
                if s == rec["tumor"]:
                    fh.write(f"\t{rec['tumor_gt']}")
                elif s == rec["normal"]:
                    fh.write(f"\t{rec['normal_gt']}")
                else:
                    fh.write("\t./.:.:.")
            fh.write("\n")
