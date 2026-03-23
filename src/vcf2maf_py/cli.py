"""Command-line interface for vcf2maf-py."""

from __future__ import annotations

import sys

import click

from vcf2maf_py import __version__


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="vcf2maf-py")
def main() -> None:
    """vcf2maf-py — Convert VCF files to MAF format and back.

    A Python reimplementation of mskcc/vcf2maf.
    """


@main.command()
@click.argument("input_vcf", type=click.Path(exists=True))
@click.option("-o", "--output-maf", type=click.Path(), default=None,
              help="Output MAF path (default: <input>.maf)")
@click.option("--tumor-id", default=None, help="Tumor sample ID (auto-detected from VCF)")
@click.option("--normal-id", default=None, help="Normal sample ID (auto-detected from VCF)")
@click.option("--ncbi-build", default="GRCh38",
              help="Reference genome build [default: GRCh38]")
@click.option("--center", default=".", help="Sequencing center name [default: .]")
@click.option("--min-hom-vaf", default=0.7, type=float,
              help="Min VAF to call homozygous when GT is absent [default: 0.7]")
@click.option("--isoform-overrides", type=click.Path(exists=True), default=None,
              help="TSV file with preferred transcript overrides (transcript_id<TAB>gene)")
@click.option("--output-columns", type=click.Choice(["core", "extended", "all"]),
              default="all", help="Which MAF columns to include [default: all]")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v info, -vv debug)")
def vcf2maf(
    input_vcf: str,
    output_maf: str | None,
    tumor_id: str | None,
    normal_id: str | None,
    ncbi_build: str,
    center: str,
    min_hom_vaf: float,
    isoform_overrides: str | None,
    output_columns: str,
    verbose: int,
) -> None:
    """Convert a VCF file to MAF format.

    Reads a VCF file (optionally annotated with VEP or SnpEff) and writes
    a MAF file suitable for downstream analysis in cBioPortal, maftools, etc.

    \b
    Examples:
      vcf2maf-py vcf2maf variants.vcf -o variants.maf
      vcf2maf-py vcf2maf somatic.vcf.gz --tumor-id TUMOR --normal-id NORMAL
      vcf2maf-py vcf2maf annotated.vcf --ncbi-build GRCh37 --center MSKCC
    """
    from vcf2maf_py.constants import ALL_MAF_COLUMNS, CORE_MAF_COLUMNS, EXTENDED_MAF_COLUMNS
    from vcf2maf_py.converter import convert_vcf_to_maf, load_isoform_overrides
    from vcf2maf_py.utils import setup_logging

    setup_logging(verbose)

    if output_maf is None:
        base = input_vcf.replace(".vcf.gz", "").replace(".vcf", "")
        output_maf = base + ".maf"

    # Choose columns
    col_map = {
        "core": CORE_MAF_COLUMNS,
        "extended": CORE_MAF_COLUMNS + EXTENDED_MAF_COLUMNS,
        "all": ALL_MAF_COLUMNS,
    }
    cols = col_map[output_columns]

    # Load isoform overrides
    overrides = None
    if isoform_overrides:
        overrides = load_isoform_overrides(isoform_overrides)

    rows = convert_vcf_to_maf(
        input_vcf=input_vcf,
        output_maf=output_maf,
        tumor_id=tumor_id,
        normal_id=normal_id,
        ncbi_build=ncbi_build,
        center=center,
        min_hom_vaf=min_hom_vaf,
        isoform_overrides=overrides,
        output_columns=cols,
    )

    click.echo(f"Converted {len(rows)} variants → {output_maf}", err=True)


@main.command()
@click.argument("input_maf", type=click.Path(exists=True))
@click.option("-o", "--output-vcf", type=click.Path(), default=None,
              help="Output VCF path (default: <input>.vcf)")
@click.option("--ref-fasta", type=click.Path(exists=True), default=None,
              help="Reference FASTA for correct indel padding bases (optional, requires pysam)")
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
def maf2vcf(
    input_maf: str,
    output_vcf: str | None,
    ref_fasta: str | None,
    verbose: int,
) -> None:
    """Convert a MAF file back to VCF format.

    Reads a MAF file and reconstructs a VCF. For correct indel representation,
    provide a reference FASTA (requires pysam). Without it, 'N' is used as
    the padding base.

    \b
    Examples:
      vcf2maf-py maf2vcf mutations.maf -o mutations.vcf
      vcf2maf-py maf2vcf mutations.maf --ref-fasta GRCh38.fa
    """
    from vcf2maf_py.maf2vcf import convert_maf_to_vcf
    from vcf2maf_py.utils import setup_logging

    setup_logging(verbose)

    if output_vcf is None:
        base = input_maf.replace(".maf.gz", "").replace(".maf", "")
        output_vcf = base + ".vcf"

    convert_maf_to_vcf(
        input_maf=input_maf,
        output_vcf=output_vcf,
        ref_fasta=ref_fasta,
    )

    click.echo(f"Converted MAF → {output_vcf}", err=True)


@main.command()
@click.argument("input_vcf", type=click.Path(exists=True))
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
def inspect(input_vcf: str, verbose: int) -> None:
    """Inspect a VCF file and report annotation status.

    Shows detected annotation source (VEP/SnpEff), sample names,
    variant count, and a summary of variant types.
    """
    from vcf2maf_py.utils import setup_logging
    from vcf2maf_py.vcf_reader import VCFReader

    setup_logging(verbose)
    reader = VCFReader(input_vcf)

    click.echo(f"File: {input_vcf}")
    click.echo(f"Samples: {', '.join(reader.header.samples) or 'none'}")

    ann = reader.header.annotation_source
    if ann == "VEP":
        fields = reader.header.csq_fields
        click.echo(f"Annotation: VEP ({len(fields or [])} CSQ fields)")
    elif ann == "SnpEff":
        fields = reader.header.ann_fields
        click.echo(f"Annotation: SnpEff ({len(fields or [])} ANN fields)")
    else:
        click.echo("Annotation: NONE — consider annotating with VEP before conversion")

    n_records = 0
    n_snps = 0
    n_indels = 0
    n_multi = 0
    for rec in reader:
        n_records += 1
        if rec.is_multiallelic:
            n_multi += 1
        for alt in rec.alt:
            if len(rec.ref) == len(alt) == 1:
                n_snps += 1
            else:
                n_indels += 1

    click.echo(f"Records: {n_records:,} ({n_snps:,} SNPs, {n_indels:,} indels, {n_multi:,} multi-allelic)")
