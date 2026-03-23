"""Tests for MAF → VCF using the standalone minimal MAF file."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2maf_py.maf2vcf import convert_maf_to_vcf


class TestMinimalMafToVcf:
    def test_converts_all_records(self, minimal_maf, tmp_vcf):
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        lines = Path(tmp_vcf).read_text().splitlines()
        data = [l for l in lines if not l.startswith("#") and l.strip()]
        assert len(data) == 6

    def test_snp_no_padding(self, minimal_maf, tmp_vcf):
        """SNPs should not have padding bases added."""
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        lines = Path(tmp_vcf).read_text().splitlines()
        data = [l for l in lines if not l.startswith("#") and l.strip()]

        # BRAF: SNP A>T at chr7:140753336
        braf = [l for l in data if "chr7" in l and "140753336" in l][0]
        fields = braf.split("\t")
        assert fields[3] == "A"   # REF
        assert fields[4] == "T"   # ALT

    def test_insertion_gets_padding(self, minimal_maf, tmp_vcf):
        """Insertions need a padding base in VCF."""
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        lines = Path(tmp_vcf).read_text().splitlines()
        data = [l for l in lines if not l.startswith("#") and l.strip()]

        # NRAS: INS at chr1:115256530-115256531, ref=-, alt=A
        nras = [l for l in data if l.startswith("chr1\t")][0]
        fields = nras.split("\t")
        ref, alt = fields[3], fields[4]
        # VCF indel: REF should be 1 base, ALT should be 2 bases (padding+A)
        assert len(ref) == 1
        assert len(alt) == 2
        assert alt.endswith("A")

    def test_deletion_gets_padding(self, minimal_maf, tmp_vcf):
        """Deletions need a padding base in VCF."""
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        lines = Path(tmp_vcf).read_text().splitlines()
        data = [l for l in lines if not l.startswith("#") and l.strip()]

        # IDH1: DEL at chr2:208248390-208248392, ref=CTG, alt=-
        idh1 = [l for l in data if "chr2" in l][0]
        fields = idh1.split("\t")
        ref, alt = fields[3], fields[4]
        # VCF: REF=padding+CTG (4 bases), ALT=padding (1 base)
        assert len(ref) == 4
        assert len(alt) == 1
        assert ref.endswith("CTG")

    def test_genotype_fields_present(self, minimal_maf, tmp_vcf):
        """VCF records should have GT:AD:DP format fields."""
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        lines = Path(tmp_vcf).read_text().splitlines()
        data = [l for l in lines if not l.startswith("#") and l.strip()]

        for line in data:
            fields = line.split("\t")
            assert fields[8] == "GT:AD:DP"
            # At least tumor sample column
            assert len(fields) >= 10

    def test_dbsnp_rs_preserved(self, minimal_maf, tmp_vcf):
        """dbSNP RS IDs should appear in VCF ID column."""
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        lines = Path(tmp_vcf).read_text().splitlines()
        data = [l for l in lines if not l.startswith("#") and l.strip()]

        braf = [l for l in data if "chr7" in l and "140753336" in l][0]
        fields = braf.split("\t")
        assert fields[2] == "rs113488022"

    def test_pairs_file(self, minimal_maf, tmp_vcf):
        """Pairs TSV should list tumor-normal pairs."""
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        pairs = Path(tmp_vcf.replace(".vcf", ".pairs.tsv")).read_text().strip()
        assert "TUMOR\tNORMAL" in pairs

    def test_vcf_header_has_build(self, minimal_maf, tmp_vcf):
        convert_maf_to_vcf(minimal_maf, tmp_vcf)
        text = Path(tmp_vcf).read_text()
        assert "##reference=GRCh38" in text
