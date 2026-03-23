"""Tests for MAF → VCF reverse conversion."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf2maf_py.converter import convert_vcf_to_maf
from vcf2maf_py.maf2vcf import convert_maf_to_vcf


class TestMaf2Vcf:
    def test_roundtrip(self, vep_vcf, tmp_maf, tmp_vcf):
        """Test VCF → MAF → VCF roundtrip preserves basic variant info."""
        # First convert VCF to MAF
        convert_vcf_to_maf(vep_vcf, tmp_maf)

        # Then convert MAF back to VCF
        convert_maf_to_vcf(tmp_maf, tmp_vcf)

        # Check VCF output exists and is valid
        vcf_path = Path(tmp_vcf)
        assert vcf_path.exists()
        assert vcf_path.stat().st_size > 0

        # Parse output VCF
        lines = vcf_path.read_text().splitlines()
        header_lines = [l for l in lines if l.startswith("##")]
        column_line = [l for l in lines if l.startswith("#CHROM")]
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        assert len(header_lines) > 0
        assert len(column_line) == 1
        assert len(data_lines) > 0

        # Check that VCF has proper format
        assert "##fileformat=VCFv4.2" in header_lines[0]
        assert "FORMAT" in column_line[0]

    def test_pairs_file_created(self, vep_vcf, tmp_maf, tmp_vcf):
        convert_vcf_to_maf(vep_vcf, tmp_maf)
        convert_maf_to_vcf(tmp_maf, tmp_vcf)

        pairs_path = Path(tmp_vcf.replace(".vcf", ".pairs.tsv"))
        assert pairs_path.exists()
        content = pairs_path.read_text().strip()
        assert "TUMOR" in content
        assert "NORMAL" in content

    def test_snp_alleles(self, vep_vcf, tmp_maf, tmp_vcf):
        convert_vcf_to_maf(vep_vcf, tmp_maf)
        convert_maf_to_vcf(tmp_maf, tmp_vcf)

        lines = Path(tmp_vcf).read_text().splitlines()
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        # Find a SNP record (BRAF: chr7)
        chr7_lines = [l for l in data_lines if l.startswith("chr7")]
        assert len(chr7_lines) >= 1
        fields = chr7_lines[0].split("\t")
        # For a SNP, REF and ALT should be single bases
        assert len(fields[3]) == 1  # REF
        assert len(fields[4]) == 1  # ALT

    def test_indel_padding(self, vep_vcf, tmp_maf, tmp_vcf):
        """Indels should get a padding base in VCF."""
        convert_vcf_to_maf(vep_vcf, tmp_maf)
        convert_maf_to_vcf(tmp_maf, tmp_vcf)

        lines = Path(tmp_vcf).read_text().splitlines()
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]

        # Check that indels have padding
        for line in data_lines:
            fields = line.split("\t")
            ref, alt = fields[3], fields[4]
            # All VCF records should have non-empty REF and ALT
            assert ref != "-"
            assert alt != "-"
            assert len(ref) > 0
            assert len(alt) > 0
