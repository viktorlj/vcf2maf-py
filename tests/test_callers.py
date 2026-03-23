"""Tests for depth extraction across different variant caller FORMAT conventions."""

from __future__ import annotations

import pytest

from vcf2maf_py.converter import convert_vcf_to_maf


class TestCallerFormats:
    """Each record in test_callers.vcf uses a different FORMAT convention."""

    @pytest.fixture(autouse=True)
    def _convert(self, callers_vcf, tmp_maf):
        self.rows = convert_vcf_to_maf(callers_vcf, tmp_maf)
        self.by_gene = {r["Hugo_Symbol"]: r for r in self.rows}

    def test_gatk_ad(self):
        """BRAF: standard GT:AD:DP format."""
        r = self.by_gene["BRAF"]
        assert r["t_depth"] == "45"
        assert r["t_ref_count"] == "30"
        assert r["t_alt_count"] == "15"
        assert r["n_depth"] == "40"
        assert r["n_ref_count"] == "40"
        assert r["n_alt_count"] == "0"

    def test_freebayes_ro_ao(self):
        """KRAS: FreeBayes GT:RO:AO:DP format."""
        r = self.by_gene["KRAS"]
        assert r["t_depth"] == "75"
        assert r["t_ref_count"] == "50"
        assert r["t_alt_count"] == "25"
        assert r["n_depth"] == "60"
        assert r["n_ref_count"] == "60"
        assert r["n_alt_count"] == "0"

    def test_strelka_snv_tiers(self):
        """TP53: Strelka SNV GT:{A,C,G,T}U format."""
        r = self.by_gene["TP53"]
        assert r["t_ref_count"] == "35"
        assert r["t_alt_count"] == "20"
        assert r["n_ref_count"] == "45"
        assert r["n_alt_count"] == "1"

    def test_dp4(self):
        """PIK3CA: GT:DP4 (ref-fwd,ref-rev,alt-fwd,alt-rev) format."""
        r = self.by_gene["PIK3CA"]
        assert r["t_ref_count"] == "40"  # 18+22
        assert r["t_alt_count"] == "18"  # 10+8
        assert r["n_ref_count"] == "55"  # 30+25
        assert r["n_alt_count"] == "1"   # 0+1

    def test_strelka_indel_tar_tir(self):
        """NRAS: Strelka indel GT:TAR:TIR:DP format."""
        r = self.by_gene["NRAS"]
        assert r["t_ref_count"] == "40"
        assert r["t_alt_count"] == "10"
        assert r["t_depth"] == "52"
        assert r["n_ref_count"] == "50"
        assert r["n_alt_count"] == "0"

    def test_nr_nv(self):
        """ALK: GT:NR:NV format."""
        r = self.by_gene["ALK"]
        assert r["t_depth"] == "60"
        assert r["t_alt_count"] == "18"
        assert r["t_ref_count"] == "42"  # 60-18
        assert r["n_depth"] == "55"
        assert r["n_alt_count"] == "1"

    def test_dv(self):
        """PTEN: GT:DV:DP format."""
        r = self.by_gene["PTEN"]
        assert r["t_depth"] == "50"
        assert r["t_alt_count"] == "14"
        assert r["t_ref_count"] == "36"  # 50-14
        assert r["n_depth"] == "45"
        assert r["n_alt_count"] == "0"

    def test_all_classifications_correct(self):
        """All variants should be classified as Missense or Frameshift."""
        for gene, expected in [
            ("BRAF", "Missense_Mutation"),
            ("KRAS", "Missense_Mutation"),
            ("TP53", "Nonsense_Mutation"),
            ("PIK3CA", "Missense_Mutation"),
            ("NRAS", "Frame_Shift_Ins"),
            ("ALK", "Missense_Mutation"),
            ("PTEN", "Missense_Mutation"),
        ]:
            assert self.by_gene[gene]["Variant_Classification"] == expected, gene
