"""Tests for edge cases in VCF → MAF conversion."""

from __future__ import annotations

import pytest

from vcf2maf_py.converter import convert_vcf_to_maf


class TestEdgeCases:
    @pytest.fixture(autouse=True)
    def _convert(self, edge_cases_vcf, tmp_maf):
        self.rows = convert_vcf_to_maf(edge_cases_vcf, tmp_maf)

    def _find(self, chrom, pos=None):
        matches = [r for r in self.rows if r["Chromosome"] == chrom]
        if pos is not None:
            matches = [r for r in matches if int(r["Start_Position"]) == pos or int(r["vcf_pos"]) == pos]
        return matches

    def test_dnp(self):
        """chr7:140753335 CA>GT is a DNP."""
        rows = self._find("chr7")
        assert len(rows) == 1
        assert rows[0]["Variant_Type"] == "DNP"
        assert rows[0]["Reference_Allele"] == "CA"
        assert rows[0]["Tumor_Seq_Allele2"] == "GT"
        assert rows[0]["Variant_Classification"] == "Missense_Mutation"

    def test_homozygous(self):
        """chr1:100 — GT=1/1 should be homozygous alt."""
        rows = self._find("chr1", 100)
        assert len(rows) == 1
        r = rows[0]
        assert r["Variant_Type"] == "DEL"
        assert r["Variant_Classification"] == "Splice_Site"
        assert r["t_depth"] == "40"
        assert r["t_alt_count"] == "40"

    def test_inframe_insertion(self):
        """chr2:200 A>ATCGA is an in-frame insertion."""
        rows = self._find("chr2")
        assert len(rows) == 1
        r = rows[0]
        assert r["Variant_Type"] == "INS"
        assert r["Reference_Allele"] == "-"
        assert r["Tumor_Seq_Allele2"] == "TCGA"
        assert r["Variant_Classification"] == "In_Frame_Ins"

    def test_onp(self):
        """chr3:300 GCGA>TCGT is an ONP (4bp substitution)."""
        rows = self._find("chr3")
        assert len(rows) == 1
        r = rows[0]
        assert r["Variant_Type"] == "ONP"
        assert r["Reference_Allele"] == "GCGA"
        assert r["Tumor_Seq_Allele2"] == "TCGT"

    def test_non_pass_filter(self):
        """chrX:100500 has LowQual FILTER — still converted."""
        rows = self._find("chrX")
        assert len(rows) == 1
        assert rows[0]["FILTER"] == "LowQual"
        assert rows[0]["Variant_Classification"] == "Silent"

    def test_intergenic(self):
        """chr5:500 is intergenic — no gene symbol."""
        rows = self._find("chr5")
        assert len(rows) == 1
        r = rows[0]
        assert r["Variant_Classification"] == "IGR"
        assert r["Hugo_Symbol"] == ""

    def test_missing_normal_genotype(self):
        """chr1:1000 has ./.:.:. for normal — depths should be empty."""
        rows = self._find("chr1", 1000)
        assert len(rows) == 1
        r = rows[0]
        assert r["t_depth"] == "45"
        assert r["t_alt_count"] == "10"
        # Normal has ./. — no depth info
        assert r["n_depth"] == ""
        assert r["n_alt_count"] == ""

    def test_star_allele_skipped(self):
        """chr6:600 has ALT=* (spanning deletion) — should be skipped."""
        rows = self._find("chr6")
        assert len(rows) == 0

    def test_no_alt_skipped(self):
        """chr7:700 has ALT=. (monomorphic) — should be skipped."""
        rows = self._find("chr7", 700)
        assert len(rows) == 0

    def test_mitochondrial(self):
        """chrM:100 — mitochondrial variant."""
        rows = self._find("chrM")
        assert len(rows) == 1
        r = rows[0]
        assert r["Chromosome"] == "chrM"
        assert r["Hugo_Symbol"] == "MT-RNR1"
        assert r["Variant_Classification"] == "Silent"
        # GT=1/1, so hom alt
        assert r["t_alt_count"] == "100"

    def test_phased_genotype(self):
        """chr2:200 has phased GT (0|1) — should still parse correctly."""
        rows = self._find("chr2")
        assert len(rows) == 1
        assert rows[0]["t_depth"] == "40"
        assert rows[0]["t_alt_count"] == "15"

    def test_total_variant_count(self):
        """Should produce 8 rows (star and no-alt skipped)."""
        assert len(self.rows) == 8


class TestTumorOnly:
    """VCF with a single sample (no matched normal)."""

    def test_single_sample(self, tumor_only_vcf, tmp_maf):
        rows = convert_vcf_to_maf(
            tumor_only_vcf, tmp_maf,
            tumor_id="SAMPLE1",
        )
        assert len(rows) == 3

        by_gene = {r["Hugo_Symbol"]: r for r in rows}
        assert by_gene["BRAF"]["Tumor_Sample_Barcode"] == "SAMPLE1"
        assert by_gene["BRAF"]["t_depth"] == "45"
        assert by_gene["BRAF"]["t_alt_count"] == "15"
        # No normal — depths should be empty
        assert by_gene["BRAF"]["n_depth"] == ""
        assert by_gene["BRAF"]["n_alt_count"] == ""

    def test_homozygous_tumor_only(self, tumor_only_vcf, tmp_maf):
        rows = convert_vcf_to_maf(tumor_only_vcf, tmp_maf)
        kras = [r for r in rows if r["Hugo_Symbol"] == "KRAS"][0]
        # GT=1/1 → homozygous alt
        assert kras["t_ref_count"] == "0"
        assert kras["t_alt_count"] == "60"
