"""Tests for transcript selection when multiple annotations exist per variant."""

from __future__ import annotations

import pytest

from vcf2maf_py.converter import convert_vcf_to_maf


class TestMultiTranscriptSelection:
    @pytest.fixture(autouse=True)
    def _convert(self, multi_transcript_vcf, tmp_maf):
        self.rows = convert_vcf_to_maf(multi_transcript_vcf, tmp_maf)
        self.by_gene = {r["Hugo_Symbol"]: r for r in self.rows}

    def test_braf_selects_canonical_coding(self):
        """BRAF has 3 transcripts: canonical protein_coding (missense),
        non-canonical protein_coding (intron), and lncRNA.
        Should pick the canonical missense."""
        r = self.by_gene["BRAF"]
        assert r["Transcript_ID"] == "ENST00000288602"
        assert r["Variant_Classification"] == "Missense_Mutation"
        assert r["HGVSp"] == "p.Val600Glu"

    def test_tp53_selects_canonical_over_non_canonical(self):
        """TP53 has canonical ENST00000269305 and non-canonical ENST00000359597
        (both missense), plus WRAP53 downstream. Should pick canonical TP53."""
        r = self.by_gene["TP53"]
        assert r["Transcript_ID"] == "ENST00000269305"
        assert r["Variant_Classification"] == "Missense_Mutation"

    def test_kras_selects_canonical_over_pseudogene(self):
        """KRAS has canonical transcript, non-canonical transcript, and a
        processed_pseudogene. Should pick canonical."""
        r = self.by_gene["KRAS"]
        assert r["Transcript_ID"] == "ENST00000256078"
        assert r["Variant_Classification"] == "Missense_Mutation"

    def test_pik3ca_selects_coding_over_noncoding(self):
        """PIK3CA has canonical protein_coding (missense), lncRNA, and
        non-canonical protein_coding (intron). Should pick canonical missense."""
        r = self.by_gene["PIK3CA"]
        assert r["Transcript_ID"] == "ENST00000263967"
        assert r["Variant_Classification"] == "Missense_Mutation"

    def test_all_effects_populated(self):
        """all_effects column should contain all transcripts, not just the best."""
        r = self.by_gene["BRAF"]
        all_eff = r["all_effects"]
        assert "ENST00000288602" in all_eff
        assert "ENST00000496384" in all_eff
        assert "ENST00000497784" in all_eff

    def test_all_effects_count(self):
        """BRAF has 3 transcripts; all_effects should have 3 entries (with trailing ;)."""
        r = self.by_gene["BRAF"]
        # Strip trailing ";" before splitting to get clean count
        entries = [e for e in r["all_effects"].split(";") if e]
        assert len(entries) == 3

    def test_tp53_all_effects_includes_wrap53(self):
        """TP53 record also annotates WRAP53 downstream — should appear in all_effects."""
        r = self.by_gene["TP53"]
        assert "WRAP53" in r["all_effects"]


class TestIsoformOverride:
    def test_override_changes_selection(self, multi_transcript_vcf, tmp_maf, tmp_path):
        """Custom isoform override should take precedence over canonical."""
        # Create an override file preferring non-canonical BRAF transcript
        override_file = tmp_path / "overrides.tsv"
        override_file.write_text("ENST00000496384\tBRAF\n")

        rows = convert_vcf_to_maf(
            multi_transcript_vcf, tmp_maf,
            isoform_overrides={"BRAF": "ENST00000496384"},
        )
        braf = [r for r in rows if r["Hugo_Symbol"] == "BRAF"][0]
        assert braf["Transcript_ID"] == "ENST00000496384"
        # This transcript has intron_variant, not missense
        assert braf["Variant_Classification"] == "Intron"
