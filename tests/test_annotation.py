"""Tests for annotation parsing and effect selection."""

from __future__ import annotations

import pytest

from vcf2maf_py.annotation import (
    TranscriptEffect,
    parse_snpeff_ann,
    parse_vep_csq,
    select_best_effect,
)
from vcf2maf_py.vcf_reader import VCFReader


class TestVEPParsing:
    def test_parse_single_effect(self):
        field_names = [
            "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
            "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON",
            "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
            "Protein_position", "Amino_acids", "Codons",
            "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
            "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "SIFT", "PolyPhen",
            "ALLELE_NUM",
        ]
        info = {
            "CSQ": (
                "T|missense_variant|MODERATE|BRAF|ENSG00000157764|"
                "Transcript|ENST00000288602|protein_coding|15/18||"
                "c.1799T>A|p.Val600Glu|1860|1799|600|V/E|gTg/gAg|"
                "rs113488022||1||HGNC|HGNC:1097|YES|deleterious(0)|"
                "probably_damaging(0.999)|1"
            )
        }
        effects = parse_vep_csq(info, field_names, alt_alleles=["T"])
        assert len(effects) == 1
        e = effects[0]
        assert e.allele == "T"
        assert e.consequence == "missense_variant"
        assert e.impact == "MODERATE"
        assert e.symbol == "BRAF"
        assert e.feature == "ENST00000288602"
        assert e.biotype == "protein_coding"
        assert e.hgvsc == "c.1799T>A"
        assert e.hgvsp == "p.Val600Glu"
        assert e.is_canonical is True

    def test_parse_multiple_effects(self):
        field_names = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
                       "Feature_type", "Feature", "BIOTYPE", "CANONICAL"]
        info = {
            "CSQ": (
                "T|missense_variant|MODERATE|GENE1|G1|Transcript|T1|protein_coding|YES,"
                "T|intron_variant|MODIFIER|GENE2|G2|Transcript|T2|lncRNA|"
            )
        }
        effects = parse_vep_csq(info, field_names)
        assert len(effects) == 2
        assert effects[0].symbol == "GENE1"
        assert effects[1].symbol == "GENE2"

    def test_allele_filtering(self):
        field_names = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
                       "Feature_type", "Feature", "BIOTYPE", "CANONICAL", "ALLELE_NUM"]
        info = {
            "CSQ": (
                "A|missense_variant|MODERATE|G1|G1|Transcript|T1|protein_coding|YES|1,"
                "T|stop_gained|HIGH|G1|G1|Transcript|T1|protein_coding|YES|2"
            )
        }
        effects = parse_vep_csq(
            info, field_names,
            alt_alleles=["A", "T"], target_allele="A",
        )
        assert len(effects) == 1
        assert effects[0].consequence == "missense_variant"


class TestSnpEffParsing:
    def test_parse_ann(self):
        info = {
            "ANN": (
                "T|missense_variant|MODERATE|BRAF|ENSG00000157764|transcript|"
                "ENST00000288602|protein_coding|15/18|c.1799T>A|p.Val600Glu|"
                "1860/8370|1799/2301|600/766||"
            )
        }
        effects = parse_snpeff_ann(info)
        assert len(effects) == 1
        e = effects[0]
        assert e.allele == "T"
        assert e.consequence == "missense_variant"
        assert e.symbol == "BRAF"
        assert e.feature == "ENST00000288602"


class TestSelectBestEffect:
    def test_prefers_protein_coding(self):
        e1 = TranscriptEffect(
            consequence="intron_variant", biotype="lncRNA", feature="T1",
        )
        e2 = TranscriptEffect(
            consequence="intron_variant", biotype="protein_coding", feature="T2",
        )
        best = select_best_effect([e1, e2])
        assert best.feature == "T2"

    def test_prefers_higher_severity(self):
        e1 = TranscriptEffect(
            consequence="synonymous_variant", biotype="protein_coding", feature="T1",
        )
        e2 = TranscriptEffect(
            consequence="missense_variant", biotype="protein_coding", feature="T2",
        )
        best = select_best_effect([e1, e2])
        assert best.feature == "T2"

    def test_prefers_canonical(self):
        e1 = TranscriptEffect(
            consequence="missense_variant", biotype="protein_coding",
            feature="T1", canonical="",
        )
        e2 = TranscriptEffect(
            consequence="missense_variant", biotype="protein_coding",
            feature="T2", canonical="YES",
        )
        best = select_best_effect([e1, e2])
        assert best.feature == "T2"

    def test_single_effect(self):
        e = TranscriptEffect(consequence="missense_variant", feature="T1")
        assert select_best_effect([e]) is e

    def test_empty_list(self):
        assert select_best_effect([]) is None

    def test_prefers_isoform_override(self):
        e1 = TranscriptEffect(
            consequence="missense_variant", biotype="protein_coding",
            symbol="BRAF", feature="ENST00000288602", canonical="YES",
        )
        e2 = TranscriptEffect(
            consequence="missense_variant", biotype="protein_coding",
            symbol="BRAF", feature="ENST00000496384", canonical="",
        )
        overrides = {"BRAF": "ENST00000496384"}
        best = select_best_effect([e1, e2], isoform_overrides=overrides)
        assert best.feature == "ENST00000496384"


class TestTranscriptEffect:
    def test_worst_consequence(self):
        e = TranscriptEffect(consequence="missense_variant&splice_region_variant")
        assert e.worst_consequence == "missense_variant"

    def test_effect_priority(self):
        e1 = TranscriptEffect(consequence="missense_variant")
        e2 = TranscriptEffect(consequence="stop_gained")
        assert e2.effect_priority < e1.effect_priority

    def test_biotype_priority(self):
        e1 = TranscriptEffect(biotype="protein_coding")
        e2 = TranscriptEffect(biotype="pseudogene")
        assert e1.biotype_priority < e2.biotype_priority
