"""Tests for the full VCF → MAF conversion pipeline."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from vcf2maf_py.converter import (
    convert_vcf_to_maf,
    extract_depths,
    get_variant_classification,
)
from vcf2maf_py.annotation import TranscriptEffect


class TestExtractDepths:
    def test_ad_field(self):
        data = {"GT": "0/1", "AD": "30,15", "DP": "45"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert depth == 45
        assert ref == 30
        assert alt == 15

    def test_ro_ao_fields(self):
        data = {"GT": "0/1", "RO": "30", "AO": "15", "DP": "45"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert ref == 30
        assert alt == 15
        assert depth == 45

    def test_dp4_field(self):
        data = {"GT": "0/1", "DP4": "10,10,8,7"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert ref == 20
        assert alt == 15
        assert depth == 35

    def test_strelka_snv(self):
        data = {"GT": "0/1", "AU": "30,32", "TU": "15,16"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert ref == 30
        assert alt == 15

    def test_strelka_indel(self):
        data = {"GT": "0/1", "TAR": "30,32", "TIR": "15,16", "DP": "50"}
        depth, ref, alt = extract_depths(data, "A", "AT")
        assert ref == 30
        assert alt == 15
        assert depth == 50

    def test_nr_nv_fields(self):
        data = {"GT": "0/1", "NR": "45", "NV": "15"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert depth == 45
        assert alt == 15
        assert ref == 30

    def test_missing_all(self):
        data = {"GT": "0/1"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert depth == -1
        assert ref == -1
        assert alt == -1

    def test_infer_depth_from_counts(self):
        data = {"AD": "30,15"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert depth == 45

    def test_dv_field(self):
        data = {"DV": "15", "DP": "45"}
        depth, ref, alt = extract_depths(data, "A", "T")
        assert alt == 15
        assert depth == 45
        assert ref == 30


class TestGetVariantClassification:
    def test_missense(self):
        e = TranscriptEffect(consequence="missense_variant")
        assert get_variant_classification(e, "SNP") == "Missense_Mutation"

    def test_nonsense(self):
        e = TranscriptEffect(consequence="stop_gained")
        assert get_variant_classification(e, "SNP") == "Nonsense_Mutation"

    def test_frameshift_del(self):
        e = TranscriptEffect(consequence="frameshift_variant")
        assert get_variant_classification(e, "DEL") == "Frame_Shift_Del"

    def test_frameshift_ins(self):
        e = TranscriptEffect(consequence="frameshift_variant")
        assert get_variant_classification(e, "INS") == "Frame_Shift_Ins"

    def test_splice_site(self):
        e = TranscriptEffect(consequence="splice_acceptor_variant")
        assert get_variant_classification(e, "SNP") == "Splice_Site"

    def test_silent(self):
        e = TranscriptEffect(consequence="synonymous_variant")
        assert get_variant_classification(e, "SNP") == "Silent"

    def test_intron(self):
        e = TranscriptEffect(consequence="intron_variant")
        assert get_variant_classification(e, "SNP") == "Intron"

    def test_compound_consequence(self):
        e = TranscriptEffect(consequence="missense_variant&splice_region_variant")
        assert get_variant_classification(e, "SNP") == "Missense_Mutation"

    def test_protein_altering_del(self):
        e = TranscriptEffect(consequence="protein_altering_variant")
        assert get_variant_classification(e, "DEL") == "In_Frame_Del"

    def test_protein_altering_ins(self):
        e = TranscriptEffect(consequence="protein_altering_variant")
        assert get_variant_classification(e, "INS") == "In_Frame_Ins"

    def test_none_effect(self):
        assert get_variant_classification(None, "SNP") == ""


class TestConvertVcfToMaf:
    def test_vep_annotated(self, vep_vcf, tmp_maf):
        rows = convert_vcf_to_maf(vep_vcf, tmp_maf, ncbi_build="GRCh38")
        assert len(rows) > 0

        # Check BRAF V600E
        braf = [r for r in rows if r["Hugo_Symbol"] == "BRAF"]
        assert len(braf) == 1
        assert braf[0]["Variant_Classification"] == "Missense_Mutation"
        assert braf[0]["Chromosome"] == "chr7"
        assert braf[0]["Variant_Type"] == "SNP"
        assert braf[0]["Reference_Allele"] == "A"
        assert braf[0]["Tumor_Seq_Allele2"] == "T"
        assert braf[0]["HGVSp"] == "p.Val600Glu"
        assert braf[0]["t_depth"] == "45"
        assert braf[0]["t_ref_count"] == "30"
        assert braf[0]["t_alt_count"] == "15"

        # Check KRAS G12V
        kras = [r for r in rows if r["Hugo_Symbol"] == "KRAS"]
        assert len(kras) == 1
        assert kras[0]["Variant_Classification"] == "Missense_Mutation"

        # Check TP53 nonsense
        tp53 = [r for r in rows if r["Hugo_Symbol"] == "TP53"]
        assert len(tp53) == 1
        assert tp53[0]["Variant_Classification"] == "Nonsense_Mutation"

        # Check NRAS frameshift insertion
        nras = [r for r in rows if r["Hugo_Symbol"] == "NRAS"]
        assert len(nras) == 1
        assert nras[0]["Variant_Classification"] == "Frame_Shift_Ins"
        assert nras[0]["Variant_Type"] == "INS"

        # Check IDH1 in-frame deletion
        idh1 = [r for r in rows if r["Hugo_Symbol"] == "IDH1"]
        assert len(idh1) == 1
        assert idh1[0]["Variant_Classification"] == "In_Frame_Del"
        assert idh1[0]["Variant_Type"] == "DEL"

    def test_multiallelic_split(self, vep_vcf, tmp_maf):
        rows = convert_vcf_to_maf(vep_vcf, tmp_maf)
        # APC has two ALT alleles — should produce two MAF rows
        apc = [r for r in rows if r["Hugo_Symbol"] == "APC"]
        assert len(apc) == 2

    def test_snpeff_annotated(self, snpeff_vcf, tmp_maf):
        rows = convert_vcf_to_maf(snpeff_vcf, tmp_maf)
        assert len(rows) == 2
        braf = [r for r in rows if r["Hugo_Symbol"] == "BRAF"]
        assert len(braf) == 1
        assert braf[0]["Variant_Classification"] == "Missense_Mutation"

    def test_unannotated_warns(self, unannotated_vcf, tmp_maf, caplog):
        import logging
        with caplog.at_level(logging.WARNING):
            rows = convert_vcf_to_maf(unannotated_vcf, tmp_maf)
        assert "No VEP" in caplog.text or "No VEP" in caplog.text
        assert len(rows) == 4  # 4 variant records

    def test_maf_file_format(self, vep_vcf, tmp_maf):
        convert_vcf_to_maf(vep_vcf, tmp_maf)
        with open(tmp_maf) as fh:
            first_line = fh.readline()
            assert first_line.startswith("#version 2.4")
            header = fh.readline()
            assert "Hugo_Symbol" in header
            assert "Variant_Classification" in header

    def test_sample_auto_detection(self, vep_vcf, tmp_maf):
        rows = convert_vcf_to_maf(vep_vcf, tmp_maf)
        assert rows[0]["Tumor_Sample_Barcode"] == "TUMOR"
        assert rows[0]["Matched_Norm_Sample_Barcode"] == "NORMAL"

    def test_custom_tumor_normal_id(self, vep_vcf, tmp_maf):
        rows = convert_vcf_to_maf(
            vep_vcf, tmp_maf,
            tumor_id="TUMOR", normal_id="NORMAL",
        )
        assert rows[0]["Tumor_Sample_Barcode"] == "TUMOR"
        assert rows[0]["Matched_Norm_Sample_Barcode"] == "NORMAL"

    def test_ncbi_build_in_output(self, vep_vcf, tmp_maf):
        rows = convert_vcf_to_maf(vep_vcf, tmp_maf, ncbi_build="GRCh37")
        assert all(r["NCBI_Build"] == "GRCh37" for r in rows)

    def test_return_without_file(self, vep_vcf):
        rows = convert_vcf_to_maf(vep_vcf, output_maf=None)
        assert len(rows) > 0

    def test_normal_depths(self, vep_vcf, tmp_maf):
        rows = convert_vcf_to_maf(vep_vcf, tmp_maf)
        braf = [r for r in rows if r["Hugo_Symbol"] == "BRAF"][0]
        assert braf["n_depth"] == "40"
        assert braf["n_ref_count"] == "40"
        assert braf["n_alt_count"] == "0"
