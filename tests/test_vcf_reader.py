"""Tests for the lightweight VCF parser."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from vcf2maf_py.vcf_reader import VCFReader


class TestVCFReader:
    def test_reads_vep_header(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        assert reader.header.annotation_source == "VEP"
        assert reader.header.samples == ["TUMOR", "NORMAL"]
        assert "CSQ" in reader.header.info_fields
        fields = reader.header.csq_fields
        assert "Allele" in fields
        assert "Consequence" in fields

    def test_reads_snpeff_header(self, snpeff_vcf):
        reader = VCFReader(snpeff_vcf)
        assert reader.header.annotation_source == "SnpEff"
        assert reader.header.ann_fields is not None

    def test_no_annotation(self, unannotated_vcf):
        reader = VCFReader(unannotated_vcf)
        assert reader.header.annotation_source is None

    def test_iterates_records(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        records = list(reader)
        assert len(records) == 6

    def test_record_fields(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        rec = next(iter(reader))
        assert rec.chrom == "chr7"
        assert rec.pos == 140753336
        assert rec.ref == "A"
        assert rec.alt == ["T"]
        assert rec.filter == "PASS"
        assert rec.id == "rs113488022"

    def test_multiallelic_detected(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        records = list(reader)
        multi = [r for r in records if r.is_multiallelic]
        assert len(multi) == 1  # APC record
        assert len(multi[0].alt) == 2

    def test_info_parsing(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        rec = next(iter(reader))
        assert "CSQ" in rec.info
        assert rec.info["CSQ"].startswith("T|")

    def test_sample_data(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        rec = next(iter(reader))
        tumor = rec.sample_data(0)
        assert tumor["GT"] == "0/1"
        assert tumor["AD"] == "30,15"
        assert tumor["DP"] == "45"

    def test_sample_data_by_name(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        rec = next(iter(reader))
        normal = rec.sample_data_by_name("NORMAL")
        assert normal["GT"] == "0/0"

    def test_single_sample_vcf(self, tumor_only_vcf):
        reader = VCFReader(tumor_only_vcf)
        assert reader.header.samples == ["SAMPLE1"]
        records = list(reader)
        assert len(records) == 3
        assert records[0].sample_data(0)["GT"] == "0/1"

    def test_gzipped_vcf(self, vep_vcf, tmp_path):
        """Test reading a gzipped VCF."""
        gz_path = tmp_path / "test.vcf.gz"
        with open(vep_vcf) as f_in, gzip.open(gz_path, "wt") as f_out:
            f_out.write(f_in.read())

        reader = VCFReader(str(gz_path))
        assert reader.header.annotation_source == "VEP"
        records = list(reader)
        assert len(records) == 6

    def test_format_fields_in_header(self, vep_vcf):
        reader = VCFReader(vep_vcf)
        assert "GT" in reader.header.format_fields
        assert "AD" in reader.header.format_fields
        assert "DP" in reader.header.format_fields

    def test_edge_case_records(self, edge_cases_vcf):
        reader = VCFReader(edge_cases_vcf)
        records = list(reader)
        # 10 data lines
        assert len(records) == 10

        # Star allele
        star = [r for r in records if r.chrom == "chr6"][0]
        assert star.alt == ["*"]

        # No-call ALT
        nocall = [r for r in records if r.chrom == "chr7" and r.pos == 700][0]
        assert nocall.alt == []
