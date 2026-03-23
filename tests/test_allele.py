"""Tests for allele trimming and coordinate conversion."""

from __future__ import annotations

import pytest

from vcf2maf_py.converter import (
    determine_variant_type,
    get_maf_coordinates,
    trim_common_bases,
)


class TestTrimCommonBases:
    def test_snp_no_trim(self):
        ref, alt, pos = trim_common_bases("A", "T", 100)
        assert ref == "A"
        assert alt == "T"
        assert pos == 100

    def test_deletion_trim_prefix(self):
        # VCF: POS=100, REF=ACGT, ALT=A
        ref, alt, pos = trim_common_bases("ACGT", "A", 100)
        assert ref == "CGT"
        assert alt == ""
        assert pos == 101

    def test_insertion_trim_prefix(self):
        # VCF: POS=100, REF=A, ALT=ATCG
        ref, alt, pos = trim_common_bases("A", "ATCG", 100)
        assert ref == ""
        assert alt == "TCG"
        assert pos == 101

    def test_complex_trim_both_ends(self):
        # REF=ATCG, ALT=AGCG → shared prefix A, shared suffix CG
        ref, alt, pos = trim_common_bases("ATCG", "AGCG", 100)
        assert ref == "T"
        assert alt == "G"
        assert pos == 101

    def test_mnv_no_trim(self):
        ref, alt, pos = trim_common_bases("AT", "GC", 100)
        assert ref == "AT"
        assert alt == "GC"
        assert pos == 100

    def test_shared_suffix_only(self):
        ref, alt, pos = trim_common_bases("TCA", "GA", 100)
        assert ref == "TC"
        assert alt == "G"
        assert pos == 100


class TestDetermineVariantType:
    def test_snp(self):
        assert determine_variant_type("A", "T") == "SNP"

    def test_dnp(self):
        assert determine_variant_type("AT", "GC") == "DNP"

    def test_tnp(self):
        assert determine_variant_type("ATC", "GCA") == "TNP"

    def test_onp(self):
        assert determine_variant_type("ATCG", "GCTA") == "ONP"

    def test_insertion(self):
        assert determine_variant_type("", "TCG") == "INS"

    def test_deletion(self):
        assert determine_variant_type("CGT", "") == "DEL"

    def test_complex_del(self):
        assert determine_variant_type("TC", "G") == "DEL"

    def test_complex_ins(self):
        assert determine_variant_type("G", "TC") == "INS"


class TestGetMafCoordinates:
    def test_snp(self):
        ref, alt, start, end, vtype = get_maf_coordinates(100, "A", "T")
        assert ref == "A"
        assert alt == "T"
        assert start == 100
        assert end == 100
        assert vtype == "SNP"

    def test_deletion(self):
        ref, alt, start, end, vtype = get_maf_coordinates(100, "ACGT", "A")
        assert ref == "CGT"
        assert alt == "-"
        assert start == 101
        assert end == 103
        assert vtype == "DEL"

    def test_insertion(self):
        ref, alt, start, end, vtype = get_maf_coordinates(100, "A", "ATCG")
        assert ref == "-"
        assert alt == "TCG"
        assert start == 100
        assert end == 101
        assert vtype == "INS"

    def test_dnp(self):
        ref, alt, start, end, vtype = get_maf_coordinates(100, "AT", "GC")
        assert ref == "AT"
        assert alt == "GC"
        assert start == 100
        assert end == 101
        assert vtype == "DNP"
