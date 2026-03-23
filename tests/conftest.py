"""Shared test fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def vep_vcf() -> str:
    return str(DATA_DIR / "test_vep.vcf")


@pytest.fixture
def snpeff_vcf() -> str:
    return str(DATA_DIR / "test_snpeff.vcf")


@pytest.fixture
def unannotated_vcf() -> str:
    return str(DATA_DIR / "test_unannotated.vcf")


@pytest.fixture
def tmp_maf(tmp_path) -> str:
    return str(tmp_path / "output.maf")


@pytest.fixture
def callers_vcf() -> str:
    return str(DATA_DIR / "test_callers.vcf")


@pytest.fixture
def edge_cases_vcf() -> str:
    return str(DATA_DIR / "test_edge_cases.vcf")


@pytest.fixture
def multi_transcript_vcf() -> str:
    return str(DATA_DIR / "test_multi_transcript.vcf")


@pytest.fixture
def tumor_only_vcf() -> str:
    return str(DATA_DIR / "test_tumor_only.vcf")


@pytest.fixture
def minimal_maf() -> str:
    return str(DATA_DIR / "test_minimal.maf")


@pytest.fixture
def tmp_maf(tmp_path) -> str:
    return str(tmp_path / "output.maf")


@pytest.fixture
def tmp_vcf(tmp_path) -> str:
    return str(tmp_path / "output.vcf")
