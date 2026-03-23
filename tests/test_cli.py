"""Tests for the CLI interface."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from vcf2maf_py.cli import main


@pytest.fixture
def runner():
    return CliRunner()


class TestVcf2MafCLI:
    def test_help(self, runner):
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "vcf2maf-py" in result.output

    def test_vcf2maf_help(self, runner):
        result = runner.invoke(main, ["vcf2maf", "--help"])
        assert result.exit_code == 0
        assert "input_vcf" in result.output.lower() or "INPUT_VCF" in result.output

    def test_vcf2maf_basic(self, runner, vep_vcf, tmp_path):
        output = str(tmp_path / "out.maf")
        result = runner.invoke(main, ["vcf2maf", vep_vcf, "-o", output])
        assert result.exit_code == 0
        assert Path(output).exists()
        assert "Converted" in result.output

    def test_vcf2maf_with_options(self, runner, vep_vcf, tmp_path):
        output = str(tmp_path / "out.maf")
        result = runner.invoke(main, [
            "vcf2maf", vep_vcf,
            "-o", output,
            "--tumor-id", "TUMOR",
            "--normal-id", "NORMAL",
            "--ncbi-build", "GRCh37",
            "--center", "TEST",
        ])
        assert result.exit_code == 0

    def test_vcf2maf_core_columns(self, runner, vep_vcf, tmp_path):
        output = str(tmp_path / "out.maf")
        result = runner.invoke(main, [
            "vcf2maf", vep_vcf,
            "-o", output,
            "--output-columns", "core",
        ])
        assert result.exit_code == 0
        with open(output) as fh:
            fh.readline()  # skip version
            header = fh.readline().strip().split("\t")
            # Core has 34 columns
            assert len(header) == 34

    def test_maf2vcf_help(self, runner):
        result = runner.invoke(main, ["maf2vcf", "--help"])
        assert result.exit_code == 0

    def test_maf2vcf_basic(self, runner, vep_vcf, tmp_path):
        maf = str(tmp_path / "test.maf")
        vcf_out = str(tmp_path / "test.vcf")
        # First create a MAF
        runner.invoke(main, ["vcf2maf", vep_vcf, "-o", maf])
        # Then convert back
        result = runner.invoke(main, ["maf2vcf", maf, "-o", vcf_out])
        assert result.exit_code == 0
        assert Path(vcf_out).exists()

    def test_inspect_command(self, runner, vep_vcf):
        result = runner.invoke(main, ["inspect", vep_vcf])
        assert result.exit_code == 0
        assert "VEP" in result.output
        assert "TUMOR" in result.output

    def test_inspect_unannotated(self, runner, unannotated_vcf):
        result = runner.invoke(main, ["inspect", unannotated_vcf])
        assert result.exit_code == 0
        assert "NONE" in result.output

    def test_version(self, runner):
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output
