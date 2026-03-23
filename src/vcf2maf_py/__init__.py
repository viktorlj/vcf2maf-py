"""vcf2maf-py — Convert VCF files to MAF format.

A Python reimplementation of mskcc/vcf2maf
(https://github.com/mskcc/vcf2maf) by Cyriac Kandoth.
"""

from __future__ import annotations

__version__ = "0.1.2"

from vcf2maf_py.converter import convert_vcf_to_maf, load_isoform_overrides
from vcf2maf_py.maf2vcf import convert_maf_to_vcf

__all__ = [
    "__version__",
    "convert_vcf_to_maf",
    "convert_maf_to_vcf",
    "load_isoform_overrides",
]
