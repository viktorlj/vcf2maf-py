# vcf2maf-py

A Python reimplementation of [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf) — convert VCF files to [MAF (Mutation Annotation Format)](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/).

## Attribution

This project is a Python port of the original Perl-based **vcf2maf** created by **Cyriac Kandoth** at [Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/). The original tool and its conversion logic, effect prioritization, and variant classification mappings are the foundation of this reimplementation.

> Cyriac Kandoth. mskcc/vcf2maf: vcf2maf v1.6. (2020). doi:[10.5281/zenodo.593251](https://doi.org/10.5281/zenodo.593251)

The original repository is available at [https://github.com/mskcc/vcf2maf](https://github.com/mskcc/vcf2maf) and is licensed under the Apache 2.0 License.

## Features

- **VCF → MAF conversion** with support for VEP and SnpEff annotations
- **MAF → VCF conversion** (reverse direction)
- **Pure Python** — no Perl, no C dependencies; installs with `pip`
- Parses **VEP CSQ** and **SnpEff ANN** annotation fields automatically
- Handles **multi-allelic sites** (splits into separate MAF rows)
- Extracts **genotype depths** from 7+ variant caller FORMAT conventions (GATK, FreeBayes, Strelka, VarScan, Delly, and more)
- **Effect prioritization** matching the original vcf2maf cascade (biotype → severity → canonical → transcript length)
- Maps **130+ Sequence Ontology terms** to MAF Variant_Classification
- Supports **custom isoform overrides** (MSKCC, MANE, or your own)
- Outputs **MAF v2.4** compatible with cBioPortal, maftools, and GDC pipelines
- **VCF inspection** command to check annotation status before conversion
- Reads gzipped VCF/MAF files transparently

## Installation

```bash
pip install vcf2maf-py
```

With optional reference FASTA support for MAF → VCF conversion:

```bash
pip install 'vcf2maf-py[pysam]'
```

## Quick start

### VCF → MAF

```bash
# Basic conversion (auto-detects annotations, samples, and genome build)
vcf2maf-py vcf2maf somatic.vcf -o somatic.maf

# With explicit options
vcf2maf-py vcf2maf somatic.vcf.gz \
    -o somatic.maf \
    --tumor-id TUMOR_SAMPLE \
    --normal-id NORMAL_SAMPLE \
    --ncbi-build GRCh38 \
    --center MSKCC

# Output only core 34 MAF columns
vcf2maf-py vcf2maf somatic.vcf -o somatic.maf --output-columns core

# With custom isoform overrides
vcf2maf-py vcf2maf somatic.vcf -o somatic.maf --isoform-overrides overrides.tsv
```

### MAF → VCF

```bash
# Basic reverse conversion
vcf2maf-py maf2vcf mutations.maf -o mutations.vcf

# With reference FASTA for correct indel padding (requires pysam)
vcf2maf-py maf2vcf mutations.maf -o mutations.vcf --ref-fasta GRCh38.fa
```

### Inspect a VCF

```bash
# Check annotation status, sample names, and variant summary
vcf2maf-py inspect somatic.vcf
```

### Python API

```python
from vcf2maf_py import convert_vcf_to_maf, convert_maf_to_vcf

# VCF → MAF
rows = convert_vcf_to_maf(
    "somatic.vcf",
    output_maf="somatic.maf",
    tumor_id="TUMOR",
    normal_id="NORMAL",
)

# MAF → VCF
convert_maf_to_vcf("mutations.maf", "mutations.vcf")
```

## Annotation support

vcf2maf-py works best with **annotated VCF files**. It auto-detects:

| Annotation tool | INFO field | Status |
|----------------|-----------|--------|
| **Ensembl VEP** | `CSQ` | Fully supported |
| **SnpEff** | `ANN` | Fully supported |
| None | — | Converts coordinates and depths; warns about missing annotations |

If your VCF is not annotated, you can annotate it with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/):

```bash
vep --input_file input.vcf --output_file annotated.vcf \
    --vcf --offline --cache \
    --symbol --canonical --biotype --sift b --polyphen b \
    --fields "Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,CANONICAL,SIFT,PolyPhen,ALLELE_NUM"
```

## How it works

For each variant in the VCF:

1. **Parse alleles** — trim shared prefix/suffix bases, convert VCF-padded indels to MAF coordinate conventions
2. **Detect annotations** — auto-detect VEP CSQ or SnpEff ANN in the VCF header
3. **Parse all transcript effects** — extract every annotated transcript for the variant
4. **Select the best effect** using a priority cascade:
   - Biotype priority (protein_coding > IG/TR genes > ncRNA > pseudogene)
   - Effect severity (~130 SO terms ranked from transcript_ablation to intergenic_variant)
   - Canonical transcript preference
   - Custom isoform overrides
5. **Map SO terms to MAF classifications** (missense_variant → Missense_Mutation, stop_gained → Nonsense_Mutation, etc.)
6. **Extract genotype depths** from FORMAT fields (AD, RO/AO, DP4, Strelka tiers, NR/NV, DV, etc.)
7. **Write MAF** with `#version 2.4` header

## Differences from the original Perl version

| Feature | Original (Perl) | This version (Python) |
|---------|-----------------|----------------------|
| Language | Perl 5 | Python 3.9+ |
| VEP dependency | Runs VEP as part of conversion | Expects pre-annotated VCF (warns if not) |
| Installation | Manual / Docker | `pip install vcf2maf-py` |
| SnpEff support | VEP field mapping only | Native ANN field parsing |
| Reference validation | Via samtools faidx | Optional (via pysam) |
| liftOver | Built-in | Not yet (planned) |
| Python API | None | Full programmatic access |

### Known output differences

Validated against the Perl vcf2maf on ~94,000 variants across 6 VCF files (GRCh37 and GRCh38). Row counts match exactly and 82/109 shared columns are identical. The remaining differences are:

| Difference | Impact | Details |
|-----------|--------|---------|
| `flanking_bps` empty | Low | Perl uses `samtools faidx` to fetch trinucleotide context from a reference FASTA. Not yet implemented. |
| `Entrez_Gene_Id` always 0 | Low | Perl has an internal gene-name-to-Entrez lookup table. Not yet implemented. |
| `FILTER` no `common_variant` | Low | Perl appends `;common_variant` to the FILTER field for variants in common databases. Not implemented. |
| `HGVSp_Short` splice sites | Low | Perl generates synthetic `p.X{position}_splice` for splice donor/acceptor variants. Not implemented (2 of 103 variants in test set). |
| Population AF `0` vs empty | Cosmetic | Python preserves VEP's `0` for absent population frequencies; Perl leaves the field empty. |
| `t_depth` edge case | Rare | For ~5% of indels, the Perl tool computes `t_depth` as `t_ref_count + t_alt_count` while Python uses the VCF `DP` field. Both are valid approaches. |

## Development

```bash
git clone https://github.com/viktorlj/vcf2maf-py
cd vcf2maf-py
uv venv && source .venv/bin/activate
uv pip install -e '.[dev]'
pytest
```

## License

Apache 2.0 — same as the [original vcf2maf](https://github.com/mskcc/vcf2maf).
