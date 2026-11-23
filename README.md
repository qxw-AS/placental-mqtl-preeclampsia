# Placental cis-mQTL and preeclampsia: SMR/HEIDI analysis code

Analysis code for:

> Li T, Qiu X. Genetically regulated placental DNA methylation and preeclampsia risk:
> a summary-data Mendelian randomization study. *BMC Genomics* (under review).

This repository contains **code only**. Figures and Supplementary Tables S1–S6 are part of
the journal submission and are not duplicated here, so that there is a single authoritative
version of each. All inputs are public summary statistics; no individual-level or
controlled-access data are redistributed.

---

## What the pipeline does

The public PACE resource releases one row per significant SNP–CpG pair. Everything
downstream works at the CpG level, so the first step collapses each CpG to its strongest
pair, and the second screens that instrument set against three outcome GWASs.

```
data/raw/pace_cis_mqtl_top_snps.tsv                       655,543 pairs
  │  01_msites.py       keep strongest pair per CpG
  ├─ 288,760 CpG-level mSites
  │  02_smr_screen.py   liftOver hg19→hg38, drop palindromic, join on position,
  │                     orient alleles, b_SMR = b_GWAS / b_mQTL
  ├─ 288,714 mapped  →  248,009 instruments  →  242,924 testable (FinnGen PE)
  ├─ 5 SMR hits
  │  03_heidi.py        BESD store from 1000G rsIDs, then SMR v1.3.1
  ├─ P_SMR-multi, P_HEIDI, n SNP per CpG (1 of the 5 is nominal by HEIDI)
  │  04_table2.py       assemble Table 2 and compare with the baseline
  └─ PASS / FAIL
```

Two implementation points that matter for anyone re-running this:

- **Harmonisation is by GRCh38 position, not by rsID.** The PACE export carries rsIDs only
  inside the composite `snp` key, so `02_smr_screen.py` joins on `(chrom, pos38)`. rsIDs are
  assigned from the 1000G `.bim` at step 3, where the SMR binary requires them.
- **The genome-wide screen is computed in Python.** `b_SMR` uses the SMR ratio estimator with
  the delta-method SE, applied to all 242,924 pairs in one pass. The SMR binary is used for
  BESD construction, HEIDI and SMR-multi. The two implementations agree on all five hits to
  within 1% in `P_SMR`; the residual difference comes from the binary writing its inputs at
  six significant digits, and Table 2 reports the binary's values.

---

## Contents

```
.
├── README.md
├── LICENSE
├── .gitignore
├── environment/
│   ├── requirements.txt          Python package versions
│   └── smr_version.txt           SMR binary provenance and checksum
├── data/
│   └── ref/
│       └── cpg_annotation.tsv    hg19 probe coordinates for the step-3 CpGs
├── scripts/
│   ├── 00_fetch_data.sh          download inputs
│   ├── 01_msites.py              PACE pairs -> CpG-level mSites
│   ├── 02_smr_screen.py          liftOver, harmonisation, genome-wide SMR screen
│   ├── 03_heidi.py               BESD inputs + HEIDI / SMR-multi
│   └── 04_table2.py              Table 2 assembly and acceptance check
└── results/
    └── table2_candidate_cpgs.tsv acceptance baseline (5 CpGs)
```

---

## Requirements

- Python 3.14.8, packages pinned in `environment/requirements.txt` (numpy, pandas, scipy,
  pyliftover)
- SMR v1.3.1 (Linux x86-64, build 2024-03-07) — **not redistributed**; download from
  https://cnsgenomics.com/software/smr/, place at `bin/smr`, `chmod +x bin/smr`
- ~20 GB free disk for the downloaded summary statistics (the PACE export dominates)

`bin/smr` is git-ignored. Record its checksum with `sha256sum bin/smr` and paste the value
into `environment/smr_version.txt` — a different build can shift `P_HEIDI` enough to change
a classification, and `cg07027493` sits at P_HEIDI = 0.088, close to the nominal boundary.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r environment/requirements.txt

bash   scripts/00_fetch_data.sh                 # see "Input data" for the manual exports
python scripts/01_msites.py
python scripts/02_smr_screen.py                 # add --strict to enforce the counts
python scripts/03_heidi.py --run                # omit --run to write inputs only
python scripts/04_table2.py
```

`02_smr_screen.py` prints every intermediate count next to the value reported in the paper.
With `--strict` a mismatch becomes a fatal error. Do not edit the expected values to make it
pass: a mismatch means an upstream file changed, and that is the finding.

## Expected output

| Step | Quantity | Value |
| --- | --- | --- |
| 01 | input pairs | 655,543 |
| 01 | CpG-level mSites | **288,760** |
| 02 | liftOver mapped (unmapped) | **288,714** (46) |
| 02 | palindromic removed | 40,705 (14.10% of mapped) |
| 02 | instruments entering the screen | **248,009** |
| 02 | testable, FinnGen PE | **242,924** |
| 02 | Bonferroni threshold | 2.016 × 10⁻⁷ |
| 02 | SMR hits | 5 |
| 03 | CpGs carried to HEIDI | 5 |
| 04 | cells compared | 50 numeric cells + HEIDI class |

The five hits and their `P_SMR`: cg20620484 (4.32e-09), cg07980026 (8.23e-09),
cg07027493 (8.47e-08), cg06204421 (9.13e-08), cg06177599 (1.57e-07).

Two intermediate counts are less stable than the others and are worth watching:

- **288,714 / 46 unmapped** depends on the release of `hg19ToHg38.over.chain`. UCSC updates
  the chain; a later release can map a handful of variants differently.
- **40,705 palindromic** depends on which pair is the strongest for each CpG, and therefore
  on the PACE export date. Recomputing from the export distributed with the manuscript gives
  40,756. The difference is 51 CpGs and does not propagate: the instrument count
  (248,009) and the testable count (242,924) are the values the pipeline is checked against.

## Acceptance test

`results/table2_candidate_cpgs.tsv` holds the values published in Table 2. `04_table2.py`
regenerates them and compares every cell against the baseline:

```bash
python scripts/04_table2.py
#   rows 5; cells 50
#   P values: |delta log10| <= 0.01
#   PASS: reproduced
```

P values are compared in log10 space, because a 3-significant-digit baseline cannot support
a relative comparison at 10⁻²⁰⁰. Effects, confidence limits and `P_HEIDI` are compared at
0.1% relative or 5 × 10⁻⁴ absolute, matching the four-decimal storage in the baseline; the
flags `--plog`, `--rtol` and `--atol` change these. `heidi_class` must match exactly.

A uniform shift in `p_smr` across all five CpGs points to a different SMR build. A change in
a single CpG points to a changed input file.

---

## Input data

`scripts/00_fetch_data.sh` downloads the machine-retrievable files. Three inputs have no
bulk endpoint and must be obtained by hand. Nothing here is committed to Git.

| Resource | Version marker | Where |
| --- | --- | --- |
| PACE placental cis-mQTL top-SNP table | 2025-12-09 | https://pace-placenta-mqtl.streamlit.app — **manual export**, save as `data/raw/pace_cis_mqtl_top_snps.tsv` |
| FinnGen R13 `O15_PREECLAMPS` | server mtime 2026-05-05 | https://r13.finngen.fi → `data/raw/` |
| FinnGen R13 `O15_PREEC_OR_FETGRO` | as above | https://r13.finngen.fi → `data/raw/` |
| GWAS Catalog `GCST90301704` | 2024-03-18 | https://www.ebi.ac.uk/gwas/studies/GCST90301704 → `data/raw/GCST90301704.h.tsv.gz` |
| UCSC `hg19ToHg38.over.chain.gz` | — | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ → `data/ref/` |
| 1000G Phase 3 EUR plink `.bim` | — | ships with the SMR reference data; per-chromosome files in `data/ref/1000G_EUR_Phase3_plink/` |
| 1000G Phase 3 EUR (489) | Zenodo 7768714 | https://zenodo.org/records/7768714 — LD reference for SMR-multi |

All resources were accessed on 16–17 November 2026 (China Standard Time). File sizes and
MD5 checksums are recorded in Table S5 of the manuscript.

**Caveat on the exposure.** The PACE resource is a live service with no bulk endpoint and
no versioned release. Two users querying it months apart can receive different tables. If
your counts differ from the expected output, check the PACE export first.

`data/ref/cpg_annotation.tsv` is the one committed file under `data/`: probes, not
observations. It holds `cpg, chr, pos_hg19` for the CpGs that step 2 carries into step 3,
taken from the Illumina EPIC manifest and identical to the `ProbeBp` column of the BESD
description. Regenerate it from the manifest if a different hit set is being analysed.

---

## Parameters fixed in the manuscript

| Parameter | Value |
| --- | --- |
| Bonferroni threshold | 0.05 / 248,009 = 2.016 × 10⁻⁷ |
| Testable pairs | 242,924 (after allele and frequency filtering) |
| HEIDI criterion | P > 0.05 not rejected; 0.05 < P < 0.10 reported as nominal |
| SMR-multi | sensitivity layer; minimum set P = 3.71 × 10⁻⁷ |
| LD reference | 1000 Genomes Phase 3 EUR, 489 individuals |

## Not reproduced by this code

- **Figure rendering** — Figs. 1–4 and the supplementary figures are built by separate
  plotting scripts that depend on the annotation and eFORGE layers, not on this pipeline.
- **eFORGE v2.0 enrichment** — a web service with no API. Parameters are in the Methods
  (top 10,000 mSites, 850k background, 200 permutations; 15-state analysis on the top
  2,000 because of server memory limits).
- **Colocalization** (`coloc`, PP.H4 > 0.8) — not evaluable: PACE releases only significant
  SNP–CpG pairs, so the full cis-region statistics needed for PP.H4 do not exist.
- **Conditional analysis** — requires GCTA and individual-level data.
- **Per-cohort mQTL sensitivity** — per-cohort summary statistics are not published.
- **Tissue specificity via blood mQTL** — the Hannon 2018 resource offers interactive
  queries only, no batch export.

---

## Licence

Code: MIT (see `LICENSE`). `results/table2_candidate_cpgs.tsv`: CC BY 4.0.
Input data remain under the terms of their original providers and are not covered by
these licences.

## Contact

Xiaowei Qiu — qxw24657@163.com — ORCID 0009-0003-8966-8073
