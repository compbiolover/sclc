# SCLC CDX single-cell RNA-seq — provenance (WS4: dynamic resistance)

Reference data for `R/emt_resistance.R` and `R/run_resistance.R`. This is the
dataset behind the WS4 question: **does EMT arise in SCLC as a function of
treatment?** — i.e. does the EMT cell-state axis (WS1) rise as a tumor evolves
from chemo-sensitive to chemo-resistant?

The expression matrices themselves are **not committed** (large; see
`.gitignore`). This file documents exactly what to obtain and how the loader
expects it laid out, so WS4 is reproducible.

## Source

> Stewart CA, Gay CM, Xi Y, Sivajothi S, Sivakamasundari V, Fujimoto J,
> Bolisetty M, Hartsfield PM, Balasubramaniyan V, Chalishazar MD, Moran C,
> Kalhor N, Stewart J, Tran H, Swisher SG, Roth JA, Zhang J, de Groot J,
> Glisson B, Oliver TG, Heymach JV, Wistuba I, Robson P, Wang J, Byers LA.
> **Single-cell analyses reveal increased intratumoral heterogeneity after the
> onset of therapy resistance in small-cell lung cancer.**
> *Nature Cancer* 2020;1(4):423–436.
> PMID **33521652** · PMC **PMC7842382** · DOI
> [10.1038/s43018-019-0020-z](https://doi.org/10.1038/s43018-019-0020-z).
> *(Citation/identifiers retrieved from PubMed.)*

Senior author L.A. Byers is the same Byers as the vendored **76GS** EMT
signature (`Data/emt_signatures/byers_76gs.tsv`, Byers et al. 2013), so this
cohort is a natural home for the EMT axis.

### Why this dataset answers the WS4 question

SCLC relapses to chemoresistance fast, but post-relapse human tissue is scarce.
The authors generated **circulating-tumor-cell–derived xenografts (CDX)** and
profiled them by single-cell RNA-seq in **paired chemo-sensitive vs
chemo-resistant** states of the *same* model. Their headline result is that
intratumoral **heterogeneity increases** with resistance, with **EMT** named
explicitly as one of the heterogeneously expressed resistance programs
(abstract). That paired design (sensitive vs resistant within a model) is what
lets WS4 test the axis with the *model*, not the cell, as the unit of
replication.

## GEO accessions

| Accession | Role | Samples | Notes |
|-----------|------|---------|-------|
| **GSE138474** | SuperSeries | 43 | Umbrella record; download via its sub-series. |
| **GSE138267** | CDX **single-cell** RNA-seq `[RNA-Seq]` | 25 | The series WS4 uses (per-cell). |
| **GSE138418** | CDX RNA-seq `[RNA-Seq]` | 18 | Companion CDX RNA-seq series. |

Sequencing platform per GEO: **Illumina HiSeq 2000 (Homo sapiens)**. Each GSM
ships a **10x-style triplet** (`barcodes.tsv` / `genes.tsv` / `matrix.mtx`)
inside `<GSM>_<title>.matrix.tar.gz`; `genes.tsv` is `<Ensembl>\t<HGNC symbol>`
on the **hg19** reference (32,738 genes, **human only**), and `matrix.mtx` is
genes × cells integer counts (~6k cells/sample).

### Sample structure (verified against GEO, and encoded in `R/fetch_cdx_data.R`)

The 25 GSMs are CDX libraries titled `<MODEL><variant>.<LIBRARY>`. **The
biological model is the `SC<number>` stem, NOT the dotted prefix** — the drug
tags (`_LY2606368`, `_Talazoparib`, `_cr`, `cis`, `-1/-2/-3`) are
treatment-derived *variants of the same model*. Keying on the dotted prefix
would split `SC4`, `SC4_LY2606368`, `SC4_Talazoparib` into three unpaired
singletons and destroy the pairing. The `treatment` characteristic gives the
state:

| `treatment` value | Condition |
|-------------------|-----------|
| `untreated`, `vehicle-treated` | **sensitive** |
| `… , relapsed` (cisplatin / prexasertib / talazoparib) | **resistant** |
| bare `cisplatin` (on-treatment, not yet relapsed) | **excluded** (NA) — neither a naive baseline nor an established resistant line |

This yields **4 paired models** (sensitive + resistant): **SC4, SC53, SC55,
SC68**, and **4 sensitive-only** models (SC16, SC39, SC49, SC75; excluded from
paired tests). One on-treatment sample (`SC55-2`, cisplatin) is excluded by
default. `cdx_classify_samples()` implements exactly this mapping and is unit-
tested against the real sheet in `tests/testthat/test-fetch_cdx_data.R`.

## How to obtain

`R/fetch_cdx_data.R` automates the whole path (needs network + GEOquery + Matrix):

```r
source("R/fetch_cdx_data.R")
fetch_cdx_data()        # GEO -> Data/cdx_scrnaseq/{cdx_counts.rds, cell_metadata.tsv}
# or step by step, to inspect/override the sensitive<->resistant mapping:
tab <- cdx_sample_table("GSE138267")      # tidy table: gsm,title,model,condition,paired,...
cdx <- load_cdx_counts(tab)               # download + read + assemble (sensitive/resistant only)
write_cdx_inputs(cdx)                      # write the two run_resistance.R inputs
```

This produces the two WS4 inputs:

1. `cdx_counts.rds` — a **genes-in-rows** (HGNC symbol) sparse counts matrix →
   scored per cell by `score_emt_singlecell()` (`R/emt_scoring.R`, WS1);
2. `cell_metadata.tsv` — per-cell `cell`, `gsm`, `library`, `model`, `condition`
   (already resolved to sensitive/resistant, so no keyword guessing downstream).

`run_resistance.R` reads both, scores, and runs the three WS4 tests.

## Layout produced by the loader / expected by the runner

```
Data/cdx_scrnaseq/
  PROVENANCE.md            <- this file (committed)
  download/                <- per-GSM tarballs + extracted triplets (gitignored)
  cdx_counts.rds           <- genes x cells sparse matrix (gitignored)
  cell_metadata.tsv        <- cell, gsm, library, model, condition (gitignored)
```

## Notes & caveats

- **Gene IDs** must be HGNC symbols for the WS1 EMT scorers; map Ensembl →
  symbol first if needed (the scorers warn on low overlap and error if < 1–3
  genes match).
- **Xenograft contamination**: the deposited matrices are aligned to a
  **human-only** hg19 reference (no `ENSMUSG` genes), so mouse stromal reads are
  not represented as mouse genes. Mouse cells could still align spuriously to
  human orthologs; the EMT axis is tumor-intrinsic, but very low-count or
  clearly non-tumor cells should be filtered (`load_cdx_counts()` keeps all
  cells — apply QC upstream if needed).
- **Unit of replication**: there are only a handful of paired CDX models, so the
  per-cell counts are large but the *paired* sample size is small. WS4 reports
  per-model deltas and sign-consistency alongside the paired Wilcoxon p-value
  for exactly this reason; do not over-read a single small-n p-value.
- **Record versions**: log the GEO download date and the GEOquery/Bioconductor
  versions in your analysis log; GEO supplementary processing can be revised.
