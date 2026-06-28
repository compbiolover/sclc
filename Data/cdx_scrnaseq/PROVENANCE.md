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

Sequencing platform per GEO: **Illumina HiSeq 2000 (Homo sapiens)**. The
single-cell capture chemistry and per-cell processing are described in the
paper's Methods; consult them before re-deriving counts (do not assume a
droplet/10x pipeline).

CDX models in the SuperSeries include **SC4, SC16, SC39, SC49, SC53, SC55,
SC68, SC75** and **TC568c**, with treatment-naive / chemo-sensitive baselines
and chemo-resistant (relapsed) counterparts. Resistant arms in the GEO sample
titles are denoted variously (e.g. a `cr` tag, or the agent used to drive
resistance such as `LY2606368` (prexasertib) or `Talazoparib`); **verify the
exact sensitive↔resistant label for each sample against the GEO sample sheet**
before scoring — pass that mapping to `prepare_resistance_emt()` via
`sensitive_labels` / `resistant_labels` rather than relying on keyword guessing.

## How to obtain

In R, the cleanest route is Bioconductor **GEOquery** (a soft dependency; not
required by the WS4 analysis functions, only to fetch raw data):

```r
# install.packages("BiocManager"); BiocManager::install("GEOquery")
library(GEOquery)
# Supplementary files (processed count matrices) for the single-cell sub-series:
getGEOSuppFiles("GSE138267", baseDir = "Data/cdx_scrnaseq")
# Sample-level metadata (titles, characteristics) to build the cell -> model /
# condition map:
gse <- getGEO("GSE138267", GSEMatrix = TRUE)
pData(gse[[1]])           # inspect sample titles & characteristics here
```

Then build the two WS4 inputs:

1. a **genes-in-rows expression matrix** (HGNC symbols) → score per cell with
   `score_emt_singlecell()` (from `R/emt_scoring.R`, WS1);
2. a **per-cell metadata** data.frame with a `model` column and a `condition`
   column (sensitive/resistant), derived from the GEO sample sheet.

`prepare_resistance_emt()` joins those two; see `R/run_resistance.R` for the
wiring.

## Layout expected by the loader / runner

```
Data/cdx_scrnaseq/
  PROVENANCE.md            <- this file (committed)
  GSE138267/               <- downloaded matrices (NOT committed; gitignored)
  cell_metadata.tsv        <- cell, model, condition (NOT committed; gitignored)
```

## Notes & caveats

- **Gene IDs** must be HGNC symbols for the WS1 EMT scorers; map Ensembl →
  symbol first if needed (the scorers warn on low overlap and error if < 1–3
  genes match).
- **Xenograft contamination**: CDX scRNA-seq can carry mouse stromal reads.
  Use the authors' human-filtered matrices, or filter to human cells, before
  scoring — mouse-derived cells would corrupt both the EMT axis and the
  per-model summaries.
- **Unit of replication**: there are only a handful of paired CDX models, so the
  per-cell counts are large but the *paired* sample size is small. WS4 reports
  per-model deltas and sign-consistency alongside the paired Wilcoxon p-value
  for exactly this reason; do not over-read a single small-n p-value.
- **Record versions**: log the GEO download date and the GEOquery/Bioconductor
  versions in your analysis log; GEO supplementary processing can be revised.
