# SCLC subtype / neuroendocrine signatures — provenance

Reference data for `R/emt_subtype_map.R`.

## Subtype markers (A / N / P)

The lineage transcription-factor subtypes are called from single marker genes,
so no gene-list file is needed:

| Subtype | Marker | Reference |
|---------|--------|-----------|
| SCLC-A | ASCL1 | Rudin et al. 2019, *Nat Rev Cancer* (PMID 30874629) |
| SCLC-N | NEUROD1 | " |
| SCLC-P | POU2F3 | " |

YAP1 is **not** treated as a subtype: YAP1-expressing SCLC lines were
reclassified as SMARCA4-deficient undifferentiated tumors (Meder et al. 2024,
*Clin Cancer Res*, PMID 38328215). `call_sclc_subtype()` instead emits a
`smarca4_ut_flag` for high-YAP1 / low-A/N/P samples as a QC signal.

## Neuroendocrine (NE) score template — `zhang_ne_50.tsv` (VENDORED)

`ne_score()` implements the 50-gene NE score of:

> Zhang W, Girard L, Zhang Y-A, et al. *Small cell lung cancer tumors and
> preclinical models display heterogeneity of neuroendocrine phenotypes.*
> Transl Lung Cancer Res 2018;7(1):32-49. PMID 29535911. DOI
> 10.21037/tlcr.2018.02.02. Open access: PMC5835590.

NE score per sample = `(cor_NE − cor_nonNE) / 2`, where `cor_NE` / `cor_nonNE`
are Pearson correlations between the sample's expression of the 50 genes and the
NE / non-NE reference vectors.

`zhang_ne_50.tsv` columns: `gene`, `ne_ref` (NE-class mean expression),
`nonne_ref` (non-NE-class mean expression), `class` (NE / non_NE). The 25 NE
and 25 non-NE genes are the top volcano-ranked genes from the paper's NE vs
non-NE comparison (Illumina human WG6 V3 platform); `ne_ref`/`nonne_ref` are the
published NE-class and non-NE-class mean expression values. Values were supplied
directly from the paper's tables (not OCR-transcribed by us).

Note: a few symbols are as-published 2018 aliases (e.g. `CYR61` = CCN1,
`TMSB15A`/`TMSB15B`); map to current HGNC at scoring time if your matrix uses
newer symbols. `ne_score()` matches case-insensitively and warns on low overlap.
