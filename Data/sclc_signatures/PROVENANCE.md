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

## Neuroendocrine (NE) score template — `zhang_ne_50.tsv` (REQUIRED, not yet vendored)

`ne_score()` implements the 50-gene NE score of:

> Zhang W, Girard L, Zhang Y-A, et al. *Small cell lung cancer tumors and
> preclinical models display heterogeneity of neuroendocrine phenotypes.*
> Transl Lung Cancer Res 2018;7(1):32-49. PMID 29535911. DOI
> 10.21037/tlcr.2018.02.02. Open access: PMC5835590.

NE score per sample = `(cor_NE − cor_nonNE) / 2`, where `cor_NE` / `cor_nonNE`
are Pearson correlations between the sample's expression of the 50 genes and the
NE / non-NE reference vectors.

**This template is not vendored automatically** because the 50 genes and their
NE/non-NE reference values live in the paper's supplementary table (PDF), and we
do not hand-transcribe gene lists (transcription/OCR error risk — see the EMT
`PROVENANCE.md` for why). To enable `ne_score()`, create
`Data/sclc_signatures/zhang_ne_50.tsv` with columns:

```
gene    ne_ref    nonne_ref
ASCL1   <ne mean> <non-NE mean>
...     ...       ...
```

from the Zhang 2018 supplement (25 NE-high + 25 non-NE genes and their
reference expression). Until then `ne_score()` errors with an actionable
message and the EMT×subtype map runs without the NE column.
