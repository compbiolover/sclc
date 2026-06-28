# Clinical / survival data — provenance

## `gse60052_survival.tsv` — GSE60052 (Jiang 2016) overall survival

Per-sample overall-survival data for the GSE60052 SCLC tumor cohort, used by the
WS3 survival arm (`R/emt_survival.R`) for the pooled George + GSE60052 analysis.

| column | meaning |
|--------|---------|
| `sample` | GEO expression sample ID (matches the columns of the GSE60052 expression matrix) |
| `os_months` | overall survival / follow-up time, months |
| `os_status` | 1 = deceased (event), 0 = alive (censored) |

**Source.** Jiang L, Huang J, Higgs BW, et al. *Genomic Landscape Survey
Identifies SRSF1 as a Key Oncodriver in Small Cell Lung Cancer.* PLoS Genetics
2016;12(4):e1005895. DOI 10.1371/journal.pgen.1005895 (open access, CC BY 4.0).
Derived from **S1 Table** ("Chinese patient clinical summary (n=99)"),
supplementary file `pgen.1005895.s012.xlsx`, columns `Months` and
`Survival Status`, restricted to patients with tumor RNA-seq (the n=48 with
matched expression + follow-up; the paper's survival KM also used n=48).

**ID crosswalk.** S1 Table patient IDs were mapped to GEO expression sample
names by normalizing both to a common key: strip a leading `B`/`T`/`N` before a
digit and a trailing `A`/`TA` (e.g. S1 `08-3503` -> GEO `B08-3503TA`; S1 `5` ->
`5A`; S1 `11-548` -> `T11-548A`). All 50 S1 RNA-seq rows mapped uniquely; 2 were
dropped for missing time/status, leaving 48. On the single key collision
(`10-5656`, present as both `B10-5656A` and `T10-5656A` in GEO) the non-`T`
(original) sample was kept.

Redistributed under the article's CC BY 4.0 licence with attribution.
