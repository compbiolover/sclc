# EMT signature gene sets — provenance

These gene sets are the canonical, published EMT signatures used by
`R/emt_scoring.R`. They are **vendored** (committed) so that EMT scoring is
fully reproducible offline. They are produced by `R/fetch_emt_signatures.R`,
which downloads them from the authoritative sources below — the lists are never
hand-typed.

To (re)generate:

```r
source("R/fetch_emt_signatures.R")
fetch_emt_signatures()          # writes the *.tsv files in this folder
```

All scores in `R/emt_scoring.R` are oriented so that **higher = more
mesenchymal**.

| File | Method | Source | Citation |
|------|--------|--------|----------|
| `byers_76gs.tsv` | 76GS | Cancer-Systems-Biology-Lab/EMT_Scoring_scRNA `Gene_signatures/76GS/EMT_signature_76GS.xlsx` (col 2) | Byers LA et al. *Clin Cancer Res* 2013;19(1):279-90. PMID 23172883. DOI 10.1158/1078-0432.CCR-12-1558 |
| `tan_ks_epithelial.tsv`, `tan_ks_mesenchymal.tsv` | KS | Cancer-Systems-Biology-Lab/EMT_Scoring_scRNA `Gene_signatures/KS/EM_gene_signature_tumor_KS.xlsx` (col1 gene, col2 Epi/Mes) | Tan TZ et al. *EMBO Mol Med* 2014;6(10):1279-93. PMID 24711451. DOI 10.15252/emmm.201404208 |
| `george_mlr.tsv` | MLR (predictor genes only) | Cancer-Systems-Biology-Lab/EMT_Scoring_scRNA `Gene_signatures/MLR/genes_for_EMT_score.txt` | George JT et al. *Cancer Res* 2017;77(22):6415-25. PMID 29038175. DOI 10.1158/0008-5472.CAN-17-1679 |
| `hallmark_emt.tsv` | Hallmark EMT | MSigDB via `msigdbr` (collection "H", `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`) | Liberzon A et al. *Cell Syst* 2015;1(6):417-25. PMID 26771021. DOI 10.1016/j.cels.2015.12.004 |

## Notes & caveats

- **Source repo**: the Jolly-lab `Cancer-Systems-Biology-Lab/EMT_Scoring_scRNA`
  repository (branch `master`) packages the original Byers/Tan/George signatures
  in a single place and is the reference implementation for these scores
  (see also Chakraborty et al., *Front Bioeng Biotechnol* 2020, PMID 32258020,
  for the comparative study of all three metrics).
- **KS variant**: by default we vendor the **tumor** signature
  (`EM_gene_signature_tumor_KS.xlsx`), appropriate for patient cohorts (George
  2015, GSE60052). For cell-line panels (DepMap/GDSC, WS3 drug response), re-run
  `fetch_emt_signatures(ks_variant = "cellLine")`.
- **MLR is incomplete by design**: `genes_for_EMT_score.txt` provides the MLR
  **predictor genes** but not the published multinomial **class coefficients**.
  `score_mlr()` therefore stays disabled (errors with an actionable message)
  until a `george_mlr.tsv` with `coef_*` columns is supplied. The consensus EMT
  axis is built from 76GS + KS + Hallmark, which is sufficient and robust.
- **Gene IDs** are HGNC symbols. When scoring a matrix that uses Ensembl IDs,
  map to symbols first (the scorers warn on low overlap and error if < 3 genes
  match).
- **Versioning**: record the `msigdbr` package version and the date of fetch in
  your analysis log; MSigDB updates the Hallmark sets between releases.
