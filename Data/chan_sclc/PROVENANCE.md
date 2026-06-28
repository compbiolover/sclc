# Chan 2021 human SCLC atlas — provenance (WS4 external validation)

Reference data for `R/fetch_chan_data.R` and the WS4 cross-sectional validation
(`emt_dispersion_groupwise`). This is the **independent cohort** used to test
whether the EMT-heterogeneity-increases-with-resistance finding from the Stewart
CDX models (`Data/cdx_scrnaseq/`) replicates in **human primary SCLC tumors**.

The expression matrices are **not committed** (large; see `.gitignore`).

## Source

> Chan JM, Quintanal-Villalonga Á, Gao VR, … Rudin CM, Pe'er D, Chan AS, et al.
> **Signatures of plasticity, metastasis, and immunosuppression in an atlas of
> human small cell lung cancer.** *Cancer Cell* 2021;39(11):1479–1496.e18.
> DOI [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008).
> Part of the Human Tumor Atlas Network (HTAN).

scRNA-seq of 155k transcriptomes from 21 human SCLC biospecimens (plus matched
LUAD/normal in the combined atlas), spanning **treatment-naive and treated**
tumors — a cross-platform (10x), between-patient contrast to the CDX models.

## How to obtain (open access)

Processed data is openly downloadable from **CZ CELLxGENE** (no DAC approval;
raw sequencing is dbGaP-controlled but not needed here):

- Collection: <https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec>
- **Combined samples** `.h5ad` (147,137 cells, ~1.46 GB) — what the loader uses:
  `https://datasets.cellxgene.cziscience.com/a9d92e38-9a6e-401b-8484-74bb15122341.h5ad`
- (A SCLC-epithelial-only `.h5ad`, ~9 GB, also exists; the combined file is
  lighter and the loader subsets SCLC tumor cells from it.)
- Authors' repo: <https://github.com/dpeerlab/SCLC_atlas-HTAN>

```r
source("R/fetch_chan_data.R")
chan <- load_chan_sclc("Data/chan_sclc/combined.h5ad", min_umi = 1000)
# -> chan$counts (symbols x cells), chan$cell_meta (cell, sample, donor, group, libsize)
```

## File structure used by the loader (verified against the .h5ad)

- `raw/X` — CSR (`csr_matrix`) **integer counts**, shape `[n_cells, n_genes]`
  (`X` is library-normalized; the loader uses `raw/X`).
- `var/feature_name` — **HGNC symbols** (`var/_index` is Ensembl).
- `obs/disease` — selects SCLC (`"small cell lung carcinoma"`).
- `obs/cell_type_general` — selects malignant cells (`"Epithelial"` ≈ 55,876
  SCLC tumor cells).
- `obs/treatment` — `"Naive"` vs platinum-containing regimens → naive/treated.
- `obs/donor_id`, `obs/HTAN_Biospecimen_ID` — patient / specimen (sample unit).
- `obs/libsize` exists but is **log10-scaled**, so the loader does **not** use it
  for QC; per-cell depth is recomputed from `raw/X` (`Matrix::colSums`) so the
  ≥1000-UMI threshold matches the CDX pipeline exactly.

## Treatment groups (SCLC arm)

`chan_treatment_group()` maps `obs/treatment` to two levels: **naive** (`"Naive"`)
vs **treated** (every platinum-doublet-containing regimen — `Platinum Doublet`,
`…,Immunotherapy`, `…,PARP inhibitor,TMZ`, etc.). The treated group therefore
mixes chemo and chemo+immunotherapy; this is a known limitation of the
validation (the immunotherapy arms affect the microenvironment, though the EMT
axis is scored on the tumor-epithelial compartment only).

## Notes & caveats

- **Design differs from the CDX discovery cohort**: Chan is cross-sectional
  (each tumor is wholly naive or treated, different patients), so the test is
  **unpaired** (`emt_dispersion_groupwise`, Mann-Whitney on per-tumor EMT
  dispersion), not the within-model paired CDX tests.
- **Sample unit**: `HTAN_Biospecimen_ID` (tumor specimen). Biospecimens from the
  same patient are not fully independent; `donor_id` is available for a more
  conservative per-patient sensitivity analysis.
- **Gene IDs**: HGNC symbols via `var/feature_name`; the EMT scorers warn on low
  overlap. Duplicate symbols are collapsed by summing.
- **Versions**: record the CELLxGENE dataset version / download date; CELLxGENE
  re-curates atlases between schema versions.
