# Data

Planned source: Breast Invasive Carcinoma (TCGA, PanCancer Atlas), via cBioPortal.

We will verify the study, expression measurements, sample identifiers, and subtype column before downloading and analyzing data.

Original downloads belong in `raw/`. Cleaned copies belong in `processed/`. Both folders are excluded from Git.

## Initial inspection (October 9, 2026)

Verified study ID: `brca_tcga_pan_can_atlas_2018`.

The `SUBTYPE` field is patient-level. Before matching to expression data, its recorded values are:

| Label | Count |
| --- | ---: |
| BRCA_LumA | 499 |
| BRCA_LumB | 197 |
| BRCA_Basal | 171 |
| BRCA_Her2 | 78 |
| BRCA_Normal | 36 |
| No subtype record | 103 |

The study contains 1,084 patients and 1,084 samples. These are not final analysis counts: expression availability and sample matching still need checking. No exclusions have been applied.

Proposed expression profile: `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna`, described by cBioPortal as batch-normalized RSEM expression from Illumina HiSeq RNASeqV2. This is already source-processed expression, not raw sequencing counts. Prefer it to the portal's precomputed cohort z-scores so our own standardization can be fitted within training folds. Source batch normalization remains a limitation to document.

To repeat the inspection from the project root:

```powershell
.\.venv\Scripts\python.exe src/inspect_dataset.py
```

This retrieves metadata and subtype records only, not the expression matrix. Untouched responses go in ignored `data/raw/`; the aggregate report is `data/dataset_inspection.json`.

Source: [cBioPortal study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018), [public API](https://www.cbioportal.org/api).
