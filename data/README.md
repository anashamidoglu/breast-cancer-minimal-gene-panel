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

The study contains 1,084 patients and 1,084 samples. These are not final analysis counts: expression availability and sample matching still need checking.

With user approval on October 9, 2026, we excluded the 36 Normal-like labels and 103 missing labels from the prepared label table, retaining 945 patients. Original responses remain unchanged. `data/processed/patient_labels.csv` contains retained labels; `data/processed/label_exclusions.csv` records each excluded identifier and reason. These files stay local. The tracked aggregate audit is `data/label_preparation.json`.

After inspection, prepare the labels with:

```powershell
.\.venv\Scripts\python.exe src/prepare_labels.py
```

Open `notebooks/01_inspect_subtype_labels.ipynb` to see the count table and chart.

Proposed expression profile: `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna`, described by cBioPortal as batch-normalized RSEM expression from Illumina HiSeq RNASeqV2. This is already source-processed expression, not raw sequencing counts. Prefer it to the portal's precomputed cohort z-scores so our own standardization can be fitted within training folds. Source batch normalization remains a limitation to document.

To repeat the inspection from the project root:

```powershell
.\.venv\Scripts\python.exe src/inspect_dataset.py
```

This retrieves metadata and subtype records only, not the expression matrix. Untouched responses go in ignored `data/raw/`; the aggregate report is `data/dataset_inspection.json`.

Source: [cBioPortal study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018), [public API](https://www.cbioportal.org/api).

## Expression download and matching

```powershell
.\.venv\Scripts\python.exe src/download_expression.py
.\.venv\Scripts\python.exe src/match_expression.py
```

The downloader uses the official [cBioPortal DataHub](https://github.com/cBioPortal/datahub/tree/master/public/brca_tcga_pan_can_atlas_2018), pins a repository commit, and verifies the expression file against its Git LFS SHA-256 checksum and byte count. Source URLs and checksums are recorded in `expression_download.json`. Re-running against a changed DataHub commit may retrieve a newer snapshot; use the recorded commit and URL to recover this exact version. The portal sample mapping is retrieved separately and its checksum is also recorded.

The source matrix has 20,531 measurement rows and 1,082 expression samples. All 945 approved patient labels match expression samples. Each matched patient has one primary solid tumor sample. No expression values are missing, nonfinite, or negative.

Twenty-four Entrez gene IDs are repeated across 50 rows with distinct measurements. We exclude all 50 ambiguous rows based solely on identifiers, leaving **20,481 unique gene features**. No measurements are averaged. Thirteen rows lack gene symbols but have unique Entrez IDs and remain included. This conservative rule may remove potentially useful genes and is recorded as a limitation.

Prepared local files:

- `data/processed/expression.parquet`: 945 sample rows by 20,481 gene columns, with sample IDs as the index. Parquet is a compact table format.
- `data/processed/sample_labels.csv`: the matching sample IDs, patient IDs, and subtype labels.
- `data/processed/gene_annotations.csv`: each feature's Entrez ID and source gene symbol.
- `data/processed/ambiguous_gene_annotations.csv`: the 50 excluded source annotations.
- `data/processed/expression_exclusions.csv`: patients without expression data (currently empty).

The aggregate audit is `expression_matching.json`. Original and processed patient-level files remain excluded from Git. No log transformation, scaling, variance filtering, feature selection, or model fitting has been performed. Source batch normalization was performed before we obtained the data and cannot be refitted within our training folds.

## Frozen split and training-only baseline

```powershell
.\.venv\Scripts\python.exe src/create_split.py
.\.venv\Scripts\python.exe src/evaluate_baseline.py
```

The split uses `train_test_split(test_size=0.20, stratify=labels['subtype'], random_state=42)` on sample labels sorted by sample ID. It does not load expression values. Each patient has one sample; uniqueness, class coverage, and patient disjointness are checked. The resulting membership lists are saved locally and the script refuses to overwrite a changed split.

| Subtype | Training | Test |
| --- | ---: | ---: |
| Luminal A | 399 | 100 |
| Luminal B | 158 | 39 |
| Basal-like | 137 | 34 |
| HER2-enriched | 62 | 16 |
| Total | 756 | 189 |

`data/processed/train_labels.csv` and `test_labels.csv` contain frozen membership and labels. `split_report.json` records aggregate counts, seed, and membership checksums. Test label counts are checked for stratification; no test predictions or performance have been inspected.

The dummy baseline always predicts the most frequent subtype in its fitting data. It uses five training-only folds from `StratifiedKFold(n_splits=5, shuffle=True, random_state=42)`. Each training sample receives one prediction from a model fitted without it. The baseline script does not open test labels or any expression matrix. The fold assignments are saved locally in `training_cv_folds.csv` for subsequent model comparisons.

Training cross-validation balanced accuracy is 0.25, macro F1 is approximately 0.173, and each one-vs-rest AUROC is 0.50. Recall is 1.0 for Luminal A and 0.0 for the other subtypes. These are baseline results, not final test results. Reports are `baseline_report.json` and `baseline_cv_folds.csv`.
