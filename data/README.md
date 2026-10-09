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

## Initial all-available-gene logistic regression

```powershell
.\.venv\Scripts\python.exe src/evaluate_logistic.py
```

This preliminary model uses fixed L2 regularization (`C=1`, `l1_ratio=0`), balanced class weights, and the saved five training folds. It reads only training expression rows using a Parquet filter. Test labels are not opened. The Pipeline performs log1p, variance filtering at `1e-8` on the log scale, and standardization within each fitting fold. The filter retains 20,156–20,170 features across folds. All five solver fits converged; a convergence warning would halt the script.

Mean training CV balanced accuracy: 0.8766, fold SD: 0.0400. Mean macro F1: 0.8694. Fold SD describes variation, not a confidence interval. Per-class recalls and one-vs-rest AUROCs are recorded in `logistic_initial_report.json`; individual fold results are in `logistic_initial_cv_folds.csv`. Sample-level predictions stay in ignored processed storage.

This fixed L2 starting reference does not replace planned tuned L1/elastic-net, random forest, and XGBoost comparisons. No panel-size analysis or final test evaluation has been performed. The explanatory notebook is `notebooks/03_first_logistic_model.ipynb`; its displayed calculations and plotting code have been verified. The confusion matrix figure explicitly refers to training cross-validation.

## Initial tree-model comparison

```powershell
.\.venv\Scripts\python.exe src/evaluate_trees.py
.\.venv\Scripts\python.exe src/plot_initial_comparison.py
```

Both tree models use the same saved training folds, log1p, and fold-local variance filtering. Standardization is unnecessary for trees. Random forest starts with 300 trees, minimum leaf size 2, square-root candidate-feature sampling, and balanced class weights. XGBoost starts with 200 trees, depth 3, learning rate 0.05, histogram splits with 64 bins, and half the features considered per tree. XGBoost uses balanced sample weights computed only from the fitting fold's labels. Seed is 42; each model uses four worker threads.

XGBoost's initial 300-tree, default-bin pilot was stopped because of runtime after one completed fold. The five reported folds all use the revised 200-tree, 64-bin, half-feature configuration. The incomplete pilot is not included in the comparison. These are exploratory fixed-setting results, not a tuned model selection or an unbiased final performance estimate.

All source gene features remain candidate inputs after fold-local variance filtering. A random feature subset per tree is not a fixed gene panel: different trees may use different genes. No test labels or test expression rows are used. Random forest and XGBoost aggregate reports and fold tables are saved alongside the logistic report. Sample-level predictions remain local and excluded from Git.

`notebooks/04_initial_model_comparison.ipynb` explains the models and displays the aggregate comparison. `figures/initial_model_comparison.png` shows mean balanced accuracy with fold standard deviations. These initial runs do not include hyperparameter search or panel-size selection; subsequent tuning is documented below.

## Nested training-only tuning

```powershell
.\.venv\Scripts\python.exe src/tune_models.py --model logistic
.\.venv\Scripts\python.exe src/tune_models.py --model random_forest
.\.venv\Scripts\python.exe src/tune_models.py --model xgboost
.\.venv\Scripts\python.exe src/plot_initial_comparison.py --stage tuned
```

The predefined bounded grids and fixed settings are recorded in `tuning_protocol.json`. For each saved outer training fold, three inner stratified folds choose the best setting by mean balanced accuracy. The selected Pipeline is then refitted on the outer fitting samples and scored on the held-aside outer samples. All preprocessing and class-balancing weights are recomputed within each fitting set, including inner folds. XGBoost's reusable balanced estimator is defined in `src/estimators.py`.

Logistic regression now uses four L1 one-versus-rest classifiers and C in [0.01, 0.1, 1.0]. Random forest uses 200 trees and minimum leaf size in [1, 3]. XGBoost uses 100 trees, learning rate 0.05, 32 histogram bins, half the features per tree, and depth in [2, 3]. Seed is 42, with inner seeds 42 plus outer fold number. This is not exhaustive optimization. Because logistic formulation and tree compute budgets differ from initial pilots, changes in scores cannot be attributed solely to tuning.

Nested scores evaluate the tuning procedure: selected settings may differ across outer folds. A separate three-fold search on all 756 training samples selects final full-feature settings and fits a local model. Its inner search score is not presented as an independent performance estimate. Fitted models and patient-level predictions stay in ignored `data/processed/`; tracked aggregate reports use the `_tuned_report.json`, `_tuned_cv_folds.csv`, and `_tuning_search.json` suffixes.

The test set remains unused. The nested tuning explanation is in `notebooks/05_nested_model_tuning.ipynb`. Model comparisons use mean balanced accuracy and fold standard deviations, not confidence intervals or claims of statistical superiority.

Completed nested results: L1 logistic regression 0.9013, random forest 0.8947, XGBoost 0.9122 mean balanced accuracy. The final training-only searches select logistic C=0.1, random-forest minimum leaf size 3, and XGBoost depth 2. Final settings are for full-feature refits and do not yet establish a gene panel. All saved outer-fold scores were verified against sample-level predictions; fitted local models reload successfully.

Before panel-size results were generated, the user approved the 2 percentage-point tolerance recorded in `panel_protocol.json`. This criterion does not establish clinical equivalence or statistical noninferiority. Gene selection must still be fitted within training folds.

## Gene-panel curves and the PAM50 gene-list benchmark

```powershell
.\.venv\Scripts\python.exe src/evaluate_panels.py --model logistic
.\.venv\Scripts\python.exe src/evaluate_panels.py --model random_forest
.\.venv\Scripts\python.exe src/evaluate_panels.py --model xgboost
.\.venv\Scripts\python.exe src/prepare_pam50.py
.\.venv\Scripts\python.exe src/evaluate_panels.py --model logistic --pam50
.\.venv\Scripts\python.exe src/evaluate_panels.py --model random_forest --pam50
.\.venv\Scripts\python.exe src/evaluate_panels.py --model xgboost --pam50
.\.venv\Scripts\python.exe src/plot_panel_curves.py
```

Automatic panels contain 1, 2, 5, 10, 20, 50, 100, or 500 gene features. ANOVA F-score selection is inside the Pipeline after fold-local variance filtering and before logistic standardization. The score measures subtype mean differences relative to within-subtype variation. It ranks genes individually and can select redundant correlated genes; it does not optimize gene interactions. All three models use the same ranking method for comparability. Every inner and outer fitting set learns its own ranking. The same inner tuning grids, seeds, and outer folds as the full-gene references are used.

The user approved the separate PAM50 gene-list benchmark during this analysis. `prepare_pam50.py` extracts all 50 genes from the published genefu `pam50` object's `centroids.map`, pins the repository commit, and records the source checksum. All 50 genes map uniquely by Entrez ID, including older symbol aliases. `pam50_gene_mapping.csv` preserves the published and expression-source symbols. Source metadata and limitations are in `pam50_reference.json`.

The PAM50 benchmark restricts each model to those 50 fixed gene features, then uses identical nested training validation and tuning grids. It does not use cohort-based gene selection. It is **not the original PAM50 nearest-centroid algorithm or the Prosigna assay**. Because target labels are PAM50-derived, these scores measure agreement with that reference, not independent clinical diagnostic validity. Sources: [original PAM50 paper](https://pubmed.ncbi.nlm.nih.gov/19204204/), [genefu documentation](https://www.bioconductor.org/packages/release/bioc/manuals/genefu/man/genefu.pdf).

| Model | Full-gene reference | Smallest tested qualifying automatic panel | Panel score | Fixed PAM50 50-gene score |
| --- | ---: | ---: | ---: | ---: |
| L1 logistic | 90.1% | 100 | 88.7% | 92.4% |
| Random forest | 89.5% | 100 | 88.6% | 93.2% |
| XGBoost | 91.2% | 500 | 89.9% | 92.1% |

All figures are mean balanced accuracy from five outer training folds. Stars select within the automatic ANOVA panel family using the preapproved 0.02 absolute tolerance; fixed PAM50 benchmarks are reported separately and pass that full-gene tolerance for all models. The 50-gene PAM50 benchmarks outperform the qualifying automatic panels in these runs. This highlights a limitation of the automatic ranking method, not a claim of clinical superiority. Gene counts are the selected input panel sizes; L1 may set additional weights to zero.

The panels vary across fitting folds: 83 of the 100 genes are shared by all five 100-gene outer refits, and 396 of the 500 genes are shared by all five 500-gene refits. This describes stability, not a final fixed list. The mean curve and tolerance rule are used for training-stage selection; they are not independent final performance estimates. A single final model/panel choice must be frozen before test evaluation.

The headline figure is `figures/panel_size_curve.png`, with a vector copy at `figures/panel_size_curve.svg`. It uses a log gene-count axis, 0–1 accuracy axis, three model curves, descriptive ±1 fold-SD bands, model-colored dashed full-gene references, a dotted dummy baseline, and practical threshold stars. Lines connect evaluated sizes; intermediate sizes have not been measured. Stars are operational plateau markers rather than geometric knees. Fixed PAM50 gene-set benchmarks are squares. Reports include `panel_size_summary.csv`, `panel_selection_report.json`, and `panel_selection_stability.csv`; per-model panel reports, tuning records, and fold gene selections are also saved. Sample-level predictions remain local in ignored storage. No test evaluation has occurred.

See `notebooks/06_gene_panel_curves.ipynb` for the plain-language explanation. Validation checked every outer-fold score, selected gene count, label alignment, probability sum, per-class recall, and AUROC, plus fixed PAM50 membership and the 2-point selection rule.
