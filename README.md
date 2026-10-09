# How few genes are enough?

A beginner-friendly investigation of how small a gene-expression panel can reproduce breast cancer molecular subtype labels while retaining performance close to a full-transcriptome classifier.

## Current stage

Python setup is verified. See [SETUP.md](SETUP.md) and start with `notebooks/00_getting_started.ipynb`. Study metadata, subtype labels, and expression availability have been checked; see [data/README.md](data/README.md). With user approval, Normal-like and missing labels were excluded. All 945 retained patients have expression data. Identifier cleanup leaves 20,481 unique gene features. Explore the counts in `notebooks/01_inspect_subtype_labels.ipynb`. No train/test split or model fitting has been performed.

## Planned approach

- Use TCGA-BRCA PanCancer Atlas expression data and verify the available subtype labels.
- Compare gene-panel sizes using balanced accuracy: the average of the fraction correctly classified within each subtype.
- Choose the panel using training data only, then evaluate it on a separate test set.
- Compare a majority-class baseline, logistic regression, random forest, and XGBoost.

The proposed panel-selection tolerance is 2 percentage points below the full-gene reference in cross-validation. This is a practical project choice, not a biological standard, and will be settled before analysis.

## Folders

- `data/`: download instructions; original and processed data stay local.
- `notebooks/`: step-by-step analyses and explanations.
- `src/`: reusable Python code.
- `figures/`: plots produced by the analysis.

## Limitations

- TCGA subtype labels were derived from PAM50, so strong performance is expected and is not a discovery.
- Results will come from one cohort, with no external validation.
- This is an educational project, not a clinical tool.
- A small computational feature set is not a validated laboratory assay.
- Source expression values were batch-normalized before download; this upstream processing cannot be refitted within our training folds.
- Gene rows with ambiguous repeated Entrez IDs were excluded, potentially removing useful measurements.

Target completion: October 31, 2026.
