# Data and reproduction

**TCGA-BRCA:** 945 retained patients; 756 development and 189 test. Source RNASeqV2 RSEM measurements from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018).

**SCAN-B GSE96058:** 3,052 four-class cases after excluding technical replicates and Normal-like labels; 2,136 development and 916 holdout. Source log2(FPKM + 0.1) measurements from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058).

Raw data, patient-level outputs, and model binaries stay local. Download and split reports track checksums. Gene lists, aggregate metrics, and frozen model specifications are public.

## Reproduce SCAN-B

Use Python 3.12 and `requirements.txt`. Download the official `GSE96058_family.soft.gz` from the GEO `soft/` directory into `data/raw/scanb/` first.

```bash
python src/download_scanb.py
python src/prepare_scanb_metadata.py
python src/prepare_scanb_expression.py
python src/develop_scanb_staged.py
python src/evaluate_scanb_holdout.py
```

Scripts refuse to overwrite completed stages. For a fresh reproduction, use a separate checkout and archive the published result/freeze files before execution.

## Other analyses

TCGA preparation runs in this order: `inspect_dataset.py`, `prepare_labels.py`, `download_expression.py`, `match_expression.py`, `create_split.py`, `evaluate_baseline.py`, and `prepare_pam50.py`. Then run `tune_models.py --model MODEL` and `evaluate_panels.py --model MODEL` for logistic, random_forest, and xgboost; add `--pam50` for the fixed-gene benchmark. `plot_panel_curves.py` precedes `freeze_final_models.py` and `evaluate_final_test.py`.

`evaluate_staged_panels.py` runs the automatic staged experiment; `--pam50-staged` runs the PAM50 follow-up. Exact design rules are in the [protocol index](../reports/README.md). SCAN-B source expression is already log-transformed; TCGA and SCAN-B fitted models are not interchangeable.
