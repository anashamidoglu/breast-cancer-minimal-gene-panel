# Data and reproduction

**TCGA-BRCA:** 945 retained patients; 756 development and 189 test. Source RNASeqV2 RSEM measurements from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018).

**SCAN-B GSE96058:** 3,052 four-class cases after excluding technical replicates and Normal-like labels; 2,136 development and 916 holdout. Source log2(FPKM + 0.1) measurements from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058).

Raw data, patient-level outputs, and model binaries stay local. Download and split reports track checksums. Gene lists, aggregate metrics, and frozen model specifications are public.

## Reproduce SCAN-B

Use Python 3.12 and install `requirements.txt` first. Create a separate workspace that contains source code, protocols, source-version pins, and the published PAM50 gene list, without the completed model/results files:

```bash
python src/create_reproduction_workspace.py
cd reproductions/fresh-run
python src/download_scanb.py
python src/prepare_scanb_metadata.py
python src/prepare_scanb_expression.py
python src/develop_scanb_staged.py
python src/evaluate_scanb_holdout.py
```

The downloader fetches both official expression and SOFT metadata, checking them against the recorded source checksums. The expression download is approximately 564 MB. Model development takes substantially longer than viewing the saved results. Use the same activated Python environment after changing directories.

The workspace is ignored by Git. An existing destination is never overwritten; choose another with `--destination reproductions/another-run`. Membership hashes in the copied split report are checked during metadata preparation. Final model serialization/checksums can vary by platform or dependency versions even when seeds and metrics match.

## Other analyses

TCGA preparation runs in this order: `inspect_dataset.py`, `prepare_labels.py`, `download_expression.py`, `match_expression.py`, `create_split.py`, `evaluate_baseline.py`, and `prepare_pam50.py`. Then run `tune_models.py --model MODEL` and `evaluate_panels.py --model MODEL` for logistic, random_forest, and xgboost; add `--pam50` for the fixed-gene benchmark. `plot_panel_curves.py` precedes `freeze_final_models.py` and `evaluate_final_test.py`.

`evaluate_staged_panels.py` runs the automatic staged experiment; `--pam50-staged` runs the PAM50 follow-up. Exact design rules are in the [protocol index](../reports/README.md). SCAN-B source expression is already log-transformed; TCGA and SCAN-B fitted models are not interchangeable.
