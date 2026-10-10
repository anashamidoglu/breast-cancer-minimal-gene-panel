# How few genes are enough?

A beginner-friendly investigation of how small a gene-expression panel can reproduce breast cancer molecular subtype labels while retaining performance close to a full-transcriptome classifier.

## Current stage

Training-only analysis finds that **100 automatically ranked genes retain performance within 2 percentage points of the full-gene reference for logistic regression and random forest; XGBoost requires 500 among the sizes tested**. The fixed published PAM50 gene set, with 50 genes, performs better than these qualifying automatic panels in all three models. This is agreement with PAM50-derived subtype labels, not independent clinical validation.

All 945 retained patients have expression data; identifier cleanup leaves 20,481 unique candidate gene features. A fixed stratified split reserves 756 samples for training and 189 for final testing. Panel evaluation uses five saved outer training folds with three inner tuning folds; gene selection and preprocessing are fitted within every fitting fold. **Final held-out testing is complete.** Random forest on the fixed PAM50 gene list was selected as the primary model before opening the test set. See [SETUP.md](SETUP.md), [data/README.md](data/README.md), and `notebooks/06_gene_panel_curves.ipynb`.

![Gene-panel performance and fixed PAM50 gene-list benchmarks](figures/panel_size_curve.png)

Nested training cross-validation compares automatic panels and full-gene references; squares show our models using the fixed PAM50 gene list, and stars mark the smallest tested automatic panels within the preapproved 2-point tolerance. Shading is fold standard deviation, not a confidence interval.

| Model | Full-gene balanced accuracy | Qualifying automatic panel | Panel balanced accuracy | Fixed PAM50 gene-list balanced accuracy |
| --- | ---: | ---: | ---: | ---: |
| L1 logistic | 90.1% | 100 genes | 88.7% | 92.4% |
| Random forest | 89.5% | 100 genes | 88.6% | 93.2% |
| XGBoost | 91.2% | 500 genes | 89.9% | 92.1% |

Scores are means across outer training folds. PAM50 gene-list models use the fixed published genes; they do not implement the original PAM50 classifier. Automatically selected genes can differ across folds; this evaluates the selection procedure rather than one final gene list. The smallest tested qualifying size is not necessarily the smallest possible panel.

![Bounded tuned model comparison on nested training cross-validation](figures/tuned_model_comparison.png)

Settings are selected within three inner folds and evaluated on five saved outer training folds; error bars show fold SD, not confidence intervals. These small searches do not establish statistical superiority. Logistic formulation and tree compute budgets differ from the initial models below.

![Initial model comparison on training cross-validation](figures/initial_model_comparison.png)

Initial models are compared on five identical training folds; error bars show fold standard deviation, not confidence intervals.


## Final held-out test

Before opening the 189 test patients, all predictors, settings, and final gene lists were frozen. The primary model was random forest on the fixed 50 PAM50 genes, selected for its highest nested training balanced accuracy and compact panel. It scored **92.7% balanced accuracy**, **91.0% macro F1**, and correctly classified **173/189 patients** (91.5% ordinary accuracy).

| Subtype | Correct / total | Recall |
| --- | ---: | ---: |
| Basal-like | 34 / 34 | 100.0% |
| HER2-enriched | 15 / 16 | 93.8% |
| Luminal A | 90 / 100 | 90.0% |
| Luminal B | 34 / 39 | 87.2% |

| Model | All 20,481 input genes | Automatic panel | Fixed PAM50 genes |
| --- | ---: | ---: | ---: |
| L1 logistic | 93.7% | 92.4% (100 genes) | 92.3% (50 genes) |
| Random forest | 88.7% | 90.6% (100 genes) | 92.7% (50 genes) |
| XGBoost | 94.4% | 95.4% (500 genes) | 92.1% (50 genes) |

All entries are held-out balanced accuracy; dummy baseline is 25.0%. Full-gene pipelines apply their fitted variance filters. The strongest test comparison was the automatic 500-gene XGBoost model; test comparisons do not change the prespecified primary model. The automatic panels retain performance within the practical 2-point margin on this holdout, but this is not statistical equivalence. Panel-size choices were made from training CV, and the curve remains a training-CV figure.

![Primary model held-out confusion matrix](figures/final_test_confusion.png)

These results measure agreement with PAM50-derived labels in one TCGA cohort. The test has only 16 HER2-enriched patients, so subtype estimates are imprecise; there is no external validation. See `data/final_model_freeze.json`, `data/final_panel_genes.csv`, `data/final_test_report.json`, and `notebooks/07_final_test.ipynb`. Patient-level predictions and model binaries stay local. Historical training reports retain `test_evaluated: false` because those reports describe the earlier training-only stage.

## Planned approach

- Use TCGA-BRCA PanCancer Atlas expression data and verify the available subtype labels.
- Compare gene-panel sizes using balanced accuracy: the average of the fraction correctly classified within each subtype.
- Choose the panel using training data only, then evaluate it on a separate test set.
- Compare a majority-class baseline, logistic regression, random forest, and XGBoost.

The user approved the panel-selection tolerance on October 9, 2026, before panel-size results were generated: select the smallest evaluated panel whose mean training cross-validation balanced accuracy is within 2 absolute percentage points of the corresponding full-gene reference. This is a practical project criterion, not proof of clinical equivalence or statistical noninferiority. See `data/panel_protocol.json`.

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
- Automatic ANOVA ranking evaluates genes individually and can select redundant measurements; its results do not rule out better small panels from other selection methods.
- The PAM50 gene-list benchmark trains our own models against PAM50-derived labels; it does not establish superiority to the PAM50 classifier or clinical assay.

Target completion: October 31, 2026.

## Gene interpretation

The automatic 100-gene panel overlaps PAM50 by 10 genes; 83 final genes were selected in all five outer training folds. The 500-gene panel overlaps by 31 genes, with 396 final genes selected in all five folds. Logistic and random forest share the same automatic 100-gene list because selection uses the same ANOVA ranking. See [gene interpretation](data/gene_interpretation.md) and `notebooks/08_gene_interpretation.ipynb` for frozen logistic weights and their limitations. This descriptive analysis uses training artifacts only and does not establish biological causes or cross-cohort reproducibility.

## Luminal error review

The frozen primary model made 12 direct Luminal A/B swaps, accounting for 12 of its 16 total held-out errors. It correctly classified 90/100 Luminal A and 34/39 Luminal B patients. Incorrect Luminal predictions generally had smaller score margins, though some were confident errors. See [the descriptive review](data/luminal_error_review.md) and `notebooks/09_luminal_error_review.ipynb`. No model or decision threshold was adjusted using this review.

## Compiled research findings

Read the [research report](reports/research_findings.md) for the research question, study design, small-panel results, held-out comparison, gene interpretation, Luminal errors, and limitations.

## Staged-panel extension

The first training-only 20-to-100-gene experiment averaged 43.5 genes and achieved 86.3% balanced accuracy, versus 88.5% for 100 genes and 93.1% for fixed PAM50 genes. The staged loss of 2.24 points narrowly missed the planned 2-point target. See the [protocol](reports/staged_panel_protocol.md) and [development findings](reports/staged_training_findings.md). No new independent evaluation has been completed.

## Independent staged-strategy replication

On 916 fresh SCAN-B holdout patients, a frozen 20-to-50 PAM50 strategy achieved 91.9% balanced accuracy using 26.4 genes on average, versus 93.2% for all 50 genes. The 1.28-point observed loss met the practical target, but its paired-bootstrap interval crossed the 2-point boundary. See the [research conclusion](reports/staged_research_conclusion.md) for all experiments and limitations. This is within-SCAN-B training and holdout evaluation, not direct TCGA model transfer.
