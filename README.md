# How few genes are enough?

A beginner-friendly investigation of how small a gene-expression panel can reproduce breast cancer molecular subtype labels while retaining performance close to a full-transcriptome classifier.

## Current stage

Training-only analysis finds that **100 automatically ranked genes retain performance within 2 percentage points of the full-gene reference for logistic regression and random forest; XGBoost requires 500 among the sizes tested**. The fixed published PAM50 gene set, with 50 genes, performs better than these qualifying automatic panels in all three models. This is agreement with PAM50-derived subtype labels, not independent clinical validation.

All 945 retained patients have expression data; identifier cleanup leaves 20,481 unique candidate gene features. A fixed stratified split reserves 756 samples for training and 189 for final testing. Panel evaluation uses five saved outer training folds with three inner tuning folds; gene selection and preprocessing are fitted within every fitting fold. **The final test set remains unused.** A fixed final panel has not yet been chosen. See [SETUP.md](SETUP.md), [data/README.md](data/README.md), and `notebooks/06_gene_panel_curves.ipynb`.

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
