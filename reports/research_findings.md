# The original gene-panel comparison

The starting question was whether a smaller selection of genes could reproduce breast cancer subtype labels about as well as thousands of genes.

TCGA provided 945 patients with four retained subtype labels. A fixed split reserved 756 for development and 189 for testing. Logistic regression, random forest, and XGBoost were compared using all 20,481 candidate genes, automatically selected subsets, and the fixed PAM50 gene list.

| Model | Full-gene test score | Automatic panel | PAM50 gene-list test score |
| --- | ---: | ---: | ---: |
| Logistic regression | 93.7% | 92.4% / 100 genes | 92.3% |
| Random forest | 88.7% | 90.6% / 100 genes | 92.7% |
| XGBoost | 94.4% | 95.4% / 500 genes | 92.1% |

Scores are balanced accuracy. Automatic sizes were selected from training validation before opening the test set. The primary model, random forest on PAM50 genes, was also chosen in advance.

![Gene count and performance](../figures/panel_size_curve.png)

The automatic panels retained strong performance with far fewer genes. PAM50 was a particularly effective compact reference during development. Most primary-model errors were between Luminal A and B: 12 of 16 errors were direct swaps.

This led to a new question: instead of choosing one panel for everyone, could additional genes be reserved for uncertain cases? See [the staged findings](staged_research_conclusion.md).

## Methods and limits

ANOVA selected genes inside each fitting fold. Five outer validation folds assessed models and three inner folds selected settings. Source log transformation, filtering, and scaling were applied within fitting subsets as appropriate. Labels were PAM50-derived, and the PAM50 comparison used the gene list rather than the original classifier. The single-cohort test does not establish clinical validity. [Detailed metrics](../data/final_test_report.json).
