# What the selected genes tell us

The fixed automatic 100-gene panels for logistic regression and random forest are identical: both use the same ANOVA ranking on all training patients. Ten genes overlap with the 50-gene PAM50 list: CCNE1, CDC20, ESR1, FOXA1, FOXC1, GPR160, MAPT, MLPH, NAT1, and PGR. That is 10% of our panel and 20% of PAM50; it is not an enrichment test.

Of these 100 final genes, 83 were also selected in every outer training fold. Of the 500 final XGBoost genes, 396 were selected in every fold, and 31 overlap with PAM50. This shows repeatability within our training cohort, not stability across other studies. Final panels contain additional genes because they were selected again using all 756 training patients.

ESR1 encodes an estrogen receptor, and PGR encodes a progesterone receptor. Their inclusion shows that the automatic panel contains genes involved in hormone signaling. These biological descriptions come from [NCBI ESR1](https://www.ncbi.nlm.nih.gov/gene/2099) and [NCBI PGR](https://www.ncbi.nlm.nih.gov/gene/5241/); the overlap and model weights are our own computed results. Gene-expression levels are not the same measurement as clinical receptor staining.

## Reading the coefficient plot

We inspect the supporting logistic model because its weights are easier to explain than tree splits. This does not change the frozen primary random forest model. Each subtype has a separate binary classifier against the remaining subtypes. A positive weight raises that classifier's score when the gene's expression increases, holding other genes fixed. A negative weight lowers it. The values apply to log-expression standardized using training means and standard deviations; they are not percentage changes in probability. Class probabilities are normalized across the four binary classifiers.

FOXC1 has the largest positive weight for Basal-like in the automatic logistic model. TPX2 has the largest positive weight for Luminal B and a large negative weight for Luminal A. These describe how this particular predictor separates labels, not causal biological effects.

FOXC1 also has a positive Luminal A weight, despite its positive Basal-like weight. Correlated genes, conditional adjustment, and separate one-versus-rest classifiers can produce such signs. Do not read a weight as saying a gene is uniquely expressed in that subtype, or that it causes that subtype. Large weights are not a complete ranking of biological importance; zero L1 weights can reflect redundant information.

![Frozen logistic gene weights](../figures/logistic_gene_weights.png)

The plot shows the eight largest absolute weights within each subtype. Panel sizes and model settings remain fixed. This analysis reads frozen training artifacts and gene selections; it neither fits models nor uses test patients. It does not perform pathway enrichment or discover new biomarkers. Source symbols are retained as in the TCGA annotation and can include historical aliases; Entrez feature IDs are the matching keys.
