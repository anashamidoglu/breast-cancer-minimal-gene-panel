# Staged gene-panel experiment: development protocol

Status: proposed development design; no staged-model results generated.

## Research question

When is a small gene panel sufficient to reproduce a breast cancer molecular subtype label, and which tumours benefit from additional gene measurements?

Primary objective: reduce the average number of distinct genes measured per patient while keeping balanced accuracy within 2 absolute percentage points of a fixed larger-panel comparator. This proposed extension carries forward the original project's practical tolerance; it is not a statistical equivalence claim.

## First experiment

Use two nested automatic panels: a 20-gene first-stage panel and a 100-gene second-stage panel. Both derive from a single ANOVA ranking fitted on the appropriate training subset, so every first-stage gene is included in the second-stage panel. The model family is class-balanced random forest, chosen before staged results because it is the original primary model family. Initial classifier settings are 200 trees, minimum leaf size 3, seed 42, and two CPU threads, held fixed for this first experiment. The 50-gene PAM50 random forest is a separate reference.

For each patient, the first-stage model predicts from 20 measurements. The routing score is the difference between its highest and second-highest predicted class probabilities. A margin below a training-selected threshold triggers the 100-gene model; otherwise the first-stage prediction is retained. Every patient receives a final prediction. This is escalation, not abstention, and the larger panel is not presumed to correct every error.

Because panels are nested, an escalated patient needs 80 additional genes, not 100 additional genes. Average gene count equals 20 + 80 × escalation fraction. For example, escalating 30% produces an average of 44 genes per patient. Reference normalization/housekeeping genes and assay overhead are outside this computational count and must be reported separately if a laboratory design is later proposed.

## Comparisons

- Fixed 20-gene model for everyone.
- Fixed 100-gene model for everyone: primary accuracy reference.
- Staged 20-to-100-gene model.
- Random escalation at the same escalation fraction, averaged over fixed random seeds. This checks whether uncertainty routing adds value beyond simply measuring more genes for some patients.
- Fixed 50-gene PAM50 model: established compact gene-list reference, distinct from the original PAM50 classifier.

The staged strategy must also be compared with the 50-gene reference for both performance and average gene count. Reporting only improvement over 100 genes could conceal a weaker tradeoff than using PAM50 for everyone.

## Training and routing rule

Development is restricted to the original 756 TCGA training patients. The 189 previously evaluated TCGA test patients will not serve as a fresh test or select the new design.

Use the saved five outer training folds. Inside each outer fitting set, generate three-fold out-of-fold predictions for both panels, repeating preprocessing and supervised gene ranking within each inner fit. Each paired inner fit uses the same ranking for the two nested panels. Select the routing threshold using those inner out-of-fold predictions only, from the fixed grid 0.00, 0.05, ..., 1.00. Route when margin < threshold. Include an explicit route-all candidate in addition to the grid so every-patient escalation is available even when margin equals 1.

Among thresholds whose combined inner balanced accuracy is at least the corresponding fixed 100-gene inner score minus 0.02, choose the threshold with the lowest escalation fraction; break ties using higher balanced accuracy, then the smaller threshold. Refit both panels on the entire outer fitting set and evaluate the chosen routing rule once on its untouched outer fold. Repeat for all five outer folds. The outer validation results assess the threshold-selection procedure, not a threshold optimized using outer outcomes.

For the final frozen candidate, repeat the same inner out-of-fold threshold selection within all 756 training patients, then refit both nested panels on all training patients. Freeze the genes, models, routing threshold, gene-count accounting, normalization, and external evaluation plan before evaluating an independent cohort.

## Outcomes and figures

Primary outcomes: balanced accuracy, average distinct genes measured, escalation fraction, and absolute performance difference from the fixed 100-gene comparator. Report macro F1, subtype recall, and escalation by recorded subtype to reveal uneven measurement burden.

Among escalated cases, report how many initially incorrect predictions become correct, how many correct predictions become incorrect, and how many stay incorrect. Compare uncertainty-based escalation with random escalation at equal burden. Keep these outcome-based assessments out of the routing rule.

Main figures: performance versus average genes measured; escalation fraction by subtype; and a transition table showing errors corrected, introduced, or unchanged after escalation. Development tradeoff curves use training validation. The independent evaluation applies one frozen rule rather than choosing a new point from its curve.

## Independent evaluation: feasibility gate

METABRIC is an independent candidate available through cBioPortal. Its expression profile uses Illumina microarrays, unlike TCGA RNA sequencing. Before design freeze, audit expression metadata, Entrez mappings, label definitions, available genes, sample identifiers, and duplicate patients. Do not inspect candidate-model performance while selecting preprocessing or the external evaluation set.

A TCGA-fitted log-RSEM model must not be applied directly to microarray values. A valid frozen cross-platform representation or an alternative compatible RNA-seq cohort must be established before external evaluation. Any learned adaptation requires a designated development subset; a separate independent evaluation subset must remain unused. Cohort-level transformations using the entire external dataset must be disclosed as transductive and cannot silently be presented as single-sample deployment.

METABRIC's PAM50-plus-claudin-low annotation must be checked before mapping labels. Claudin-low, Normal-like, and missing labels must not be silently reassigned to the four target classes. Missing selected genes must not be replaced with zero or dropped after seeing results. External eligibility and handling rules must be frozen in advance.

## Interpretation and limits

A useful result could be selective measurement savings with preserved accuracy, or evidence that adding genes rarely resolves the cases selected as uncertain. Model uncertainty must not be equated with biological ambiguity without independent evidence. The strongest scientific extension would test whether reproducibly difficult tumours relate to proliferation, tumour composition, or clinical-marker disagreement using prespecified analyses.

Probability margins are uncalibrated model scores. The staged approach is a computational simulation using already-measured expression, not a demonstrated staged laboratory assay. Measurement count is not validated monetary cost, turnaround time, or clinical utility. Labels remain PAM50-derived. The first design uses a chosen 20/100 pair and fixed model settings; further variants would be development analyses requiring nested evaluation and an independent final test. Novelty is not established by this protocol and requires a focused literature review.

## Dataset references

[cBioPortal METABRIC study](https://www.cbioportal.org/study?id=brca_metabric) and [official DataHub study files](https://github.com/cBioPortal/datahub/tree/master/public/brca_metabric).
