# SCAN-B staged strategy: independent-cohort replication protocol

## Objective and evaluation

Evaluate whether uncertainty-based escalation from a training-selected 20-gene subset of PAM50 to all 50 PAM50 genes can preserve balanced accuracy within 2 absolute percentage points while reducing average distinct genes measured. This is independent-cohort replication after within-cohort fitting, not direct transfer of the TCGA models.

## Frozen split

Official GSE96058 SOFT metadata defines 3,273 primary cases plus 136 technical replicates. Retain primary titles only; verify each replicate title maps to an original by removing `repl`. BioSample accessions differ between primary and replicate records, so canonical original titles are the grouping key. Exclude 221 Normal-like primary cases. The remaining 3,052 four-class cases are split into 2,136 development and 916 final-holdout cases, stratified with seed 20261010. Saved membership hashes must match before every modeling stage. Final labels may be used to establish stratification and counts, not to tune or select a model.

## Expression feasibility

The official processed data is log2(FPKM + 0.1), with transcripts collapsed by gene symbol. Do not apply another log transform. Match the 50 published PAM50 genes using mapped current/source symbols and predeclared aliases. Every required gene must have one unambiguous expression row; stop if mapping is ambiguous or incomplete. No gene is dropped or imputed based on results. Exclude technical replicates before modeling. Retain only required gene rows for this experiment.

## Candidate fixed before independent predictions

Use random forest with 200 trees, balanced class weights, minimum leaf size 3, seed 42, and two threads. ANOVA ranks the 50 genes within each fitting set; select the top 20. The larger model uses all 50. Source-transformed expression is used directly. Gene order is fixed consistently for final fitting and prediction. No model family, panel size, or hyperparameter grid is selected from SCAN-B holdout results.

A more conservative routing rule is specified because the two TCGA development experiments narrowly missed the outer 2-point target: select the smallest escalation fraction meeting a **1-point inner-development margin** against the full 50-gene model, while the **final evaluation target remains 2 points**. This is a new hypothesis recorded before SCAN-B performance inspection. Use three stratified inner folds for out-of-fold predictions, margin = largest minus second-largest class score, thresholds 0.00 to 1.00 in increments of 0.05 plus route-all, and the original tie-breaking rule. All gene selection is refitted within each inner fit. Five outer development folds assess this procedure. No alternate threshold is selected from those outer results.

Then repeat the same three-fold threshold procedure on all 2,136 development cases, fit both final models there, and freeze the selected genes, threshold, classifier artifacts, and all candidate comparisons before evaluating the 916 holdout cases once. The complete design is applied even if outer development results fail the target. No retries to obtain a preferred outcome.

## Comparisons and outcomes

Compare fixed 20 genes, fixed 50 genes, uncertainty routing, and random escalation at the exact same holdout escalation count (100 fixed draws). Report balanced accuracy, macro F1, subtype recall, average genes = 20 + 30 × escalation fraction, and errors corrected/introduced/remaining. Primary descriptive success requires staged balanced accuracy >= fixed-50 score minus 0.02 and average genes < 50. Report paired stratified-bootstrap uncertainty for the accuracy difference and state whether its interval supports the practical retention claim. A point-estimate pass alone is not statistical equivalence.

## Limits

This design replicates the strategy in a separate cohort; it does not establish model transport across processing pipelines, clinical validity, a laboratory cost benefit, or novelty. All prior unsuccessful designs remain part of the research record. The conservative rule was motivated by TCGA results, not selected using the fresh SCAN-B holdout.
