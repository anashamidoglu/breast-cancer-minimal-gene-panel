# PAM50-restricted staged development protocol

This follow-up is motivated by the completed automatic 20-to-100 experiment. Before generating follow-up results, the design is fixed as follows: select 20 genes by ANOVA from the published PAM50 list within each fitting fold, escalating to all 50 PAM50 genes. The larger panel contains every first-stage gene. Average distinct gene count is 20 + 30 × escalation fraction.

The model settings, original 756 training patients, saved five outer folds, three inner folds, margin routing score, threshold grid, tie-breaking rule, and 2-point tolerance remain as specified in the original staged protocol. The previous TCGA test set is not used. The same fixed 50-gene comparator is evaluated with and without the ANOVA ordering; random forest feature sampling can depend on column ordering, so numerical equality of predictions is not assumed. Both configurations will be reported rather than silently selecting the better one.

This is a hypothesis-driven follow-up development analysis, not an independent confirmation. No further panel sizes or rules will be searched in this experiment. External feasibility and validation remain required before the desired staged-measurement conclusion can be established.
