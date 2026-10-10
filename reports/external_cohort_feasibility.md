# Independent cohort feasibility: initial audit

SCAN-B GSE96058 is a promising independent RNA-seq candidate. The official GEO record reports 3,273 cases plus 136 technical replicates (3,409 profiles), PAM50 subtyping, and a 564.3 MB compressed transformed gene-expression file. Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058

This improves platform comparability relative to METABRIC microarrays, but does not resolve normalization compatibility. Original TCGA RSEM and SCAN-B FPKM values cannot be assumed interchangeable. Before external prediction, verify the exact transformation and gene mappings, remove technical replicate duplication with a predefined rule, inspect label metadata, and freeze preprocessing. No sample-level external data or model performance have been inspected yet.

A separate development/calibration portion may be necessary if adaptation is learned. Its membership must be frozen before adaptation, with a patient-disjoint final portion left unopened. If the design is retrained on SCAN-B, that tests independent-cohort replication, not direct transport of the TCGA-fitted predictor. These questions must be stated separately.
