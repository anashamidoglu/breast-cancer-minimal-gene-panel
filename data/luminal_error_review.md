# TCGA Luminal errors

The primary 50-gene random forest correctly classified 90/100 Luminal A and 34/39 Luminal B test patients. Eight A cases were predicted as B and four B cases as A: 12 of the model's 16 errors were direct swaps. Three other Luminal cases were predicted as HER2-enriched.

Incorrect predictions generally had smaller score gaps, though some errors had a strong model preference. This descriptive review did not change any model settings. Detailed counts are in `luminal_error_summary.json`.
