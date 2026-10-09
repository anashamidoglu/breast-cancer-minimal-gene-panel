# Understanding the Luminal mistakes

This review describes saved predictions from the frozen primary random forest on the 50 PAM50 genes. It does not refit the model, change its genes, or select a confidence threshold.

| Recorded subtype | Correct | Predicted as the other Luminal subtype | Predicted HER2-enriched |
| --- | ---: | ---: | ---: |
| Luminal A (100 patients) | 90 | 8 | 2 |
| Luminal B (39 patients) | 34 | 4 | 1 |

The model correctly classified 124 of 139 Luminal patients. Twelve of the primary model's 16 total errors were direct A/B swaps (75%). Although there were more A-to-B errors in count, their rates were 8/100 = 8.0% and 4/39 = 10.3% for B-to-A. Counts alone would obscure the smaller Luminal B group. These are descriptive estimates from small groups, not a significant difference claim.

## Why the distinction can be difficult

Luminal A and B share hormone-related biology. NCI describes both as hormone receptor-positive, with Luminal B often showing more cell division activity: [Luminal A](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/luminal-a-breast-cancer), [Luminal B](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/luminal-b-breast-cancer). This is context for why similar expression patterns may be difficult to separate; our error review does not establish the biological reason for individual mistakes. Our recorded labels are gene-expression subtypes, and are not interchangeable with clinical receptor-test categories.

## Were the errors close calls?

We calculate the difference between the model's largest and second-largest probability estimates, called the margin. A small margin means its two leading scores are close. For example, 0.46 versus 0.44 gives a margin of 0.02; 0.90 versus 0.08 gives 0.82. These estimates have not been shown to be calibrated, so they are not reliable percentages of diagnostic certainty.

For recorded Luminal A, the median margin was 0.874 for correct predictions and 0.212 for incorrect ones. For recorded Luminal B, it was 0.611 for correct predictions and 0.133 for incorrect ones. Errors generally had closer competing scores in this sample. Five of the twelve direct A/B swaps had margins below 0.05, but one had a margin of 0.621: some mistakes occurred despite a strong preference. The 0.05 value is only an illustrative description, not a chosen decision threshold.

The figure shows all 139 patients. On the right, positive values mean the model assigned a higher score to A than B; negative values mean B was higher. A third class can still win, so this difference does not by itself determine the final prediction. Red crosses mark every incorrect prediction, including three HER2-enriched predictions. Vertical offsets separate points for readability and have no biological meaning.

![Luminal error review](../figures/luminal_error_review.png)

These observations can guide discussion and future independent studies. They cannot justify tuning this model using the already-opened test patients. Aggregate results are public; the patient-level review stays in ignored local processed data.
