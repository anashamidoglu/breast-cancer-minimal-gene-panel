# Staging within PAM50

The second experiment selected 20 first-stage genes from PAM50 and expanded uncertain cases to all 50.

| Strategy | Average genes | Balanced accuracy |
| --- | ---: | ---: |
| Fixed 20-gene subset | 20.0 | 89.7% |
| Random escalation | 22.2 | 89.8% |
| Uncertainty escalation | 22.2 | 90.8% |
| Fixed 50 genes | 50.0 | 93.1% |

Better panel content improved accuracy while using fewer measurements than the first experiment. However, the 2.33-point loss still missed the planned target. Escalation corrected 12 predictions and introduced three errors.

## Methods and limits

This was a training-only follow-up motivated by the first experiment. The model family and validation design remained fixed. The result suggested using a more conservative routing rule, which was recorded before development in a separate cohort. See [the independent replication](staged_research_conclusion.md). [Protocol](pam50_staged_protocol.md).
