# Automatic staged panels: first experiment

Start with 20 automatically selected genes, then expand to 100 when the model's leading subtype scores are close.

| Strategy | Average genes | Balanced accuracy |
| --- | ---: | ---: |
| Fixed 20 genes | 20.0 | 81.7% |
| Random escalation | 43.5 | 83.9% |
| Uncertainty escalation | 43.5 | 86.3% |
| Fixed 100 genes | 100.0 | 88.5% |
| Fixed PAM50 genes | 50.0 | 93.1% |

Uncertainty routing worked better than random routing, but its 2.24-point loss missed the planned 2-point target. Additional genes corrected 69 predictions and introduced 21 errors. Luminal B cases were escalated most frequently.

## Methods and limits

These are nested validation results from the 756 TCGA development patients. Gene ranking and routing thresholds were selected within fitting subsets. No fresh test was used. PAM50 offered a stronger accuracy-versus-gene-count tradeoff, motivating [the next experiment](pam50_staged_findings.md). [Protocol](staged_panel_protocol.md).
