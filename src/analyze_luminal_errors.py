"""Describe saved primary-model Luminal errors; never fit or tune a model."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tune_models import ROOT


def main():
    plan = json.loads((ROOT / 'data/final_model_freeze.json').read_text())
    report = json.loads((ROOT / 'data/final_test_report.json').read_text())
    assert plan['primary_case'] == report['primary_case'] == 'random_forest_pam50'
    predictions = pd.read_csv(ROOT / 'data/processed/final_test_random_forest_pam50_predictions.csv')
    primary = next(c for c in report['cases'] if c['case_id'] == plan['primary_case'])
    table = pd.crosstab(predictions.subtype, predictions.predicted_subtype).reindex(index=report['classes'], columns=report['classes'], fill_value=0)
    assert table.to_numpy().tolist() == primary['confusion_matrix']
    luminal = predictions[predictions.subtype.isin(['Luminal A', 'Luminal B'])].copy()
    probabilities = luminal[[f'probability_{label}' for label in report['classes']]].to_numpy()
    ordered = np.sort(probabilities, axis=1)
    luminal['top_two_margin'] = ordered[:, -1] - ordered[:, -2]
    luminal['correct'] = luminal.subtype == luminal.predicted_subtype
    luminal['a_minus_b'] = luminal['probability_Luminal A'] - luminal['probability_Luminal B']
    luminal.to_csv(ROOT / 'data/processed/luminal_error_review.csv', index=False)
    rows = []
    for subtype in ['Luminal A', 'Luminal B']:
        subset = luminal[luminal.subtype == subtype]
        counts = subset.predicted_subtype.value_counts()
        rows.append(dict(recorded_subtype=subtype, total=len(subset), prediction_counts={label: int(counts.get(label, 0)) for label in report['classes']}, median_top_two_margin_correct=float(subset.loc[subset.correct, 'top_two_margin'].median()), median_top_two_margin_incorrect=float(subset.loc[~subset.correct, 'top_two_margin'].median())))
    swaps = luminal[(~luminal.correct) & luminal.predicted_subtype.isin(['Luminal A', 'Luminal B'])]
    result = dict(primary_case=plan['primary_case'], evaluation='Descriptive post-test review of saved predictions; no refitting or selection.', luminal_patients=len(luminal), correct=int(luminal.correct.sum()), direct_a_b_swaps=len(swaps), all_primary_errors=int((predictions.subtype != predictions.predicted_subtype).sum()), per_subtype=rows, direct_swap_top_two_margins=sorted(swaps.top_two_margin.tolist()), note='Random forest probability estimates are not demonstrated to be calibrated. Margins describe model scores, not certainty about the true diagnosis. No threshold selected from these results.')
    (ROOT / 'data/luminal_error_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    labels = ['Basal-like', 'HER2-enriched', 'Luminal A', 'Luminal B']
    matrix = table.loc[['Luminal A', 'Luminal B'], labels].to_numpy()
    axes[0].imshow(matrix, cmap='Blues', aspect='auto')
    for i in range(2):
        for j in range(4):
            axes[0].text(j, i, str(matrix[i,j]), ha='center', va='center', color='white' if matrix[i,j] > matrix.max()/2 else 'black')
    axes[0].set(xticks=range(4), xticklabels=labels, yticks=range(2), yticklabels=['Luminal A (100)', 'Luminal B (39)'], xlabel='Predicted subtype', ylabel='Recorded subtype', title='Where the 139 Luminal patients went')
    plt.setp(axes[0].get_xticklabels(), rotation=25, ha='right')
    for i, subtype in enumerate(['Luminal A', 'Luminal B']):
        subset = luminal[luminal.subtype == subtype]
        for correct, marker, color in [(True, 'o', '#4477AA'), (False, 'x', '#CC6677')]:
            values = subset[subset.correct == correct].a_minus_b.to_numpy()
            offsets = np.linspace(-.16, .16, len(values)) if len(values) else []
            axes[1].scatter(values, i + np.array(offsets), marker=marker, color=color, alpha=.7, label=('Correct' if correct else 'Incorrect') if i == 0 else None)
    axes[1].axvline(0, color='black', linewidth=.8)
    axes[1].set(xlim=(-1.05,1.05), yticks=[0,1], yticklabels=['Recorded Luminal A','Recorded Luminal B'], xlabel='Model score difference: P(Luminal A) − P(Luminal B)', title='Right favors A; left favors B')
    axes[1].legend()
    fig.suptitle('Frozen primary random forest on PAM50 genes: held-out Luminal review')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/luminal_error_review.png', dpi=180)
    fig.savefig(ROOT / 'figures/luminal_error_review.svg')
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
