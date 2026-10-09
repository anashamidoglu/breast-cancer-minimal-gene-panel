"""Describe frozen panels and coefficients without fitting or using test data."""
import json
import sys

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tune_models import ROOT
from freeze_final_models import sha


def main():
    plan = json.loads((ROOT / 'data/final_model_freeze.json').read_text())
    panels = pd.read_csv(ROOT / 'data/final_panel_genes.csv')
    pam = set(pd.read_csv(ROOT / 'data/pam50_gene_mapping.csv').feature_id)
    summaries, memberships = [], []
    for name in ['logistic', 'random_forest', 'xgboost']:
        case_id = name + '_automatic'
        final = panels[panels.case_id == case_id].copy()
        folds = pd.read_csv(ROOT / f'data/{name}_panel_gene_selections.csv')
        folds = folds[folds.panel_size == len(final)]
        counts = folds.groupby('feature_id').fold.nunique()
        final['outer_folds_selected'] = final.feature_id.map(counts).fillna(0).astype(int)
        final['in_pam50'] = final.feature_id.isin(pam)
        memberships.extend(final.to_dict('records'))
        summaries.append(dict(case_id=case_id, genes=len(final), pam50_overlap=int(final.in_pam50.sum()), pam50_total=50, final_genes_selected_in_all_five_folds=int((final.outer_folds_selected == 5).sum()), overlap_symbols=sorted(final.loc[final.in_pam50, 'gene_symbol'].dropna().tolist())))
    pd.DataFrame(memberships).to_csv(ROOT / 'data/gene_interpretation_membership.csv', index=False)
    (ROOT / 'data/gene_interpretation_summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
    rows = []
    for family in ['automatic', 'pam50']:
        case_id = 'logistic_' + family
        case = next(c for c in plan['cases'] if c['case_id'] == case_id)
        path = ROOT / case['artifact']
        assert sha(path) == case['artifact_sha256']
        saved = joblib.load(path)
        symbols = panels[panels.case_id == case_id].set_index('feature_id').gene_symbol
        for label, estimator in zip(saved['class_names'], saved['pipeline']['model'].estimators_):
            for feature, weight in zip(saved['feature_ids'], estimator.coef_[0]):
                rows.append(dict(case_id=case_id, subtype=label, feature_id=feature, gene_symbol=symbols.loc[feature], coefficient=float(weight)))
    coefficients = pd.DataFrame(rows)
    coefficients.to_csv(ROOT / 'data/logistic_gene_coefficients.csv', index=False)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    selected = coefficients[coefficients.case_id == 'logistic_automatic']
    for ax, label in zip(axes.flat, plan['classes']):
        table = selected[selected.subtype == label].copy()
        table = table.loc[table.coefficient.abs().nlargest(8).index].sort_values('coefficient')
        ax.barh(table.gene_symbol.fillna(table.feature_id), table.coefficient, color=['#4477AA' if w > 0 else '#CC6677' for w in table.coefficient])
        ax.axvline(0, color='black', linewidth=.7)
        ax.set(title=label, xlabel='Weight per training-standardized log-expression unit')
    fig.suptitle('Frozen automatic 100-gene logistic model: strongest weights\nPositive favors this subtype versus the rest; associations, not biological causes', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(ROOT / 'figures/logistic_gene_weights.png', dpi=180)
    fig.savefig(ROOT / 'figures/logistic_gene_weights.svg')
    plt.close(fig)
    print(json.dumps(summaries, indent=2))
    for label in plan['classes']:
        table = selected[selected.subtype == label]
        print(label, table.nlargest(4, 'coefficient')[['gene_symbol', 'coefficient']].to_dict('records'))


if __name__ == '__main__':
    main()
