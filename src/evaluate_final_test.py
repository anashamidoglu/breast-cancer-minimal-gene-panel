"""Evaluate frozen predictors once on the held-out patients; never fit here."""
import json
from datetime import datetime, timezone

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score, roc_auc_score
from threadpoolctl import threadpool_limits

from freeze_final_models import ROOT, sha


def main():
    output = ROOT / 'data/final_test_report.json'
    if output.exists():
        raise RuntimeError('Test already evaluated; inspect the saved report.')
    plan_path = ROOT / 'data/final_model_freeze.json'
    plan = json.loads(plan_path.read_text())
    processed = ROOT / 'data/processed'
    for case in plan['cases']:
        assert sha(ROOT / case['artifact']) == case['artifact_sha256']
    for filename, digest in plan['membership_sha256'].items():
        assert sha(processed / filename) == digest
    train = pd.read_csv(processed / 'train_labels.csv')
    test = pd.read_csv(processed / 'test_labels.csv')
    assert not set(train.patientId) & set(test.patientId)
    assert len(test) == 189 and test.sampleId.is_unique
    x = pd.read_parquet(processed / 'expression.parquet', filters=[('sampleId', 'in', test.sampleId.tolist())]).loc[test.sampleId]
    classes = plan['classes']
    y = np.array([classes.index(label) for label in test.subtype])
    reports, rows = [], []
    for case in plan['cases']:
        saved = joblib.load(ROOT / case['artifact'])
        assert saved['class_names'] == classes
        assert saved['training_membership_sha256'] == plan['membership_sha256']['train_labels.csv']
        inputs = x.loc[:, saved['feature_ids']] if saved['feature_ids'] else np.zeros((len(test), 1))
        with threadpool_limits(limits=2):
            probability = saved['pipeline'].predict_proba(inputs)
            prediction = saved['pipeline'].predict(inputs).astype(int)
        assert np.isfinite(probability).all()
        np.testing.assert_allclose(probability.sum(axis=1), 1, atol=1e-6)
        matrix = confusion_matrix(y, prediction, labels=np.arange(len(classes)))
        score = float(balanced_accuracy_score(y, prediction))
        macro = float(f1_score(y, prediction, average='macro', zero_division=0))
        report = dict(case_id=case['case_id'], genes=case['genes'], balanced_accuracy=score, macro_f1=macro, accuracy=float(np.mean(y == prediction)), confusion_matrix=matrix.tolist(), per_class={label: dict(correct=int(matrix[i, i]), total=int(matrix[i].sum()), recall=float(matrix[i, i] / matrix[i].sum()), one_vs_rest_auroc=float(roc_auc_score(y == i, probability[:, i]))) for i, label in enumerate(classes)})
        reports.append(report)
        rows.append(dict(case_id=case['case_id'], genes=case['genes'], primary=case['case_id'] == plan['primary_case'], training_cv_balanced_accuracy=case['training_cv_balanced_accuracy'], test_balanced_accuracy=score, test_macro_f1=macro))
        predictions = test.assign(predicted_subtype=[classes[i] for i in prediction])
        for i, label in enumerate(classes):
            predictions[f'probability_{label}'] = probability[:, i]
        predictions.to_csv(processed / f"final_test_{case['case_id']}_predictions.csv", index=False)
        print(f"{case['case_id']}: balanced accuracy {score:.4f}, macro F1 {macro:.4f}", flush=True)
    result = dict(evaluated_at_utc=datetime.now(timezone.utc).isoformat(), frozen_plan_sha256=sha(plan_path), primary_case=plan['primary_case'], test_samples=len(test), classes=classes, confusion_matrix_orientation='Recorded labels in rows; predictions in columns.', cases=reports, limitations=['PAM50-derived labels: agreement with those labels, not independent clinical diagnosis.', 'Single-cohort holdout; external validation remains necessary.', 'Training CV selected the candidates; supporting test comparisons do not select a new winner.', 'The 2 percentage-point panel rule was a training selection rule, not statistical equivalence.'])
    output.write_text(json.dumps(result, indent=2) + '\n')
    pd.DataFrame(rows).to_csv(ROOT / 'data/final_test_comparison.csv', index=False)
    primary = next(r for r in reports if r['case_id'] == plan['primary_case'])
    fig, ax = plt.subplots(figsize=(7.5, 6))
    matrix = np.array(primary['confusion_matrix'])
    ax.imshow(matrix, cmap='Blues')
    for i in range(4):
        for j in range(4):
            ax.text(j, i, str(matrix[i, j]), ha='center', va='center', color='white' if matrix[i, j] > matrix.max()/2 else 'black')
    ax.set(xticks=range(4), yticks=range(4), xticklabels=classes, yticklabels=classes, xlabel='Predicted subtype', ylabel='Recorded subtype', title=f"Frozen primary: random forest on PAM50 genes\n189 held-out patients; balanced accuracy {primary['balanced_accuracy']:.1%}")
    plt.setp(ax.get_xticklabels(), rotation=25, ha='right')
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/final_test_confusion.png', dpi=180)
    fig.savefig(ROOT / 'figures/final_test_confusion.svg')
    plt.close(fig)


if __name__ == '__main__':
    main()
