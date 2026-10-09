"""Freeze final candidates using training data only, before opening the test set."""
import hashlib
import json
import warnings
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, FunctionTransformer
from threadpoolctl import threadpool_limits

from evaluate_panels import panel_specification
from tune_models import ROOT, SEED, specification, search_rows


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    plan_path = ROOT / 'data/final_model_freeze.json'
    if plan_path.exists() or (ROOT / 'data/final_test_report.json').exists():
        raise RuntimeError('Final candidates already frozen; refusing to refit.')
    processed = ROOT / 'data/processed'
    split = json.loads((ROOT / 'data/split_report.json').read_text())
    train_path = processed / 'train_labels.csv'
    assert sha(train_path) == split['membership_sha256'][train_path.name]
    train = pd.read_csv(train_path)
    x = pd.read_parquet(processed / 'expression.parquet', filters=[('sampleId', 'in', train.sampleId.tolist())]).loc[train.sampleId]
    encoder = LabelEncoder().fit(train.subtype)
    y = encoder.transform(train.subtype)
    annotations = pd.read_csv(processed / 'gene_annotations.csv').set_index('feature_id')
    pam = pd.read_csv(ROOT / 'data/pam50_gene_mapping.csv').feature_id.tolist()
    selections = json.loads((ROOT / 'data/panel_selection_report.json').read_text())['selection']
    cases, genes, searches = [], [], []
    for selection in selections:
        name = selection['model']
        for family in ['pam50', 'automatic']:
            case_id = f'{name}_{family}'
            k = 50 if family == 'pam50' else selection['smallest_evaluated_panel_within_tolerance']
            inputs = x.loc[:, pam] if family == 'pam50' else x
            pipeline, grid = specification(name) if family == 'pam50' else panel_specification(name, k)
            search = GridSearchCV(pipeline, grid, scoring='balanced_accuracy', cv=StratifiedKFold(3, shuffle=True, random_state=SEED), n_jobs=1, error_score='raise')
            print(f'Fitting {case_id} on 756 training patients only...', flush=True)
            with warnings.catch_warnings(), threadpool_limits(limits=2):
                warnings.simplefilter('error', ConvergenceWarning)
                search.fit(inputs, y)
            fitted = search.best_estimator_
            chosen = inputs.columns[fitted['variance'].get_support()]
            if family == 'automatic':
                chosen = chosen[fitted['genes'].get_support()]
            assert len(chosen) == k and chosen.is_unique
            # The final predictor needs only the selected genes. Preserve the fitted
            # scaling/classifier, and verify equivalence before opening the test set.
            compact = Pipeline([('log', FunctionTransformer(np.log1p))] + ([('scale', fitted['scale'])] if name == 'logistic' else []) + [('model', fitted['model'])])
            with threadpool_limits(limits=2):
                np.testing.assert_allclose(compact.predict_proba(x.loc[:, chosen]), fitted.predict_proba(inputs), atol=1e-12)
            artifact = processed / f'{case_id}_final_model.joblib'
            joblib.dump(dict(pipeline=compact, feature_ids=chosen.tolist(), class_names=encoder.classes_.tolist(), training_membership_sha256=sha(train_path)), artifact)
            cases.append(dict(case_id=case_id, model=name, family=family, genes=k, artifact=str(artifact.relative_to(ROOT)), artifact_sha256=sha(artifact), parameters=search.best_params_, training_cv_balanced_accuracy=selection['pam50_gene_set_cv_mean' if family == 'pam50' else 'selected_panel_cv_mean']))
            searches.extend(search_rows(search, case_id))
            genes.extend(dict(case_id=case_id, feature_id=f, gene_symbol=None if pd.isna(annotations.loc[f, 'Hugo_Symbol']) else annotations.loc[f, 'Hugo_Symbol']) for f in chosen)
        artifact = processed / f'{name}_tuned_full_feature_model.joblib'
        saved = joblib.load(artifact)
        assert saved['training_membership_sha256'] == sha(train_path)
        cases.append(dict(case_id=f'{name}_full', model=name, family='full', genes=len(saved['feature_ids']), artifact=str(artifact.relative_to(ROOT)), artifact_sha256=sha(artifact), training_cv_balanced_accuracy=selection['full_gene_reference']))
    dummy = DummyClassifier(strategy='most_frequent').fit(np.zeros((len(y), 1)), y)
    artifact = processed / 'dummy_final_model.joblib'
    joblib.dump(dict(pipeline=dummy, feature_ids=[], class_names=encoder.classes_.tolist(), training_membership_sha256=sha(train_path)), artifact)
    cases.append(dict(case_id='dummy', model='dummy', family='dummy', genes=0, artifact=str(artifact.relative_to(ROOT)), artifact_sha256=sha(artifact), training_cv_balanced_accuracy=0.25))
    pd.DataFrame(genes).to_csv(ROOT / 'data/final_panel_genes.csv', index=False)
    (ROOT / 'data/final_training_search.json').write_text(json.dumps(searches, indent=2) + '\n')
    plan = dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(), primary_case='random_forest_pam50', primary_reason='Highest nested training balanced accuracy, with a fixed 50-gene panel.', comparison_policy='All comparisons are descriptive; test results cannot change the primary model or settings.', training_samples=len(train), membership_sha256=split['membership_sha256'], classes=encoder.classes_.tolist(), cases=cases, test_opened=False)
    plan_path.write_text(json.dumps(plan, indent=2) + '\n')
    print('All ten candidates frozen. Test data have not been opened.', flush=True)


if __name__ == '__main__':
    main()

