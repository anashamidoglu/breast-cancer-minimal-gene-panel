"""Nested training-only assessment of the prespecified 20-to-100 gene strategy."""
import argparse
import json
import hashlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import balanced_accuracy_score, f1_score, recall_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from threadpoolctl import threadpool_limits

from evaluate_panels import panel_specification
from tune_models import ROOT, specification


SMALL_COUNT = 20
LARGE_COUNT = 100
PAM_ONLY = False

def margin(probability):
    ordered = np.sort(probability, axis=1)
    return ordered[:, -1] - ordered[:, -2]


def fit_pair(x, y):
    small, _ = panel_specification('random_forest', SMALL_COUNT)
    large, _ = panel_specification('random_forest', LARGE_COUNT)
    for pipeline in [small, large]:
        pipeline.set_params(model__min_samples_leaf=3).fit(x, y)
    survivors = x.columns[large['variance'].get_support()]
    chosen_small = set(survivors[small['genes'].get_support()])
    chosen_large = set(survivors[large['genes'].get_support()])
    assert chosen_small <= chosen_large and len(chosen_small) == SMALL_COUNT and len(chosen_large) == LARGE_COUNT
    np.testing.assert_allclose(small['genes'].scores_, large['genes'].scores_)
    return small, large


def choose_threshold(x, y, seed):
    p_small = np.zeros((len(y), 4))
    p_large = np.zeros((len(y), 4))
    for fit, validation in StratifiedKFold(3, shuffle=True, random_state=seed).split(x, y):
        small, large = fit_pair(x.iloc[fit], y[fit])
        p_small[validation] = small.predict_proba(x.iloc[validation])
        p_large[validation] = large.predict_proba(x.iloc[validation])
    a, b = p_small.argmax(axis=1), p_large.argmax(axis=1)
    reference = balanced_accuracy_score(y, b)
    rows = []
    for threshold in [round(i * .05, 2) for i in range(21)] + ['route_all']:
        route = np.ones(len(y), dtype=bool) if threshold == 'route_all' else margin(p_small) < threshold
        score = balanced_accuracy_score(y, np.where(route, b, a))
        rows.append(dict(threshold=threshold, escalation_fraction=float(route.mean()), balanced_accuracy=float(score), eligible=bool(score >= reference - .02)))
    eligible = [r for r in rows if r['eligible']]
    best = min(eligible, key=lambda r: (r['escalation_fraction'], -r['balanced_accuracy'], 2 if r['threshold'] == 'route_all' else r['threshold']))
    return best, rows, float(reference)


def main():
    global LARGE_COUNT, PAM_ONLY
    parser = argparse.ArgumentParser()
    parser.add_argument('--pam50-staged', action='store_true')
    args = parser.parse_args()
    PAM_ONLY = args.pam50_staged
    LARGE_COUNT = 50 if PAM_ONLY else 100
    prefix = 'pam50_staged' if PAM_ONLY else 'staged'
    output = ROOT / f'data/{prefix}_training_report.json'
    if output.exists():
        raise RuntimeError('Results already exist; refusing to repeat or overwrite.')
    processed = ROOT / 'data/processed'
    split = json.loads((ROOT / 'data/split_report.json').read_text())
    train_path = processed / 'train_labels.csv'
    assert hashlib.sha256(train_path.read_bytes()).hexdigest() == split['membership_sha256'][train_path.name]
    train = pd.read_csv(train_path)
    folds = pd.read_csv(processed / 'training_cv_folds.csv')
    assert folds[['sampleId', 'patientId', 'subtype']].equals(train)
    x = pd.read_parquet(processed / 'expression.parquet', filters=[('sampleId', 'in', train.sampleId.tolist())]).loc[train.sampleId]
    encoder = LabelEncoder().fit(train.subtype)
    y = encoder.transform(train.subtype)
    pam = pd.read_csv(ROOT / 'data/pam50_gene_mapping.csv').feature_id.tolist()
    if PAM_ONLY:
        assert len(pam) == 50 and len(set(pam)) == 50
        x = x.loc[:, pam]
    oof = train.copy()
    scores, searches, curves = [], [], []
    for name in ['small', 'large', 'staged', 'pam50']:
        oof[name + '_prediction'] = -1
    oof['escalated'] = False
    oof['small_margin'] = 0.
    for fold in sorted(folds.cv_fold.unique()):
        fit = np.flatnonzero(folds.cv_fold.to_numpy() != fold)
        validation = np.flatnonzero(folds.cv_fold.to_numpy() == fold)
        print(f'Outer fold {fold}/5: training-only threshold selection...', flush=True)
        with threadpool_limits(limits=2):
            selected, inner, reference = choose_threshold(x.iloc[fit], y[fit], 42 + int(fold))
            small, large = fit_pair(x.iloc[fit], y[fit])
            pam_model, _ = specification('random_forest')
            pam_model.set_params(model__min_samples_leaf=3).fit(x.iloc[fit].loc[:, pam], y[fit])
            p = small.predict_proba(x.iloc[validation])
            a = p.argmax(axis=1)
            b = large.predict(x.iloc[validation]).astype(int)
            c = pam_model.predict(x.iloc[validation].loc[:, pam]).astype(int)
        threshold = selected['threshold']
        route = np.ones(len(validation), dtype=bool) if threshold == 'route_all' else margin(p) < threshold
        staged = np.where(route, b, a)
        predictions = dict(small=a, large=b, staged=staged, pam50=c)
        for name, prediction in predictions.items():
            oof.loc[validation, name + '_prediction'] = prediction
        oof.loc[validation, 'escalated'] = route
        oof.loc[validation, 'small_margin'] = margin(p)
        count = int(route.sum())
        random_scores = []
        for seed in range(100):
            rng = np.random.default_rng(42000 + int(fold)*100 + seed)
            random_route = np.zeros(len(validation), dtype=bool)
            random_route[rng.choice(len(validation), count, replace=False)] = True
            random_scores.append(balanced_accuracy_score(y[validation], np.where(random_route, b, a)))
        row = dict(fold=int(fold), selected_threshold=threshold, inner_large_reference=reference, inner_staged_score=selected['balanced_accuracy'], escalation_fraction=float(route.mean()), average_genes=float(SMALL_COUNT + (LARGE_COUNT-SMALL_COUNT)*route.mean()), random_escalation_balanced_accuracy_mean=float(np.mean(random_scores)), random_escalation_balanced_accuracy_sd=float(np.std(random_scores, ddof=1)))
        for name, prediction in predictions.items():
            row[name + '_balanced_accuracy'] = float(balanced_accuracy_score(y[validation], prediction))
            row[name + '_macro_f1'] = float(f1_score(y[validation], prediction, average='macro'))
        row['errors_corrected'] = int((route & (a != y[validation]) & (b == y[validation])).sum())
        row['errors_introduced'] = int((route & (a == y[validation]) & (b != y[validation])).sum())
        row['errors_remaining'] = int((route & (a != y[validation]) & (b != y[validation])).sum())
        scores.append(row)
        searches.extend(dict(fold=int(fold), **r) for r in inner)
        # A descriptive outer curve, not used to select thresholds or a new model.
        for t in [round(i*.05, 2) for i in range(21)] + ['route_all']:
            mask = np.ones(len(validation), dtype=bool) if t == 'route_all' else margin(p) < t
            curves.append(dict(fold=int(fold), threshold=t, average_genes=float(SMALL_COUNT+(LARGE_COUNT-SMALL_COUNT)*mask.mean()), balanced_accuracy=float(balanced_accuracy_score(y[validation], np.where(mask,b,a)))))
        print(f"Fold {fold}: threshold={threshold}; genes={row['average_genes']:.1f}; staged={row['staged_balanced_accuracy']:.3f}; large={row['large_balanced_accuracy']:.3f}", flush=True)
    frame = pd.DataFrame(scores)
    assert not (oof.filter(like='_prediction') < 0).any().any()
    oof.to_csv(processed / f'{prefix}_training_oof_predictions.csv', index=False)
    frame.to_csv(ROOT / f'data/{prefix}_training_cv_folds.csv', index=False)
    pd.DataFrame(curves).to_csv(ROOT / f'data/{prefix}_training_tradeoff.csv', index=False)
    (ROOT / f'data/{prefix}_inner_threshold_search.json').write_text(json.dumps(searches, indent=2)+'\n')
    models = {}
    for name in ['small','large','staged','pam50']:
        models[name] = dict(balanced_accuracy_cv_mean=float(frame[name+'_balanced_accuracy'].mean()), balanced_accuracy_cv_sd=float(frame[name+'_balanced_accuracy'].std(ddof=1)), macro_f1_cv_mean=float(frame[name+'_macro_f1'].mean()), pooled_recall=dict(zip(encoder.classes_, recall_score(y,oof[name+'_prediction'],average=None).tolist())))
    result = dict(panel_family='PAM50 subset' if PAM_ONLY else 'All-gene ANOVA', small_genes=SMALL_COUNT, large_genes=LARGE_COUNT, evaluation='Five outer training folds with three inner folds selecting routing thresholds; fixed classifier settings.', training_samples=len(train), test_used=False, external_evaluation_completed=False, models=models, average_genes=float(SMALL_COUNT+(LARGE_COUNT-SMALL_COUNT)*oof.escalated.mean()), escalation_fraction=float(oof.escalated.mean()), random_escalation_balanced_accuracy_cv_mean=float(frame.random_escalation_balanced_accuracy_mean.mean()), errors_corrected=int(frame.errors_corrected.sum()), errors_introduced=int(frame.errors_introduced.sum()), errors_remaining=int(frame.errors_remaining.sum()), escalation_by_subtype=oof.groupby('subtype').escalated.mean().to_dict(), tolerance=.02, difference_from_large=models['staged']['balanced_accuracy_cv_mean']-models['large']['balanced_accuracy_cv_mean'], notes=['No fresh external validation yet.', 'Outer tradeoff curve is descriptive and must not select a new threshold.', 'Random-routing SD measures variation among routing draws, not a confidence interval.', 'PAM50 comparator uses fixed leaf size 3; earlier PAM50 results used nested setting searches.'])
    output.write_text(json.dumps(result, indent=2)+'\n')
    fig, ax = plt.subplots(figsize=(8,5))
    curve = pd.DataFrame(curves).groupby('threshold',sort=False)[['average_genes','balanced_accuracy']].mean().sort_values('average_genes')
    ax.plot(curve.average_genes,curve.balanced_accuracy,color='#999999',label='Descriptive threshold curve')
    for name,count,color in [('small',SMALL_COUNT,'#4477AA'),('large',LARGE_COUNT,'#228833'),('pam50',50,'#AA3377'),('staged',result['average_genes'],'#CC6677')]:
        ax.scatter(count,models[name]['balanced_accuracy_cv_mean'],color=color,s=70,label=name)
    ax.scatter(result['average_genes'],result['random_escalation_balanced_accuracy_cv_mean'],marker='x',color='black',label='Random escalation, equal burden')
    ax.set(xlim=(15,LARGE_COUNT+5),ylim=(0,1),xlabel='Average distinct genes per patient',ylabel='Mean outer-fold balanced accuracy',title='Training-only staged panel evaluation')
    ax.legend(loc='lower right');fig.tight_layout();fig.savefig(ROOT/f'figures/{prefix}_training_tradeoff.png',dpi=180);fig.savefig(ROOT/f'figures/{prefix}_training_tradeoff.svg');plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    main()
