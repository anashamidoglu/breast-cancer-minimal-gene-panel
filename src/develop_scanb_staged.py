"""Develop and freeze prespecified staged predictors without reading holdout labels."""
import hashlib
import json
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder
from threadpoolctl import threadpool_limits

from tune_models import ROOT
from evaluate_staged_panels import margin


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair(x, y):
    kwargs = dict(n_estimators=200, class_weight='balanced', min_samples_leaf=3, random_state=42, n_jobs=2)
    small = Pipeline([('genes', SelectKBest(f_classif, k=20)), ('model', RandomForestClassifier(**kwargs))]).fit(x, y)
    large = RandomForestClassifier(**kwargs).fit(x, y)
    assert small['genes'].get_support().sum() == 20
    assert np.isfinite(small['genes'].scores_).all()
    return small, large


def select_rule(x, y, seed):
    a, b, margins = np.zeros(len(y), dtype=int), np.zeros(len(y), dtype=int), np.zeros(len(y))
    for fit, validation in StratifiedKFold(3, shuffle=True, random_state=seed).split(x, y):
        small, large = pair(x.iloc[fit], y[fit])
        probability = small.predict_proba(x.iloc[validation])
        a[validation] = probability.argmax(axis=1)
        b[validation] = large.predict(x.iloc[validation])
        margins[validation] = margin(probability)
    reference = balanced_accuracy_score(y, b)
    rows = []
    for t in [round(i*.05,2) for i in range(21)] + ['route_all']:
        route = np.ones(len(y), dtype=bool) if t == 'route_all' else margins < t
        score = balanced_accuracy_score(y, np.where(route,b,a))
        rows.append(dict(threshold=t, escalation_fraction=float(route.mean()), balanced_accuracy=float(score), eligible=bool(score >= reference-.01)))
    best = min([r for r in rows if r['eligible']], key=lambda r: (r['escalation_fraction'], -r['balanced_accuracy'], 2 if r['threshold']=='route_all' else r['threshold']))
    return best, rows


def main():
    freeze_path = ROOT/'data/scanb_model_freeze.json'
    if freeze_path.exists():
        raise RuntimeError('Models already frozen; no repeat fitting.')
    directory = ROOT/'data/processed/scanb'
    split = json.loads((ROOT/'data/scanb_split_report.json').read_text())
    path = directory/'development_labels.csv'
    assert sha(path) == split['membership_sha256'][path.name]
    labels = pd.read_csv(path)
    x = pd.read_parquet(directory/'pam50_expression.parquet', filters=[('sampleId','in',labels.sampleId.tolist())]).loc[labels.sampleId]
    encoder = LabelEncoder().fit(labels.subtype);y=encoder.transform(labels.subtype)
    folds = StratifiedKFold(5, shuffle=True, random_state=20261010)
    predictions = labels.copy()
    rows, searches = [], []
    for fold,(fit,validation) in enumerate(folds.split(x,y),1):
        print(f'SCAN-B development outer fold {fold}/5...',flush=True)
        with threadpool_limits(limits=2):
            rule, candidates = select_rule(x.iloc[fit],y[fit],42+fold)
            small,large = pair(x.iloc[fit],y[fit])
            p=small.predict_proba(x.iloc[validation]);a=p.argmax(axis=1);b=large.predict(x.iloc[validation])
        route=np.ones(len(validation),dtype=bool) if rule['threshold']=='route_all' else margin(p)<rule['threshold']
        staged=np.where(route,b,a)
        for name,value in [('small',a),('large',b),('staged',staged)]:predictions.loc[validation,name+'_prediction']=value
        predictions.loc[validation,'escalated']=route
        predictions.loc[validation,'fold']=fold
        row=dict(fold=fold,threshold=rule['threshold'],average_genes=float(20+30*route.mean()),escalation_fraction=float(route.mean()))
        for name,value in [('small',a),('large',b),('staged',staged)]:row[name+'_balanced_accuracy']=float(balanced_accuracy_score(y[validation],value))
        rows.append(row);searches.extend(dict(stage=f'outer_{fold}',**r) for r in candidates)
        print(f"Fold {fold}: {row}",flush=True)
    print('Selecting final rule on development patients only...',flush=True)
    with threadpool_limits(limits=2):
        rule,candidates=select_rule(x,y,42)
        small,large=pair(x,y)
    searches.extend(dict(stage='final_development',**r) for r in candidates)
    chosen=x.columns[small['genes'].get_support()].tolist()
    artifact=directory/'staged_models.joblib'
    joblib.dump(dict(small=small,large=large,feature_ids=x.columns.tolist(),small_feature_ids=chosen,class_names=encoder.classes_.tolist()),artifact)
    annotations=pd.read_csv(ROOT/'data/scanb_pam50_mapping.csv').set_index('feature_id')
    pd.DataFrame([dict(feature_id=f,scanb_symbol=annotations.loc[f,'scanb_symbol']) for f in chosen]).to_csv(ROOT/'data/scanb_final_small_panel.csv',index=False)
    predictions.to_csv(directory/'development_oof_predictions.csv',index=False)
    frame=pd.DataFrame(rows);frame.to_csv(ROOT/'data/scanb_development_cv_folds.csv',index=False)
    (ROOT/'data/scanb_threshold_search.json').write_text(json.dumps(searches,indent=2)+'\n')
    report=dict(development_samples=len(labels),models={name:dict(balanced_accuracy_cv_mean=float(frame[name+'_balanced_accuracy'].mean()),balanced_accuracy_cv_sd=float(frame[name+'_balanced_accuracy'].std(ddof=1))) for name in ['small','large','staged']},mean_fold_average_genes=float(frame.average_genes.mean()),final_rule=rule,holdout_used=False)
    (ROOT/'data/scanb_development_report.json').write_text(json.dumps(report,indent=2)+'\n')
    freeze=dict(frozen_at_utc=datetime.now(timezone.utc).isoformat(),threshold=rule['threshold'],artifact=str(artifact.relative_to(ROOT)),artifact_sha256=sha(artifact),small_genes=20,large_genes=50,inner_margin=.01,final_margin=.02,membership_sha256=split['membership_sha256'],expression_mapping_sha256=sha(ROOT/'data/scanb_pam50_mapping.csv'),expression_artifact_sha256=sha(directory/'pam50_expression.parquet'),protocol_sha256=sha(ROOT/'reports/scanb_replication_protocol.md'),holdout_evaluated=False)
    freeze_path.write_text(json.dumps(freeze,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    main()
