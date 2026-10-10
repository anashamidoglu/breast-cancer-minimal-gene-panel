"""One-time independent-cohort holdout assessment of frozen staged predictors."""
import json
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score, recall_score
from threadpoolctl import threadpool_limits

from develop_scanb_staged import ROOT, sha
from evaluate_staged_panels import margin


def main():
    output=ROOT/'data/scanb_holdout_report.json'
    if output.exists():raise RuntimeError('Holdout already evaluated; inspect saved results.')
    plan_path=ROOT/'data/scanb_model_freeze.json';plan=json.loads(plan_path.read_text())
    directory=ROOT/'data/processed/scanb'
    assert sha(ROOT/plan['artifact'])==plan['artifact_sha256']
    assert sha(directory/'pam50_expression.parquet')==plan['expression_artifact_sha256']
    assert sha(ROOT/'reports/scanb_replication_protocol.md')==plan['protocol_sha256']
    for name,digest in plan['membership_sha256'].items():assert sha(directory/name)==digest
    labels=pd.read_csv(directory/'holdout_labels.csv');development=pd.read_csv(directory/'development_labels.csv')
    assert len(labels)==916 and not set(labels.patientId)&set(development.patientId)
    saved=joblib.load(ROOT/plan['artifact']);classes=saved['class_names'];y=np.array([classes.index(v) for v in labels.subtype])
    x=pd.read_parquet(directory/'pam50_expression.parquet',filters=[('sampleId','in',labels.sampleId.tolist())]).loc[labels.sampleId,saved['feature_ids']]
    with threadpool_limits(limits=2):
        p=saved['small'].predict_proba(x);a=p.argmax(axis=1);b=saved['large'].predict(x).astype(int)
    np.testing.assert_allclose(p.sum(axis=1),1)
    threshold=plan['threshold'];route=np.ones(len(y),dtype=bool) if threshold=='route_all' else margin(p)<threshold
    staged=np.where(route,b,a)
    models={}
    for name,prediction in [('small',a),('large',b),('staged',staged)]:
        models[name]=dict(balanced_accuracy=float(balanced_accuracy_score(y,prediction)),macro_f1=float(f1_score(y,prediction,average='macro')),recall=dict(zip(classes,recall_score(y,prediction,average=None).tolist())),correct=int((y==prediction).sum()))
    random_scores=[]
    for seed in range(100):
        rng=np.random.default_rng(42000+seed);mask=np.zeros(len(y),dtype=bool);mask[rng.choice(len(y),int(route.sum()),replace=False)]=True
        random_scores.append(balanced_accuracy_score(y,np.where(mask,b,a)))
    # Resample the same patients for both strategies within each recorded class.
    rng=np.random.default_rng(20261010);groups=[np.flatnonzero(y==i) for i in range(4)];differences=[]
    for _ in range(5000):
        differences.append(float(np.mean([np.mean(staged[idx]==y[idx])-np.mean(b[idx]==y[idx]) for group in groups for idx in [rng.choice(group,len(group),replace=True)]])))
    low,high=np.quantile(differences,[.025,.975]);difference=models['staged']['balanced_accuracy']-models['large']['balanced_accuracy']
    report=dict(evaluation='Fresh patient-disjoint SCAN-B holdout after SCAN-B development; independent-cohort strategy replication, not TCGA model transfer.',evaluated_at_utc=datetime.now(timezone.utc).isoformat(),freeze_sha256=sha(plan_path),holdout_samples=len(y),threshold=threshold,models=models,average_genes=float(20+30*route.mean()),escalation_fraction=float(route.mean()),escalated_cases=int(route.sum()),difference_from_large=float(difference),paired_stratified_bootstrap_95_interval=[float(low),float(high)],bootstrap_repetitions=5000,point_estimate_target_met=bool(difference>=-.02 and route.mean()<1),interval_lower_bound_above_minus_2_points=bool(low>-.02),random_escalation_balanced_accuracy_mean=float(np.mean(random_scores)),random_escalation_balanced_accuracy_sd=float(np.std(random_scores,ddof=1)),errors_corrected=int((route&(a!=y)&(b==y)).sum()),errors_introduced=int((route&(a==y)&(b!=y)).sum()),errors_remaining=int((route&(a!=y)&(b!=y)).sum()),escalation_by_subtype={label:float(route[y==i].mean()) for i,label in enumerate(classes)})
    predictions=labels.assign(small_prediction=a,large_prediction=b,staged_prediction=staged,escalated=route,small_margin=margin(p))
    predictions.to_csv(directory/'holdout_predictions.csv',index=False)
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
