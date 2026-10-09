"""First tree-model comparison on saved training folds; fixed settings, no test use."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score, recall_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, LabelEncoder
from sklearn.utils.class_weight import compute_sample_weight
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parents[1]


def make_pipeline(name):
    if name == "random_forest":
        classifier = RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                            class_weight="balanced", random_state=42, n_jobs=4)
    else:
        classifier = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                                   max_bin=64, colsample_bytree=0.5,
                                   tree_method="hist", objective="multi:softprob",
                                   eval_metric="mlogloss", random_state=42, n_jobs=4)
    # Trees do not need standardized units. Keep log/variance processing inside folds.
    return Pipeline([
        ("log", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
        ("variance", VarianceThreshold(threshold=1e-8)),
        ("model", classifier),
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", nargs="+", choices=["random_forest", "xgboost"],
                        default=["random_forest", "xgboost"])
    args = parser.parse_args()
    processed = ROOT / "data" / "processed"
    split = json.loads((ROOT / "data" / "split_report.json").read_text())
    path = processed / "train_labels.csv"
    if hashlib.sha256(path.read_bytes()).hexdigest() != split["membership_sha256"][path.name]:
        raise ValueError("Frozen training membership changed.")
    train = pd.read_csv(path)
    membership = pd.read_csv(processed / "training_cv_folds.csv")
    if not membership[["sampleId", "patientId", "subtype"]].equals(train):
        raise ValueError("Saved folds differ from training membership.")
    x = pd.read_parquet(processed / "expression.parquet",
                        filters=[("sampleId", "in", train["sampleId"].tolist())]).loc[train["sampleId"]]
    # XGBoost expects integer class codes. Encoding class names learns no expression patterns.
    encoder = LabelEncoder()
    y = encoder.fit_transform(train["subtype"])
    classes = encoder.classes_
    for name in args.models:
        predictions = np.zeros(len(y), dtype=int)
        probabilities = np.zeros((len(y), len(classes)))
        rows = []
        for fold in sorted(membership["cv_fold"].unique()):
            validation = np.flatnonzero(membership["cv_fold"].to_numpy() == fold)
            fitting = np.flatnonzero(membership["cv_fold"].to_numpy() != fold)
            model = make_pipeline(name)
            fit_params = {}
            if name == "xgboost":
                # XGBoost has no multiclass class_weight option: compute equivalent
                # per-sample weights from this fitting fold's labels only.
                fit_params["model__sample_weight"] = compute_sample_weight("balanced", y[fitting])
            print(f"{name}: fitting training fold {fold}/5...", flush=True)
            model.fit(x.iloc[fitting], y[fitting], **fit_params)
            prediction = model.predict(x.iloc[validation]).astype(int)
            probability = model.predict_proba(x.iloc[validation])
            if not np.array_equal(model.classes_, np.arange(len(classes))):
                raise ValueError("Unexpected class order.")
            predictions[validation] = prediction
            probabilities[validation] = probability
            score = balanced_accuracy_score(y[validation], prediction)
            rows.append({"fold": int(fold), "fit_samples": len(fitting), "validation_samples": len(validation),
                         "gene_features_after_variance_filter": int(model["variance"].get_support().sum()),
                         "balanced_accuracy": score,
                         "macro_f1": f1_score(y[validation], prediction, average="macro", zero_division=0)})
            print(f"{name} fold {fold} balanced accuracy: {score:.3f}", flush=True)
        folds = pd.DataFrame(rows)
        folds.to_csv(ROOT / "data" / f"{name}_initial_cv_folds.csv", index=False)
        recall = recall_score(y, predictions, labels=np.arange(len(classes)), average=None, zero_division=0)
        report = {
            "model": name, "evaluation": "Fixed initial settings on five saved training folds; no tuning yet",
            "training_samples": len(y), "input_gene_features": x.shape[1], "classes": classes.tolist(),
            "parameters": {key: value for key, value in make_pipeline(name)["model"].get_params().items()
                           if key in ["n_estimators", "max_depth", "min_samples_leaf", "max_features", "learning_rate", "class_weight", "random_state", "tree_method", "max_bin", "colsample_bytree"]},
            "preprocessing": "log1p and variance threshold 1e-8 inside each fold; no standardization needed for trees",
            "imbalance_handling": "balanced class weights" if name == "random_forest" else "balanced sample weights computed separately from each fitting fold",
            "balanced_accuracy_cv_mean": float(folds["balanced_accuracy"].mean()),
            "balanced_accuracy_cv_sd": float(folds["balanced_accuracy"].std(ddof=1)),
            "macro_f1_cv_mean": float(folds["macro_f1"].mean()),
            "macro_f1_cv_sd": float(folds["macro_f1"].std(ddof=1)),
            "out_of_fold_per_class_recall": dict(zip(classes, recall.tolist())),
            "out_of_fold_confusion_matrix": confusion_matrix(y, predictions).tolist(),
            "confusion_matrix_orientation": "Recorded labels in rows, predictions in columns, class order as listed",
            "out_of_fold_one_vs_rest_auroc": {label: roc_auc_score((y == i).astype(int), probabilities[:, i])
                                             for i, label in enumerate(classes)},
            "test_evaluated": False,
            "note": "Preliminary full-feature model comparison, not tuned or small-panel results. Fold SD is not a confidence interval.",
        }
        (ROOT / "data" / f"{name}_initial_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        output = train.assign(predicted_subtype=encoder.inverse_transform(predictions))
        for i, label in enumerate(classes):
            output[f"probability_{label}"] = probabilities[:, i]
        output.to_csv(processed / f"{name}_initial_oof_predictions.csv", index=False)
        print(f"{name} mean balanced accuracy: {report['balanced_accuracy_cv_mean']:.3f}", flush=True)


if __name__ == "__main__":
    main()
