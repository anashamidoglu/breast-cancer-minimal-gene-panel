"""Nested, training-only model tuning with a small predefined search budget."""

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score, recall_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.multiclass import OneVsRestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, LabelEncoder, StandardScaler
from threadpoolctl import threadpool_limits
from estimators import BalancedXGBClassifier

ROOT = Path(__file__).resolve().parents[1]
SEED = 42


def specification(name):
    steps = [("log", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
             ("variance", VarianceThreshold(threshold=1e-8))]
    if name == "logistic":
        steps.append(("scale", StandardScaler()))
        # liblinear solves four binary L1 models efficiently for this wide table.
        model = OneVsRestClassifier(LogisticRegression(solver="liblinear", l1_ratio=1.0,
                                   class_weight="balanced", max_iter=3000, random_state=SEED))
        grid = {"model__estimator__C": [0.01, 0.1, 1.0]}
    elif name == "random_forest":
        model = RandomForestClassifier(n_estimators=200, class_weight="balanced", random_state=SEED, n_jobs=2)
        grid = {"model__min_samples_leaf": [1, 3]}
    else:
        model = BalancedXGBClassifier(n_estimators=100, learning_rate=0.05, max_bin=32,
                                     colsample_bytree=0.5, tree_method="hist", objective="multi:softprob",
                                     eval_metric="mlogloss", random_state=SEED, n_jobs=2)
        grid = {"model__max_depth": [2, 3]}
    steps.append(("model", model))
    return Pipeline(steps), grid


def search_model(name, x, y, seed):
    pipeline, grid = specification(name)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    search = GridSearchCV(pipeline, grid, scoring="balanced_accuracy", cv=inner,
                          n_jobs=1, refit=True, error_score="raise", return_train_score=False)
    with warnings.catch_warnings(), threadpool_limits(limits=2):
        warnings.simplefilter("error", ConvergenceWarning)
        search.fit(x, y)
    return search


def search_rows(search, stage):
    return [{"stage": stage, "parameters": params,
             "mean_inner_balanced_accuracy": float(mean), "sd_inner_balanced_accuracy": float(sd),
             "rank": int(rank)} for params, mean, sd, rank in zip(
                 search.cv_results_["params"], search.cv_results_["mean_test_score"],
                 search.cv_results_["std_test_score"], search.cv_results_["rank_test_score"])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, choices=["logistic", "random_forest", "xgboost"])
    args = parser.parse_args()
    name = args.model
    processed = ROOT / "data" / "processed"
    split = json.loads((ROOT / "data" / "split_report.json").read_text())
    train_path = processed / "train_labels.csv"
    if hashlib.sha256(train_path.read_bytes()).hexdigest() != split["membership_sha256"][train_path.name]:
        raise ValueError("Training membership changed.")
    train = pd.read_csv(train_path)
    membership = pd.read_csv(processed / "training_cv_folds.csv")
    if not membership[["sampleId", "patientId", "subtype"]].equals(train):
        raise ValueError("Saved folds and training labels differ.")
    x = pd.read_parquet(processed / "expression.parquet",
                        filters=[("sampleId", "in", train["sampleId"].tolist())]).loc[train["sampleId"]]
    encoder = LabelEncoder()
    y = encoder.fit_transform(train["subtype"])
    classes = encoder.classes_
    predictions = np.zeros(len(y), dtype=int)
    probabilities = np.zeros((len(y), len(classes)))
    rows, searches = [], []
    for fold in sorted(membership["cv_fold"].unique()):
        validation = np.flatnonzero(membership["cv_fold"].to_numpy() == fold)
        fitting = np.flatnonzero(membership["cv_fold"].to_numpy() != fold)
        print(f"{name}: tuning within outer training fold {fold}/5...", flush=True)
        search = search_model(name, x.iloc[fitting], y[fitting], SEED + int(fold))
        prediction = search.predict(x.iloc[validation]).astype(int)
        probability = search.predict_proba(x.iloc[validation])
        if not np.array_equal(search.classes_, np.arange(len(classes))):
            raise ValueError("Class ordering changed.")
        predictions[validation] = prediction
        probabilities[validation] = probability
        score = balanced_accuracy_score(y[validation], prediction)
        rows.append({"fold": int(fold), "fit_samples": len(fitting), "validation_samples": len(validation),
                     "balanced_accuracy": score,
                     "macro_f1": f1_score(y[validation], prediction, average="macro", zero_division=0),
                     "best_parameters": json.dumps(search.best_params_, sort_keys=True),
                     "features_after_variance_filter": int(search.best_estimator_["variance"].get_support().sum())})
        searches.extend(search_rows(search, f"outer_fold_{fold}"))
        print(f"{name} outer fold {fold}: {score:.3f}; selected {search.best_params_}", flush=True)
        # Checkpoint aggregate progress, so completed folds are reviewable.
        pd.DataFrame(rows).to_csv(ROOT / "data" / f"{name}_tuned_cv_folds.csv", index=False)
    print(f"{name}: selecting deployment settings using all training samples only...", flush=True)
    final_search = search_model(name, x, y, SEED)
    searches.extend(search_rows(final_search, "all_training_setting_selection"))
    pipeline, grid = specification(name)
    folds = pd.DataFrame(rows)
    recall = recall_score(y, predictions, labels=np.arange(len(classes)), average=None, zero_division=0)
    report = {
        "model": name,
        "evaluation": "Five saved outer training folds; three inner stratified folds select settings separately within each outer fit set",
        "search_grid": grid, "fixed_parameters": pipeline["model"].get_params(deep=False),
        "final_training_selected_parameters": final_search.best_params_,
        "training_samples": len(y), "input_gene_features": x.shape[1], "classes": classes.tolist(),
        "balanced_accuracy_cv_mean": float(folds["balanced_accuracy"].mean()),
        "balanced_accuracy_cv_sd": float(folds["balanced_accuracy"].std(ddof=1)),
        "macro_f1_cv_mean": float(folds["macro_f1"].mean()),
        "macro_f1_cv_sd": float(folds["macro_f1"].std(ddof=1)),
        "out_of_fold_per_class_recall": dict(zip(classes, recall.tolist())),
        "out_of_fold_confusion_matrix": confusion_matrix(y, predictions).tolist(),
        "confusion_matrix_orientation": "Recorded labels in rows; predictions in columns; class order as listed",
        "out_of_fold_one_vs_rest_auroc": {label: roc_auc_score((y == i).astype(int), probabilities[:, i])
                                         for i, label in enumerate(classes)},
        "test_evaluated": False,
        "notes": ["Bounded search, not exhaustive optimization.",
                  "Logistic is L1 one-vs-rest, differing from the initial multinomial L2 model.",
                  "Tree count/bin settings are smaller than initial pilots for nested-CV runtime.",
                  "Preprocessing and class-balancing weights are fitted within every inner/outer fit.",
                  "Nested CV assesses the tuning procedure; final settings are selected separately on all training data.",
                  "Fold SD is not a confidence interval; final test and small-panel analysis remain pending."],
    }
    # Estimator objects in fixed parameters are described as text in aggregate JSON.
    (ROOT / "data" / f"{name}_tuned_report.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    (ROOT / "data" / f"{name}_tuning_search.json").write_text(json.dumps(searches, indent=2) + "\n", encoding="utf-8")
    output = train.assign(predicted_subtype=encoder.inverse_transform(predictions))
    for i, label in enumerate(classes):
        output[f"probability_{label}"] = probabilities[:, i]
    output.to_csv(processed / f"{name}_tuned_oof_predictions.csv", index=False)
    joblib.dump({"pipeline": final_search.best_estimator_, "class_names": classes.tolist(),
                 "feature_ids": x.columns.tolist(), "training_membership_sha256": split["membership_sha256"][train_path.name]},
                processed / f"{name}_tuned_full_feature_model.joblib")
    print(f"{name} nested-CV balanced accuracy: {report['balanced_accuracy_cv_mean']:.4f}", flush=True)


if __name__ == "__main__":
    main()
