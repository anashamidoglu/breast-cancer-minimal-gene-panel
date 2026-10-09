"""Evaluate a fixed first logistic model on training folds; never use the test set."""

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score, recall_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]


def make_pipeline():
    # Every learned preprocessing step is fitted separately within each fit fold.
    # log1p is log(1 + expression), which handles zero expression safely.
    return Pipeline([
        ("log", FunctionTransformer(np.log1p, feature_names_out="one-to-one")),
        ("variance", VarianceThreshold(threshold=1e-8)),
        ("scale", StandardScaler()),
        # l1_ratio=0 gives L2 regularization in this installed sklearn version.
        # This fixed initial reference precedes tuned L1/elastic-net comparisons.
        ("model", LogisticRegression(C=1.0, l1_ratio=0.0, solver="lbfgs",
                                     class_weight="balanced", max_iter=2000, random_state=42)),
    ])


def main():
    processed = ROOT / "data" / "processed"
    split = json.loads((ROOT / "data" / "split_report.json").read_text())
    train_path = processed / "train_labels.csv"
    if hashlib.sha256(train_path.read_bytes()).hexdigest() != split["membership_sha256"][train_path.name]:
        raise ValueError("Frozen training membership changed.")
    train = pd.read_csv(train_path)
    membership = pd.read_csv(processed / "training_cv_folds.csv")
    if not membership[["sampleId", "patientId", "subtype"]].equals(train):
        raise ValueError("Saved folds do not match the frozen training labels.")
    # Parquet row filtering retrieves only training samples, then aligns their order.
    x = pd.read_parquet(processed / "expression.parquet",
                        filters=[("sampleId", "in", train["sampleId"].tolist())])
    x = x.loc[train["sampleId"]]
    if len(x) != len(train) or not np.isfinite(x.to_numpy()).all() or (x.to_numpy() < 0).any():
        raise ValueError("Training expression failed alignment/value checks.")
    y = train["subtype"].to_numpy()
    classes = np.unique(y)
    predictions = np.empty(len(y), dtype=object)
    probabilities = np.zeros((len(y), len(classes)))
    rows = []
    for fold in sorted(membership["cv_fold"].unique()):
        validation = np.flatnonzero(membership["cv_fold"].to_numpy() == fold)
        fitting = np.flatnonzero(membership["cv_fold"].to_numpy() != fold)
        model = make_pipeline()
        print(f"Fitting training fold {fold}/5...", flush=True)
        # A convergence failure should stop the run, not produce a misleading score.
        with warnings.catch_warnings(), threadpool_limits(limits=4):
            warnings.simplefilter("error", ConvergenceWarning)
            model.fit(x.iloc[fitting], y[fitting])
            prediction = model.predict(x.iloc[validation])
            probability = model.predict_proba(x.iloc[validation])
        if not np.array_equal(model.classes_, classes):
            raise ValueError("Unexpected class order.")
        predictions[validation] = prediction
        probabilities[validation] = probability
        score = balanced_accuracy_score(y[validation], prediction)
        rows.append({"fold": int(fold), "fit_samples": len(fitting), "validation_samples": len(validation),
                     "gene_features_after_variance_filter": int(model["variance"].get_support().sum()),
                     "balanced_accuracy": score,
                     "macro_f1": f1_score(y[validation], prediction, average="macro", zero_division=0),
                     "max_solver_iterations": int(model["model"].n_iter_.max())})
        print(f"Fold {fold} balanced accuracy: {score:.3f}", flush=True)
    folds = pd.DataFrame(rows)
    folds.to_csv(ROOT / "data" / "logistic_initial_cv_folds.csv", index=False)
    recall = recall_score(y, predictions, labels=classes, average=None, zero_division=0)
    report = {
        "model": "Initial all-available-gene logistic regression, fixed L2 regularization",
        "evaluation": "Five saved stratified training folds; no hyperparameter search yet",
        "parameters": {"C": 1.0, "l1_ratio": 0.0, "solver": "lbfgs", "class_weight": "balanced", "max_iter": 2000},
        "preprocessing": "log1p, variance threshold 1e-8 on log scale, StandardScaler, all inside Pipeline",
        "training_samples": len(y), "input_gene_features": x.shape[1], "classes": classes.tolist(),
        "balanced_accuracy_cv_mean": float(folds["balanced_accuracy"].mean()),
        "balanced_accuracy_cv_sd": float(folds["balanced_accuracy"].std(ddof=1)),
        "macro_f1_cv_mean": float(folds["macro_f1"].mean()),
        "macro_f1_cv_sd": float(folds["macro_f1"].std(ddof=1)),
        "out_of_fold_per_class_recall": dict(zip(classes, recall.tolist())),
        "out_of_fold_confusion_matrix": confusion_matrix(y, predictions, labels=classes).tolist(),
        "confusion_matrix_orientation": "Recorded labels in rows, predictions in columns, class order as listed",
        "out_of_fold_one_vs_rest_auroc": {label: roc_auc_score((y == label).astype(int), probabilities[:, i])
                                         for i, label in enumerate(classes)},
        "test_evaluated": False,
        "note": "Preliminary training CV result, not a tuned model or final test estimate. Fold SD is not a confidence interval. Source batch normalization predates the split.",
    }
    (ROOT / "data" / "logistic_initial_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    # Patient-level outputs stay in ignored storage.
    output = train.assign(predicted_subtype=predictions)
    for i, label in enumerate(classes):
        output[f"probability_{label}"] = probabilities[:, i]
    output.to_csv(processed / "logistic_initial_oof_predictions.csv", index=False)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
