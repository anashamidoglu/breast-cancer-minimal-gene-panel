"""Evaluate a majority-class baseline with five folds of training data only."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score, recall_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
SEED = 42


def main():
    processed = ROOT / "data" / "processed"
    split = json.loads((ROOT / "data" / "split_report.json").read_text())
    path = processed / "train_labels.csv"
    if hashlib.sha256(path.read_bytes()).hexdigest() != split["membership_sha256"][path.name]:
        raise ValueError("Training membership changed after the split was frozen.")
    train = pd.read_csv(path)
    y = train["subtype"].to_numpy()
    classes = np.unique(y)
    # The dummy ignores gene values, so placeholder inputs suffice. No expression
    # file or test labels are opened by this script.
    x = np.zeros((len(y), 1))
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    predictions = np.empty(len(y), dtype=object)
    probabilities = np.zeros((len(y), len(classes)))
    fold_assignment = np.zeros(len(y), dtype=int)
    folds = []
    for number, (fit_rows, validation_rows) in enumerate(cv.split(x, y), start=1):
        model = DummyClassifier(strategy="most_frequent")
        model.fit(x[fit_rows], y[fit_rows])
        prediction = model.predict(x[validation_rows])
        probability = model.predict_proba(x[validation_rows])
        if not np.array_equal(model.classes_, classes):
            raise ValueError("Fold class ordering changed.")
        predictions[validation_rows] = prediction
        probabilities[validation_rows] = probability
        fold_assignment[validation_rows] = number
        folds.append({
            "fold": number, "fit_samples": len(fit_rows), "validation_samples": len(validation_rows),
            "predicted_subtype": str(prediction[0]),
            "balanced_accuracy": balanced_accuracy_score(y[validation_rows], prediction),
            "macro_f1": f1_score(y[validation_rows], prediction, average="macro", zero_division=0),
        })
    fold_table = pd.DataFrame(folds)
    fold_table.to_csv(ROOT / "data" / "baseline_cv_folds.csv", index=False)
    # Persist training-only fold membership for subsequent model comparisons.
    train.assign(cv_fold=fold_assignment).to_csv(processed / "training_cv_folds.csv", index=False)
    recall = recall_score(y, predictions, labels=classes, average=None, zero_division=0)
    auc = {label: roc_auc_score((y == label).astype(int), probabilities[:, column])
           for column, label in enumerate(classes)}
    report = {
        "model": "DummyClassifier(strategy='most_frequent')",
        "evaluation": "Five-fold stratified cross-validation on the frozen training set only",
        "random_seed": SEED, "training_samples": len(y), "classes": classes.tolist(),
        "balanced_accuracy_cv_mean": float(fold_table["balanced_accuracy"].mean()),
        "balanced_accuracy_cv_sd": float(fold_table["balanced_accuracy"].std(ddof=1)),
        "macro_f1_cv_mean": float(fold_table["macro_f1"].mean()),
        "macro_f1_cv_sd": float(fold_table["macro_f1"].std(ddof=1)),
        "out_of_fold_per_class_recall": dict(zip(classes, recall.tolist())),
        "out_of_fold_confusion_matrix": confusion_matrix(y, predictions, labels=classes).tolist(),
        "confusion_matrix_orientation": "Rows are recorded labels; columns are predicted labels; class order as above.",
        "out_of_fold_one_vs_rest_auroc": auc,
        "test_evaluated": False,
        "explanation": "Each validation sample is predicted by a dummy fitted without that sample. Fold SD is descriptive variation, not a confidence interval.",
    }
    (ROOT / "data" / "baseline_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
