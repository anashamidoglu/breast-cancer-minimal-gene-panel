"""Nested training-CV panel curves with gene selection inside every fit."""

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score, recall_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder
from threadpoolctl import threadpool_limits

from tune_models import ROOT, SEED, specification

PANEL_SIZES = [1, 2, 5, 10, 20, 50, 100, 500]


def panel_specification(name, k):
    base, grid = specification(name)
    # ANOVA ranks genes by differences among subtype means relative to variation
    # within subtypes. It is supervised: it must be refitted inside every fold.
    steps = base.steps[:2] + [("genes", SelectKBest(score_func=f_classif, k=k))] + base.steps[2:]
    return Pipeline(steps), grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, choices=["logistic", "random_forest", "xgboost"])
    parser.add_argument("--sizes", nargs="+", type=int, default=PANEL_SIZES)
    parser.add_argument("--pam50", action="store_true", help="Use the fixed published PAM50 gene list instead of selecting genes")
    args = parser.parse_args()
    name = args.model
    prefix = f"{name}_pam50" if args.pam50 else f"{name}_panel"
    if args.pam50:
        args.sizes = [50]
    if not args.sizes or any(k < 1 or k > 500 for k in args.sizes) or len(set(args.sizes)) != len(args.sizes):
        raise ValueError("Use distinct panel sizes between 1 and 500.")
    processed = ROOT / "data" / "processed"
    split = json.loads((ROOT / "data" / "split_report.json").read_text())
    path = processed / "train_labels.csv"
    if hashlib.sha256(path.read_bytes()).hexdigest() != split["membership_sha256"][path.name]:
        raise ValueError("Training membership changed.")
    train = pd.read_csv(path)
    membership = pd.read_csv(processed / "training_cv_folds.csv")
    if not membership[["sampleId", "patientId", "subtype"]].equals(train):
        raise ValueError("Saved folds and labels differ.")
    x = pd.read_parquet(processed / "expression.parquet",
                        filters=[("sampleId", "in", train["sampleId"].tolist())]).loc[train["sampleId"]]
    encoder = LabelEncoder()
    y = encoder.fit_transform(train["subtype"])
    classes = encoder.classes_
    annotations = pd.read_csv(processed / "gene_annotations.csv").set_index("feature_id")
    if args.pam50:
        pam50 = pd.read_csv(ROOT / "data" / "pam50_gene_mapping.csv")
        if len(pam50) != 50 or not pam50["feature_id"].is_unique:
            raise ValueError("PAM50 mapping must have 50 unique gene IDs.")
        x = x.loc[:, pam50["feature_id"].tolist()]
    output_dir = processed / (f"{name}_pam50" if args.pam50 else f"{name}_panels")
    output_dir.mkdir(exist_ok=True)
    results, selections, reports, tuning_rows = [], [], [], []
    for k in sorted(args.sizes):
        predictions = np.zeros(len(y), dtype=int)
        probabilities = np.zeros((len(y), len(classes)))
        size_rows = []
        for fold in sorted(membership["cv_fold"].unique()):
            validation = np.flatnonzero(membership["cv_fold"].to_numpy() == fold)
            fitting = np.flatnonzero(membership["cv_fold"].to_numpy() != fold)
            pipeline, grid = specification(name) if args.pam50 else panel_specification(name, k)
            cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED + int(fold))
            search = GridSearchCV(pipeline, grid, scoring="balanced_accuracy", cv=cv,
                                  n_jobs=1, refit=True, error_score="raise")
            print(f"{prefix}: {k} genes, outer fold {fold}/5...", flush=True)
            with warnings.catch_warnings(), threadpool_limits(limits=2):
                warnings.simplefilter("error", ConvergenceWarning)
                search.fit(x.iloc[fitting], y[fitting])
                prediction = search.predict(x.iloc[validation]).astype(int)
                probability = search.predict_proba(x.iloc[validation])
            if not np.array_equal(search.classes_, np.arange(len(classes))):
                raise ValueError("Unexpected class order.")
            fitted = search.best_estimator_
            survivors = x.columns[fitted["variance"].get_support()]
            chosen = survivors if args.pam50 else survivors[fitted["genes"].get_support()]
            if len(chosen) != k or not chosen.is_unique:
                raise ValueError("Panel does not contain exactly k distinct gene IDs.")
            for feature_id in chosen:
                symbol = annotations.loc[feature_id, "Hugo_Symbol"]
                selections.append({"model": name, "panel_size": k, "fold": int(fold),
                                   "feature_id": feature_id, "gene_symbol": None if pd.isna(symbol) else str(symbol)})
            for params, mean, sd in zip(search.cv_results_["params"], search.cv_results_["mean_test_score"],
                                        search.cv_results_["std_test_score"]):
                tuning_rows.append({"model": name, "panel_size": k, "fold": int(fold), "parameters": params,
                                    "mean_inner_balanced_accuracy": float(mean), "sd_inner_balanced_accuracy": float(sd)})
            predictions[validation] = prediction
            probabilities[validation] = probability
            score = balanced_accuracy_score(y[validation], prediction)
            row = {"model": name, "panel_size": k, "fold": int(fold),
                   "balanced_accuracy": score,
                   "macro_f1": f1_score(y[validation], prediction, average="macro", zero_division=0),
                   "best_parameters": json.dumps(search.best_params_, sort_keys=True)}
            results.append(row)
            size_rows.append(row)
            pd.DataFrame(results).to_csv(ROOT / "data" / f"{prefix}_cv_folds.csv", index=False)
        folds = pd.DataFrame(size_rows)
        recall = recall_score(y, predictions, labels=np.arange(len(classes)), average=None, zero_division=0)
        reports.append({"model": name, "panel_size": k, "classes": classes.tolist(),
                        "gene_set": "Fixed published PAM50 genes" if args.pam50 else "ANOVA selection within every fitting fold",
                        "balanced_accuracy_cv_mean": float(folds["balanced_accuracy"].mean()),
                        "balanced_accuracy_cv_sd": float(folds["balanced_accuracy"].std(ddof=1)),
                        "macro_f1_cv_mean": float(folds["macro_f1"].mean()),
                        "out_of_fold_per_class_recall": dict(zip(classes, recall.tolist())),
                        "out_of_fold_one_vs_rest_auroc": {label: roc_auc_score((y == i).astype(int), probabilities[:, i])
                                                         for i, label in enumerate(classes)},
                        "out_of_fold_confusion_matrix": confusion_matrix(y, predictions).tolist(),
                        "test_evaluated": False})
        output = train.assign(predicted_subtype=encoder.inverse_transform(predictions))
        for i, label in enumerate(classes):
            output[f"probability_{label}"] = probabilities[:, i]
        output.to_csv(output_dir / f"{k}_genes_oof_predictions.csv", index=False)
        pd.DataFrame(selections).to_csv(ROOT / "data" / f"{prefix}_gene_selections.csv", index=False)
        (ROOT / "data" / f"{prefix}_reports.json").write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
        (ROOT / "data" / f"{prefix}_tuning_search.json").write_text(json.dumps(tuning_rows, indent=2) + "\n", encoding="utf-8")
        print(f"{name}: {k} genes mean balanced accuracy = {reports[-1]['balanced_accuracy_cv_mean']:.4f}", flush=True)


if __name__ == "__main__":
    main()
