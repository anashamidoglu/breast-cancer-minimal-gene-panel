"""Plot aggregate training-CV results without opening patient-level data."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["initial", "tuned"], default="initial")
    args = parser.parse_args()
    rows = []
    for prefix, name in [("baseline", "Dummy"), ("logistic", "Logistic regression"),
                         ("random_forest", "Random forest"), ("xgboost", "XGBoost")]:
        suffix = "_report.json" if prefix == "baseline" else f"_{args.stage}_report.json"
        report = json.loads((ROOT / "data" / f"{prefix}{suffix}").read_text())
        if report["test_evaluated"]:
            raise ValueError("Expected training-only comparison.")
        rows.append({"model": name, "balanced_accuracy_mean": report["balanced_accuracy_cv_mean"],
                     "balanced_accuracy_fold_sd": report["balanced_accuracy_cv_sd"],
                     "macro_f1_mean": report["macro_f1_cv_mean"]})
    results = pd.DataFrame(rows)
    results.to_csv(ROOT / "data" / f"{args.stage}_model_comparison.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(results["model"], results["balanced_accuracy_mean"],
                  yerr=results["balanced_accuracy_fold_sd"], capsize=5,
                  color=["#999999", "#4477AA", "#228833", "#CCBB44"])
    for bar, score in zip(bars, results["balanced_accuracy_mean"]):
        ax.text(bar.get_x() + bar.get_width() / 2, score - 0.07, f"{score:.1%}", ha="center")
    title = "Initial model comparison: training cross-validation" if args.stage == "initial" else "Tuned model comparison: nested training cross-validation"
    caption = "Fixed settings; five identical training folds." if args.stage == "initial" else "Three inner folds tune settings; five saved outer training folds evaluate them."
    ax.set(ylim=(0, 1), ylabel="Balanced accuracy", title=title)
    fig.text(0.5, 0.015, caption + "\nError bars: fold SD, not confidence intervals.",
             ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(ROOT / "figures" / f"{args.stage}_model_comparison.png", dpi=160)
    plt.close(fig)
    print(results.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
