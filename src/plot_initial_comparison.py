"""Plot aggregate training-CV results without opening patient-level data."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    rows = []
    for prefix, name in [("baseline", "Dummy"), ("logistic", "Logistic regression"),
                         ("random_forest", "Random forest"), ("xgboost", "XGBoost")]:
        suffix = "_report.json" if prefix == "baseline" else "_initial_report.json"
        report = json.loads((ROOT / "data" / f"{prefix}{suffix}").read_text())
        if report["test_evaluated"]:
            raise ValueError("Expected training-only comparison.")
        rows.append({"model": name, "balanced_accuracy_mean": report["balanced_accuracy_cv_mean"],
                     "balanced_accuracy_fold_sd": report["balanced_accuracy_cv_sd"],
                     "macro_f1_mean": report["macro_f1_cv_mean"]})
    results = pd.DataFrame(rows)
    results.to_csv(ROOT / "data" / "initial_model_comparison.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(results["model"], results["balanced_accuracy_mean"],
                  yerr=results["balanced_accuracy_fold_sd"], capsize=5,
                  color=["#999999", "#4477AA", "#228833", "#CCBB44"])
    for bar, score in zip(bars, results["balanced_accuracy_mean"]):
        ax.text(bar.get_x() + bar.get_width() / 2, score - 0.07, f"{score:.1%}", ha="center")
    ax.set(ylim=(0, 1), ylabel="Balanced accuracy", title="Initial model comparison: training cross-validation")
    fig.text(0.5, 0.015, "Fixed settings; five identical training folds. Error bars: fold SD, not confidence intervals.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(ROOT / "figures" / "initial_model_comparison.png", dpi=160)
    plt.close(fig)
    print(results.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
