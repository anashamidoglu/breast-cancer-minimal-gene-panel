"""Summarize nested panel CV, apply the fixed tolerance, and plot PAM50 benchmarks."""

import json
from itertools import combinations

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from tune_models import ROOT
from evaluate_panels import PANEL_SIZES

MODELS = {"logistic": ("L1 logistic", "#4477AA"),
          "random_forest": ("Random forest", "#228833"),
          "xgboost": ("XGBoost", "#AA3377")}


def main():
    protocol = json.loads((ROOT / "data" / "panel_protocol.json").read_text())
    tolerance = protocol["maximum_absolute_drop"]
    summaries, choices, stability_rows = [], [], []
    fig, ax = plt.subplots(figsize=(11, 6.5))
    full_count = None
    for name, (label, color) in MODELS.items():
        reports = json.loads((ROOT / "data" / f"{name}_panel_reports.json").read_text())
        if sorted(r["panel_size"] for r in reports) != PANEL_SIZES or any(r["test_evaluated"] for r in reports):
            raise ValueError(f"Incomplete or test-based panel results for {name}.")
        full = json.loads((ROOT / "data" / f"{name}_tuned_report.json").read_text())
        full_count = full["input_gene_features"]
        reference = full["balanced_accuracy_cv_mean"]
        frame = pd.DataFrame([{"model": name, "panel_size": r["panel_size"],
                               "balanced_accuracy_mean": r["balanced_accuracy_cv_mean"],
                               "balanced_accuracy_fold_sd": r["balanced_accuracy_cv_sd"],
                               "macro_f1_mean": r["macro_f1_cv_mean"], "gene_set": "Fold-selected ANOVA"}
                              for r in reports]).sort_values("panel_size")
        frame["drop_from_full_reference"] = reference - frame["balanced_accuracy_mean"]
        frame["within_tolerance"] = frame["drop_from_full_reference"] <= tolerance
        summaries.extend(frame.to_dict("records"))
        # Include the same model's previously evaluated full-feature point.
        x = np.append(frame["panel_size"].to_numpy(), full_count)
        mean = np.append(frame["balanced_accuracy_mean"].to_numpy(), reference)
        sd = np.append(frame["balanced_accuracy_fold_sd"].to_numpy(), full["balanced_accuracy_cv_sd"])
        ax.plot(x, mean, color=color, marker="o", markersize=5, linewidth=2, label=label)
        ax.fill_between(x, np.clip(mean - sd, 0, 1), np.clip(mean + sd, 0, 1), color=color, alpha=0.10)
        ax.axhline(reference, color=color, linestyle="--", linewidth=1, alpha=0.7)
        eligible = frame.loc[frame["within_tolerance"]]
        k = None if eligible.empty else int(eligible.iloc[0]["panel_size"])
        selected_score = None if k is None else float(eligible.iloc[0]["balanced_accuracy_mean"])
        if k is not None:
            ax.scatter([k], [selected_score], marker="*", s=190, color=color, edgecolor="white", zorder=6)
        pam = json.loads((ROOT / "data" / f"{name}_pam50_reports.json").read_text())
        if len(pam) != 1 or pam[0]["panel_size"] != 50 or pam[0]["test_evaluated"]:
            raise ValueError("Expected one training-only fixed PAM50 benchmark.")
        pam_score = pam[0]["balanced_accuracy_cv_mean"]
        ax.scatter([50], [pam_score], marker="s", s=75, facecolor="white", edgecolor=color, linewidth=2, zorder=7)
        summaries.append({"model": name, "panel_size": 50, "balanced_accuracy_mean": pam_score,
                          "balanced_accuracy_fold_sd": pam[0]["balanced_accuracy_cv_sd"],
                          "macro_f1_mean": pam[0]["macro_f1_cv_mean"], "gene_set": "Fixed published PAM50"})
        choices.append({"model": name, "full_gene_reference": reference,
                        "threshold": reference - tolerance, "smallest_evaluated_panel_within_tolerance": k,
                        "selected_panel_cv_mean": selected_score, "pam50_gene_set_cv_mean": pam_score,
                        "selected_panel_family": "Automatically selected ANOVA panels only; PAM50 benchmark considered separately",
                        "pam50_within_full_gene_tolerance": bool(pam_score >= reference - tolerance),
                        "interpretation": "Fold-specific panels and nested training validation; not one fixed final panel or clinical equivalence."})
        # Stability describes repeated selection, without defining a final fixed panel.
        selections = pd.read_csv(ROOT / "data" / f"{name}_panel_gene_selections.csv")
        if k is not None:
            selected = selections.loc[selections["panel_size"] == k]
            counts = selected.groupby("feature_id")["fold"].nunique()
            gene_sets = [set(selected.loc[selected["fold"] == fold, "feature_id"]) for fold in range(1, 6)]
            overlap = [len(a & b) / len(a | b) for a, b in combinations(gene_sets, 2)]
            choices[-1]["selected_in_all_five_outer_folds"] = int((counts == 5).sum())
            choices[-1]["mean_pairwise_selection_jaccard"] = float(np.mean(overlap))
            for feature_id, count in counts.items():
                symbol = selected.loc[selected["feature_id"] == feature_id, "gene_symbol"].iloc[0]
                stability_rows.append({"model": name, "panel_size": k, "feature_id": feature_id,
                                       "gene_symbol": None if pd.isna(symbol) else symbol,
                                       "selected_in_outer_folds": int(count), "selection_fraction": float(count / 5)})
    ax.axhline(0.25, color="#666666", linestyle=":", linewidth=1.5)
    ax.set(xscale="log", ylim=(0, 1), xlabel="Number of gene features (log scale)",
           ylabel="Balanced accuracy", title="How few genes retain breast cancer subtype performance?")
    ax.set_xticks(PANEL_SIZES + [full_count], [str(k) for k in PANEL_SIZES] + [f"All\n{full_count:,}"])
    ax.grid(axis="y", alpha=0.15)
    handles, labels = ax.get_legend_handles_labels()
    handles += [Line2D([0], [0], color="#666666", linestyle="--", label="Full-gene references (model colors)"),
                Line2D([0], [0], color="#666666", linestyle=":", label="Dummy baseline"),
                Line2D([0], [0], marker="*", markersize=12, color="#666666", linestyle="None", label="Smallest ANOVA panel within 2 points"),
                Line2D([0], [0], marker="s", markersize=7, markerfacecolor="white", color="#666666", linestyle="None", label="Fixed PAM50 gene-set benchmark")]
    ax.legend(handles=handles, loc="lower right", fontsize=9, framealpha=0.95)
    fig.text(0.5, 0.02, "Five outer training folds; selection and tuning inside fitting folds. Shading: mean ± fold SD, not confidence intervals.\n"
             "Squares use our models on the PAM50 gene list; they are not the original PAM50 classifier. Test set unused.\n"
             "Smallest qualifying automatic panels: L1 logistic 100 genes; random forest 100 genes; XGBoost 500 genes.",
             ha="center", va="bottom", fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(ROOT / "figures" / "panel_size_curve.png", dpi=180)
    fig.savefig(ROOT / "figures" / "panel_size_curve.svg")
    plt.close(fig)
    pd.DataFrame(summaries).to_csv(ROOT / "data" / "panel_size_summary.csv", index=False)
    pd.DataFrame(stability_rows).to_csv(ROOT / "data" / "panel_selection_stability.csv", index=False)
    report = {"absolute_tolerance": tolerance, "selection_method": "ANOVA F-score within each inner/outer fit",
              "selection": choices, "test_evaluated": False,
              "note": "Stars are the prespecified practical elbow markers, not geometric knees or statistical noninferiority claims. Panel size was selected using training CV; independent final test evaluation remains pending."}
    (ROOT / "data" / "panel_selection_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
