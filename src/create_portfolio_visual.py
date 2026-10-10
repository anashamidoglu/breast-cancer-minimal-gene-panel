"""Plot public-facing comparisons from the frozen SCAN-B results."""
import json
import matplotlib.pyplot as plt
from tune_models import ROOT


def main():
    r = json.loads((ROOT / "data/scanb_holdout_report.json").read_text())
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9, 5.2))
    strategies = [
        ("20 genes for everyone", 20, r["models"]["small"]["balanced_accuracy"], "#64748B"),
        ("Random escalation", r["average_genes"], r["random_escalation_balanced_accuracy_mean"], "#D97706"),
        ("Escalate uncertain cases", r["average_genes"], r["models"]["staged"]["balanced_accuracy"], "#0F766E"),
        ("50 genes for everyone", 50, r["models"]["large"]["balanced_accuracy"], "#2563EB"),
    ]
    for label, genes, score, color in strategies:
        ax.scatter(genes, score * 100, s=105, color=color, zorder=3)
        offset = (12, -8) if genes < 40 else (-12, -8)
        ax.annotate(f"{label}\n{score:.1%}", (genes, score * 100), xytext=offset, textcoords="offset points", ha="left" if genes < 40 else "right", fontsize=10)
    ax.set(xlim=(15, 55), ylim=(85, 96), xlabel="Average distinct genes measured per patient", ylabel="Balanced accuracy (%)", title="Fewer measurements, a small accuracy tradeoff")
    ax.grid(axis="y", color="#E2E8F0"); ax.set_axisbelow(True)
    fig.text(.13, .02, "916 reserved SCAN-B patients. Y-axis starts at 85% to show differences; points are evaluated strategies.", fontsize=9, color="#475569")
    fig.tight_layout(rect=(0, .06, 1, 1))
    for extension in ["png", "svg"]:
        fig.savefig(ROOT / f"figures/measurement_tradeoff.{extension}", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    rates = r["escalation_by_subtype"]
    labels = ["Basal-like", "Luminal A", "Luminal B", "HER2-enriched"]
    values = [rates[label] * 100 for label in labels]
    bars = ax.barh(labels, values, color="#0F766E", height=.55)
    ax.invert_yaxis()
    ax.bar_label(bars, labels=[f"{v:.1f}%" for v in values], padding=8)
    ax.set(xlim=(0, 40), xlabel="Patients escalated to the 50-gene panel (%)", title="Which subtypes needed more measurements?")
    ax.grid(axis="x", color="#E2E8F0"); ax.set_axisbelow(True)
    fig.text(.13, .02, "Escalation uses model uncertainty. Recorded subtypes are used only to summarize results afterwards.", fontsize=9, color="#475569")
    fig.tight_layout(rect=(0, .06, 1, 1))
    for extension in ["png", "svg"]:
        fig.savefig(ROOT / f"figures/escalation_by_subtype.{extension}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
