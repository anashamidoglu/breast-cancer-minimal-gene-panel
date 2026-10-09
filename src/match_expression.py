"""Match expression samples to approved labels and audit structural cleanup."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    raw = ROOT / "data" / "raw"
    processed = ROOT / "data" / "processed"
    expression = pd.read_csv(raw / "data_mrna_seq_v2_rsem.txt", sep="\t")
    if expression.columns[:2].tolist() != ["Hugo_Symbol", "Entrez_Gene_Id"]:
        raise ValueError("Unexpected gene annotation columns.")
    samples = pd.DataFrame(json.loads((raw / "samples.json").read_text()))
    labels = pd.read_csv(processed / "patient_labels.csv")
    if samples["sampleId"].duplicated().any() or samples["patientId"].duplicated().any():
        raise ValueError("Repeated samples or patients need an explicit resolution rule.")
    if set(samples["sampleType"]) != {"Primary Solid Tumor"}:
        raise ValueError("Unexpected sample types require review.")
    expression_samples = expression.columns[2:].tolist()
    if len(set(expression_samples)) != len(expression_samples):
        raise ValueError("Duplicate expression sample identifiers.")
    unknown = set(expression_samples) - set(samples["sampleId"])
    if unknown:
        raise ValueError(f"Expression samples absent from portal mapping: {len(unknown)}")
    joined = samples[["sampleId", "patientId"]].merge(labels, on="patientId", validate="one_to_one")
    matched = joined.loc[joined["sampleId"].isin(expression_samples)].sort_values("sampleId")
    missing = joined.loc[~joined["sampleId"].isin(expression_samples)].copy()
    missing["reason"] = "No expression column in source matrix"
    missing.to_csv(processed / "expression_exclusions.csv", index=False)

    # Repeated Entrez IDs have distinct measurements, not exact duplicate rows.
    # Exclude every ambiguous ID rather than arbitrarily choosing or averaging rows.
    # This rule uses identifiers only, never subtype or expression statistics.
    ambiguous = expression["Entrez_Gene_Id"].duplicated(keep=False)
    expression.loc[ambiguous, ["Hugo_Symbol", "Entrez_Gene_Id"]].to_csv(
        processed / "ambiguous_gene_annotations.csv", index=False
    )
    genes = expression.loc[~ambiguous, ["Hugo_Symbol", "Entrez_Gene_Id"]].copy()
    genes["feature_id"] = "ENTREZ_" + genes["Entrez_Gene_Id"].astype(str)
    genes.to_csv(processed / "gene_annotations.csv", index=False)
    values = expression.loc[~ambiguous, matched["sampleId"].tolist()].to_numpy(dtype=np.float64).T
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Missing, nonfinite, or negative expression values require investigation.")
    # Rows become samples; columns become genes. Keep identifiers alongside labels.
    matrix = pd.DataFrame(values, index=matched["sampleId"], columns=genes["feature_id"])
    matrix.index.name = "sampleId"
    matrix.to_parquet(processed / "expression.parquet")
    matched.to_csv(processed / "sample_labels.csv", index=False)
    report = {
        "source_gene_rows": len(expression),
        "source_expression_samples": len(expression_samples),
        "approved_label_patients": len(labels),
        "matched_samples_and_patients": len(matched),
        "patients_without_expression": len(missing),
        "retained_subtype_counts": matched["subtype"].value_counts().to_dict(),
        "ambiguous_entrez_ids": int(expression.loc[ambiguous, "Entrez_Gene_Id"].nunique()),
        "ambiguous_gene_rows_excluded": int(ambiguous.sum()),
        "gene_features_retained": len(genes),
        "missing_gene_symbols_retained_with_entrez_id": int(genes["Hugo_Symbol"].isna().sum()),
        "missing_nonfinite_or_negative_values": 0,
        "sample_type": "Primary Solid Tumor",
        "repeated_patients": 0,
        "gene_rule": "Exclude all rows with ambiguous repeated Entrez IDs; retain unique IDs even when symbols are absent.",
        "note": "Structural matching only. No transformation, variance filtering, feature selection, splitting, or model fitting.",
    }
    (ROOT / "data" / "expression_matching.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
