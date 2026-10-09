"""Apply the agreed four-subtype definition; preserve an exclusion audit."""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LABEL_NAMES = {
    "BRCA_LumA": "Luminal A",
    "BRCA_LumB": "Luminal B",
    "BRCA_Her2": "HER2-enriched",
    "BRCA_Basal": "Basal-like",
}


def main():
    raw = ROOT / "data" / "raw"
    patients = pd.DataFrame(json.loads((raw / "patients.json").read_text()))
    labels = pd.DataFrame(json.loads((raw / "subtype_labels.json").read_text()))
    if patients["patientId"].duplicated().any() or labels["patientId"].duplicated().any():
        raise ValueError("Duplicate patient identifiers require investigation.")
    unexpected = set(labels["value"].dropna()) - set(LABEL_NAMES) - {"BRCA_Normal"}
    if unexpected:
        raise ValueError(f"Unexpected subtype values: {unexpected}")

    # A left join keeps every patient, including those without a subtype label.
    table = patients[["patientId"]].merge(
        labels[["patientId", "value"]], on="patientId", how="left", validate="one_to_one"
    )
    table["subtype"] = table["value"].map(LABEL_NAMES)
    retained = table.loc[table["subtype"].notna(), ["patientId", "subtype"]]
    excluded = table.loc[table["subtype"].isna(), ["patientId", "value"]].copy()
    excluded["reason"] = excluded["value"].map({"BRCA_Normal": "Normal-like outside four-subtype scope"})
    excluded.loc[excluded["value"].isna(), "reason"] = "Missing subtype label"

    processed = ROOT / "data" / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    retained.to_csv(processed / "patient_labels.csv", index=False)
    excluded.to_csv(processed / "label_exclusions.csv", index=False)
    report = {
        "decision": "User approved excluding Normal-like and missing subtype labels on October 9, 2026.",
        "patients_before": len(table),
        "patients_retained": len(retained),
        "excluded_by_reason": excluded["reason"].value_counts().to_dict(),
        "retained_subtype_counts": retained["subtype"].value_counts().to_dict(),
        "note": "Patient labels only; expression availability has not yet been matched. No train/test split or model fitting.",
    }
    (ROOT / "data" / "label_preparation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
