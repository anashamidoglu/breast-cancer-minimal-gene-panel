"""Read public study metadata and subtype labels without downloading expression data."""

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import requests

STUDY_ID = "brca_tcga_pan_can_atlas_2018"
API = "https://www.cbioportal.org/api"
ROOT = Path(__file__).resolve().parents[1]


def get_records(endpoint, **params):
    """Request a public table and stop if the server reports an error."""
    response = requests.get(f"{API}/{endpoint}", params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def main():
    study = get_records(f"studies/{STUDY_ID}")
    profiles = get_records(f"studies/{STUDY_ID}/molecular-profiles")
    attributes = get_records(f"studies/{STUDY_ID}/clinical-attributes")
    # Verify the column before requesting its values: do not guess its meaning.
    subtype_attributes = [a for a in attributes if a["clinicalAttributeId"] == "SUBTYPE"]
    if len(subtype_attributes) != 1 or not subtype_attributes[0]["patientAttribute"]:
        raise ValueError("Expected patient-level SUBTYPE attribute; inspect study metadata.")
    labels = get_records(
        f"studies/{STUDY_ID}/clinical-data",
        clinicalDataType="PATIENT", attributeId="SUBTYPE", projection="DETAILED",
    )
    patients = get_records(f"studies/{STUDY_ID}/patients", pageSize=10000)
    counts = Counter(row["value"] for row in labels)
    missing = len({p["patientId"] for p in patients} - {r["patientId"] for r in labels})

    # Save untouched responses locally for provenance; these files are ignored by Git.
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    for name, records in [("study", study), ("molecular_profiles", profiles),
                          ("clinical_attributes", attributes), ("subtype_labels", labels),
                          ("patients", patients)]:
        (raw / f"{name}.json").write_text(json.dumps(records, indent=2), encoding="utf-8")

    # Only aggregate counts and public metadata go into the tracked report.
    report = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": STUDY_ID,
        "study_name": study["name"],
        "study_sample_count": study["allSampleCount"],
        "patient_count": len(patients),
        "label_column": "SUBTYPE",
        "label_level": "patient",
        "subtype_counts_before_expression_matching": dict(sorted(counts.items())),
        "patients_without_subtype_record": missing,
        "expression_profiles": [p for p in profiles if p["molecularAlterationType"] == "MRNA_EXPRESSION"],
        "api_base": API,
        "note": "Counts are before matching labels to expression samples. No exclusions applied.",
    }
    (ROOT / "data" / "dataset_inspection.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
