"""Create a fixed patient-disjoint split, and refuse to overwrite a different one."""

import hashlib
import json
from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
SEED = 42
TEST_FRACTION = 0.20


def main():
    processed = ROOT / "data" / "processed"
    labels = pd.read_csv(processed / "sample_labels.csv").sort_values("sampleId").reset_index(drop=True)
    if labels["sampleId"].duplicated().any() or labels["patientId"].duplicated().any():
        raise ValueError("Repeated samples/patients need a grouped split.")
    # Stratification preserves the relative sizes of the four subtypes.
    train, test = train_test_split(labels, test_size=TEST_FRACTION, stratify=labels["subtype"], random_state=SEED)
    train = train.sort_values("sampleId")
    test = test.sort_values("sampleId")
    if set(train["patientId"]) & set(test["patientId"]):
        raise ValueError("A patient appears in both sets.")
    if len(train) + len(test) != len(labels) or set(train["subtype"]) != set(test["subtype"]):
        raise ValueError("Incomplete split or absent class.")
    # Save labels and membership, not expression copies. The test expression stays unused.
    for name, table in [("train_labels.csv", train), ("test_labels.csv", test)]:
        path = processed / name
        content = table.to_csv(index=False, lineterminator="\n")
        if path.exists() and path.read_text(encoding="utf-8") != content:
            raise ValueError(f"Refusing to replace frozen {name}. Investigate the data change.")
        path.write_text(content, encoding="utf-8")
    report = {
        "random_seed": SEED, "test_fraction": TEST_FRACTION,
        "training_samples": len(train), "test_samples": len(test),
        "training_subtype_counts": train["subtype"].value_counts().to_dict(),
        "test_subtype_counts": test["subtype"].value_counts().to_dict(),
        "patient_overlap": 0,
        "membership_sha256": {name: hashlib.sha256((processed / name).read_bytes()).hexdigest()
                              for name in ["train_labels.csv", "test_labels.csv"]},
        "policy": "Test membership and aggregate label counts may be checked; test expression and prediction performance remain unused until final evaluation.",
    }
    (ROOT / "data" / "split_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
