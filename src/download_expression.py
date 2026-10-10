"""Download the verified expression file from official cBioPortal DataHub."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

from inspect_dataset import API, ROOT, STUDY_ID


def main():
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / "data/expression_download.json"
    if manifest_path.exists():
        commit = json.loads(manifest_path.read_text())["datahub_commit"]
    else:
        response = requests.get("https://api.github.com/repos/cBioPortal/datahub/commits/master", timeout=60)
        response.raise_for_status()
        commit = response.json()["sha"]
    folder = f"public/{STUDY_ID}"
    source = f"https://raw.githubusercontent.com/cBioPortal/datahub/{commit}/{folder}"
    metadata = requests.get(f"{source}/meta_mrna_seq_v2_rsem.txt", timeout=60)
    metadata.raise_for_status()
    if "datatype: CONTINUOUS" not in metadata.text or "stable_id: rna_seq_v2_mrna" not in metadata.text:
        raise ValueError("Unexpected expression metadata.")
    (raw / "meta_mrna_seq_v2_rsem.txt").write_bytes(metadata.content)
    pointer = requests.get(f"{source}/data_mrna_seq_v2_rsem.txt", timeout=60)
    pointer.raise_for_status()
    if not pointer.text.startswith("version https://git-lfs.github.com/spec/v1"):
        raise ValueError("Expected a Git LFS pointer; investigate DataHub format.")
    expected_hash = pointer.text.split("oid sha256:")[1].splitlines()[0].strip()
    expected_size = int(pointer.text.split("size ")[1].strip())
    url = f"https://media.githubusercontent.com/media/cBioPortal/datahub/{commit}/{folder}/data_mrna_seq_v2_rsem.txt"
    target = raw / "data_mrna_seq_v2_rsem.txt"
    # Reuse only a complete, verified download.
    if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() == expected_hash:
        print("Using verified existing expression download.", flush=True)
    else:
        partial = target.with_suffix(".download")
        print(f"Downloading {expected_size / 1_000_000:.1f} MB of expression data...", flush=True)
        with requests.get(url, stream=True, timeout=(30, 120)) as download:
            download.raise_for_status()
            with partial.open("wb") as output:
                for chunk in download.iter_content(chunk_size=1024 * 1024):
                    output.write(chunk)
        if partial.stat().st_size != expected_size or hashlib.sha256(partial.read_bytes()).hexdigest() != expected_hash:
            raise ValueError("Downloaded expression file failed size or checksum verification.")
        partial.replace(target)

    # Use the portal's explicit mapping; never infer patient IDs by truncating sample IDs.
    response = requests.get(f"{API}/studies/{STUDY_ID}/samples", params={"projection": "DETAILED", "pageSize": 10000}, timeout=60)
    response.raise_for_status()
    samples = response.json()
    (raw / "samples.json").write_text(json.dumps(samples, indent=2), encoding="utf-8")
    manifest = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": STUDY_ID, "datahub_commit": commit,
        "expression_url": url, "expression_sha256": expected_hash,
        "expression_bytes": expected_size,
        "sample_mapping_url": response.url,
        "sample_mapping_sha256": hashlib.sha256((raw / "samples.json").read_bytes()).hexdigest(),
        "profile": "rna_seq_v2_mrna", "source_processing": "Batch-normalized RSEM; not raw sequencing counts",
    }
    (ROOT / "data" / "expression_download.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("Download and checksum verification complete.", flush=True)


if __name__ == "__main__":
    main()
