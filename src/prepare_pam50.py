"""Extract a version-pinned published PAM50 gene list and map by Entrez ID."""

import hashlib
import json

import pandas as pd
import rdata
import requests

from inspect_dataset import ROOT


def main():
    raw = ROOT / "data" / "raw"
    response = requests.get("https://api.github.com/repos/bhklab/genefu", timeout=60)
    response.raise_for_status()
    branch = response.json()["default_branch"]
    response = requests.get(f"https://api.github.com/repos/bhklab/genefu/commits/{branch}", timeout=60)
    response.raise_for_status()
    commit = response.json()["sha"]
    url = f"https://raw.githubusercontent.com/bhklab/genefu/{commit}/data/pam50.rda"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    path = raw / "pam50_reference.rda"
    path.write_bytes(response.content)
    reference = rdata.read_rda(path)["pam50"]
    mapping = reference["centroids.map"].copy()
    if len(mapping) != 50 or mapping["EntrezGene.ID"].duplicated().any():
        raise ValueError("Expected 50 unique published PAM50 gene IDs.")
    genes = pd.DataFrame({"pam50_source_symbol": mapping["probe.centroids"].astype(str),
                          "Entrez_Gene_Id": mapping["EntrezGene.ID"].astype(int)})
    annotations = pd.read_csv(ROOT / "data" / "processed" / "gene_annotations.csv")
    matched = genes.merge(annotations, on="Entrez_Gene_Id", how="left", validate="one_to_one")
    if matched["feature_id"].isna().any() or not matched["feature_id"].is_unique:
        raise ValueError("Some PAM50 genes are unavailable or ambiguous; investigate before benchmarking.")
    matched.to_csv(ROOT / "data" / "pam50_gene_mapping.csv", index=False)
    report = {
        "reference": "genefu pam50 centroids.map; gene list only, not centroid classification",
        "repository_commit": commit, "source_url": url,
        "source_sha256": hashlib.sha256(response.content).hexdigest(),
        "original_paper": "https://pubmed.ncbi.nlm.nih.gov/19204204/",
        "package_documentation": "https://www.bioconductor.org/packages/release/bioc/manuals/genefu/man/genefu.pdf",
        "published_genes": 50, "matched_unique_features": len(matched),
        "mapping_rule": "Use published Entrez IDs, preserving original and source-expression symbols; no data-driven selection.",
        "benchmark": "Our three models fitted on the fixed PAM50 gene set using the same nested training validation as selected panels.",
        "limitation": "Not the original PAM50 algorithm or Prosigna assay; agreement with PAM50-derived target labels is not independent diagnostic validation.",
    }
    (ROOT / "data" / "pam50_reference.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
