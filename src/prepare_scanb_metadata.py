"""Freeze patient-disjoint SCAN-B replication splits before expression modeling."""
import gzip
import hashlib
import json

import pandas as pd
from sklearn.model_selection import train_test_split

from tune_models import ROOT


def main():
    source = ROOT / 'data/raw/scanb/GSE96058_family.soft.gz'
    text = gzip.decompress(source.read_bytes()).decode()
    rows = []
    for block in text.split('^SAMPLE = ')[1:]:
        row = {'geo_accession': block.splitlines()[0].strip()}
        for line in block.splitlines():
            if line.startswith('!Sample_title = '):
                row['sampleId'] = line.split(' = ', 1)[1]
            if line.startswith('!Sample_characteristics_ch1 = '):
                key, value = line.split(' = ', 1)[1].split(': ', 1)
                row[key] = value
            if line.startswith('!Sample_relation = BioSample: '):
                row['biosample'] = line.rsplit('/', 1)[1]
        rows.append(row)
    metadata = pd.DataFrame(rows)
    assert len(metadata) == 3409 and metadata.sampleId.is_unique
    metadata['patientId'] = metadata.sampleId.str.removesuffix('repl')
    metadata['technical_replicate'] = metadata.sampleId.str.endswith('repl')
    primary = metadata[~metadata.technical_replicate].copy()
    assert len(primary) == 3273 and primary.patientId.is_unique and primary.biosample.is_unique
    for row in metadata[metadata.technical_replicate].itertuples():
        original = primary[primary.sampleId == row.sampleId.removesuffix('repl')]
        assert len(original) == 1 and original.iloc[0].patientId == row.patientId
    mapping = {'LumA': 'Luminal A', 'LumB': 'Luminal B', 'Basal': 'Basal-like', 'Her2': 'HER2-enriched'}
    assert set(primary['pam50 subtype']) == set(mapping) | {'Normal'}
    primary['subtype'] = primary['pam50 subtype'].map(mapping)
    retained = primary[primary.subtype.notna()][['sampleId', 'patientId', 'subtype']].sort_values('sampleId').reset_index(drop=True)
    assert len(retained) == 3052
    development, holdout = train_test_split(retained, test_size=.30, stratify=retained.subtype, random_state=20261010)
    directory = ROOT / 'data/processed/scanb'
    directory.mkdir(exist_ok=True)
    checksums = {}
    for name, frame in [('metadata.csv', metadata), ('development_labels.csv', development.sort_values('sampleId')), ('holdout_labels.csv', holdout.sort_values('sampleId'))]:
        path = directory / name
        data = frame.to_csv(index=False).encode()
        if path.exists() and path.read_bytes() != data:
            raise RuntimeError('Refusing to change frozen SCAN-B metadata or split.')
        path.write_bytes(data)
        checksums[name] = hashlib.sha256(data).hexdigest()
    assert not set(development.patientId) & set(holdout.patientId)
    report = dict(source='GSE96058 official SOFT metadata', patient_key='Canonical primary sample title; technical replicate suffix removed. BioSample IDs differ for replicate records.', source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), profiles=3409, technical_replicates_excluded=136, primary_cases=3273, normal_like_excluded=221, retained_cases=3052, development_cases=len(development), fresh_holdout_cases=len(holdout), development_counts=development.subtype.value_counts().to_dict(), holdout_counts=holdout.subtype.value_counts().to_dict(), seed=20261010, holdout_fraction=.30, membership_sha256=checksums, model_performance_inspected=False, evaluation_scope='Independent-cohort replication after within-cohort training, not direct TCGA model transport.')
    (ROOT / 'data/scanb_split_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
