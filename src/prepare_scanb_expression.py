"""Extract unambiguous PAM50 rows from source-transformed SCAN-B expression."""
import csv
import gzip
import hashlib
import json

import numpy as np
import pandas as pd

from tune_models import ROOT


def main():
    source = ROOT / 'data/raw/scanb/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz'
    reference = pd.read_csv(ROOT / 'data/pam50_gene_mapping.csv')
    aliases = {'ORC6L': 'ORC6'}
    candidates = {row.feature_id: set([row.Hugo_Symbol, row.pam50_source_symbol, aliases.get(row.Hugo_Symbol, row.Hugo_Symbol)]) for row in reference.itertuples()}
    wanted = set.union(*candidates.values())
    rows = {}
    with gzip.open(source, 'rt') as stream:
        header = next(csv.reader([next(stream)]))
        assert len(header) == 3410 and len(set(header[1:])) == 3409
        for i, line in enumerate(stream):
            symbol = next(csv.reader([line.split(',', 1)[0]]))[0]
            if symbol in wanted:
                if symbol in rows:
                    raise ValueError('Duplicate required gene symbol: ' + symbol)
                values = next(csv.reader([line]))
                assert len(values) == len(header)
                rows[symbol] = np.array(values[1:], dtype=float)
            if i % 5000 == 0:
                print(f'Checked {i} gene rows; found {len(rows)} candidate PAM50 symbols.', flush=True)
    selected, mapping = [], []
    for row in reference.itertuples():
        matches = candidates[row.feature_id] & set(rows)
        if len(matches) != 1:
            raise ValueError(f'Unresolved mapping {row.feature_id}: {matches}')
        symbol = next(iter(matches))
        selected.append(rows[symbol])
        mapping.append(dict(feature_id=row.feature_id, pam50_symbol=row.pam50_source_symbol, tcga_symbol=row.Hugo_Symbol, scanb_symbol=symbol))
    matrix = pd.DataFrame(np.array(selected).T, index=pd.Index(header[1:], name='sampleId'), columns=reference.feature_id.tolist())
    assert np.isfinite(matrix.to_numpy()).all() and matrix.shape == (3409, 50)
    directory = ROOT / 'data/processed/scanb'
    primary = pd.read_csv(directory / 'metadata.csv')
    primary = primary[~primary.technical_replicate]
    matrix = matrix.loc[primary.sampleId]
    matrix.to_parquet(directory / 'pam50_expression.parquet')
    pd.DataFrame(mapping).to_csv(ROOT / 'data/scanb_pam50_mapping.csv', index=False)
    report = dict(required_genes=50, mapped_genes=50, primary_samples=len(matrix), source_rows=i+1, expression_transform='Source log2(FPKM + 0.1); no further log transform.', aliases=aliases, source_sha256=json.loads((ROOT/'data/scanb_download.json').read_text())['sha256'], finite=True, model_performance_inspected=False)
    (ROOT/'data/scanb_expression_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print('All 50 genes mapped uniquely; 3,273 primary expression profiles prepared.', flush=True)


if __name__ == '__main__':
    main()
