"""Create a separate analysis workspace without overwriting published outputs."""
import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', default='reproductions/fresh-run')
    args = parser.parse_args()
    destination = Path(args.destination).resolve()
    if any(destination.is_relative_to(ROOT / name) for name in ['src', 'reports', 'data', 'figures']):
        raise ValueError('Choose a destination outside the source and recorded-output directories.')
    # Existing destinations are never merged, reset, or removed.
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copytree(ROOT / 'src', destination / 'src', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    shutil.copytree(ROOT / 'reports', destination / 'reports')
    shutil.copy2(ROOT / 'requirements.txt', destination / 'requirements.txt')
    (destination / 'data').mkdir()
    (destination / 'figures').mkdir()
    # Source pins and the published gene list are inputs, not fitted outputs.
    for name in ['pam50_gene_mapping.csv', 'pam50_reference.json',
                 'expression_download.json', 'scanb_download.json', 'scanb_split_report.json']:
        shutil.copy2(ROOT / 'data' / name, destination / 'data' / name)
    print(f'Fresh workspace: {destination}')
    print('Run analysis commands from this directory; raw data will be downloaded again.')


if __name__ == '__main__':
    main()
