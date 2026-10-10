"""Download official SCAN-B processed expression with checkpoints and provenance."""
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAME = 'GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz'
URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/' + NAME


def download_metadata(directory):
    """Fetch SOFT metadata and check against the published snapshot when available."""
    target = directory / 'GSE96058_family.soft.gz'
    record = ROOT / 'data/scanb_split_report.json'
    expected = json.loads(record.read_text())['source_sha256'] if record.exists() else None
    if not target.exists():
        url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/soft/GSE96058_family.soft.gz'
        partial = target.with_suffix('.partial')
        with urllib.request.urlopen(url, timeout=60) as response, partial.open('wb') as stream:
            while chunk := response.read(1024 * 1024):
                stream.write(chunk)
        with partial.open('rb') as stream:
            if expected and hashlib.file_digest(stream, 'sha256').hexdigest() != expected:
                raise ValueError('SOFT metadata differs from the published snapshot.')
        partial.replace(target)
    with target.open('rb') as stream:
        if expected and hashlib.file_digest(stream, 'sha256').hexdigest() != expected:
            raise ValueError('Existing SOFT metadata failed checksum verification.')


def main():
    directory = ROOT / 'data/raw/scanb'
    directory.mkdir(parents=True, exist_ok=True)
    download_metadata(directory)
    record = ROOT / "data/scanb_download.json"
    recorded = json.loads(record.read_text()) if record.exists() else None
    target = directory / NAME
    partial = directory / (NAME + '.partial')
    if target.exists():
        with target.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        if recorded and (digest != recorded['sha256'] or target.stat().st_size != recorded['bytes']):
            raise ValueError('Existing expression download failed checksum verification.')
        if not recorded:
            raise ValueError('Existing download lacks provenance; retain it separately and download again.')
        print('Verified existing expression download retained.', flush=True)
        return
    offset = partial.stat().st_size if partial.exists() else 0
    request = urllib.request.Request(URL, headers={'Range': f'bytes={offset}-'} if offset else {})
    with urllib.request.urlopen(request, timeout=60) as response:
        resumed = response.status == 206
        if resumed and not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
            raise ValueError("Unexpected resume range; partial file retained.")
        if offset and not resumed:
            offset = 0
        expected = offset + int(response.headers['Content-Length']) if response.headers.get('Content-Length') else None
        total, checkpoint = offset, offset
        with partial.open('ab' if resumed else 'wb') as stream:
            while chunk := response.read(1024 * 1024):
                stream.write(chunk)
                total += len(chunk)
                if total - checkpoint >= 50 * 1024 * 1024:
                    print(f'Downloaded {total / 1024**2:.0f} MiB' + (f' / {expected / 1024**2:.0f} MiB' if expected else ''), flush=True)
                    checkpoint = total
        if expected is not None and total != expected:
            raise ValueError('Incomplete download; partial file retained for resume.')
    with partial.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    if recorded and (digest != recorded['sha256'] or total != recorded['bytes']):
        raise ValueError('Downloaded expression differs from the published snapshot.')
    partial.replace(target)
    if recorded:
        print('Download verified against the published snapshot.', flush=True)
        return
    (ROOT / 'data/scanb_download.json').write_text(json.dumps(dict(url=URL, bytes=total, sha256=digest, downloaded_at_utc=datetime.now(timezone.utc).isoformat(), source='NCBI GEO GSE96058', expression_transform='log2(FPKM + 0.1), per official sample processing metadata'), indent=2)+'\n')
    print('Download complete; SHA256 recorded.', flush=True)


if __name__ == '__main__':
    main()
