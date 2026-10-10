"""Download official SCAN-B processed expression with checkpoints and provenance."""
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAME = 'GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz'
URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/' + NAME


def main():
    directory = ROOT / 'data/raw/scanb'
    directory.mkdir(exist_ok=True)
    target = directory / NAME
    partial = directory / (NAME + '.partial')
    if target.exists():
        print('Existing complete download retained; not overwritten.', flush=True)
        return
    offset = partial.stat().st_size if partial.exists() else 0
    request = urllib.request.Request(URL, headers={'Range': f'bytes={offset}-'} if offset else {})
    with urllib.request.urlopen(request, timeout=60) as response:
        resumed = response.status == 206
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
    digest = hashlib.file_digest(partial.open('rb'), 'sha256').hexdigest()
    partial.replace(target)
    (ROOT / 'data/scanb_download.json').write_text(json.dumps(dict(url=URL, bytes=total, sha256=digest, downloaded_at_utc=datetime.now(timezone.utc).isoformat(), source='NCBI GEO GSE96058', expression_transform='log2(FPKM + 0.1), per official sample processing metadata'), indent=2)+'\n')
    print('Download complete; SHA256 recorded.', flush=True)


if __name__ == '__main__':
    main()
