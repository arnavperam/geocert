"""Download and verify the public 10x breast two-FOV format-test source."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import urllib.request
import zipfile

URL = 'https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip'
BYTES = 379018013
SHA256 = 'cc1e987b06aa748a6b24d3d6f51fc0d6765daa4836c483f14ec4f1bd18b1779b'
MEMBERS = ['cells.zarr.zip', 'transcripts.parquet', 'cells.parquet', 'experiment.xenium',
           'gene_panel.json', 'cell_boundaries.parquet', 'nucleus_boundaries.parquet']


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    if not __debug__:
        raise RuntimeError('Run without Python -O; all source integrity checks are required')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = args.output_dir
    root.mkdir(parents=True, exist_ok=True)
    archive = root / URL.rsplit('/', 1)[-1]
    elapsed = None
    if not archive.exists():
        partial = archive.with_suffix('.zip.partial')
        request = urllib.request.Request(URL, headers={'User-Agent': 'Mozilla/5.0', 'Range': f'bytes=0-{BYTES-1}'})
        started = time.perf_counter()
        total = 0
        with urllib.request.urlopen(request, timeout=90) as response, partial.open('wb') as output:
            if response.status == 206:
                assert response.headers.get('Content-Range') == f'bytes 0-{BYTES-1}/{BYTES}'
            elif response.status != 200:
                raise ValueError(f'Unexpected response status: {response.status}')
            while block := response.read(8 * 1024 * 1024):
                total += len(block)
                if total > BYTES:
                    raise ValueError('Response exceeds the fixed source archive size')
                output.write(block)
        assert total == BYTES and digest(partial) == SHA256
        partial.rename(archive)
        elapsed = time.perf_counter() - started
    assert archive.stat().st_size == BYTES and digest(archive) == SHA256
    source = root / 'breast_source'
    source.mkdir(exist_ok=True)
    provenance = {'url': URL, 'bytes': BYTES, 'sha256': SHA256, 'download_seconds': elapsed,
                  'source': '10x Genomics Xenium 2.0.0 breast 2-FOV format-test dataset', 'members': []}
    with zipfile.ZipFile(archive) as bundle:
        for name in MEMBERS:
            bundle.extract(name, source)
            path = source / name
            provenance['members'].append({'name': name, 'bytes': path.stat().st_size,
                                           'sha256': digest(path), 'zip_crc32': bundle.getinfo(name).CRC})
    (root / 'breast_archive_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps(provenance), flush=True)


if __name__ == '__main__':
    main()
