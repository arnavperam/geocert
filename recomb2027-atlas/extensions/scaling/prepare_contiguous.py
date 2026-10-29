"""Prepare central contiguous squares from the original breast format-test data.

Selection depends only on the raster extent. Original transcript coordinates and
IDs are retained; reference owners are the corresponding cell-mask labels. The
vendor owner is retained separately and every disagreement is counted.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import time

import numpy as np
import pyarrow.parquet as pq
import zarr

ARCHIVE_SHA256 = 'cc1e987b06aa748a6b24d3d6f51fc0d6765daa4836c483f14ec4f1bd18b1779b'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if not __debug__:
        raise RuntimeError('Run without Python -O; all source checks are required')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--source-provenance', type=Path, required=True)
    parser.add_argument('--selection-protocol', type=Path, required=True)
    parser.add_argument('--frozen-code-root', type=Path, required=True)
    parser.add_argument('--side-um', type=int, choices=[400, 600], required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    dest = args.output_dir
    if dest.exists() and any(dest.iterdir()):
        raise ValueError('A fresh output directory is required')
    provenance = json.loads(args.source_provenance.read_text())
    assert provenance['sha256'] == ARCHIVE_SHA256
    expected = {row['name']: row['sha256'] for row in provenance['members']}
    names = ['cells.zarr.zip', 'transcripts.parquet', 'cells.parquet', 'gene_panel.json']
    for name in names:
        assert sha(args.source_dir / name) == expected[name], name
    protocol = json.loads(args.selection_protocol.read_text())
    assert protocol['source_archive_sha256'] == ARCHIVE_SHA256
    region = next(r for r in protocol['regions'] if r['nominal_side_um'] == args.side_um)
    store = zarr.ZipStore(str(args.source_dir / 'cells.zarr.zip'), mode='r')
    group = zarr.open_group(store=store, mode='r')
    transform = group['masks/homogeneous_transform'][:].tolist()
    pixel = 1 / Fraction(float(transform[0][0]))
    shape = group['masks/1'].shape
    assert list(shape) == protocol['raster_shape_yx']
    n = math.floor(Fraction(args.side_um) / pixel)
    x0, y0 = (shape[1] - n) // 2, (shape[0] - n) // 2
    assert region['bbox_pixels'] == [x0, y0, x0 + n, y0 + n]
    cell = group['masks/1'][y0:y0+n, x0:x0+n]
    nucleus = group['masks/0'][y0:y0+n, x0:x0+n]
    ids = pq.read_table(args.source_dir / 'cells.parquet', columns=['cell_id']).column(0).to_pylist()
    encoded = group['cell_id'][:]
    decoded = [''.join(chr(ord('a') + int(c, 16)) for c in f'{int(a):08x}') + f'-{int(b)}' for a,b in encoded]
    assert ids == decoded
    store.close()
    label_map = {str(int(k)): ids[int(k)-1] for k in np.unique(cell) if k > 0}
    panel = json.loads((args.source_dir / 'gene_panel.json').read_text())
    genes = {t['type']['data']['name'] for t in panel['payload']['targets'] if t['type']['descriptor'] == 'gene'}
    table = pq.read_table(args.source_dir / 'transcripts.parquet', columns=[
        'transcript_id', 'cell_id', 'feature_name', 'x_location', 'y_location', 'qv'])
    columns = {name: table[name].to_numpy(zero_copy_only=False) for name in table.column_names}
    xp = np.floor(columns['x_location'].astype(np.float64) / float(pixel)).astype(np.int64)
    yp = np.floor(columns['y_location'].astype(np.float64) / float(pixel)).astype(np.int64)
    selected = np.flatnonzero((xp >= x0) & (xp < x0+n) & (yp >= y0) & (yp < y0+n)
                              & (columns['qv'] >= 20) & np.isin(columns['feature_name'], list(genes)))
    roi = region['roi']
    crop = dest / 'data/source_crops' / roi
    crop.mkdir(parents=True)
    (dest / 'results').mkdir()
    (dest / 'docs').mkdir()
    shutil.copyfile(args.frozen_code_root / 'docs/ATLAS_PROTOCOL.json', dest / 'docs/ATLAS_PROTOCOL.json')
    np.savez_compressed(crop / 'masks.npz', cell_mask=cell, nucleus_mask=nucleus)
    counts = defaultdict(Counter)
    records, csv_rows = [], []
    disagreements = 0
    for i in selected:
        x, y = float(columns['x_location'][i]), float(columns['y_location'][i])
        px, py = math.floor(Fraction(x) / pixel)-x0, math.floor(Fraction(y) / pixel)-y0
        assert px == int(xp[i])-x0 and py == int(yp[i])-y0
        label = int(cell[py, px])
        reference = 'UNASSIGNED' if label == 0 else label_map[str(label)]
        vendor = str(columns['cell_id'][i])
        disagreements += reference != vendor
        ident, gene = str(int(columns['transcript_id'][i])), str(columns['feature_name'][i])
        records.append({'id': ident, 'gene': gene, 'cell': label, 'pixel_row': py, 'pixel_col': px,
                        'x_um': repr(x), 'y_um': repr(y)})
        csv_rows.append({'transcript_id': ident, 'cell_id': reference, 'source_cell_id': vendor,
                         'feature_name': gene, 'x_location': repr(x), 'y_location': repr(y),
                         'qv': repr(float(columns['qv'][i])), 'mask_label': label})
        counts[label][gene] += 1
    assert len({r['id'] for r in records}) == len(records)
    with (crop / 'transcripts.csv.gz').open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
                writer.writeheader(); writer.writerows(csv_rows)
    border = set(map(int, np.unique(np.r_[cell[0], cell[-1], cell[:,0], cell[:,-1]])))
    metadata = {'roi_id': roi, 'pixel_size_um': float(pixel), 'bbox_pixels': region['bbox_pixels'],
                'bounds_um': [float(pixel*v) for v in region['bbox_pixels']], 'transform': transform,
                'selection': {'method': protocol['selection'], 'nominal_side_um': args.side_um,
                              'actual_side_um': float(n*pixel), 'protocol_sha256': sha(args.selection_protocol)},
                'q20_gene_rows': len(records), 'source_assignment_disagreements': disagreements,
                'ownership_semantics': 'Reference cell-mask owner; original vendor cell_id retained as source_cell_id',
                'complete_source_labels': [int(k) for k in np.unique(cell) if k > 0 and k not in border],
                'label_to_cell_id': label_map, 'source_archive_sha256': ARCHIVE_SHA256,
                'source_archive_url': provenance['url'], 'source_members_sha256': {name: expected[name] for name in names},
                'scope': '10x breast format-test data; computation and geometry evaluation, not biological validation'}
    (crop / 'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    hashes = {name: sha(crop / name) for name in ['masks.npz', 'metadata.json', 'transcripts.csv.gz']}
    observed = {'scope': metadata['scope'], 'checked_source_records': len(records), 'rois': {
        roi: {'counts': dict(counts), 'records': records, 'source_hashes': hashes}}}
    (dest / 'data/observed_crop_counts.json').write_text(json.dumps(observed, separators=(',', ':'))+'\n')
    summary = {'roi': roi, 'retained_transcripts': len(records), 'source_assignment_disagreements': disagreements,
               'source_hashes': hashes, 'preparation_seconds': time.perf_counter()-started,
               'selection_protocol_sha256': sha(args.selection_protocol)}
    (dest / 'preparation_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
