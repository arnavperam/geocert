"""Check prepared contiguous crops directly against the original source archive.

This checker imports neither the preparation adapter nor any GeoCert module.
It reconstructs the complete Q20 measured-gene census and exact rational
coordinate-to-pixel map, then compares every prepared observation and mask.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import zipfile

import numpy as np
import pyarrow.parquet as pq
import zarr

ARCHIVE_SHA256 = 'cc1e987b06aa748a6b24d3d6f51fc0d6765daa4836c483f14ec4f1bd18b1779b'
MEMBERS = ('cells.zarr.zip', 'transcripts.parquet', 'cells.parquet',
           'experiment.xenium', 'gene_panel.json', 'cell_boundaries.parquet',
           'nucleus_boundaries.parquet')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def unique_index(rows, key, description):
    result = {}
    for row in rows:
        value = row[key]
        require(value not in result, 'duplicate ' + description)
        result[value] = row
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--selection-protocol', type=Path, required=True)
    parser.add_argument('--prepared-root', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    require(sha(args.archive) == ARCHIVE_SHA256, 'original archive checksum differs')
    source_hashes = {}
    with zipfile.ZipFile(args.archive) as archive:
        for name in MEMBERS:
            content = (args.source_dir / name).read_bytes()
            require(content == archive.read(name), 'extracted member differs: ' + name)
            source_hashes[name] = hashlib.sha256(content).hexdigest()

    plan = json.loads(args.selection_protocol.read_text())
    plan_hash = sha(args.selection_protocol)
    require(plan['source_archive_sha256'] == ARCHIVE_SHA256, 'selection source differs')
    store = zarr.ZipStore(str(args.source_dir / 'cells.zarr.zip'), mode='r')
    try:
        group = zarr.open_group(store=store, mode='r')
        transform = group['masks/homogeneous_transform'][:]
        require(transform.shape == (4, 4), 'unexpected coordinate transform shape')
        require(transform[0, 0] > 0 and transform[0, 0] == transform[1, 1]
                and all(transform[i, j] == 0 for i in (0, 1) for j in range(4) if i != j),
                'source transform is not isotropic and unshifted')
        scale = Fraction(float(transform[0, 0]))
        cell_mask, nuclear_mask = group['masks/1'][:], group['masks/0'][:]
        encoded_ids = group['cell_id'][:]
    finally:
        store.close()
    require(list(cell_mask.shape) == plan['raster_shape_yx'], 'selection raster extent differs')
    cell_ids = pq.read_table(args.source_dir / 'cells.parquet', columns=['cell_id']).column(0).to_pylist()
    decoded = [''.join(chr(97 + int(c, 16)) for c in f'{int(a):08x}') + f'-{int(b)}'
               for a, b in encoded_ids]
    require(cell_ids == decoded, 'cell-label and original cell-ID mapping differs')
    panel = json.loads((args.source_dir / 'gene_panel.json').read_text())
    genes = {entry['type']['data']['name'] for entry in panel['payload']['targets']
             if entry['type']['descriptor'] == 'gene'}
    columns = pq.read_table(args.source_dir / 'transcripts.parquet', columns=[
        'transcript_id', 'cell_id', 'feature_name', 'x_location', 'y_location', 'qv']).to_pydict()
    selected = [i for i, (qv, gene) in enumerate(zip(columns['qv'], columns['feature_name']))
                if qv >= 20 and gene in genes]
    # Convert every qualifying original coordinate exactly, including rows outside
    # the crops. A floating-point prefilter cannot silently discard a boundary row.
    pixels = {i: (math.floor(Fraction(float(columns['x_location'][i])) * scale),
                  math.floor(Fraction(float(columns['y_location'][i])) * scale)) for i in selected}
    print(f'Checked archive, cell-ID map and {len(selected)} Q20 measured-gene rows', flush=True)

    regions = []
    seen_roots = set()
    for root in args.prepared_root:
        root = root.resolve()
        require(root not in seen_roots, 'duplicate prepared root')
        seen_roots.add(root)
        observed_path = root / 'data/observed_crop_counts.json'
        observed = json.loads(observed_path.read_text())
        require(len(observed['rois']) == 1, 'one contiguous crop per prepared root required')
        roi, supplied = next(iter(observed['rois'].items()))
        region = next(item for item in plan['regions'] if item['roi'] == roi)
        crop = root / 'data/source_crops' / roi
        metadata = json.loads((crop / 'metadata.json').read_text())
        require(metadata['source_archive_sha256'] == ARCHIVE_SHA256, 'prepared archive binding differs')
        require(metadata['selection']['protocol_sha256'] == plan_hash, 'prepared selection binding differs')
        require(metadata['transform'] == transform.tolist(), 'prepared transform differs')
        require(metadata['pixel_size_um'] == float(1 / scale), 'prepared pixel size differs')
        require(all(source_hashes[name] == value for name, value in metadata['source_members_sha256'].items()),
                'prepared source-member binding differs')
        side = math.floor(Fraction(region['nominal_side_um']) * scale)
        x0, y0 = (cell_mask.shape[1] - side) // 2, (cell_mask.shape[0] - side) // 2
        bbox = [x0, y0, x0 + side, y0 + side]
        require(region['bbox_pixels'] == metadata['bbox_pixels'] == bbox, 'crop is not the specified central square')
        x1, y1 = bbox[2:]
        window = cell_mask[y0:y1, x0:x1]
        labels = [int(k) for k in np.unique(window) if k > 0]
        require(metadata['label_to_cell_id'] == {str(k): cell_ids[k - 1] for k in labels},
                'prepared cell-label map differs')
        border = set(map(int, np.unique(np.concatenate((window[0], window[-1], window[:, 0], window[:, -1])))))
        require(metadata['complete_source_labels'] == [k for k in labels if k not in border],
                'prepared crop-contained label list differs')
        with np.load(crop / 'masks.npz') as prepared_masks:
            require(np.array_equal(prepared_masks['cell_mask'], cell_mask[y0:y1, x0:x1]), 'cell mask differs')
            require(np.array_equal(prepared_masks['nucleus_mask'], nuclear_mask[y0:y1, x0:x1]), 'nuclear mask differs')
        expected = {}
        for i, (x, y) in pixels.items():
            if x0 <= x < x1 and y0 <= y < y1:
                ident = str(int(columns['transcript_id'][i]))
                require(ident not in expected, 'duplicate retained source transcript ID')
                expected[ident] = i
        records = unique_index(supplied['records'], 'id', 'prepared observation ID')
        with gzip.open(crop / 'transcripts.csv.gz', 'rt', newline='') as stream:
            exported = unique_index(csv.DictReader(stream), 'transcript_id', 'prepared CSV transcript ID')
        require(set(records) == set(exported) == set(expected), 'retained transcript census differs')
        counts = defaultdict(Counter)
        disagreements = 0
        for ident, i in expected.items():
            x, y = pixels[i]
            label, gene = int(cell_mask[y, x]), columns['feature_name'][i]
            source_owner = columns['cell_id'][i]
            raster_owner = 'UNASSIGNED' if label == 0 else cell_ids[label - 1]
            disagreements += source_owner != raster_owner
            record, row = records[ident], exported[ident]
            require(record['cell'] == label and record['gene'] == gene
                    and record['pixel_col'] == x - x0 and record['pixel_row'] == y - y0,
                    'prepared observation gene, owner or pixel differs')
            for dimension in ('x', 'y'):
                original = float(columns[dimension + '_location'][i])
                require(float(record[dimension + '_um']) == original
                        and float(row[dimension + '_location']) == original, 'prepared coordinate differs')
            require(row['cell_id'] == raster_owner and row['source_cell_id'] == source_owner
                    and int(row['mask_label']) == label and row['feature_name'] == gene
                    and float(row['qv']) == float(columns['qv'][i]), 'prepared CSV metadata differs')
            counts[str(label)][gene] += 1
        require({cell: dict(values) for cell, values in counts.items()} == supplied['counts'], 'prepared source counts differ')
        require(disagreements == metadata['source_assignment_disagreements'], 'owner-disagreement count differs')
        require(len(expected) == metadata['q20_gene_rows'], 'metadata transcript census differs')
        files = {'data/observed_crop_counts.json': sha(observed_path)}
        for name in ('masks.npz', 'metadata.json', 'transcripts.csv.gz'):
            actual = sha(crop / name)
            require(supplied['source_hashes'][name] == actual, 'prepared file hash binding differs: ' + name)
            files[f'data/source_crops/{roi}/{name}'] = actual
        regions.append({'roi': roi, 'nominal_side_um': region['nominal_side_um'], 'bbox_pixels': bbox,
                        'retained_transcripts': len(expected), 'source_owner_disagreements': disagreements,
                        'prepared_files_sha256': files, 'status': 'passed'})
        print(f'{roi}: all {len(expected)} retained observations and both masks match', flush=True)
    output = {'status': 'passed', 'scope': 'Original-source identity, complete Q20 measured-gene census, exact rational raster pixels, masks, coordinates, genes, owners and source counts; no geometric optimization is performed.',
              'independence': 'Imports no preparation adapter or GeoCert module.',
              'archive_sha256': ARCHIVE_SHA256, 'archive_bytes': args.archive.stat().st_size,
              'source_members_sha256': source_hashes, 'selection_protocol_sha256': plan_hash,
              'original_transcript_rows': len(columns['qv']), 'q20_measured_gene_rows': len(selected),
              'regions': regions, 'checker_sha256': sha(Path(__file__)),
              'python': sys.version, 'seconds': time.perf_counter() - start}
    require(not any(args.output.resolve().is_relative_to(root) for root in seen_roots),
            'write the audit report outside prepared inputs')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
