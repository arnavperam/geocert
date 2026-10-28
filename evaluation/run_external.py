"""Prepare the protocol's two crops; this stage performs no optimization."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import numpy as np
import pandas as pd
import zarr

root = Path('run/external')
out = Path('evaluation/results')
schema = json.loads((out / 'schema.json').read_text())
store = zarr.ZipStore(str(root / 'cells.zarr.zip'), mode='r')
group = zarr.open_group(store, mode='r')
transform = np.asarray(group['masks/homogeneous_transform'])
expected = np.diag([transform[0, 0], transform[1, 1], 1, 1])
if not np.array_equal(transform, expected) or transform[0, 0] != transform[1, 1]:
    raise ValueError('The locked raster adapter requires an isotropic axis-aligned transform.')
scale = float(transform[0, 0])
ids = np.asarray(group['cell_id'])
translation = str.maketrans('0123456789abcdef', 'abcdefghijklmnop')
label_to_id = {i + 1: f'{int(prefix):08x}'.translate(translation) + '-' + str(int(suffix))
               for i, (prefix, suffix) in enumerate(ids)}
id_to_label = {v: k for k, v in label_to_id.items()}
panel = json.loads((root / 'gene_panel.json').read_text())
genes = {target['type']['data']['name'] for target in panel['payload']['targets']
         if target['type']['descriptor'] == 'gene'}
all_rna = pd.read_parquet(root / 'transcripts.parquet')
all_rna = all_rna[(all_rna.qv >= 20) & all_rna.feature_name.isin(genes)].copy()
report = {'protocol_sha256': hashlib.sha256(Path('evaluation/PROTOCOL.md').read_bytes()).hexdigest(),
          'dataset_sha256': schema['sha256'], 'q20_measured_gene_RNA': len(all_rna),
          'panel_gene_count': len(genes), 'crops': []}
for index, selected in enumerate(schema['locked_crop_centers'], 1):
    cx, cy = selected['center_um']
    physical = [cx - 100, cy - 100, cx + 100, cy + 100]
    x0, y0 = [math.floor(v * scale) for v in physical[:2]]
    x1, y1 = [math.ceil(v * scale) for v in physical[2:]]
    shape = group['masks/1'].shape
    if not (0 <= x0 < x1 <= shape[1] and 0 <= y0 < y1 <= shape[0]):
        raise ValueError('A locked crop leaves the raster: retain the failure, do not substitute.')
    cell = np.asarray(group['masks/1'][y0:y1, x0:x1])
    nucleus = np.asarray(group['masks/0'][y0:y1, x0:x1])
    rows = all_rna[(all_rna.x_location >= physical[0]) &
                   (all_rna.x_location < physical[2]) &
                   (all_rna.y_location >= physical[1]) &
                   (all_rna.y_location < physical[3])].copy()
    # The exact same binary-rational transform is used for each measured coordinate.
    s = Fraction(scale)
    px = np.array([math.floor(Fraction(float(v)) * s) - x0 for v in rows.x_location])
    py = np.array([math.floor(Fraction(float(v)) * s) - y0 for v in rows.y_location])
    rows['mask_label'] = cell[py, px]
    supplied = rows.cell_id.map(id_to_label).fillna(0).astype('int64')
    mismatches = rows.loc[supplied.to_numpy() != rows.mask_label.to_numpy(),
                          ['transcript_id', 'cell_id', 'mask_label', 'x_location', 'y_location']]
    roi = f'breast_{index:02d}'
    folder = root / roi
    folder.mkdir(exist_ok=True)
    np.savez_compressed(folder / 'masks.npz', cell_mask=cell, nucleus_mask=nucleus)
    # Seventeen significant digits retain each original float32 coordinate exactly.
    rows.to_csv(folder / 'transcripts.csv.gz', index=False, float_format='%.17g')
    complete = sorted(set(map(int, np.unique(cell))) -
                      set(map(int, np.concatenate([cell[0], cell[-1], cell[:, 0], cell[:, -1]]))) - {0})
    meta = {'roi_id': roi, 'pixel_size_um': 1 / scale, 'bbox_pixels': [x0, y0, x1, y1],
            'bounds_um': physical, 'transform': transform.tolist(), 'selection': selected,
            'q20_gene_rows': len(rows), 'source_assignment_disagreements': len(mismatches),
            'complete_source_labels': complete,
            'label_to_cell_id': {str(k): label_to_id[k] for k in np.unique(cell) if k}}
    (folder / 'metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
    mismatches.to_csv(folder / 'assignment_disagreements.csv', index=False, float_format='%.17g')
    report['crops'].append(meta)
store.close()
(out / 'crops.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'prepared_crops': len(report['crops']),
                 'rna_records': [r['q20_gene_rows'] for r in report['crops']],
                 'assignment_disagreements': [r['source_assignment_disagreements'] for r in report['crops']]}, indent=2))
