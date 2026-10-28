"""Inspect the fixed external dataset without selecting molecular outcomes."""
from pathlib import Path
import hashlib
import json
import zipfile
import urllib.request

import numpy as np
import pandas as pd
import zarr

URL = 'https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_human_Breast_2fov/Xenium_V1_human_Breast_2fov_outs.zip'
EXPECTED_MD5 = '4632914eca973a1d532231ea646e10cc'
ROOT = Path('run/external')
ROOT.mkdir(parents=True, exist_ok=True)
ARCHIVE = Path('/tmp/geocert-breast-2fov.zip')
if not ARCHIVE.exists():
    with urllib.request.urlopen(URL, timeout=90) as source, ARCHIVE.open('wb') as target:
        count = 0
        while chunk := source.read(1024 * 1024):
            count += len(chunk)
            if count > 450_000_000:
                raise RuntimeError('Download exceeds fixed size allowance.')
            target.write(chunk)
md5 = hashlib.md5(ARCHIVE.read_bytes()).hexdigest()
if md5 != EXPECTED_MD5:
    raise RuntimeError(f'Archive checksum mismatch: {md5}')
wanted = ('cells.zarr.zip', 'transcripts.parquet', 'cells.parquet',
          'gene_panel.json', 'experiment.xenium', 'features.tsv.gz')
selected = {}
with zipfile.ZipFile(ARCHIVE) as source:
    members = [{'name': item.filename, 'bytes': item.file_size} for item in source.infolist()]
    for item in source.infolist():
        name = Path(item.filename).name
        if name in wanted:
            if name in selected:
                raise RuntimeError(f'Ambiguous archive member: {name}')
            destination = ROOT / name
            destination.write_bytes(source.read(item))
            selected[name] = str(destination)

report = {'url': URL, 'md5': md5,
          'sha256': hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
          'archive_bytes': ARCHIVE.stat().st_size, 'members': members,
          'selected': selected}
for name in ('cells.parquet', 'transcripts.parquet'):
    if name in selected:
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(selected[name])
        report[name] = {'schema': str(pf.schema), 'rows': pf.metadata.num_rows,
                        'row_groups': pf.metadata.num_row_groups}
if 'cells.parquet' in selected:
    cells = pd.read_parquet(selected['cells.parquet'])
    report['cell_columns'] = list(cells.columns)
    chosen = []
    for quantile in (0.25, 0.75):
        point = cells[['x_centroid', 'y_centroid']].quantile(quantile)
        distances = ((cells[['x_centroid', 'y_centroid']] - point) ** 2).sum(axis=1)
        tied = cells.loc[distances == distances.min()].sort_values('cell_id')
        cell = tied.iloc[0]
        chosen.append({'quantile': quantile, 'cell_id': str(cell.cell_id),
                       'center_um': [float(cell.x_centroid), float(cell.y_centroid)]})
    report['locked_crop_centers'] = chosen
if 'cells.zarr.zip' in selected:
    store = zarr.ZipStore(selected['cells.zarr.zip'], mode='r')
    group = zarr.open_group(store, mode='r')
    report['zarr_tree'] = str(group.tree())
    report['zarr_root_attributes'] = dict(group.attrs)
    details = {}
    def inspect(name, obj):
        entry = {'attributes': dict(obj.attrs)}
        if isinstance(obj, zarr.Array):
            entry.update(shape=list(obj.shape), dtype=str(obj.dtype), chunks=list(obj.chunks))
            if obj.size <= 100:
                entry['small_values'] = np.asarray(obj).tolist()
        details[name] = entry
    group.visititems(inspect)
    report['zarr_details'] = details
    store.close()
if 'experiment.xenium' in selected:
    experiment = json.loads(Path(selected['experiment.xenium']).read_text())
    report['experiment_format'] = {k: experiment.get(k) for k in
        ('major_version', 'minor_version', 'pixel_size', 'analysis_sw_version',
         'panel_name', 'panel_tissue_type', 'segmentation_stain')}
if 'gene_panel.json' in selected:
    panel = json.loads(Path(selected['gene_panel.json']).read_text())
    report['panel_top_keys'] = list(panel)
    report['panel_payload_keys'] = list(panel.get('payload', {}))
    report['panel_example_target'] = panel.get('payload', {}).get('targets', [])[:1]
output = Path('evaluation/results/schema.json')
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(report, indent=2, default=str) + '\n')
print(json.dumps({'download_verified': True, 'selected': selected,
                  'locked_crop_centers': report.get('locked_crop_centers'),
                  'zarr_tree': report.get('zarr_tree')}, indent=2))
