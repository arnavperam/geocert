"""Re-ingest supplied transcript crops and verify prepared observations/counts.

The restricted fallback reader is inherited, not an independent ingestion engine.
The breast crops use ordinary CSV/gzip, not the Parquet decoder.
"""
from pathlib import Path
import sys, json, gzip, csv, math
from collections import Counter
from fractions import Fraction as F
import numpy as np
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / 'analysis/vendor'))
try:
    import pyarrow.parquet as pq
    def read_parquet(path):return pq.read_table(path).to_pandas()
    reader = 'PyArrow'
except ImportError:
    from parquet_minimal import read_flat_parquet as read_parquet
    reader = 'Inherited strict flat-Parquet decoder; shared ingestion dependency'
results = []
for base in [R, R / 'external']:
    supplied = json.loads((base/'data/observed_crop_counts.json').read_text())
    for roi, prepared in sorted(supplied['rois'].items()):
        d = base/'data/source_crops'/roi
        mask = np.load(d/'masks.npz')['cell_mask']
        meta = json.loads((d/'metadata.json').read_text())
        expected = {v['id']: v for v in prepared['records']}
        counts = Counter();seen = set()
        if (d/'transcripts.parquet').is_file():
            raw = read_parquet(d/'transcripts.parquet')
            records = ({'id': str(int(v.transcript_id)), 'gene': str(v.feature_name), 'cell': int(v.mask_label),
                        'pixel_row': int(v.pixel_row), 'pixel_col': int(v.pixel_col),
                        'x_um': str(float(v.x_location)), 'y_um': str(float(v.y_location))} for v in raw.itertuples(index=False))
        else:
            p = 1/F(float(meta['transform'][0][0]));x0,y0 = meta['bbox_pixels'][:2]
            with gzip.open(d/'transcripts.csv.gz','rt') as f:
                records = []
                for row in csv.DictReader(f):
                    x,y = F(float(row['x_location'])),F(float(row['y_location']))
                    records.append({'id': row['transcript_id'], 'gene': row['feature_name'], 'cell': int(row['mask_label']),
                                    'pixel_row': math.floor(y/p)-y0, 'pixel_col': math.floor(x/p)-x0,
                                    'x_um': row['x_location'], 'y_um': row['y_location']})
        for v in records:
            ident = v['id']
            if ident in seen:raise AssertionError('duplicate source transcript')
            seen.add(ident)
            if v != expected[ident]:raise AssertionError(('prepared record differs',roi,ident))
            if int(mask[v['pixel_row'],v['pixel_col']]) != v['cell']:raise AssertionError('source owner differs')
            counts[(str(v['cell']),v['gene'])] += 1
        target = Counter({(str(c),g):n for c,gg in prepared['counts'].items() for g,n in gg.items()})
        if counts != target or seen != set(expected):raise AssertionError('source counts or membership differ')
        results.append({'roi':roi,'source_transcripts':len(seen),'prepared_records_match':True,'count_keys':len(counts),'mask_owners_match':True})
result = {'reader_lung':reader,'reader_breast':'Python CSV/gzip','results':results,'total_source_records':sum(v['source_transcripts'] for v in results),'failures':0}
(R/'results/source_input_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
