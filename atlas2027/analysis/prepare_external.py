"""Import immutable cropped observations from a prior GitHub acquisition artifact.

No segmentation is inferred. Recheck all measured records against the masks.
Vendor describes this trimmed dataset as a file-format test, so this evaluation
is a transfer test, not biological ground truth or a donor calibration study.
"""
from pathlib import Path
import sys,json,csv,gzip,hashlib,shutil,math
from collections import defaultdict,Counter
from fractions import Fraction as F
import numpy as np
BASE=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--artifact-directory',type=Path,required=True);args=ap.parse_args()
 src=args.artifact_directory;dest=BASE/'external';(dest/'data/source_crops').mkdir(parents=True,exist_ok=True);(dest/'results').mkdir(exist_ok=True);(dest/'docs').mkdir(exist_ok=True)
 schema=json.loads((src/'evaluation/results/schema.json').read_text());out={'scope':'All Q20 measured-gene records retained by the immutable two-FOV breast crop artifact', 'reader':'Python csv/gzip; no Parquet decoding', 'rois':{}}
 total=0
 for orig in sorted((src/'analysis/external_inputs').iterdir()):
  d=dest/'data/source_crops'/orig.name;d.mkdir(exist_ok=True)
  for name in ['masks.npz','transcripts.csv.gz'] :shutil.copy2(orig/name,d/name)
  shutil.copy2(orig/'metadata.json',d/'source_metadata.json')
  meta=json.loads((orig/'metadata.json').read_text());meta['source_archive_sha256']=schema['sha256'];meta['source_archive_url']=schema['url']
  (d/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
  z=np.load(d/'masks.npz');mask=z['cell_mask'];x0,y0=meta['bbox_pixels'][:2];p=1/F(float(meta['transform'][0][0]));records=[];counts=defaultdict(Counter);seen=set()
  with gzip.open(d/'transcripts.csv.gz','rt') as f:
   for row in csv.DictReader(f):
    ident=row['transcript_id'];assert ident not in seen;seen.add(ident)
    x,y=F(float(row['x_location'])),F(float(row['y_location']));px=math.floor(x/p)-x0;py=math.floor(y/p)-y0
    label=int(row['mask_label']);assert 0<=px<mask.shape[1] and 0<=py<mask.shape[0]
    assert int(mask[py,px])==label and float(row['qv'])>=20
    expected='UNASSIGNED' if label==0 else meta['label_to_cell_id'][str(label)]
    assert row['cell_id']==expected
    records.append({'id':ident,'gene':row['feature_name'],'cell':label,'pixel_row':py,'pixel_col':px,'x_um':row['x_location'],'y_um':row['y_location']})
    counts[label][row['feature_name']]+=1
  hashes={name:sha(d/name) for name in ['masks.npz','metadata.json','source_metadata.json','transcripts.csv.gz']}
  out['rois'][orig.name]={'counts':dict(counts),'records':records,'source_hashes':hashes}
  total+=len(records);assert len(records)==meta['q20_gene_rows']
 out['checked_source_records']=total
 (dest/'data/observed_crop_counts.json').write_text(json.dumps(out,indent=2)+'\n')
 print({'rois':len(out['rois']),'checked_RNA':total,'assignment_discrepancies':0,'source_dataset_sha256':schema['sha256']})
if __name__=='__main__':main()
