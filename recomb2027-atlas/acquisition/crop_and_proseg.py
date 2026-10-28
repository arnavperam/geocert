#!/usr/bin/env python3
"""Frozen prospective crop selection and five Proseg executions per crop.

Not executed in the reported study. Keep output labels and exclusions; subsequent
Atlas ingestion needs a separately verified mask/coordinate adapter.
"""
from pathlib import Path
import argparse,json,hashlib,subprocess,time,shutil
import numpy as np
import pandas as pd

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--dataset',choices=['breast','lung'],required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 cells=pd.read_parquet(a.input/'cells.parquet');rng=np.random.default_rng(20261028);xy=cells[['x_centroid','y_centroid']].to_numpy();chosen=[]
 exclusions=[]
 if a.dataset=='lung':
  root=Path(__file__).resolve().parents[1]
  for f in sorted((root/'data/source_crops').glob('*/metadata.json')):exclusions.append(json.loads(f.read_text())['selection']['bounds_um'])
 for i in rng.permutation(len(cells)):
  x,y=map(float,xy[i]);box=[x-100,y-100,x+100,y+100]
  if any(not (box[2]+100<q[0] or q[2]+100<box[0] or box[3]+100<q[1] or q[3]+100<box[1]) for q in exclusions+[q['bounds_um'] for q in chosen]):continue
  count=int(((abs(xy[:,0]-x)<=100)&(abs(xy[:,1]-y)<=100)).sum())
  if count<30:continue
  chosen.append({'cell_id':str(cells.iloc[i]['cell_id']),'bounds_um':box,'centroids_in_core':count})
  if len(chosen)==2:break
 if len(chosen)!=2:raise ValueError('Fewer than two crops satisfy the frozen rule')
 # Selection is persisted before reading any transcript assignments or Proseg output.
 (a.out/'crop_selection.json').write_text(json.dumps({'seed':20261028,'core_width_um':200,'halo_um':100,'crops':chosen},indent=2)+'\n')
 version=subprocess.check_output(['proseg','--version'],text=True).strip();assert '3.2.0' in version
 report={'dataset':a.dataset,'proseg_version':version,'runs':[],'selected_crops':chosen}
 for j,q in enumerate(chosen,1):
  root=a.out/f'crop_{j:02d}';root.mkdir(exist_ok=True);x0,y0,x1,y1=q['bounds_um']
  # PyArrow filters reduce memory and preserve the unmodified transcript schema.
  t=pd.read_parquet(a.input/'transcripts.parquet',filters=[('x_location','>=',x0-100),('x_location','<=',x1+100),('y_location','>=',y0-100),('y_location','<=',y1+100)])
  t.to_csv(root/'input.csv.gz',index=False,compression='gzip');q['input_rows']=len(t)
  for rep in range(5):
   dest=root/f'run_{rep+1:02d}';dest.mkdir(exist_ok=True)
   cmd=['proseg','--xenium',str((root/'input.csv.gz').resolve()),'--nthreads','2','--output-path',str(dest.resolve()),'--output-transcript-metadata','transcript-metadata.csv.gz','--output-cell-metadata','cell-metadata.csv.gz','--output-gene-metadata','gene-metadata.csv.gz','--output-counts','counts.mtx.gz']
   start=time.time()
   with (dest/'execution.log').open('w') as f:
    try:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=3600);status=r.returncode
    except subprocess.TimeoutExpired:status='time_limit'
   rec={'crop':j,'replicate':rep+1,'command':cmd,'returncode':status,'seconds':time.time()-start,'random_seed':'upstream does not expose one; stochastic run'}
   report['runs'].append(rec);(dest/'run.json').write_text(json.dumps(rec,indent=2)+'\n')
   # No undisclosed fallback: failure remains a failure, with all logs retained.
   (a.out/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
