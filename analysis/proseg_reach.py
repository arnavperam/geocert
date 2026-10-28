"""Verify imported Proseg records and measure their reach in each atlas.

The alternative is one recorded stochastic crop-level run, with no factorization.
Only common observed molecules enter count comparisons. Excluded molecules remain
geometric constraints; comparisons with an excluded editable gene are withheld.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import csv,gzip,json,hashlib,math,sys
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
BUDGETS=list(map(F,'0,.01,.025,.05,.075,.1,.15,.2,.25,.35,.5,.75,1,2,4,8,16'.split(',')))
def rows(p):
 with gzip.open(p,'rt',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
 out={};audits=[]
 for roi in ('roi_01','roi_02'):
  d=R/'data/proseg'/roi;run=json.loads((d/'run.json').read_text())
  for fn,h in run['output_sha256'].items():assert sha(d/fn)==h,(roi,fn)
  assert sha(R/'data/source_crops'/roi/'transcripts.parquet')==run['input_parquet_sha256']
  meta=json.loads((R/'data/source_crops'/roi/'metadata.json').read_text());name2label={v:int(k) for k,v in meta['label_to_cell_id'].items()};name2label['UNASSIGNED']=0
  src={r['transcript_id']:r for r in rows(d/'input.csv.gz')};cells=rows(d/'cell-metadata.csv.gz');cellmap={r['cell']:name2label[r['original_cell_id']] for r in cells}
  assert len(src)==run['comparison']['source_transcripts'];assert len(cellmap)==len(cells)
  obs=json.loads((R/'data/observed_crop_counts.json').read_text())['rois'][roi]['records'];lookup={r['id']:r for r in obs}
  assert set(lookup)==set(src)
  for tid,r in src.items():assert r['feature_name']==lookup[tid]['gene'] and int(r['mask_label'])==lookup[tid]['cell']
  target={};origcounts=Counter();altcounts=Counter();seen=set();maxxy=0
  for q in rows(d/'transcript-metadata.csv.gz'):
   tid=q['transcript_id'];assert tid not in seen and tid in src;seen.add(tid);r=src[tid]
   assert q['gene']==r['feature_name']
   for axis in ['x','y']:
    err=abs(float(q['observed_'+axis])-float(r[axis+'_location']));maxxy=max(maxxy,err);assert err<=.002
   k=cellmap[q['assignment']] if q['assignment'] else 0;target[tid]=k
   origcounts[(lookup[tid]['cell'],q['gene'])]+=1;altcounts[(k,q['gene'])]+=1
  excluded=set(src)-seen;assert excluded==set(run['comparison']['excluded_source_ids'])
  audits.append({'roi':roi,'input_records':len(src),'output_records':len(target),'excluded_records':len(excluded),'verified_output_hashes':len(run['output_sha256']),
   'changed_owners':sum(lookup[t]['cell']!=k for t,k in target.items()),'background':sum(k==0 for k in target.values()),'max_observed_xy_printing_difference_um':maxxy,
   'scope':'Raw saved runs checked; no Proseg rerun, truth, halo, seed replication, or native-accuracy inference.'})
  out[roi]=(target,excluded,origcounts,altcounts,lookup)
 return out,audits

def main():
 alternatives,audits=prepare();summaries=[];allrows=[]
 for name in ['cap4','cap8','cap16','cap20']:
  base=R/'reach_study'/name/'lung';specs={(q['roi'],q['cell'],q['gene']):q for q in json.loads((base/'results/cell_spectra.json').read_text())}
  for ap in sorted((base/'data').glob('atlas_*.json')):
   a=json.loads(ap.read_text());roi=a['roi'];alt,excluded,source,target,lookup=alternatives[roi];complete=set(a['complete_cells'])
   excluded_src=Counter((lookup[t]['cell'],lookup[t]['gene']) for t in excluded)
   excluded_mutable={(k,r['gene']) for p in a['patches'] for r in p['points'] if r['transcript_id'] in excluded for k in p['cell_labels'] if k>0}
   keys={k for k in source.keys()|target.keys() if k[0] in complete}|{(q['cell'],q['gene']) for q in specs.values() if q['roi']==roi}
   rr=[]
   for k,g in sorted(keys):
    if (k,g) in excluded_mutable:continue
    q=specs.get((roi,k,g));removed=excluded_src[k,g]
    n=a['source_counts'].get(str(k),{}).get(g,0);s={int(v)-removed:F(c) for v,c in q['minima'].items()} if q else {n-removed:F(0)}
    old=source[k,g];new=target[k,g];assert old==n-removed and s.get(old)==0
    minimum=s.get(new)
    r={'width':name,'atlas_id':a['atlas_id'],'roi':roi,'cell':k,'gene':g,'source_count_common':old,'proseg_count_common':new,
      'changed':old!=new,'attainable_at_some_budget':minimum is not None,'minimum_area_um2':str(minimum) if minimum is not None else 'inf'}
    rr.append(r);allrows.append(r)
   changed=[r for r in rr if r['changed']];summary={'width':name,'roi':roi,'complete_cells':len(complete),'count_queries':len(rr),'excluded_editable_gene_queries':len(excluded_mutable),
    'changed_count_queries':len(changed),'changed_count_attainable_at_some_budget':sum(r['attainable_at_some_budget'] for r in changed),'by_budget':[]}
   for b in BUDGETS:
    summary['by_budget'].append({'budget_um2':str(b),'attainable_changed_counts':sum(r['attainable_at_some_budget'] and F(r['minimum_area_um2'])<=b for r in changed),'changed_counts':len(changed)})
   summaries.append(summary)
 (R/'results/proseg_import_audit.json').write_text(json.dumps(audits,indent=2)+'\n')
 (R/'results/proseg_reach.json').write_text(json.dumps({'source':'Imported single recorded Proseg 3.2.0 run per clipped lung crop, --no-factorization','audit':audits,'summaries':summaries},indent=2)+'\n')
 with (R/'results/proseg_reach_counts.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
