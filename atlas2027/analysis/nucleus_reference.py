"""Deterministic image-only reference and incompatibility audit; not Proseg."""
from pathlib import Path
import sys,json,csv,time,os,hashlib
from fractions import Fraction as F
from collections import defaultdict,Counter
import numpy as np
from scipy.ndimage import distance_transform_edt
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));D=Path(os.environ.get('GEOCERT_DATA_ROOT',R))
from atlas.model import make_patch
protocol=json.loads((R/'docs/PAIRED_REFERENCE_PROTOCOL.json').read_text());radius=protocol['radius_um'];allrows=[];summaries=[];start=time.time()
known={(q['roi'],q['cell'],q['gene']):q for q in json.loads((D/'results/cell_spectra.json').read_text())}
for ap in sorted((D/'data').glob('atlas_*.json')):
 a=json.loads(ap.read_text());roi=a['roi'];p=float(F(a['pixel_size_um']));z=np.load(D/'data/source_crops'/roi/'masks.npz');m=z['cell_mask'];n=z['nucleus_mask'];h,w=m.shape
 nuclei=np.where(n>0,m,0);dist,inds=distance_transform_edt(nuclei==0,sampling=p,return_indices=True)
 alt=nuclei[inds[0],inds[1]].copy();alt[dist>radius]=0
 source=defaultdict(Counter);target=defaultdict(Counter);removed=defaultdict(Counter);newowners={};eligible={}
 for r in a['records']:
  row,col=r['pixel_row'],r['pixel_col'];k=r['cell'];gene=r['gene'];edge=min(row+.5,col+.5,h-row-.5,w-col-.5)*p
  if edge<radius:removed[k][gene]+=1;continue
  source[k][gene]+=1;o=int(alt[row,col]);target[o][gene]+=1;newowners[r['id']]=o;eligible[r['id']]=r
 patchbyRNA={r['transcript_id']:(q,i) for q in a['patches'] for i,r in enumerate(q['points'])}
 # Exclusion must not change any editable RNA in a certified target cell.
 # If it does, that gene is unavailable rather than silently altering a spectrum.
 excluded_mutable={(k,r['gene']) for q in a['patches'] for r in q['points'] if r['transcript_id'] not in eligible for k in q['cell_labels'] if k>0}
 keys={(k,g) for k in a['complete_cells'] for g in set(source[k])|set(target[k])}|{(q['cell'],q['gene']) for q in known.values() if q['roi']==roi}
 cnt=Counter();badoutside=badinside=0;totalcost=F(0);infeasiblepatch=0
 for k,g in sorted(keys):
  if (k,g) in excluded_mutable:cnt['excluded_boundary_gene']+=1;continue
  rec=known.get((roi,k,g));basecount=a['source_counts'].get(str(k),{}).get(g,0)
  spec={int(v)-removed[k][g]:F(d) for v,d in rec['minima'].items()} if rec else {basecount-removed[k][g]:F(0)}
  observed=source[k][g];wanted=target[k][g]
  assert observed in spec,(roi,k,g,observed,spec)
  minimum=spec.get(wanted);cnt['queries']+=1;cnt['same_count']+=wanted==observed;cnt['different_count']+=wanted!=observed
  cnt['target_count_attainable']+=minimum is not None;cnt['changed_target_attainable']+=(minimum is not None and wanted!=observed)
  allrows.append({'roi':roi,'cell':k,'gene':g,'source_count':observed,'reference_count':wanted,'attainable':minimum is not None,'minimum_area_um2':str(minimum) if minimum is not None else 'inf'})
 for rid,r in eligible.items():
  if r['cell']==newowners[rid]:continue
  if rid not in patchbyRNA:badoutside+=1
  else:
   q,i=patchbyRNA[rid]
   if newowners[rid] not in q['cell_labels']:badinside+=1
 for q in a['patches']:
  pins={}
  for i,r in enumerate(q['points']):
   rid=r['transcript_id']
   if rid not in eligible:continue
   if newowners[rid] not in q['cell_labels']:continue
   pins[i]=q['cell_labels'].index(newowners[rid])
  if not pins:continue
  e=make_patch(q);ans=e.clipped(pins)
  if ans['feasible']:totalcost+=F(ans['area'])*F(q['area_scale_um2'])
  else:infeasiblepatch+=1
 rows=[r for r in allrows if r['roi']==roi]
 cnt.update({'retained_comparison_RNA':len(eligible),'changed_RNA':sum(r['cell']!=newowners[rid] for rid,r in eligible.items()),'changed_fixed_exterior_RNA':badoutside,'unavailable_patch_destination_RNA':badinside,'infeasible_pinned_patches':infeasiblepatch})
 summaries.append({'roi':roi,**cnt,'retained_reference_vector_attainable':not(badoutside or badinside or infeasiblepatch),
 'budget_coverage':[{'budget_um2':str(b),'covered_count_entries':sum(r['attainable'] and F(r['minimum_area_um2'])<=b for r in rows),'count_entries':len(rows)} for b in map(F,['0','.05','.1','.25','.5','1','2','4'])]})
 print(summaries[-1],flush=True)
(D/'results/nucleus_reference_summary.json').write_text(json.dumps({'protocol':protocol,'rois':summaries,'seconds':time.time()-start},indent=2)+'\n')
with (D/'results/nucleus_reference_counts.csv').open('w',newline='') as f:
 out=csv.DictWriter(f,fieldnames=list(allrows[0]));out.writeheader();out.writerows(allrows)
