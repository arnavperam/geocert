"""Recompute all-budget denominators and external-reference reach on a fixed key set."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import json,csv,math,gzip,hashlib
import numpy as np
from scipy.ndimage import distance_transform_edt
R=Path(__file__).resolve().parents[2]
BUDGETS=list(map(F,json.loads((R/'results/reach/plan_identity.json').read_text())['budgets_um2']))

def readcsv(p):
 with p.open() as f:return list(csv.DictReader(f))
def models(dataset):
 base=R if dataset=='lung' else R/'external'
 return [('baseline',base)]+[(str(r),R/'reach'/dataset/('r'+str(r).replace('.','p'))) for r in ['0.5','1','2','3']]
def loadset(base):
 atlas={a['roi']:a for f in sorted((base/'data').glob('atlas_*.json')) for a in [json.loads(f.read_text())]}
 spectra={(s['roi'],s['cell'],s['gene']):s for s in json.loads((base/'results/cell_spectra.json').read_text())}
 return atlas,spectra

def main():
 summaries=[]; budget_rows=[];reference_rows=[];exclusions=[]; saturation=[]
 for dataset in ('lung','breast'):
  base=R if dataset=='lung' else R/'external';loaded={label:loadset(root) for label,root in models(dataset)}
  # Reference construction precedes its reach evaluation. All widths use the same
  # eligible cell/gene key set and the same input observations.
  references={};exclusionkeys=defaultdict(set)
  baseline=loaded['baseline'][0]
  for roi,a in baseline.items():
   p=float(F(a['pixel_size_um']));z=np.load(base/'data/source_crops'/roi/'masks.npz');mask=z['cell_mask'];nuc=z['nucleus_mask'];h,w=mask.shape
   labelled=np.where(nuc>0,mask,0);dist,inds=distance_transform_edt(labelled==0,sampling=p,return_indices=True);alt=labelled[inds[0],inds[1]];alt[dist>5]=0
   records={r['id']:r for r in a['records']};nucowner={rid:int(alt[r['pixel_row'],r['pixel_col']]) for rid,r in records.items() if min(r['pixel_row']+.5,r['pixel_col']+.5,h-r['pixel_row']-.5,w-r['pixel_col']-.5)*p>=5}
   refpairs=[('nuclear',nucowner)]
   if dataset=='lung':refpairs.append(('proseg',json.loads((R/'results/reach'/f'proseg_join_{roi}.json').read_text())['target_owners']))
   for ref,owners in refpairs:
    source=Counter((r['cell'],r['gene']) for rid,r in records.items() if rid in owners)
    target=Counter((owners[rid],records[rid]['gene']) for rid in owners)
    missing=set(records)-set(owners);excluded={(records[rid]['cell'],records[rid]['gene']) for rid in missing}
    for atl,_ in loaded.values():
     for patch in atl[roi]['patches']:
      for point in patch['points']:
       if point['transcript_id'] in missing:excluded.update((cell,point['gene']) for cell in patch['cell_labels'])
    keys={k for k in source.keys()|target.keys() if k[0] in set(a['complete_cells']) and k not in excluded}
    references[(roi,ref)]={'owners':owners,'source':source,'target':target,'keys':keys,'records':records}
    exclusions.append({'dataset':dataset,'roi':roi,'reference':ref,'retained_transcripts':len(owners),'excluded_transcripts':len(missing),'excluded_complete_cell_gene_keys':len({k for k in excluded if k[0] in set(a['complete_cells'])}),'comparison_keys':len(keys),'changed_count_keys':sum(source[k]!=target[k] for k in keys)})
  for label,root in models(dataset):
   atlas,spectra=loaded[label];rr=readcsv(R/'results/reach'/f'baseline_{dataset}_fine_intervals.csv' if label=='baseline' else root/'results/cell_expression_intervals.csv')
   # This sensitivity table uses the common 17-budget grid, separately from the original primary table.
   grouped=defaultdict(list)
   eligible=set();tight=set();affected=set()
   for row in rr:
    grouped[F(row['budget_um2'])].append(row);key=(row['roi'],int(row['cell']),row['gene'])
    if int(row['independent_upper'])-int(row['independent_lower'])>=2:eligible.add(key)
    if (row['lower'],row['upper'])!=(row['independent_lower'],row['independent_upper']):tight.add(key);affected.add(key[:2])
   all_edges=0;selected=0;patches=0;editable=0;displacements=[];radiuses=[]
   for roi,a in atlas.items():
    m=np.load(base/'data/source_crops'/roi/'masks.npz')['cell_mask'];cells=set(a['complete_cells'])
    for left,right in [(m[:,:-1],m[:,1:]),(m[:-1,:],m[1:,:])]:
     q=left!=right;u=left[q];v=right[q]
     all_edges+=sum(int(x) in cells or int(y) in cells for x,y in zip(u,v))
    for pat in a['patches']:
     p=float(F(a['pixel_size_um']));radiuses.append(p*float(F(pat['halfwidth_pixels']))/math.sqrt(2));displacements.append(p*float(F(pat['band_halfwidth_pixels']))/math.sqrt(2))
    selected+=a['summary']['selected_source_edges'];patches+=len(a['patches']);editable+=a['summary']['unique_editable_RNA']
   sat=[]
   for s in spectra.values():
    d={int(k):F(v) for k,v in s['minima'].items()};sat.append(max(d[min(d)],d[max(d)]))
   summary={'dataset':dataset,'setting':label,'patches':patches,'selected_edges':selected,'incident_edges':all_edges,'edge_coverage':selected/all_edges,'complete_cells':sum(len(a['complete_cells']) for a in atlas.values()),'editable_RNA':editable,'count_questions':len(spectra),'maximum_support_um':max(radiuses,default=0),'maximum_displacement_cap_um':max(displacements,default=0),'marginal_width_ge2_questions':len(eligible),'tightened_questions':len(tight),'cells_with_tightening':len(affected),'exact_rows':len(rr),'last_endpoint_saturation_budget_um2':str(max(sat,default=F(0))),'last_endpoint_saturation_float':float(max(sat,default=F(0)))}
   summaries.append(summary)
   previous={}
   for budget,rows in sorted(grouped.items()):
    now={(q['roi'],q['cell'],q['gene']):(int(q['lower']),int(q['upper'])) for q in rows}
    qualifying=[q for q in rows if int(q['independent_upper'])-int(q['independent_lower'])>=2]
    budget_rows.append({'dataset':dataset,'setting':label,'budget_um2':str(budget),'questions':len(rows),'non_singleton':sum(l!=u for l,u in now.values()),'summed_interval_width':sum(u-l for l,u in now.values()),'endpoint_changes_from_previous_budget':sum(now[k]!=previous[k] for k in now) if previous else 0,'eligible_marginal_width_ge2':len(qualifying),'tightened_questions':sum((q['lower'],q['upper'])!=(q['independent_lower'],q['independent_upper']) for q in rows)})
    previous=now
   # Same reference keys, same missing-record policy, across every selected atlas.
   for ref in ('nuclear','proseg'):
    if dataset=='breast' and ref=='proseg':continue
    changed=[];allq=[];fixedchanged=0;wrongdest=0
    for roi,a in atlas.items():
     r=references[(roi,ref)];point_patch={p['transcript_id']:pat for pat in a['patches'] for p in pat['points']}
     for rid,owner in r['owners'].items():
      if owner==r['records'][rid]['cell']:continue
      if rid not in point_patch:fixedchanged+=1
      elif owner not in point_patch[rid]['cell_labels']:wrongdest+=1
     for cell,gene in sorted(r['keys']):
      source=r['source'][(cell,gene)];want=r['target'][(cell,gene)]
      s=spectra.get((roi,cell,gene));spec={int(v):F(d) for v,d in s['minima'].items()} if s else {a['source_counts'].get(str(cell),{}).get(gene,0):F(0)}
      assert source in spec and spec[source]==0,(dataset,roi,cell,gene)
      minimum=spec.get(want);row={'dataset':dataset,'setting':label,'roi':roi,'reference':ref,'cell':cell,'gene':gene,'source_count':source,'reference_count':want,'changed':source!=want,'minimum_area_um2':str(minimum) if minimum is not None else 'inf'}
      allq.append(row)
      if source!=want:changed.append(row)
    summary[f'{ref}_changed_count_keys']=len(changed);summary[f'{ref}_attainable_changed_counts']=sum(r['minimum_area_um2']!='inf' for r in changed)
    summary[f'{ref}_asymptotic_reach']=sum(r['minimum_area_um2']!='inf' for r in changed)/len(changed) if changed else None
    summary[f'{ref}_changed_fixed_RNA']=fixedchanged;summary[f'{ref}_unavailable_destination_RNA']=wrongdest
    summary[f'{ref}_full_vector_ruled_out_by_fixed_or_destination']=bool(fixedchanged or wrongdest)
    reference_rows.extend(allq)
    for budget in BUDGETS:
     saturation.append({'dataset':dataset,'setting':label,'reference':ref,'budget_um2':str(budget),'changed_count_keys':len(changed),'attainable_changed_counts':sum(q['minimum_area_um2']!='inf' and F(q['minimum_area_um2'])<=budget for q in changed)})
  print(dataset,'complete',flush=True)
 out=R/'results/reach'
 for name,items in [('summary',summaries),('budget_summary',budget_rows),('reference_budget_summary',saturation),('reference_exclusions',exclusions)]:
  (out/(name+'.json')).write_text(json.dumps(items,indent=2)+'\n')
  with (out/(name+'.csv')).open('w',newline='') as f:
   fields=list(dict.fromkeys(k for row in items for k in row));w=csv.DictWriter(f,fields);w.writeheader();w.writerows(items)
 with gzip.open(out/'reference_counts.csv.gz','wt',newline='') as f:
  w=csv.DictWriter(f,list(reference_rows[0]));w.writeheader();w.writerows(reference_rows)
 print('All settings:',len(summaries),'reference rows:',len(reference_rows))
if __name__=='__main__':main()
