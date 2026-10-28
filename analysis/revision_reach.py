"""Frozen lung-only ribbon-width/reference-reach sensitivity study.

Original observed coordinates and all geometric RNA constraints are retained.
Reference counts are scored only on a shared observation population. This is
model reach, not assignment accuracy, calibration, or an independent donor test.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
from types import SimpleNamespace
import argparse,contextlib,csv,gzip,hashlib,json,math,os,shutil,sys,time
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
OUT=R/'results/revision_reach'
PROSEG=R/'results/proseg'
from atlas.model import make_patch,verify_witness,compose,serialized
from atlas.chain_spectrum import compile_chain,spectrum as chain_spectrum

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(path,obj):Path(path).write_text(json.dumps(serialized(obj),indent=2)+'\n')
def load_csv(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else Path(p).open()) as f:return list(csv.DictReader(f))
def protocol():return json.loads((R/'docs/REVISION_REACH_PROTOCOL.json').read_text())
def instantiate(cap):
 dest=OUT/f'cap_{cap}';(dest/'data').mkdir(parents=True,exist_ok=True);(dest/'results').mkdir(exist_ok=True)
 link=dest/'data/source_crops'
 if not link.exists():link.symlink_to(R/'data/source_crops',target_is_directory=True)
 if cap==4:
  for p in (R/'data').glob('atlas_roi_*.json'):shutil.copyfile(p,dest/'data'/p.name)
  shutil.copyfile(R/'results/atlas_extraction.json',dest/'results/atlas_extraction.json')
 else:
  # Parameterize a frozen builder without changing the production source file.
  source=(R/'analysis/build_atlas.py').read_text()
  old='for width in (F(4),F(2),F(1),F(1,2)):'
  assert source.count(old)==1
  source=source.replace(old,'for width in STUDY_WIDTHS:')
  source=source.replace("(ROOT/f'data/atlas_{d.name}.json').write_text", "(DEST/f'data/atlas_{d.name}.json').write_text")
  source=source.replace("(ROOT/('results/atlas_extraction_'", "(DEST/('results/atlas_extraction_'")
  ns={'__file__':str(R/'analysis/build_atlas.py'),'__name__':'revision_parameterized_builder'}
  oldenv=os.environ.pop('GEOCERT_DATA_ROOT',None)
  try:
   exec(compile(source,str(R/'analysis/build_atlas.py'),'exec'),ns)
   ns.update(DEST=dest,STUDY_WIDTHS=list(map(F,[str(cap),'4','2','1','1/2'])))
   with (dest/'extraction.log').open('w') as log,contextlib.redirect_stdout(log):ns['main']()
  finally:
   if oldenv is not None:os.environ['GEOCERT_DATA_ROOT']=oldenv
 # Repeat the separate rational geometric acceptance check on every width.
 import subprocess
 env=dict(os.environ,GEOCERT_DATA_ROOT=str(dest))
 with (dest/'verification.log').open('w') as log:
  subprocess.run([sys.executable,str(R/'analysis/verify_atlas.py')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 return dest

def references(a):
 import numpy as np
 from scipy.ndimage import distance_transform_edt
 roi=a['roi'];p=float(F(a['pixel_size_um']));records={r['id']:r for r in a['records']}
 raw=PROSEG/roi
 source=load_csv(raw/'input.csv.gz');transcripts=load_csv(raw/'transcript-metadata.csv.gz');cells=load_csv(raw/'cell-metadata.csv.gz')
 ids={};labelmap={}
 for s in source:
  tid=s['transcript_id'];assert tid not in ids;ids[tid]=s
  if s['cell_id'] in labelmap:assert labelmap[s['cell_id']]==int(s['mask_label'])
  labelmap[s['cell_id']]=int(s['mask_label'])
 assert set(ids)==set(records)
 cmap={s['cell']:labelmap[s['original_cell_id']] for s in cells}
 owners={}
 for q in transcripts:
  tid=q['transcript_id'];assert tid not in owners and tid in records
  r=records[tid];assert r['gene']==q['gene']
  for axis in ['x','y']:assert abs(float(r[axis+'_um'])-float(q['observed_'+axis]))<=.002
  owners[tid]=cmap[q['assignment']] if q['assignment'] else 0
  assert int(ids[tid]['mask_label'])==r['cell']
 # Nuclear expansion: exact same original observed positions, with fixed 5um border guard.
 z=np.load(R/'data/source_crops'/roi/'masks.npz');mask=z['cell_mask'];nuc=z['nucleus_mask'];h,w=mask.shape
 nuclei=np.where(nuc>0,mask,0);distance,indices=distance_transform_edt(nuclei==0,sampling=p,return_indices=True)
 alt=nuclei[indices[0],indices[1]].copy();alt[distance>5]=0
 nuclear={}
 for tid,r in records.items():
  row,col=r['pixel_row'],r['pixel_col']
  if min(row+.5,col+.5,h-row-.5,w-col-.5)*p>=5:nuclear[tid]=int(alt[row,col])
 return {'proseg':owners,'nuclear':nuclear}

def evaluate(a,reference,name,cap,dest):
 budgets=list(map(F,protocol()['budgets_um2']));pop=set(reference);source=Counter();target=Counter()
 records={r['id']:r for r in a['records']}
 for tid,o in reference.items():
  r=records[tid];source[(r['cell'],r['gene'])]+=1;target[(o,r['gene'])]+=1
 parts=defaultdict(list);singles=defaultdict(list);local=[];chain_queries=0;enums=0;verification=0
 for j,p in enumerate(a['patches']):
  selected=[i for i,r in enumerate(p['points']) if r['transcript_id'] in pop]
  if not selected:continue
  e=make_patch(p);scale=F(p['area_scale_um2']);genehist=Counter(p['points'][i]['gene'] for i in selected)
  for i in selected:
   r=p['points'][i];ans=e.clipped({i:1-r['source']})
   if ans['feasible']:
    cost=F(ans['area'])*scale;verify_witness(p,ans['owners'],ans['boundaries'],cost);verification+=1
    for k in p['cell_labels']:
     if k>0:singles[(k,r['gene'])].append((p['cell_labels'][r['source']],cost))
  mesh=compile_chain(e.problem,e.pins,seconds=120) if max(genehist.values())>1 else None
  for gene,n in sorted(genehist.items()):
   scores=[[0,0] for _ in range(e.problem.n)]
   for i in selected:
    if p['points'][i]['gene']==gene:scores[i]=[1,0]
   if n>1:
    ans=chain_spectrum(mesh,scores,seconds=120);assert ans['complete'];sp=SimpleNamespace(minima=ans['minima'],witnesses=ans['witnesses']);chain_queries+=1
   else:
    sp=e.spectrum(scores);assert sp.complete;enums+=sp.examined
   lsource=sum(p['points'][i]['gene']==gene and p['points'][i]['source']==0 for i in selected)
   entry={'patch':p['id'],'gene':gene,'minima':{},'witnesses':{}}
   for value,cost in sp.minima.items():
    w=sp.witnesses[value];physical=F(cost)*scale;verify_witness(p,w['owners'],w['boundaries'],physical);verification+=1
    entry['minima'][value]=physical;entry['witnesses'][value]={'area_um2':physical,'owners':w['owners'],'boundaries':w['boundaries']}
   local.append(entry)
   for side,k in enumerate(p['cell_labels']):
    if k>0:parts[k,gene].append((p['id'],{(v-lsource)*(1 if side==0 else -1):c*scale for v,c in sp.minima.items()}))
  if j%100==0:print('spectrum',cap,a['roi'],name,j,flush=True)
 keys={k for k in source.keys()|target.keys()|parts.keys() if k[0] in a['complete_cells']}
 queries=[];rows=[]
 for k,g in sorted(keys):
  state=compose(parts[(k,g)],source[k,g]);wanted=target[k,g];cost=state.get(wanted)
  queries.append({'cap':cap,'roi':a['roi'],'reference':name,'cell':k,'gene':g,'source_count':source[k,g],
   'reference_count':wanted,'changed':source[k,g]!=wanted,'attainable':cost is not None,'minimum_area_um2':str(cost[0]) if cost else 'inf',
   'attaining_trace':cost[1] if cost else {},'minima':{str(v):str(c) for v,(c,_) in state.items()},'witness_traces':{str(v):t for v,(_,t) in state.items()}})
  for b in budgets:
   values=[v for v,(c,_) in state.items() if c<=b];assert values
   lo,hi=min(values),max(values);ml=source[k,g]-sum(owner==k and c<=b for owner,c in singles[k,g]);mh=source[k,g]+sum(owner!=k and c<=b for owner,c in singles[k,g])
   assert ml<=lo<=hi<=mh
   rows.append({'cap':cap,'roi':a['roi'],'reference':name,'cell':k,'gene':g,'budget_um2':str(b),
      'source_count':source[k,g],'reference_count':wanted,'lower':lo,'upper':hi,'marginal_lower':ml,'marginal_upper':mh,
      'reference_count_attainable':cost is not None and cost[0]<=b})
 with gzip.open(dest/f'{a["roi"]}_{name}_witnesses.json.gz','wt') as f:json.dump(serialized(local),f)
 write(dest/f'{a["roi"]}_{name}_queries.json',queries)
 with (dest/f'{a["roi"]}_{name}_intervals.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 changed=[q for q in queries if q['changed']]
 return {'roi':a['roi'],'reference':name,'cap_chart_units':cap,'shared_transcripts':len(pop),'excluded_from_count_population':len(records)-len(pop),
  'geometric_transcripts_retained':len(records),'queries':len(queries),'changed_queries':len(changed),'changed_attainable_unlimited':sum(q['attainable'] for q in changed),
  'chain_queries':chain_queries,'enumerated_singleton_patterns':enums,'witness_checks':verification,
  'by_budget':[{'budget_um2':str(b),'changed_attainable':sum(q['attainable'] and F(q['minimum_area_um2'])<=b for q in changed),
   'non_singleton':sum(r['lower']!=r['upper'] for r in rows if F(r['budget_um2'])==b),
   'marginal_width_ge2':sum(r['marginal_upper']-r['marginal_lower']>=2 for r in rows if F(r['budget_um2'])==b),
   'tightened':sum((r['lower'],r['upper'])!=(r['marginal_lower'],r['marginal_upper']) for r in rows if F(r['budget_um2'])==b)} for b in budgets]}

def run(cap):
 t=time.time();dest=instantiate(cap);summaries=[];geo=[]
 for af in sorted((dest/'data').glob('atlas_roi_*.json')):
  a=json.loads(af.read_text());p=float(F(a['pixel_size_um']));s=a['summary']
  geo.append({'roi':a['roi'],'atlas_id':a['atlas_id'],'proposed_cap_chart_units':cap,'proposed_displacement_cap_um':3*cap*p/(4*math.sqrt(2)),
    'actual_max_displacement_cap_um':max(float(F(q['band_halfwidth_pixels']))*p/math.sqrt(2) for q in a['patches']),
    'accepted_patches':len(a['patches']),'candidates':s['candidates'],'selected_edges':s['selected_source_edges'],'eligible_edges':s['eligible_shared_edges'],
    'edge_coverage':s['selected_source_edges']/s['eligible_shared_edges'],'patch_acceptance':len(a['patches'])/s['candidates'],
    'editable_transcripts':s['unique_editable_RNA'],'width_histogram':s['width_histogram']})
  for name,ref in references(a).items():summaries.append(evaluate(a,ref,name,cap,dest/'results'))
 write(dest/'summary.json',{'protocol_sha256':sha(R/'docs/REVISION_REACH_PROTOCOL.json'),'geometry':geo,'references':summaries,'seconds':time.time()-t})
 (dest/'data/source_crops').unlink();print('COMPLETE',cap,round(time.time()-t,2),flush=True)

def summarize():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 data=[json.loads((OUT/f'cap_{c}/summary.json').read_text()) for c in protocol()['caps_chart_units']]
 # The primary reach denominator must be identical across the width models.
 for name in ['nuclear','proseg']:
  populations=[]
  for cap in protocol()['caps_chart_units']:
   pop={}
   for path in (OUT/f'cap_{cap}/results').glob(f'*_{name}_queries.json'):
    for q in json.loads(path.read_text()):
     if q['changed']:pop[(q['roi'],q['cell'],q['gene'])]=(q['source_count'],q['reference_count'])
   populations.append(pop)
  assert all(pop==populations[0] for pop in populations), 'reference comparison population changed across widths'
 result=[]
 for d in data:
  geo=d['geometry'];cap=geo[0]['proposed_cap_chart_units']
  denominators={q['roi']:q['all_incident_unique_edges'] for q in json.loads((R/'results/coverage_summary.json').read_text())}
  for g in geo:g['all_incident_edges']=denominators[g['roi']];g['all_incident_edge_coverage']=g['selected_edges']/g['all_incident_edges']
  coverage=sum(x['selected_edges'] for x in geo)/sum(x['all_incident_edges'] for x in geo)
  for name in ['nuclear','proseg']:
   refs=[q for q in d['references'] if q['reference']==name];n=sum(q['changed_queries'] for q in refs)
   result.append({'cap':cap,'proposed_max_displacement_um':max(g['proposed_displacement_cap_um'] for g in geo),'actual_max_displacement_um':max(g['actual_max_displacement_cap_um'] for g in geo),
    'accepted_patches':sum(g['accepted_patches'] for g in geo),'candidates':sum(g['candidates'] for g in geo),'boundary_coverage':coverage,'editable_transcripts':sum(g['editable_transcripts'] for g in geo),
    'reference':name,'changed_queries':n,'attainable_unlimited':sum(q['changed_attainable_unlimited'] for q in refs),
    'by_budget':[{'budget_um2':b,'attainable':sum(next(v['changed_attainable'] for v in q['by_budget'] if v['budget_um2']==b) for q in refs)} for b in protocol()['budgets_um2']]})
 write(OUT/'summary.json',{'protocol':protocol(),'aggregate':result,'interpretation':'Conditional count-level reach on fixed shared observation populations, not joint-vector reach, accuracy or confidence. Width-cap atlases are not nested.'})
 plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,'legend.fontsize':10})
 fig,axes=plt.subplots(1,2,figsize=(6.3,2.55),constrained_layout=True)
 for name,color,marker in [('nuclear','#236a89','o'),('proseg','#ab4e25','s')]:
  rr=[r for r in result if r['reference']==name];x=[r['proposed_max_displacement_um'] for r in rr]
  axes[0].plot(x,[100*r['attainable_unlimited']/r['changed_queries'] for r in rr],marker=marker,color=color,label=name.capitalize())
 rr=[r for r in result if r['reference']=='nuclear']
 axes[1].plot([r['proposed_max_displacement_um'] for r in rr],[100*r['boundary_coverage'] for r in rr],color='#5d6265',marker='^')
 axes[0].set(xlabel='Displacement cap (µm)',ylabel='Changed counts reached (%)',title='A. Exact count reach',ylim=(0,30));axes[0].legend(loc='upper left')
 axes[1].set(xlabel='Displacement cap (µm)',ylabel='Editable boundary (%)',title='B. Boundary coverage',ylim=(0,100))
 for ax in axes:ax.grid(alpha=.15);ax.set_xticks([.45,1.35,3.16],['0.45','1.35','3.16'])
 for ext in ['pdf','png']:fig.savefig(R/'manuscript/figures'/f'reach_width.{ext}',dpi=220)
 plt.close(fig)
 fig,ax=plt.subplots(figsize=(6.3,3.1),constrained_layout=True)
 for name,color in [('nuclear','#236a89'),('proseg','#ab4e25')]:
  for r in [q for q in result if q['reference']==name]:
   ax.plot([float(F(z['budget_um2'])) for z in r['by_budget']],[100*z['attainable']/r['changed_queries'] for z in r['by_budget']],color=color,linestyle={4:':',12:'--',28:'-'}[r['cap']],label=f"{name.capitalize()}, cap {r['proposed_max_displacement_um']:.2f} µm")
 ax.set(xlabel='Shared area budget (µm²)',ylabel='Changed counts reached (%)',title='Exact count reach by budget',ylim=(0,30));ax.legend(loc='lower right',ncol=2,fontsize=8);ax.grid(alpha=.15)
 for ext in ['pdf','png']:fig.savefig(R/'manuscript/figures'/f'reach_budget.{ext}',dpi=220)
 plt.close(fig)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--cap',type=int,choices=[4,12,28]);ap.add_argument('--summarize',action='store_true');args=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 if args.cap:run(args.cap)
 if args.summarize:summarize()
