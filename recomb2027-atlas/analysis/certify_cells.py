"""All complete-reference-cell count spectra under one common atlas per crop."""
from pathlib import Path
import sys,json,csv,hashlib,time
from fractions import Fraction as F
from collections import defaultdict,Counter
from types import SimpleNamespace
import os
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
from atlas.chain_spectrum import compile_chain,spectrum as chain_spectrum
from atlas.binding import verify_prepared_atlas
from atlas.model import make_patch,count_scores,serialized,verify_witness,compose,interval
BUDGETS=list(map(F,['0','0.05','0.1','0.25','0.5','1','2','4']))

def main():
 t0=time.time();allrows=[];summaries=[];allstates=[];allpatches=[]
 method=os.environ.get('GEOCERT_COUNT_METHOD','enumeration')
 if method not in ('chain','enumeration'):raise ValueError('unknown count method')
 destination=Path(os.environ.get('GEOCERT_RESULT_DIR',ROOT/'results'));destination.mkdir(parents=True,exist_ok=True)
 audit_path=ROOT/'results/atlas_geometric_validation.json'
 if not audit_path.is_file():raise ValueError('run analysis/verify_atlas.py before certifying a new atlas')
 audit=json.loads(audit_path.read_text())
 atlases=[json.loads(path.read_text()) for path in sorted((ROOT/'data').glob('atlas_*.json'))]
 if not atlases:raise ValueError('no prepared atlases supplied')
 for a in atlases:verify_prepared_atlas(a,ROOT/'data/source_crops'/a['roi'],audit)
 for a in atlases:
  roi=a['roi'];parts=defaultdict(list);singles=defaultdict(list);local_records=[]
  counts=a['source_counts']; examined=0;calls=0;maxlocal=0;nontrivlocal=0;chainqueries=0;chainstates=0
  for j,p in enumerate(a['patches']):
   if not p['points']:continue
   e=make_patch(p);scale=F(p['area_scale_um2']);source=e.problem.source
   # Singleton minima are a strong marginal envelope; they do not compose a budget.
   for i,r in enumerate(p['points']):
    ans=e.clipped({i:1-r['source']});calls+=1
    if ans['feasible']:
     cost=F(ans['area'])*scale
     verify_witness(p,ans['owners'],ans['boundaries'],cost)
     for k in p['cell_labels']:
      if k>0:singles[(k,r['gene'])].append({'id':r['transcript_id'],'source':p['cell_labels'][r['source']],'cost':str(cost)})
   genehist=Counter(r['gene'] for r in p['points']);maxlocal=max(maxlocal,max(genehist.values()))
   mesh=compile_chain(e.problem,e.pins,seconds=60) if method=='chain' and max(genehist.values())>1 else None
   for gene,n in sorted(genehist.items()):
    # A two-cell patch conserves its local gene total, so one spectrum serves both sides.
    scores=count_scores(p,e,p['cell_labels'][0],gene)
    if mesh is not None and n>1:
     ans=chain_spectrum(mesh,scores,seconds=60)
     for w in ans['witnesses'].values():w['pins']={}
     sp=SimpleNamespace(complete=ans['complete'],minima=ans['minima'],witnesses=ans['witnesses'],examined=0)
     chainqueries+=1;chainstates+=ans.get('DP_states',0)
    else:sp=e.spectrum(scores,pattern_limit=65536)
    examined+=sp.examined
    if not sp.complete:raise RuntimeError('exact target spectrum incomplete; no exact output written')
    lsource=sum(r['gene']==gene and r['source']==0 for r in p['points'])
    entry={'patch':p['id'],'gene':gene,'copies':n,'source_left':lsource,'examined':sp.examined,
           'minima':{str(k):str(v*scale) for k,v in sorted(sp.minima.items())},'witnesses':{}}
    nontrivlocal+=len(sp.minima)>1
    for value,w in sp.witnesses.items():
     cost=F(w['area'])*scale
     verify_witness(p,w['owners'],w['boundaries'],cost)
     entry['witnesses'][str(value)]={'physical_area_um2':str(cost),'owners':w['owners'],'boundaries':w['boundaries'],'pins':w['pins']}
    local_records.append(entry)
    for side,k in enumerate(p['cell_labels']):
     if k==0:continue
     delta={str((v-lsource)*(1 if side==0 else -1)):str(cost*scale) for v,cost in sp.minima.items()}
     parts[(k,gene)].append((p['id'],delta))
   if j%50==0:print(roi,j,'patches',len(local_records),'spectra',flush=True)
  positive={(int(k),g) for k,v in counts.items() if int(k) in a['complete_cells'] for g,n in v.items() if n>0}
  keys=positive|set(parts);states=[];rows=[]
  for k,g in sorted(keys):
   n=counts.get(str(k),{}).get(g,0);state=compose(parts[(k,g)],n)
   rec={'atlas_id':a['atlas_id'],'roi':roi,'cell':k,'gene':g,'source_count':n,'patches':[p for p,v in parts[(k,g)]],
        'minima':{str(v):str(d) for v,(d,_) in sorted(state.items())},'witness_traces':{str(v):trace for v,(d,trace) in sorted(state.items())}}
   states.append(rec)
   for budget in BUDGETS:
    q=interval(state,budget)
    low=n-sum(r['source']==k and F(r['cost'])<=budget for r in singles[(k,g)])
    high=n+sum(r['source']!=k and F(r['cost'])<=budget for r in singles[(k,g)])
    assert low<=q['lower']<=q['upper']<=high
    row={'atlas_id':a['atlas_id'],'roi':roi,'cell':k,'gene':g,'budget_um2':str(budget),
         'source_count':n,'lower':q['lower'],'upper':q['upper'],'status':q['status'],
         'mutable_patches':len(parts[(k,g)]),'independent_lower':low,'independent_upper':high,
         'lower_cost_um2':q['lower_cost'],'upper_cost_um2':q['upper_cost']}
    rows.append(row)
  summary={'roi':roi,'complete_cells':len(a['complete_cells']),'positive_or_editable_cell_gene_keys':len(keys),
           'editable_cell_gene_keys':len(parts.keys() & keys)-sum(not parts[k] for k in keys),
           'local_gene_spectra':len(local_records),'local_variable_gene_spectra':nontrivlocal,
           'count_method':method,'chain_queries':chainqueries,'chain_DP_states':chainstates,'local_target_assignments':examined,'singleton_questions':calls,'maximum_gene_copies_per_patch':maxlocal,
           'maximum_patches_per_cell_gene':max(map(lambda x:len(x['patches']),states)),
           'maximum_cell_gene_observed_count':max(r['source_count'] for r in states),
           'interval_rows':len(rows),'by_budget':[]}
  # Do not count dictionary entries created merely by querying a fixed gene as active.
  summary['editable_cell_gene_keys']=sum(bool(r['patches']) for r in states)
  for b in BUDGETS:
   rr=[r for r in rows if r['budget_um2']==str(b)]
   summary['by_budget'].append({'budget_um2':str(b),'rows':len(rr),'non_singleton':sum(r['lower']!=r['upper'] for r in rr),
         'strictly_tighter_than_independent':sum((r['lower'],r['upper'])!=(r['independent_lower'],r['independent_upper']) for r in rr),
         'max_width':max(r['upper']-r['lower'] for r in rr),'source_detection_can_disappear':sum(r['source_count']>0 and r['lower']==0 for r in rr),
         'source_zero_can_appear':sum(r['source_count']==0 and r['upper']>0 for r in rr)})
  allrows.extend(rows);allstates.extend(states);allpatches.extend({'roi':roi,**r} for r in local_records);summaries.append(summary)
  print(summary,flush=True)
 (destination/'cell_spectra.json').write_text(json.dumps(serialized(allstates),indent=2)+'\n')
 (destination/'patch_spectra.json').write_text(json.dumps(serialized(allpatches),indent=2)+'\n')
 with (destination/'cell_expression_intervals.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 if len({(r['roi'],r['cell'],r['gene'],r['budget_um2']) for r in allrows})!=len(allrows):raise AssertionError('duplicate cell-level key')
 result={'rois':summaries,'seconds':time.time()-t0,'unique_key_check':True,'source_scope':'All retained observed counts for complete reference cells; one geometry-only atlas per crop; omitted boundary portions fixed.',
    'record_semantics':'Marginal extrema under one common model and one global physical-area budget. Different endpoint rows need not coexist.'}
 (destination/'cell_expression_summary.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 if not __debug__:raise RuntimeError('Scientific checks require assertions: do not run Python with -O')
 main()
