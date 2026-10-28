"""Predeclared co-detection predicates under the same global atlas.

A marker conjunction is an observed-count query, not a false-double-positive
label, a cell annotation, differential expression, or evidence of signaling.
"""
from pathlib import Path
import sys,json,csv,time
from fractions import Fraction as F
from collections import Counter
import os
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
from atlas.model import make_patch,count_scores,verify_witness,serialized
from geocert_expression.joint import outcome_spectrum
BUDGETS=list(map(F,['0','0.05','0.1','0.25','0.5','1','2','4']))
PAIRS=[('EPCAM','PTPRC'),('EPCAM','CD3D'),('EPCAM','PECAM1')]

def convolve(state,local,pid):
 out={}
 for x,(d,tr) in state.items():
  for y,(e,idx) in local.items():
   key=tuple(a+b for a,b in zip(x,y));z=d+e
   if key not in out or z<out[key][0]:out[key]=(z,{**tr,pid:idx})
 return out

def main():
 t=time.time();rows=[];records=[];checks=0;local_queries=0
 for path in sorted((ROOT/'data').glob('atlas_*.json')):
  a=json.loads(path.read_text());counts=a['source_counts'];cache={}
  for cell in a['complete_cells']:
   incident=[p for p in a['patches'] if cell in p['cell_labels'] and p['points']]
   localgenes={r['gene'] for p in incident for r in p['points']}
   for pair in PAIRS:
    source=tuple(counts.get(str(cell),{}).get(g,0) for g in pair)
    if not any(source) and not localgenes.intersection(pair):continue
    state={source:(F(0),{})};components=[]
    for p in incident:
     if not any(r['gene'] in pair for r in p['points']):continue
     engine=cache.setdefault(p['id'],None)
     if engine is None:engine=make_patch(p);cache[p['id']]=engine
     ch=[count_scores(p,engine,cell,g) for g in pair]
     sp=outcome_spectrum(engine,ch,pattern_limit=65536);local_queries+=1
     if not sp['complete']:raise RuntimeError('incomplete joint spectrum')
     n=tuple(sum(r['gene']==g and p['cell_labels'][r['source']]==cell for r in p['points']) for g in pair)
     local={}
     for j,w in enumerate(sp['outcomes']):
      d=F(w['area'])*F(p['area_scale_um2']);verify_witness(p,w['owners'],w['boundaries'],d);checks+=1
      delta=tuple(v-u for v,u in zip(w['counts'],n));local[delta]=(d,j)
     state=convolve(state,local,p['id']);components.append({'patch':p['id'],'source':n,'spectrum':sp})
    records.append({'atlas_id':a['atlas_id'],'roi':a['roi'],'cell':cell,'genes':pair,'source':source,
                    'components':components,'outcomes':[{'counts':k,'physical_area_um2':str(d),'trace':tr} for k,(d,tr) in sorted(state.items())]})
    for b in BUDGETS:
     vals=[x for x,(d,_) in state.items() if d<=b];truth=[all(z>0 for z in x) for x in vals]
     status='always' if all(truth) else ('variable' if any(truth) else 'never')
     marginal_can_both=all(max(x[j] for x in vals)>0 for j in range(2))
     marginal_forces_both=all(min(x[j] for x in vals)>0 for j in range(2))
     rows.append({'atlas_id':a['atlas_id'],'roi':a['roi'],'cell':cell,'gene1':pair[0],'gene2':pair[1],
                  'source1':source[0],'source2':source[1],'budget_um2':str(b),'predicate':'both_counts_positive',
                  'status':status,'marginal_upper_allows_both':marginal_can_both,'marginal_lower_forces_both':marginal_forces_both,
                  'joint_exclusion_beyond_marginals':marginal_can_both and status=='never'})
 out={'cell_pair_questions':len(records),'budgeted_predicates':len(rows),'local_joint_queries':local_queries,'checked_witnesses':checks,
      'by_pair':{},'seconds':time.time()-t,'interpretation':'Joint observed co-detection under a declared global atlas; no ground-truth misassignment labels.'}
 for pair in PAIRS:
  rr=[r for r in rows if (r['gene1'],r['gene2'])==pair]
  out['by_pair']['/'.join(pair)]={'cell_questions':len(rr)//len(BUDGETS),'by_budget':[
       {'budget_um2':str(b),'statuses':dict(Counter(r['status'] for r in rr if r['budget_um2']==str(b))),
        'source_double_positive_can_be_removed':sum(r['source1']>0 and r['source2']>0 and r['status']!='always' for r in rr if r['budget_um2']==str(b)),
        'joint_exclusion_beyond_marginals':sum(r['joint_exclusion_beyond_marginals'] for r in rr if r['budget_um2']==str(b))} for b in BUDGETS]}
 (ROOT/'results/joint_marker_spectra.json').write_text(json.dumps(serialized(records),indent=2)+'\n')
 (ROOT/'results/joint_marker_summary.json').write_text(json.dumps(out,indent=2)+'\n')
 with (ROOT/'results/joint_marker_predicates.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print(out)
if __name__=='__main__':main()
