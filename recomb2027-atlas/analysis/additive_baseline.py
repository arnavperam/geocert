"""Feasibility-matched additive singleton cost for every atlas count question.

This is an INNER budget family. It is an ablation of shared physical cost, not
an alternative exact solver for the same objective or a safe exclusion bound.
The same feasible target assignments and collateral ownership are retained.
"""
from pathlib import Path
import sys,os,json,csv,time
from fractions import Fraction as F
from collections import defaultdict,Counter
from itertools import product
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
from atlas.scalar_reference import prepare,query
from atlas.model import compose,interval

def main():
 t=time.time();localout=[];parts=defaultdict(list);pattern_checks=0
 localref={(r['roi'],r['patch'],r['gene']):r for r in json.loads((ROOT/'results/patch_spectra.json').read_text())}
 for af in sorted((ROOT/'data').glob('atlas_*.json')):
  a=json.loads(af.read_text())
  for p in a['patches']:
   if not p['points']:continue
   pre=prepare(p);singles={}
   for i,r in enumerate(p['points']):
    q=query(p,{i:1-r['source']},pre)
    if q is not None:singles[i]=q[1]
   for gene,n in sorted(Counter(r['gene'] for r in p['points']).items()):
    ids=[i for i,r in enumerate(p['points']) if r['gene']==gene];amin={};dmin={};source=sum(p['points'][i]['source']==0 for i in ids)
    for target in product([0,1],repeat=len(ids)):
     answer=query(p,dict(zip(ids,target)),pre)
     if answer is None:continue
     owners,cost=answer;add=sum((singles[i] for i,o in enumerate(owners) if o!=p['points'][i]['source']),F(0));value=sum(owners[i]==0 for i in ids)
     if add<cost:raise AssertionError('additive upper bound failed')
     amin[value]=min(amin.get(value,add),add);dmin[value]=min(dmin.get(value,cost),cost);pattern_checks+=1
    truth={int(k):F(v) for k,v in localref[(a['roi'],p['id'],gene)]['minima'].items()}
    if dmin!=truth:raise AssertionError('independent exact reference mismatch')
    localout.append({'roi':a['roi'],'patch':p['id'],'gene':gene,'source_left':source,
                     'exact_minima':{str(k):str(v) for k,v in dmin.items()},'additive_minima':{str(k):str(v) for k,v in amin.items()}})
    for side,cell in enumerate(p['cell_labels']):
     if cell>0:parts[(a['roi'],cell,gene)].append((p['id'],{str((v-source)*(1 if side==0 else -1)):d for v,d in amin.items()}))
 states={}
 for r in json.loads((ROOT/'results/cell_spectra.json').read_text()):
  key=(r['roi'],r['cell'],r['gene']);states[key]=compose(parts[key],r['source_count'])
 rows=[];by_budget=defaultdict(Counter)
 for row in csv.DictReader((ROOT/'results/cell_expression_intervals.csv').open()):
  key=(row['roi'],int(row['cell']),row['gene']);b=F(row['budget_um2']);pred=interval(states[key],b);lo,hi=int(row['lower']),int(row['upper'])
  assert lo<=pred['lower']<=pred['upper']<=hi
  missed=lo!=pred['lower'] or hi!=pred['upper']
  r={**row,'additive_inner_lower':pred['lower'],'additive_inner_upper':pred['upper'],'additive_misses_endpoint':missed};rows.append(r)
  by_budget[row['budget_um2']]['queries']+=1;by_budget[row['budget_um2']]['missed_endpoints']+=missed
  by_budget[row['budget_um2']]['strict_outer_tightening']+=(lo,hi)!=(int(row['independent_lower']),int(row['independent_upper']))
 out={'feasible_pattern_inequality_checks':pattern_checks,'local_count_spectra':len(localout),
      'cell_gene_questions':len(states),'budgeted_questions':len(rows),'rows_with_additive_missed_endpoint':sum(r['additive_misses_endpoint'] for r in rows),
      'affected_cell_gene_questions':len({(r['roi'],r['cell'],r['gene']) for r in rows if r['additive_misses_endpoint']}),
      'by_budget':[{'budget_um2':b,**by_budget[b]} for b in sorted(by_budget,key=F)],'seconds':time.time()-t,
      'interpretation':'Feasibility-matched INNER approximation replaces shared area by summed individually minimal costs for every changed RNA, including collateral changes. No deployed-segmenter or solver-performance claim.'}
 (ROOT/'results/additive_local_spectra.json').write_text(json.dumps(localout,indent=2)+'\n')
 (ROOT/'results/additive_summary.json').write_text(json.dumps(out,indent=2)+'\n')
 with (ROOT/'results/count_baseline_comparison.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print(out)
if __name__=='__main__':main()
