"""Verify external witness accounting and a hash-defined exact height subset."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'core'))
from data_api import variant,reassignment_costs
from mixed import make_problem,verify_mixed
from layered import realize
from scc_solver import prepare,solve_compiled

if not __debug__:raise RuntimeError('Assertions required.')
charts={c['id']:c for c in json.loads((ROOT/'results/external/charts.json').read_text())}
d=json.loads((ROOT/'results/external/analysis.json').read_text());assert d['complete']
models={};cache={};witness_records=0;height_checks=0
for cid,c in charts.items():
 p,pins,n=make_problem(c)
 for setting,q in [('baseline',p),('maximal',variant(p))]:models[cid,setting]=(q,pins,n)

def witness(cid,setting,owners,area):
 global witness_records
 p,pins,n=models[cid,setting];labels=list(p.source);labels[:n]=owners
 assert len(owners)==n
 for g,o in pins.items():assert labels[g]==o
 key=(cid,setting,tuple(labels))
 if key not in cache:
  r=realize(p,labels);assert r['feasible'];verify_mixed(p,labels,r['boundaries'],r['area']);cache[key]=r['area']
 assert cache[key]==F(area);witness_records+=1
 return labels

hc={}
for row in d['optional']:
 cid,setting=row['chart'],row['setting'];p,pins,n=models[cid,setting]
 for prefix in ['joint','additive']:
  owners=witness(cid,setting,row[prefix+'_owners'],row[prefix+'_area'])
  gain=sum(owners[g]!=p.source[g] for g in range(n));assert gain==row[prefix+'_changed']
  objective=F(row['penalty'])*F(row[prefix+'_area'])-gain
  if prefix=='joint':assert objective==F(row['joint_objective'])
  else:assert objective-F(row['joint_objective'])==F(row['regret'])
 if int(hashlib.sha256(cid.encode()).hexdigest(),16)%5==0:
  if (cid,setting) not in hc:hc[cid,setting]=prepare(p,pins,'height',True)
  c,plan,_=hc[cid,setting]
  a=solve_compiled(c,'height',reassignment_costs(p,n),F(row['penalty']),plan)
  assert a['feasible'] and a['objective']==F(row['joint_objective']);height_checks+=1
for row in d['readouts']:
 ch=charts[row['chart']];k=row['cell_index'];gene=row['gene']
 assert ch['cell_labels'][k]==row['cell_label']
 assert row['source_count']==row['fixed_count']+sum(z['source']==k and z['gene']==gene for z in ch['points'])
 for endpoint in ['minimum','maximum']:
  a=row[endpoint];owners=witness(row['chart'],row['setting'],a['owners'],a['area'])
  assert a['count']==row['fixed_count']+sum(owners[g]==k and z['gene']==gene for g,z in enumerate(ch['points']))
 assert row['minimum']['count']<=row['source_count']<=row['maximum']['count']
summary={'optional_records':len(d['optional']),'readouts':len(d['readouts']),
         'witness_records':witness_records,'unique_model_witnesses':len(cache),
         'hash_selected_height_objective_checks':height_checks,'failures':0}
(ROOT/'results/external/verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
