from pathlib import Path
import json,gzip,sys
from fractions import Fraction as F
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from atlas.model import verify_witness
from analysis.revision_reach import references
report=[]
for cap in [4,12,28]:
 D=R/f'results/revision_reach/cap_{cap}'
 if not (D/'summary.json').exists():continue
 count=0;localchecks=0
 for af in (D/'data').glob('atlas_*.json'):
  a=json.loads(af.read_text());ps={p['id']:p for p in a['patches']}
  for name,ref in references(a).items():
   with gzip.open(D/f'results/{a["roi"]}_{name}_witnesses.json.gz','rt') as f:locals=json.load(f)
   lm={(r['patch'],r['gene']):r for r in locals}
   for r in locals:
    p=ps[r['patch']]
    for v,w in r['witnesses'].items():
     verify_witness(p,w['owners'],w['boundaries'],w['area_um2'])
     got=sum(q['gene']==r['gene'] and q['transcript_id'] in ref and w['owners'][i]==0 for i,q in enumerate(p['points']))
     assert got==int(v);localchecks+=1
   qs=json.loads((D/f'results/{a["roi"]}_{name}_queries.json').read_text())
   for q in qs:
    for v,trace in q['witness_traces'].items():
     got=q['source_count'];cost=F(0)
     for pid,delta in trace.items():
      p=ps[pid];side=p['cell_labels'].index(q['cell']);ls=sum(r['gene']==q['gene'] and r['transcript_id'] in ref and r['source']==0 for r in p['points']);lv=ls+int(delta)*(1 if side==0 else -1)
      w=lm[pid,q['gene']]['witnesses'][str(lv)];cost+=F(w['area_um2']);got+=int(delta)
     assert got==int(v) and cost==F(q['minima'][v]);count+=1
 report.append({'cap':cap,'local_witness_counts_and_geometry':localchecks,'global_state_count_area_traces':count,'failures':0})
print(json.dumps(report,indent=2))
(R/'results/revision_reach/witness_validation.json').write_text(json.dumps(report,indent=2)+'\n')
