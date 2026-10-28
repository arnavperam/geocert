"""Independent-spectrum samples and full global/interval checks of the width study."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import csv,json,time,sys,os
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from atlas.scalar_reference import spectrum as scalar_truth
from atlas.model import verify_witness

def check(d):
 t=time.time(); atlases={q['roi']:q for q in [json.loads(p.read_text()) for p in sorted((d/'data').glob('atlas_*.json'))]}
 patches={(roi,p['id']):p for roi,a in atlases.items() for p in a['patches']};local=json.loads((d/'results/patch_spectra.json').read_text())
 chosen=local if d.name=='lung' else sorted([q for q in local if q['copies']>1],key=lambda q:(q['roi'],q['patch'],q['gene']))[:64]+sorted([q for q in local if q['copies']==1],key=lambda q:(q['roi'],q['patch'],q['gene']))[:64]
 patterns=0
 for q in chosen:
  p=patches[q['roi'],q['patch']];true,n=scalar_truth(p,q['gene']);patterns+=n
  assert true=={int(k):F(v) for k,v in q['minima'].items()},(d,q['patch'],q['gene'])
 # Every saved local witness is checked against original model inequalities.
 lookup={};wcount=0
 for q in local:
  p=patches[q['roi'],q['patch']];lookup[q['roi'],q['patch'],q['gene']]=q
  for count,w in q['witnesses'].items():
   cost=verify_witness(p,w['owners'],w['boundaries'],w['physical_area_um2']);assert cost==F(q['minima'][count])
   assert sum(r['gene']==q['gene'] and o==0 for r,o in zip(p['points'],w['owners']))==int(count);wcount+=1
 # Independent dictionary convolution; no production compose() call.
 globals_=json.loads((d/'results/cell_spectra.json').read_text());out={};states=0
 for q in globals_:
  counts=atlases[q['roi']]['source_counts'];assert q['source_count']==counts.get(str(q['cell']),{}).get(q['gene'],0)
  state={q['source_count']:F(0)}
  for pid in q['patches']:
   p=patches[q['roi'],pid];x=lookup[q['roi'],pid,q['gene']];sgn=1 if q['cell']==p['cell_labels'][0] else -1
   new={}
   for n,cost in state.items():
    for value,area in x['minima'].items():
     v=n+sgn*(int(value)-x['source_left']);c=cost+F(area)
     if v not in new or c<new[v]:new[v]=c
   state=new;states+=len(new)
  assert state=={int(v):F(c) for v,c in q['minima'].items()};out[(q['roi'],str(q['cell']),q['gene'])]=state
 nr=0;seen=set()
 with (d/'results/cell_expression_intervals.csv').open() as f:
  for row in csv.DictReader(f):
   k=(row['roi'],row['cell'],row['gene']);budget=F(row['budget_um2']);ident=(*k,budget);assert ident not in seen;seen.add(ident)
   values=[v for v,c in out[k].items() if c<=budget]
   assert (min(values),max(values))==(int(row['lower']),int(row['upper'])) and row['status']=='exact';nr+=1
 record={'directory':str(d.relative_to(R)),'local_scalar_spectra':len(chosen),'scalar_target_patterns':patterns,'local_boundary_witnesses':wcount,'global_spectra':len(globals_),
  'independent_convolution_states':states,'interval_rows':nr,'failures':0,'seconds':time.time()-t,
  'scope':'All local witnesses and global minima/intervals; independent scalar truth on all lung local spectra and predetermined first 64 multicopy + 64 single-copy breast spectra. Shared prepared observations and geometric witness checker.'}
 (d/'results/reach_validation.json').write_text(json.dumps(record,indent=2)+'\n');print(record,flush=True);return record

def main():
 records=[]
 for width in ['cap4','cap8','cap16','cap20']:
  for ds in ['lung','breast']:records.append(check(R/'reach_study'/width/ds))
 (R/'results/reach_validation.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
