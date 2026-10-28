"""Check stored width spectra, complete composition and interval semantics.

All global compositions and rows are checked. A deterministic hash-selected subset
of at most 32 multicopy local questions per setting (2--6 target copies) is
independently enumerated; untested local minima are not called independent truth.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import csv,hashlib,json,sys,time,signal
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from atlas.scalar_reference import spectrum as scalar

def main():
 out=[]
 for ds in ('lung','breast'):
  for radius in ('0p5','1','2','3'):
   base=R/'reach'/ds/('r'+radius);t=time.time()
   aa={a['roi']:a for f in sorted((base/'data').glob('atlas_*.json')) for a in [json.loads(f.read_text())]}
   pp={(roi,p['id']):p for roi,a in aa.items() for p in a['patches']}
   local=json.loads((base/'results/patch_spectra.json').read_text());li={(r['roi'],r['patch'],r['gene']):r for r in local}
   assert len(li)==len(local)
   sample=sorted([r for r in local if 2<=r['copies']<=6],key=lambda r:hashlib.sha256(f"{r['roi']}:{r['patch']}:{r['gene']}".encode()).hexdigest())[:32]
   checked=patterns=limits=0
   def alarm(*_):raise TimeoutError('independent scalar check exceeded 15 seconds')
   old=signal.signal(signal.SIGALRM,alarm)
   for r in sample:
    signal.alarm(15)
    try:
     truth,n=scalar(pp[(r['roi'],r['patch'])],r['gene']);assert truth=={int(v):F(d) for v,d in r['minima'].items()}
     checked+=1;patterns+=n
    except TimeoutError:limits+=1
    finally:signal.alarm(0)
   signal.signal(signal.SIGALRM,old)
   gs=json.loads((base/'results/cell_spectra.json').read_text());index={};globalchecked=tracechecked=0
   for r in gs:
    key=r['roi'],r['cell'],r['gene'];assert key not in index;index[key]=r
    # A separately written pairwise min-plus recurrence uses recorded local spectra.
    values={r['source_count']:F(0)}
    for pid in r['patches']:
     q=li[(r['roi'],pid,r['gene'])];p=pp[(r['roi'],pid)];sgn=1 if p['cell_labels'][0]==r['cell'] else -1
     nxt={}
     for k,v in q['minima'].items():
      delta=sgn*(int(k)-q['source_left'])
      for n,cost in values.items():
       u=n+delta;d=cost+F(v)
       if u not in nxt or d<nxt[u]:nxt[u]=d
     values=nxt
    assert values=={int(v):F(d) for v,d in r['minima'].items()};globalchecked+=1
    assert aa[r['roi']]['source_counts'].get(str(r['cell']),{}).get(r['gene'],0)==r['source_count']
    for value,trace in r['witness_traces'].items():
     assert int(value)==r['source_count']+sum(trace.values())
     cost=F(0)
     for pid,delta in trace.items():
      q=li[(r['roi'],pid,r['gene'])];p=pp[(r['roi'],pid)];sgn=1 if p['cell_labels'][0]==r['cell'] else -1
      lc=q['source_left']+sgn*delta;w=q['witnesses'][str(lc)]
      assert sum(x['gene']==r['gene'] and owner==0 for x,owner in zip(p['points'],w['owners']))==lc
      assert F(w['physical_area_um2'])==F(q['minima'][str(lc)])
      cost+=F(w['physical_area_um2'])
     assert cost==F(r['minima'][value]);tracechecked+=1
   seen=set();rowcount=0
   with (base/'results/cell_expression_intervals.csv').open() as f:
    for r in csv.DictReader(f):
     key=r['roi'],int(r['cell']),r['gene'];d=F(r['budget_um2']);ident=(*key,d);assert ident not in seen;seen.add(ident)
     s=index[key];v=[int(k) for k,c in s['minima'].items() if F(c)<=d]
     assert r['status']=='exact' and (int(r['lower']),int(r['upper']))==(min(v),max(v))
     assert int(r['independent_lower'])<=min(v)<=max(v)<=int(r['independent_upper']);rowcount+=1
   z={'dataset':ds,'radius':radius.replace('p','.'),'all_local_spectra':len(local),'independent_local_sample':checked,'independent_target_patterns':patterns,'sample_resource_limits':limits,'global_spectra_checked':globalchecked,'global_witness_traces_checked':tracechecked,'interval_rows_checked':rowcount,'failures':0,'seconds':time.time()-t}
   out.append(z);print(z,flush=True)
 (R/'results/reach/output_validation.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
