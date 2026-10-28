"""All multi-copy local gene queries by the polynomial count-chain algorithm.

Prepared physical models are unchanged. Single-copy queries retain their cheaper
pin-closure path. This check compares complete exact spectra, not selected bounds.
"""
from pathlib import Path
import sys,json,time,os,hashlib
from collections import defaultdict
from fractions import Fraction as F
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
D=Path(os.environ.get('GEOCERT_DATA_ROOT',R))
from atlas.chain_spectrum import compile_chain,spectrum
from atlas.model import make_patch,count_scores
want=defaultdict(list)
for q in json.loads((D/'results/patch_spectra.json').read_text()):
 if q['copies']>1:want[(q['roi'],q['patch'])].append(q)
records=[];start=time.time();matched=fail=0
for path in sorted((D/'data').glob('atlas_*.json')):
 atlas=json.loads(path.read_text())
 for p in atlas['patches']:
  targets=want.get((atlas['roi'],p['id']))
  if not targets:continue
  engine=make_patch(p);rec={'roi':atlas['roi'],'patch':p['id'],'groups':len(p['points']),'queries':[]}
  try:
   mesh=compile_chain(engine.problem,engine.pins,seconds=30,candidate_cap=200000)
   rec['mesh']=mesh.stats
   for q in targets:
    sc=count_scores(p,engine,p['cell_labels'][0],q['gene'])
    ans=spectrum(mesh,sc,seconds=30,state_cap=500000)
    physical={str(k):str(d*F(p['area_scale_um2'])) for k,d in ans['minima'].items()}
    assert physical==q['minima'],(p['id'],q['gene'],physical,q['minima'])
    rec['queries'].append({'gene':q['gene'],'copies':q['copies'],'count_levels':len(physical),'exact_match':True,'DP_states':ans['DP_states'],'solve_seconds':ans['seconds']})
    matched+=1
   rec['status']='complete'
  except (MemoryError,TimeoutError) as e:
   rec['status']='resource_limit';rec['message']=str(e);fail+=1
  records.append(rec)
  if len(records)%25==0:print(len(records),matched,'seconds',time.time()-start,flush=True)
(D/'results/chain_native_records.json').write_text(json.dumps(records,indent=2)+'\n')
summary={'eligible_multi_copy_queries':sum(map(len,want.values())),'eligible_patches':len(want),'completed_queries':matched,
 'resource_limited_patches':fail,'failures':0,'seconds':time.time()-start,'records':len(records),
 'scope':'Each completed full local spectrum matches separately enumerated native truth; no timing superiority claimed.'}
(D/'results/chain_native_validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
