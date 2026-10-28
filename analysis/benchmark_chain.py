"""Exact full count spectra on a declared equal-tent family; no owner enumeration."""
from pathlib import Path
import sys,json,argparse,time,resource
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from atlas.chain_spectrum import compile_chain,spectrum
from atlas.model import serialized
from mixed import MixedLayers,TaggedBox
ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
n=args.n;rho=F(1,1000);h=F(1,4);d=F(1,10);T=2+(n-1)*d
p=MixedLayers([[(0,0),(T,0)]],[TaggedBox(h-rho,h-rho,1+i*d,1+i*d,i,False) for i in range(n)],
 [1]*n,[1],[],[],clearance=rho,bands=[(-F(3,10),F(3,10))],mode='nominal')
t0=time.perf_counter()
try:
 mesh=compile_chain(p,seconds=60,candidate_cap=300000)
 out=spectrum(mesh,[[1,0] for _ in range(n)],seconds=90,state_cap=1000000)
 expected={0:F(0),**{k:h*h+(k-1)*(h*d-d*d/4) for k in range(1,n+1)}}
 assert out['minima']==expected,(n,out['minima'],expected)
 result={'n':n,'status':'exact','count_levels':len(out['minima']),'all_analytical_costs_match':True,
 'compile_seconds':mesh.stats['compile_seconds'],'solve_and_check_seconds':out['seconds'],
 'total_seconds':time.perf_counter()-t0,'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
 'mesh':mesh.stats,'DP_states':out['DP_states'],'transition_evaluations':out['feasible_transition_evaluations'],
 'minimum_area_all_changed':str(out['minima'][n]),'witnesses_checked':len(out['witnesses'])}
except (MemoryError,TimeoutError) as e:result={'n':n,'status':'resource_limit','message':str(e),'seconds':time.perf_counter()-t0}
Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(result)
