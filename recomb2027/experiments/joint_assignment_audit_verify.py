"""Independent min-cut replay of every saved pair, using integer flow certificates.

The audit obtains minimum area by implication closure and envelope reconstruction.
This replay instead compiles the event graph and minimizes area with pair pins.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,csv,time,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from mixed import make_problem
from data_api import variant
from pinned_compiler import solve_pinned
from check_flow import check_certificate

def main():
    d=ROOT/'results/joint_assignment_audit';rows=list(csv.DictReader((d/'pairs.csv').open()))
    charts={c['id']:c for c in json.loads((ROOT/'data/charts.json').read_text())}
    checks=[];t=time.perf_counter()
    for cid,setting in sorted({(r['corridor'],r['setting']) for r in rows}):
        p,pins,n=make_problem(charts[cid]);p=variant(p) if setting=='maximal' else p
        en=ExpressionEngine(p,pins);comp=en.compiled
        count=0
        for r in rows:
            if (r['corridor'],r['setting'])!=(cid,setting):continue
            pp={int(r['group_1']):int(r['owner_1']),int(r['group_2']):int(r['owner_2'])}
            ans=solve_pinned(comp,[[0]*(p.m+1) for _ in range(p.n)],pins=pp,area_weight=1,certificate=True,reconstruct=False)
            check_certificate(ans['certificate'])
            assert ans['feasible']==(r['joint_feasible']=='True')
            if ans['feasible']:assert F(ans['area'])==F(r['joint_area'])
            checks.append({'corridor':cid,'setting':setting,'pins':pp,'feasible':ans['feasible'],'area':str(ans.get('area','')),'certificate':ans['certificate']})
            count+=1
        print(cid,setting,count,flush=True)
    out={'method':'independent event-compiled integer-capacity min-cut replay with max-flow/min-cut equality and primal/dual certificate checks',
         'pairs_sha256':hashlib.sha256((d/'pairs.csv').read_bytes()).hexdigest(),
         'checked':len(checks),'failures':0,'seconds':time.perf_counter()-t}
    (d/'cut_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    (d/'cut_certificates.json').write_text(json.dumps(checks,default=str,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
