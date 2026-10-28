"""Exact cell-count sensitivity with fixed out-of-corridor contributions.

No segmentation or classification accuracy is inferred. Run with assertions enabled.
The exact compiler is unchanged. All optimization inputs are published rational charts.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, sys, time, hashlib
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'core'))
from data_api import load_cases, encode
from scc_solver import prepare, solve_compiled
from readout_inputs import build_inputs

def run(out):
    if not __debug__: raise RuntimeError('Assertions must remain enabled.')
    inputs=build_inputs()
    bychart={}
    for r in inputs['readouts']:bychart.setdefault(r['chart'],[]).append(r)
    results=[];penalized=[];reference_checks=[];t0=time.perf_counter()
    for ix,ch,setting,p,pins,n,expected in load_cases():
        if ch['id'] not in bychart: continue
        compiled,plan,stages=prepare(p,pins,'direct',True)
        hcache=None
        # Exact lexicographic scalarization: first integer count, then D.
        area_bound=(p.hi-p.lo)*sum(hi-lo for lo,hi in p.bands)
        dominant=int(area_bound)+2
        for r in bychart[ch['id']]:
            target=r['cell_index'];gene=r['gene'];ends=[]
            gene_inds=[g for g,z in enumerate(ch['points']) if z['gene']==gene]
            for sign in [1,-1]:
                costs=[[0]*(p.m+1) for _ in range(p.n)]
                for g in gene_inds:costs[g][target]=sign*dominant
                a=solve_compiled(compiled,'direct',costs,F(1),plan)
                assert a['feasible']
                count=r['fixed_count']+sum(a['owners'][g]==target for g in gene_inds)
                assert 0<=a['area']<=area_bound
                ends.append({'count':count,'area':str(a['area']),'owners':a['owners'][:n]})
            assert ends[0]['count']<=r['source_count']<=ends[1]['count']
            record={**r,'setting':setting,'minimum':ends[0],'maximum':ends[1]}
            results.append(record)
            # Outcome-blind, hash-defined exact competing representation check.
            key=f"{ch['id']}|{target}|{gene}"
            if int(hashlib.sha256(key.encode()).hexdigest(),16)%25==0:
                if hcache is None:hcache=prepare(p,pins,'height',True)
                hc,hp,_=hcache
                for direction,sign in enumerate([1,-1]):
                    costs=[[0]*(p.m+1) for _ in range(p.n)]
                    for g in gene_inds:costs[g][target]=sign*dominant
                    a=solve_compiled(hc,'height',costs,F(1),hp)
                    count=r['fixed_count']+sum(a['owners'][g]==target for g in gene_inds)
                    assert a['feasible'] and count==ends[direction]['count'] and str(a['area'])==ends[direction]['area']
                    reference_checks.append({'key':key,'setting':setting,'direction':direction,'count':count,'area':str(a['area'])})
            if setting=='baseline' and ends[0]['count']<ends[1]['count']:
                for lam in [F(1,4),F(1),F(4)]:
                    for sign in [1,-1]:
                        costs=[[0]*(p.m+1) for _ in range(p.n)]
                        for g in gene_inds:costs[g][target]=sign
                        a=solve_compiled(compiled,'direct',costs,lam,plan)
                        count=r['fixed_count']+sum(a['owners'][g]==target for g in gene_inds)
                        assert ends[0]['count']<=count<=ends[1]['count']
                        penalized.append({'chart':ch['id'],'cell':r['cell_label'],'gene':gene,'source_count':r['source_count'],'lambda':str(lam),'direction':'decrease' if sign==1 else 'increase','count':count,'area':str(a['area']),'objective_with_fixed':str(a['objective']+sign*r['fixed_count']),'owners':a['owners'][:n]})
        print(ix,setting,ch['id'],len(results),round(time.perf_counter()-t0,2),flush=True)
    payload={'protocol_sha256':inputs['protocol_sha256'],'rows':results,'penalized':penalized,'height_checks':reference_checks,'elapsed_seconds':time.perf_counter()-t0,'scope':'Existing development inputs; full-source-cell accounting with local alternatives and fixed out-of-corridor counts.'}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2,default=encode)+'\n')
    print('COMPLETE',len(results),len(penalized),len(reference_checks),flush=True)
if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=BASE/'run/readout_results.json');ar=pa.parse_args();run(ar.output)
