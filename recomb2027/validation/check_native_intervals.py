"""Recompute every native score spectrum using the older cut compiler.

No ExpressionEngine, PinClosure, or target-spectrum implementation is imported.
Target owners are exhaustively enumerated; all collateral owners are optimized
by an independent min-cut path. Geometry primitives and ingestion remain shared.
"""
from pathlib import Path
import sys,json,itertools,time
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'core')]
from mixed import make_problem,verify_mixed
from data_api import variant
from pinned_compiler import compile_pinned,solve_pinned
from layered import realize
from check_flow import check_certificate

def main():
    t0=time.perf_counter();charts={c['id']:c for c in json.loads((ROOT/'data/charts.json').read_text())}
    source=json.loads((ROOT/'data/observed_crop_counts.json').read_text())
    queries=cuts=patterns=bounds=boundaries=proofs=0
    for path in sorted((ROOT/'results/native').glob('*.json')):
        rec=json.loads(path.read_text());ch=charts[rec['input_chart']];p,pins,n=make_problem(ch)
        if rec['summary']['setting']=='maximal':p=variant(p)
        c=compile_pinned(p,pins);cache={};wc={}
        for q in rec['queries']:
            target=q['local_indices'];k=q['cell_index'];fixed=q['fixed_outside_count'];truth={}
            total=source['rois'][ch['roi']]['counts'].get(str(q['cell_id']),{}).get(q['gene'],0)
            assert fixed==total-sum(ch['points'][i]['source']==k for i in target)
            for own in itertools.product(range(p.m+1),repeat=len(target)):
                pp=tuple(zip(target,own));patterns+=1
                if pp not in cache:
                    cache[pp]=solve_pinned(c,[[0]*(p.m+1) for _ in range(p.n)],area_weight=1,pins=dict(pp),reconstruct=False,certificate=False);cuts+=1
                a=cache[pp]
                if not a['feasible']:continue
                val=fixed+sum(o==k for o in own);truth[val]=min(a['area'],truth.get(val,a['area']))
                owners=tuple(a['owners'])
                if owners not in wc:
                    w=realize(p,list(owners));assert w['feasible'] and w['area']==a['area'];verify_mixed(p,list(owners),w['boundaries'],w['area']);wc[owners]=w;boundaries+=1
            assert truth=={int(k):F(v) for k,v in q['spectrum']['minima'].items()},(ch['id'],q['gene'],k)
            for v,w in q['spectrum']['witnesses'].items():
                assert fixed+sum(w['owners'][i]==k for i in target)==int(v)
                verify_mixed(p,w['owners'],[[tuple(map(F,s)) for s in bb] for bb in w['boundaries']],F(w['area']))
                assert F(w['area'])==truth[int(v)]
            for a in q['dual_queries']:
                if 'certificate' in a:check_certificate(a['certificate']);proofs+=1
            for r in rec['rows']:
                if r['gene']!=q['gene'] or r['cell_index']!=k:continue
                vals=[v for v,d in truth.items() if d<=F(r['budget_um2'])]
                assert (min(vals),max(vals))==(r['lower'],r['upper'])
                assert r['dual_lower']<=min(vals)<=max(vals)<=r['dual_upper'];bounds+=1
            queries+=1
    out={'models':84,'independently_recomputed_score_spectra':queries,'full_target_patterns':patterns,
         'unique_pinned_cut_queries':cuts,'verified_interval_rows':bounds,'unique_minimum_boundaries_reconstructed':boundaries,
         'stored_root_flow_proofs_replayed':proofs,'failures':0,'seconds':time.perf_counter()-t0,
         'independence':'different optimizer and target enumeration; shared geometry primitives and source ingestion; internal code-path validation'}
    (ROOT/'results/native_interval_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
