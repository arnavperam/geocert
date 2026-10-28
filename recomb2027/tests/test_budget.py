"""Compare count certificates with exhaustive full-owner geometry, plus a real gap."""
from pathlib import Path
import sys,json,itertools,random,time
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine,convolve_score_tables,ScoreSpectrum
from layered import Layers,Box,realize
from mixed import verify_mixed

def model(seed):
    r=random.Random(seed);m=1+seed%2;n=6 if m==1 else 4
    boxes=[];source=[];choices=[F(-2,5),F(2,5)] if m==1 else [F(-2,5),F(2,5),F(8,5),F(12,5)]
    for g in range(n):
        x=r.choice(choices);t=F(r.randrange(8,73),10)
        boxes.append(Box(x,x,t,t,g));source.append(sum(x>2*a for a in range(m)))
    p=Layers([[(0,2*a),(8,2*a)] for a in range(m)],boxes,source,[F(1),F(3,4)][:m],
             [F(1,2)]*(m-1),[F(7,2)]*(m-1),clearance=F(1,20),
             bands=[(2*a-1,2*a+1) for a in range(m)],mode='nominal')
    return p,r

def main():
    t0=time.perf_counter();owner_checks=bounds=costs=cuts=incomplete=0
    for seed in range(1000,1060):
        p,r=model(seed);engine=ExpressionEngine(p);rows=[[0]*(p.m+1) for _ in range(p.n)]
        for g in r.sample(range(p.n),2):rows[g]=[r.randrange(-2,4) for _ in range(p.m+1)]
        const=3;truth={}
        for owners in itertools.product(range(p.m+1),repeat=p.n):
            owner_checks+=1;ans=realize(p,list(owners))
            if ans['feasible']:
                v=const+sum(row[o] for row,o in zip(rows,owners));d=ans['area'];truth[v]=min(d,truth.get(v,d))
        sp=engine.spectrum(rows,const);assert sp.complete and sp.minima==truth
        costs+=len(truth);qs=engine.dual_queries(rows,[0,F(1,4),1,4,16]);cuts+=len(qs)
        for budget in sorted(set([F(0),F(1,100),F(1,4),F(1),F(4)]+list(truth.values()))):
            vals=[v for v,d in truth.items() if d<=budget];exact=sp.interval(budget)
            enc=engine.enclosure(rows,budget,const,queries=qs)
            assert enc['lower']<=min(vals)<=max(vals)<=enc['upper']
            assert (exact['lower'],exact['upper'])==(min(vals),max(vals))
            assert enc['minimum_witness_value']>=min(vals) and enc['maximum_witness_value']<=max(vals)
            for key in ['minimum_witness','maximum_witness']:
                w=enc[key];verify_mixed(p,w['owners'],w['boundaries'],F(w['area']));assert F(w['area'])<=budget
            bounds+=1
        part=engine.spectrum(rows,const,pattern_limit=1)
        if part.total_patterns>1:
            assert not part.complete;incomplete+=1
            try:part.interval(1)
            except ValueError:pass
            else:raise AssertionError('partial spectrum accepted')
    boxes=[Box(F(19,10),F(19,10),5,5,0),Box(F(29,10),F(29,10),15,15,1),Box(F(29,10),F(29,10),25,25,2)]
    p=Layers([[(0,0),(30,0)]],boxes,[1,1,1],1,[],clearance=F(1,10),bands=[(-4,4)],mode='nominal')
    e=ExpressionEngine(p);scores=[[2,0],[3,0],[3,0]];sp=e.spectrum(scores)
    enc=e.enclosure(scores,10,multipliers=[0,F(1,3),F(1,2),1,4])
    assert sp.interval(10)['upper']==3 and enc['upper']==4
    exact=e.enclosure(scores,10,multipliers=[0,F(1,3)],refine=True);assert exact['upper']==3
    (ROOT/'results/dual_gap.json').write_text(json.dumps({'spectrum':sp.minima,'budget':'10','root_outer':[enc['lower'],enc['upper']],
        'exact_interval':[exact['lower'],exact['upper']],'queries':e.dual_queries(scores,[0,F(1,3)])},default=str,indent=2)+'\n')
    parts=[ScoreSpectrum({0:F(0),1:F(3,5),2:F(7,8)},{},True,0,0,0),ScoreSpectrum({0:F(0),1:F(1,2),3:F(5,4)},{},True,0,0,0)]
    combined=convolve_score_tables(parts);cart={}
    for (v,d),(w,a) in itertools.product(parts[0].minima.items(),parts[1].minima.items()):cart[v+w]=min(cart.get(v+w,d+a),d+a)
    assert combined==cart
    out={'models':60,'full_owner_vectors_checked':owner_checks,'exact_score_minima_checked':costs,'budget_enclosures_checked':bounds,
         'penalized_exact_cuts':cuts,'incomplete_spectra_rejected':incomplete,'deliberate_gap':{'exact_maximum':3,'penalty_bound':4},
         'composition_check':True,'failures':0,'seconds':time.perf_counter()-t0}
    (ROOT/'results/budget_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
