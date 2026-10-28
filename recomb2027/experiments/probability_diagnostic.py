"""Profile-Gibbs and perturb-and-MAP diagnostic on a fully enumerated toy model.

No posterior claim is made for a continuous membrane prior. Unary perturbations
are approximate; finite-precision Gumbel draws are not rational probability proofs.
"""
from pathlib import Path
import sys,json,math,random,collections,itertools,time
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from layered import Layers,Box,realize
from pinned_compiler import solve_pinned

def gumbel(r):return -math.log(-math.log(max(1e-15,min(1-1e-15,r.random()))))

def main():
    t0=time.perf_counter()
    p=Layers([[(0,0),(3,0)]],[Box(F(49,100),F(49,100),1,1,0),Box(F(49,100),F(49,100),F(26,25),F(26,25),1)],
             [1,1],1,[],clearance=F(1,100),bands=[(-F(3,4),F(3,4))],mode='nominal')
    engine=ExpressionEngine(p);lam=F(20);reward=F(337,125)
    states=[];energies=[]
    for owners in itertools.product(range(2),repeat=2):
        a=realize(p,list(owners));assert a['feasible'];states.append(list(owners));energies.append(lam*a['area']-reward*sum(o==0 for o in owners))
    weights=[math.exp(-float(e)) for e in energies];z=sum(weights);prob=[w/z for w in weights]
    runs=[];N=2500
    for seed in [1701,1702,1703]:
        rng=random.Random(seed);full=collections.Counter();unary=collections.Counter()
        for _ in range(N):
            j=min(range(4),key=lambda j:float(energies[j])-gumbel(rng));full[j]+=1
            costs=[[(-reward if o==0 else F(0))-F(f'{gumbel(rng):.12f}') for o in range(2)] for g in range(2)]
            ans=solve_pinned(engine.compiled,costs,area_weight=lam,certificate=False,reconstruct=False)
            unary[states.index(ans['owners'])]+=1
        f=[full[i]/N for i in range(4)];u=[unary[i]/N for i in range(4)]
        runs.append({'seed':seed,'draws':N,'configuration_Gumbel_frequency':f,'unary_Gumbel_cut_frequency':u,
                     'configuration_TV':sum(abs(a-b) for a,b in zip(f,prob))/2,'unary_TV':sum(abs(a-b) for a,b in zip(u,prob))/2})
    out={'model':'finite ownership profile Gibbs, probabilities proportional to exp(-lambda*minimum_deformation + unary reward)',
         'not_continuous_boundary_posterior':True,'states':states,'energies':list(map(str,energies)),
         'enumerated_profile_probabilities':prob,'runs':runs,'seconds':time.perf_counter()-t0,
         'interpretation':'Configuration-level Gumbel is the standard exact-distribution identity in ideal arithmetic; unary Gumbel perturbations generally have a different distribution. Observed frequencies are numerical Monte Carlo diagnostics.'}
    (ROOT/'results/probability_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
