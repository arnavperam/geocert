"""Locked synthetic calibration/test split, with known generated owner truth.

The generated curves are paired measurements only in simulation; this does not
replace membrane annotation, real segmenters, independent donors or optical noise.
"""
from pathlib import Path
import sys,json,random,time,statistics,collections,math
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from geocert_expression.calibration import threshold,wilson
from layered import Layers,Box,realize,value

def sample(rng):
    base=[(F(0),F(0)),(F(10),F(0))]
    truth=[(F(0),F(0))]+[(F(t),F(rng.randrange(-70,71),100)) for t in [2,5,8]]+[(F(10),F(0))]
    pts=[];source=[];owners=[]
    while len(pts)<6:
        t=F(rng.randrange(50,951),100);x=F(rng.randrange(-100,101),100)
        if abs(x)<F(2,100) or abs(x-value(truth,t))<F(2,100):continue
        pts.append((x,t));source.append(int(x>0));owners.append(int(x>value(truth,t)))
    return {'base':base,'truth':truth,'points':pts,'source':source,'owners':owners}

def build(row,slope,band):
    boxes=[Box(x,x,t,t,i) for i,(x,t) in enumerate(row['points'])]
    return Layers([row['base']],boxes,row['source'],slope,[],clearance=F(1,100),bands=[(-band,band)],mode='nominal')

def maxslope(curve):return max(abs((v-u)/(b-a)) for (a,u),(b,v) in zip(curve,curve[1:]))

def auc(scores,labels):
    positive=[s for s,y in zip(scores,labels) if y];negative=[s for s,y in zip(scores,labels) if not y]
    if not positive or not negative:return None
    return sum((a>b)+F(1,2)*(a==b) for a in positive for b in negative)/len(positive)/len(negative)

def main():
    t0=time.perf_counter();protocol=json.loads((ROOT/'docs/EXPERIMENT_PROTOCOL.json').read_text())['synthetic_calibration']
    rng=random.Random(protocol['training_seed']);train=[sample(rng) for _ in range(protocol['training_models'])]
    # Fit family from training curves only. Rank calibration then separately
    # handles family misses: no failed calibration unit is dropped.
    slope=max(maxslope(r['truth']) for r in train)
    band=max(abs(h) for r in train for _,h in r['truth'])
    calrng=random.Random(protocol['calibration_seed']);cal=[]
    for _ in range(protocol['calibration_models']):
        row=sample(calrng);p=build(row,slope,band);ans=realize(p,row['owners'])
        cal.append(ans['area'] if ans['feasible'] else None)
    q=threshold(cal,F(protocol['alpha']));delta=None if q['threshold'] is None else F(q['threshold'])
    rng=random.Random(protocol['test_seed']);records=[];failures=0;ownercover=countcover=sourcecover=0;widths=[]
    allflags=[];allfragility=[];alldistance=[]
    for unit in range(protocol['test_models']):
        row=sample(rng);p=build(row,slope,band);truth=realize(p,row['owners']);engine=ExpressionEngine(p)
        contains=delta is None or (truth['feasible'] and truth['area']<=delta)
        ownercover+=contains;failures+=not truth['feasible'];intervals=[];allinside=True;sourcecorrect=True
        for gene in range(3):
            inds=[i for i in range(6) if i%3==gene]
            for cell in range(2):
                scores=[[int(g in inds and c==cell) for c in range(2)] for g in range(6)]
                sp=engine.spectrum(scores);d=max(sp.minima.values()) if delta is None else delta
                ex=sp.interval(d) if delta is not None else {'lower':0,'upper':len(inds)};v=sum(row['owners'][i]==cell for i in inds);v0=sum(row['source'][i]==cell for i in inds)
                allinside &= ex['lower']<=v<=ex['upper'];sourcecorrect &= v==v0;widths.append(ex['upper']-ex['lower'])
                intervals.append({'gene':gene,'cell':cell,'lower':ex['lower'],'upper':ex['upper'],'truth':v,'source':v0})
        assert not contains or allinside
        countcover+=allinside;sourcecover+=sourcecorrect
        for g,(x,t) in enumerate(row['points']):
            w=engine.clipped({g:1-row['source'][g]})
            fragility=-w['area'] if w['feasible'] else F(-10**6)
            allflags.append(row['source'][g]!=row['owners'][g]);allfragility.append(fragility);alldistance.append(-abs(x))
        records.append({'unit':unit,'owner_truth_in_family':truth['feasible'],'truth_minimum_deformation':str(truth['area']) if truth['feasible'] else None,
                        'owner_truth_within_budget':contains,'all_counts_covered':allinside,'counts':intervals})
    # Explicit misspecification guard test: too many incompatible annotations
    # yield an unbounded calibration threshold instead of a false finite radius.
    rejected=threshold([None]*20+[F(0)]*79,F(1,10));assert rejected['status']=='unbounded'
    out={'type':'simulation-only paired-boundary calibration','training_units':len(train),'family_fit':{'slope':str(slope),'band':str(band)},
         'calibration':q,'calibration_scores':[str(x) if x is not None else None for x in cal],
         'test_units':len(records),'test_family_failures':failures,'owner_coverage':ownercover,'simultaneous_count_coverage':countcover,
         'count_coverage_fraction':countcover/len(records),'count_coverage_wilson95':wilson(countcover,len(records)),
         'source_point_matrix_correct':sourcecover,'mean_interval_width':statistics.mean(widths),
         'pooled_misassignment_auc':{'exact_fragility':float(auc(allfragility,allflags)),'distance_to_flat_source':float(auc(alldistance,allflags)),
                                   'interpretation':'descriptive correlated transcript-level scores in simulation, not clinical/real-platform calibration'},
         'family_miss_guard':rejected,'records':records,'seconds':time.perf_counter()-t0}
    (ROOT/'results/synthetic_calibration.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['records','calibration_scores']},indent=2))
if __name__=='__main__':main()
