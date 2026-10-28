"""Current input interface; source-derived exact JSON, not a raw-data reader."""
from pathlib import Path
from fractions import Fraction as F
import json
from mixed import make_problem,MixedLayers
ROOT=Path(__file__).resolve().parents[1]

def variant(p,end=F(1,2),slope=2,band=2):
    extension=F(end)*(p.hi-p.lo)
    baselines=[([(b[0][0]-extension,b[0][1])]+b+[(b[-1][0]+extension,b[-1][1])]) if extension else list(b) for b in p.baselines]
    bands=[]
    for b,(lo,hi) in zip(p.baselines,p.bands):
        mn=min(v for _,v in b);mx=max(v for _,v in b)
        bands.append((mn-F(band)*(mn-lo),mx+F(band)*(hi-mx)))
    return MixedLayers(baselines,p.boxes,p.source,[F(slope)*s for s in p.slopes],p.minimum_gaps,p.maximum_gaps,clearance=p.clearance,bands=bands,mode=p.mode)

def load_cases():
    charts=json.loads((ROOT/'data/charts.json').read_text());expected=json.loads((ROOT/'data/expected_optional.json').read_text());out=[]
    for i,ch in enumerate(charts):
        p,pins,n=make_problem(ch)
        for j,(tag,q) in enumerate([('baseline',p),('maximal',variant(p))]):
            e=expected[2*i+j]
            if (e['id'],e['setting'])!=(ch['id'],tag):raise ValueError('input manifest order mismatch')
            out.append((i,ch,tag,q,pins,n,e))
    return out

def reassignment_costs(p,n):
    return [[0 if g>=n or c==p.source[g] else -1 for c in range(p.m+1)] for g in range(p.n)]

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x))
