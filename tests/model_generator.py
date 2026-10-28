import sys,json,random
from pathlib import Path
from fractions import Fraction as F
from layered import Layers,Box

def dump(path,obj):
    def default(x):
        if isinstance(x,F):return str(x)
        raise TypeError(type(x))
    Path(path).write_text(json.dumps(obj,indent=2,default=default))

def make(seed,n=None,m=None):
    rng=random.Random(seed);m=m or rng.randint(1,3);n=n or rng.randint(1,4)
    bases=[]
    for k in range(m):
        base=2*k
        bases.append([(0,F(base)),(2,F(base)+F(rng.randint(-1,1),10)),(4,F(base)+F(rng.randint(-1,1),10)),(6,F(base))])
    gaps=[rng.choice([F(1,3),F(1),F(3,2)]) for _ in range(m-1)]
    maxg=[rng.choice([None,F(5,2),F(3)]) for _ in range(m-1)]
    src=[];boxes=[]
    for g in range(n):
        cell=rng.randrange(m+1);src.append(cell)
        x=F(-4,5) if cell==0 else F(2*(m-1))+F(4,5) if cell==m else F(2*(cell-1))+rng.choice([F(2,5),F(1),F(8,5)])
        t=F(rng.randint(5,55),10);ex=rng.choice([F(0),F(1,20)]);ey=rng.choice([F(0),F(1,10),F(1,5)])
        boxes.append(Box(x-ex,x+ex,t-ey,t+ey,g))
        if rng.random()<.2:boxes.append(Box(x-ex,x+ex,t-ey,t+ey+F(1,20),g))
    slopes=[rng.choice([F(1,4),F(1),F(2)]) for _ in range(m)] if seed%2 else F(1)
    p=Layers(bases,boxes,src,slopes,gaps,maxg,clearance=F(1,20),bands=[(-3,2*m+1)]*m)
    return p,rng
