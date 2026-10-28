"""Nonflat sources, coincident RNA times, mandatory regions and exact score spectra."""
from pathlib import Path
import sys,json,random,time,itertools
from fractions import Fraction as F
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from atlas.chain_spectrum import compile_chain,spectrum
from mixed import MixedLayers,TaggedBox
from layered import realize,value
from atlas.model import serialized
rng=random.Random(4217);patterns=checks=witnesses=0;start=time.time()
for seed in range(80):
 n=1+seed%6;baseline=[(F(0),F(0)),(F(1),F(rng.randint(-3,3),10)),(F(2),F(rng.randint(-3,3),10)),(F(3),F(0))]
 rho=F(1,100);boxes=[];source=[]
 for g in range(n):
  t=F(rng.randint(1,11),4);x=value(baseline,t)+F(rng.choice([-1,1])*rng.randint(1,5),10)
  boxes.append(TaggedBox(x,x,t,t,g,False));source.append(int(x>value(baseline,t)))
 pins={}
 if seed%4==0:
  g=len(source);source.append(int(-F(1,5)>value(baseline,F(1))));pins[g]=0
  # May force positive source deformation; enumeration does not use pin clipping.
  boxes.append(TaggedBox(-F(1,5),-F(1,5),F(3,4),F(5,4),g,True))
 p=MixedLayers([baseline],boxes,source,[1],[],[],clearance=rho,bands=[(-1,1)],mode='nominal')
 scores=[[int(g%2==0),0] for g in range(n)]+[[0,0] for g in pins]
 truth={}
 for v in itertools.product((0,1),repeat=n):
  owners=list(v)+[0 for _ in pins];patterns+=1
  ans=realize(p,owners)
  if not ans['feasible']:continue
  area=ans['area'];q=sum(row[o] for row,o in zip(scores,owners))
  if q not in truth or area<truth[q]:truth[q]=area
 mesh=compile_chain(p,pins);got=spectrum(mesh,scores)
 assert got['minima']==truth,(seed,got['minima'],truth)
 checks+=1;witnesses+=len(got['witnesses'])
out={'models':checks,'exhaustive_owners':patterns,'returned_witnesses':witnesses,'failures':0,'seconds':time.time()-start,
'independence':'Fixed-owner analytic reconstruction is separate from critical-height construction and count-chain optimization; geometry definitions shared.'}
(R/'results/chain_nonflat_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
