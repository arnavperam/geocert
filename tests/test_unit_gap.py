"""A genuine hard-budget gap with eight distinct unit-count point RNAs."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,itertools
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from atlas.chain_spectrum import compile_chain,spectrum
from atlas.model import serialized
from mixed import MixedLayers,TaggedBox
from layered import realize
rho=F(1,10);epsilon=F(1,100);boxes=[];sites=[]
for center,height,copies in [(5,2,2),(15,3,3),(25,3,3)]:
 for j in range(copies):
  t=F(center)+j*epsilon;x=F(height)-j*epsilon-rho;g=len(boxes)
  boxes.append(TaggedBox(x,x,t,t,g,False));sites.append([str(t),str(x)])
p=MixedLayers([[(0,0),(30,0)]],boxes,[1]*8,[1],[],[],clearance=rho,bands=[(-4,4)],mode='nominal')
a=spectrum(compile_chain(p),[[1,0]]*8);expected={0:F(0),2:F(4),3:F(9),5:F(13),6:F(18),8:F(22)}
assert a['minima']==expected
truth={};feasible=0
for owners in itertools.product((0,1),repeat=8):
 z=realize(p,list(owners))
 if z['feasible']:
  feasible+=1;q=owners.count(0);truth[q]=min(truth.get(q,z['area']),z['area'])
assert truth==expected
assert max(q for q,d in expected.items() if d<=10)==3
mu=F(1,3);b=mu*10-min(mu*d-q for q,d in expected.items());assert b==4
out={'distinct_points':len(set(map(tuple,sites))),'unit_RNAs':8,'feasible_owner_vectors':feasible,'all_owner_vectors_checked':256,
 'spectrum':{str(q):str(d) for q,d in expected.items()},'budget':'10','exact_maximum':3,'best_penalty_upper_bound':'4',
 'attaining_multiplier':'1/3','sites':sites,'witness_count_three':serialized(a['witnesses'][3]),
 'proof':'For mu<=1/3 the count-eight witness gives B(mu)>=8-12mu>=4. For mu>=1/3 the count-two witness gives B(mu)>=2+6mu>=4. At 1/3 equality holds.',
 'interpretation':'Distinct explicit unit transcripts. No binary-encoded RNA multiplicity or free bundle variable.'}
(R/'results/unit_molecule_gap.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='witness_count_three'})
