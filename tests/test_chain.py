"""Test the exact count chain against independently enumerated scalar truth."""
from pathlib import Path
import sys,json,random,time
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from atlas.chain_spectrum import compile_chain,spectrum
from atlas.model import make_patch,count_scores,serialized,verify_witness
from atlas.scalar_reference import spectrum as independent_spectrum
from geocert_expression.budget import ExpressionEngine
from mixed import MixedLayers,TaggedBox
rng=random.Random(632);checks=states=witnesses=0;t0=time.time()
for seed in range(100):
 n=rng.randrange(1,8);T=F(8);rho=F(1,100);points=[]
 for i in range(n):
  x=F(rng.choice([-1,1])*rng.randrange(1,10),10);t=F(rng.randrange(1,80),10)
  points.append({'transcript_id':str(i),'gene':'A' if i%3 else 'B','source':int(x>0),'t':str(t),'x':str(x)})
 p={'points':points,'nuclear_constraints':[],'clearance_pixels':str(rho),'domain_pixels':[0,str(T)],
    'band_halfwidth_pixels':'3/4','halfwidth_pixels':'1','slope':'1','area_scale_um2':'1','cell_labels':[1,2]}
 if seed%3==0:p['nuclear_constraints']=[[1,'1/2','2','3']]
 e=make_patch(p);cm=compile_chain(e.problem,e.pins)
 for gene in ['A','B']:
  sc=count_scores(p,e,1,gene);ans=spectrum(cm,sc);truth,nn=independent_spectrum(p,gene)
  assert ans['minima']==truth,(seed,gene,ans['minima'],truth);checks+=1;states+=ans['DP_states']
  for w in ans['witnesses'].values():verify_witness(p,w['owners'],w['boundaries'],w['area']);witnesses+=1
# Geometrically admissible source, but a pinned owner forces positive area.
p=MixedLayers([[(0,0),(2,0)]],[TaggedBox(F(1,5),F(1,5),1,1,0,False)],[1],[1],[],[],clearance=F(1,20),bands=[(-1,1)],mode='nominal')
cm=compile_chain(p,{0:0});ans=spectrum(cm,[[1,0]]);assert ans['minima']=={1:F(1,16)};checks+=1
# A free regional group is deliberately outside the finite-local-count contract.
p=MixedLayers([[(0,0),(2,0)]],[TaggedBox(F(1,5),F(1,5),F(1,2),F(3,2),0,True)],[1],[1],[],[],clearance=F(1,20),bands=[(-1,1)],mode='universal')
try:compile_chain(p)
except ValueError:pass
else:raise AssertionError('accepted free interval owner without a memory state')
out={'count_spectrum_comparisons':checks,'DP_states':states,'checked_boundary_witnesses':witnesses,
     'positive_mandatory_area_test':True,'unsupported_region_rejected':True,'failures':0,'seconds':time.time()-t0}
(ROOT/'results/chain_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
