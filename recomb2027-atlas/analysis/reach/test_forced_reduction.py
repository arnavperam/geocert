"""Compare hard-closure candidate pruning with the unpruned exact count chain."""
from pathlib import Path
from fractions import Fraction as F
import json,random,sys,time
R=Path(__file__).resolve().parents[2];sys.path.insert(0,str(R))
from atlas.model import make_patch,count_scores,verify_witness
from atlas.chain_spectrum import compile_chain,spectrum
from atlas.scalar_reference import spectrum as scalar
rng=random.Random(20261028);compared=owners=witnesses=0;t=time.time()
for seed in range(80):
 points=[]
 for i in range(rng.randrange(1,8)):
  x=F(rng.choice([-1,1])*rng.randrange(1,15),10)
  points.append({'transcript_id':str(i),'gene':'A' if i%3 else 'B','source':int(x>0),'t':str(F(rng.randrange(1,80),10)),'x':str(x)})
 p={'points':points,'nuclear_constraints':[],'clearance_pixels':'1/100','domain_pixels':[0,'8'],'band_halfwidth_pixels':'3/4','halfwidth_pixels':'1','slope':'1','area_scale_um2':'1','cell_labels':[1,2]}
 if seed%2==0:p['nuclear_constraints']=[[1,'1/2','2','3']]
 e=make_patch(p);cl=e.source_engine.closure({})
 forced={**{g:0 for g in cl['forced_ones']},**{g:1 for g in cl['forced_zeros']}}
 raw=compile_chain(e.problem,e.pins);pruned=compile_chain(e.problem,{**e.pins,**forced})
 for gene in ('A','B'):
  scores=count_scores(p,e,1,gene);a=spectrum(raw,scores);b=spectrum(pruned,scores);truth,n=scalar(p,gene)
  assert a['minima']==b['minima']==truth,(seed,gene)
  compared+=1;owners+=n
  for w in b['witnesses'].values():verify_witness(p,w['owners'],w['boundaries'],w['area']);witnesses+=1
out={'models':80,'exact_count_spectra_compared':compared,'independently_enumerated_target_patterns':owners,'checked_pruned_witnesses':witnesses,'failures':0,'seconds':time.time()-t,'scope':'New fixed-bit pruning compared with unpruned chain and separately integrated scalar truth; same prepared geometry.'}
(R/'results/reach/forced_reduction_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
