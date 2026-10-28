"""Adopted point-contract/binding guards and equivalent forced-owner dispatch."""
from pathlib import Path
import sys, json, copy, random
from fractions import Fraction as F
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R))
from atlas.chain_spectrum import compile_chain,spectrum
from atlas.binding import verify_prepared_atlas,atlas_digest
from atlas.model import make_patch,count_scores
from mixed import MixedLayers,TaggedBox
checks=0
for mode in ['nominal','universal','possible_outer']:
 for x0,x1,t0,t1 in [(F(1,5),F(2,5),1,1),(F(1,5),F(1,5),F(1,2),F(3,2))]:
  p=MixedLayers([[(0,0),(2,0)]],[TaggedBox(x0,x1,t0,t1,0,True)],[1],[1],[],[],clearance=F(1,20),bands=[(-1,1)],mode=mode)
  try:compile_chain(p)
  except ValueError:checks+=1
  else:raise AssertionError('region input accepted by point-only chain')
a=json.loads((R/'data/atlas_roi_01.json').read_text());audit=json.loads((R/'results/atlas_geometric_validation.json').read_text())
verify_prepared_atlas(a,R/'data/source_crops/roi_01',audit)
for field in ['pixel_size_um','source_counts','points','duplicate']:
 q=copy.deepcopy(a)
 if field=='pixel_size_um':q[field]='1'
 elif field=='source_counts':
  k=next(iter(q[field]));g=next(iter(q[field][k]));q[field][k][g]+=1
 elif field=='points':next(p for p in q['patches'] if p['points'])['points'][0]['gene']='ALTERED'
 else:q['patches'].append(copy.deepcopy(q['patches'][0]))
 try:verify_prepared_atlas(q,R/'data/source_crops/roi_01',audit)
 except ValueError:checks+=1
 else:raise AssertionError('unbound altered atlas accepted')
 # Rehashing the input still cannot bind it to the old checked geometry.
 q['atlas_id']=atlas_digest(q)
 try:verify_prepared_atlas(q,R/'data/source_crops/roi_01',audit)
 except ValueError:checks+=1
 else:raise AssertionError('altered atlas with changed digest accepted')
rng=random.Random(110019);dispatch=0
for repeat in range(30):
 pts=[]
 for i in range(7):
  x=F(rng.choice([-12,-9,-4,2,5,11]),10);pts.append({'transcript_id':str(i),'gene':'G' if i%2 else 'H','source':int(x>0),'t':str(F(rng.randrange(1,70),10)),'x':str(x)})
 p={'points':pts,'nuclear_constraints':[], 'clearance_pixels':'1/100','domain_pixels':[0,'8'],'band_halfwidth_pixels':'3/4','halfwidth_pixels':'2','slope':'1','area_scale_um2':'1','cell_labels':[1,2]}
 e=make_patch(p);pins=dict(e.pins)
 for i,r in enumerate(pts):
  if not e.source_engine.closure({i:1-r['source']})['feasible']:pins[i]=r['source']
 raw=compile_chain(e.problem,e.pins);pruned=compile_chain(e.problem,pins)
 for gene in ['G','H']:
  scores=count_scores(p,e,1,gene);u=spectrum(raw,scores);v=spectrum(pruned,scores);w=e.spectrum(scores)
  assert u['minima']==v['minima']==w.minima;dispatch+=1
out={'rejected_invalid_contract_or_binding_cases':checks,'dispatch_equivalence_spectra':dispatch,'failures':0,'scope':'Adopted public input guard and atlas binding; exhaustive target projection versus pinned/unpinned critical chains.'}
(R/'results/adopted_guards_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
