"""Reject location uncertainty that destroys unique nodal RNA ownership."""
from pathlib import Path
from fractions import Fraction as F
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from atlas.chain_spectrum import compile_chain,spectrum
from mixed import MixedLayers,TaggedBox
from layered import realize

def problem(box,mode):
    return MixedLayers([[(0,0),(2,0)]],[box],[int(box.xlo+box.xhi>0)],[1],[],[],clearance=F(1,10),bands=[(-1,1)],mode=mode)

# The inherited possible-location relaxation has overlapping owner obligations.
# It is a valid core model, but outside the point-count-chain theorem.
box=TaggedBox(-F(1,5),F(1,5),F(1),F(1),0,False)
p=problem(box,'possible_outer')
assert realize(p,[0])['area']==realize(p,[1])['area']==0
rejections=0
for box in [box,TaggedBox(F(1,5),F(1,5),F(1,2),F(3,2),0,False)]:
    for mode in ('universal','nominal','possible_outer'):
        try:compile_chain(problem(box,mode))
        except ValueError as e:
            assert 'point-valued' in str(e);rejections+=1
        else:raise AssertionError('unsupported free location region accepted')
multi=MixedLayers([[(0,0),(2,0)]],
    [TaggedBox(F(1,5),F(1,5),F(1,2),F(1,2),0,False),
     TaggedBox(F(1,5),F(1,5),F(3,2),F(3,2),0,False)],
    [1],[1],[],[],clearance=F(1,10),bands=[(-1,1)],mode='nominal')
try:compile_chain(multi)
except ValueError as e:
    assert 'one point observation' in str(e);rejections+=1
else:raise AssertionError('free shared owner at multiple locations accepted')
checks=0
point=TaggedBox(F(1,5),F(1,5),F(1),F(1),0,False)
for mode in ('universal','nominal','possible_outer'):
    p=problem(point,mode);got=spectrum(compile_chain(p),[[1,0]])['minima']
    truth={1-o:realize(p,[o])['area'] for o in (0,1)}
    assert got==truth;checks+=1
out={'unsupported_free_location_or_multisite_groups_rejected':rejections,
     'supported_point_mode_spectra_verified':checks,'failures':0,
     'counterexample':'With possible_outer horizontal location uncertainty, both owners have zero deformation. The old chain silently chose one nodal owner and overcharged the other; revised chain rejects the unsupported input.'}
(ROOT/'results/chain_contract_tests.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
