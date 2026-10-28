from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.evaluation import evaluate_aligned,fwer_gate
from geocert_expression.calibration import threshold
r=[{'region':'a','id':1,'gene':'G','truth_owner':1,'predicted_owner':1,'fragility':0},
   {'region':'a','id':2,'gene':'G','truth_owner':1,'predicted_owner':2,'fragility':1},
   {'region':'b','id':3,'gene':'G','truth_owner':2,'predicted_owner':2,'fragility':0}]
ii=[{'region':'a','cell':1,'gene':'G','lower':1,'upper':2},
    {'region':'a','cell':2,'gene':'G','lower':0,'upper':1},
    {'region':'b','cell':2,'gene':'G','lower':1,'upper':1}]
o=evaluate_aligned(r,ii);assert o['simultaneous_region_coverage']==2 and o['pooled_AUROC']==1
f=fwer_gate(['1/1000','1/10'],family_coverage_error='1/100',familywise_error='1/20');assert f['reject']==[True,False]
rejected=0
for bad in [r+[r[0]]]:
 try:evaluate_aligned(bad,ii)
 except ValueError:rejected+=1
try:evaluate_aligned(r,ii[:1])
except ValueError:rejected+=1
try:fwer_gate([0],family_coverage_error='1/10',familywise_error='1/20')
except ValueError:rejected+=1
assert rejected==3
out={'aligned_example':o,'FWER_bridge':f,'malformed_or_inapplicable_inputs_rejected':rejected,'failures':0}
(ROOT/'results/evaluation_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
