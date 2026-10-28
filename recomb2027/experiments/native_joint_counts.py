"""Translate the previously identified native thickness obstruction into counts.

Post-hoc illustration; the named corridor and target genes were supplied in prior
review, not searched to maximize a new success fraction.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from geocert_expression.joint import outcome_spectrum,predicate_status
from mixed import make_problem,MixedLayers
from data_api import variant
ch=next(c for c in json.loads((ROOT/'data/charts.json').read_text()) if c['id']=='c8f8188430d092')
p0,pins,n=make_problem(ch);p=variant(p0)
obs=json.loads((ROOT/'data/observed_crop_counts.json').read_text())['rois'][ch['roi']]['counts']
channels=[];offsets=[];source=[];queries=[]
for gene,k in [('MDM2',2),('TFPI',0)]:
 ids=[i for i,z in enumerate(ch['points']) if z['gene']==gene]
 total=obs.get(str(ch['cell_labels'][k]),{}).get(gene,0);fixed=total-sum(ch['points'][i]['source']==k for i in ids)
 channels.append([[int(g in ids and c==k) for c in range(p.m+1)] for g in range(p.n)])
 offsets.append(fixed);source.append(total);queries.append({'gene':gene,'cell_id':ch['cell_labels'][k],'source_count':total,'fixed_outside':fixed,'indices':ids})
def run(model):
 e=ExpressionEngine(model,pins);sp=outcome_spectrum(e,channels,offsets)
 single=[e.spectrum(s,c).interval(F(4)) for s,c in zip(channels,offsets)]
 predicate=predicate_status(sp,F(4),lambda z:all(v>=a+1 for v,a in zip(z,source)))
 return {'individual_intervals':single,'joint_outcomes':sp,'both_increased_by_one':predicate}
answer=run(p)
# This is the already specified minimum-only ablation, not a newly tuned gap.
gaps=[max(F(1,10**9),p.bands[k+1][0]-p.bands[k][1]) for k in range(p.m-1)]
relaxed=MixedLayers(p.baselines,p.boxes,p.source,p.slopes,gaps,p.maximum_gaps,clearance=p.clearance,bands=p.bands,mode=p.mode)
other=run(relaxed)
assert answer['both_increased_by_one']['status']=='never'
assert other['both_increased_by_one']['status']=='variable'
out={'corridor':ch['id'],'setting':'maximal','budget_um2':'4','queries':queries,'original_minimum_gaps':list(map(str,p.minimum_gaps)),
     'relaxed_minimum_gaps':list(map(str,gaps)),'original':answer,'minimum_only_relaxation':other,
     'interpretation':'conditional simultaneous count increases; not actual cell-state, corrected-owner or signaling evidence'}
(ROOT/'results/native_joint_counts.json').write_text(json.dumps(out,default=str,indent=2)+'\n')
print(json.dumps({'queries':queries,'intervals':[[r['lower'],r['upper']] for r in answer['individual_intervals']],
                  'joint':answer['both_increased_by_one']['status'],'minimum_relaxed':other['both_increased_by_one']['status']},indent=2))
