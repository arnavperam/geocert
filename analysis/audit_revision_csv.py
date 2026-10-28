"""Read-only exact census of the supplied Atlas interval tables."""
from pathlib import Path
import csv,json,hashlib
from collections import defaultdict,Counter
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'
FILES={'lung':ROOT/'results/cell_expression_intervals.csv','breast':ROOT/'external/results/cell_expression_intervals.csv'}
def census(path):
 rows=list(csv.DictReader(path.open()));byq=defaultdict(list);byb=defaultdict(list)
 for r in rows:
  for k in ['source_count','lower','upper','independent_lower','independent_upper','mutable_patches']:r[k]=int(r[k])
  r['budget']=F(r['budget_um2']);r['tight']=r['lower']>r['independent_lower'] or r['upper']<r['independent_upper']
  q=(r['atlas_id'],r['roi'],r['cell'],r['gene']);byq[q].append(r);byb[r['budget']].append(r)
 assert len({(q,r['budget']) for q,rr in byq.items() for r in rr})==len(rows)
 budgets=sorted(byb);keys=sorted(byq)
 table=[]
 for b in budgets:
  rr=byb[b];table.append({'budget':str(b),'rows':len(rr),'nonsingleton':sum(r['upper']>r['lower'] for r in rr),'tightened':sum(r['tight'] for r in rr),'width_sum':sum(r['upper']-r['lower'] for r in rr),'max_width':max(r['upper']-r['lower'] for r in rr),'marginal_width_at_least_2':sum(r['independent_upper']-r['independent_lower']>=2 for r in rr),'mean_width':sum(r['upper']-r['lower'] for r in rr)/len(rr)})
 changed=[]
 for a,b in zip(budgets,budgets[1:]):
  dif=[]
  for q,rr in byq.items():
   ra=next(r for r in rr if r['budget']==a);rb=next(r for r in rr if r['budget']==b)
   if (ra['lower'],ra['upper'])!=(rb['lower'],rb['upper']):dif.append({'roi':q[1],'cell':q[2],'gene':q[3],'from':[ra['lower'],ra['upper']],'to':[rb['lower'],rb['upper']]})
  changed.append({'from_budget':str(a),'to_budget':str(b),'questions_changed':len(dif),'examples':dif[:5]})
 effectq=[q for q,rr in byq.items() if any(r['tight'] for r in rr)]
 multiq=[q for q,rr in byq.items() if max(r['independent_upper']-r['independent_lower'] for r in rr)>=2]
 positive_count_q=[q for q,rr in byq.items() if rr[0]['source_count']>0]
 values={'csv_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(rows),'questions':len(byq),'statuses':dict(Counter(r['status'] for r in rows)),
 'cells_in_csv':len({q[:3] for q in byq}),'positive_source_questions':len(positive_count_q),'tightened_rows':sum(r['tight'] for r in rows),'questions_ever_tightened':len(effectq),'questions_with_max_marginal_width_at_least_2':len(multiq),'tightened_questions_in_multimove_denominator':len(set(effectq)&set(multiq)),
 'cells_with_ever_tightened_question':len({q[:3] for q in effectq}),'cells_with_max_marginal_width_at_least_2':len({q[:3] for q in multiq}),
 'by_budget':table,'endpoint_changes_between_budgets':changed}
 examples=[r for r in rows if r['roi']=='breast_01' and r['cell']=='81' and r['gene']=='LDHB']
 values['ldhb_cell81']=[{k:v for k,v in r.items() if k!='budget'} for r in examples]
 return values
result={k:census(p) for k,p in FILES.items()}
(OUT/'comparison_csv_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
