"""Regenerate supplementary budget tables from exact recorded CSVs/JSONs."""
from pathlib import Path
from fractions import Fraction as F
import json
R=Path(__file__).resolve().parents[2]
summary=json.loads((R/'results/reach/summary.json').read_text()); budgets=json.loads((R/'results/reach/budget_summary.json').read_text())
settings=['baseline','0.5','1','2','3'];labels=['Adaptive','0.5','1','2','3']
t=[]
for ds,name in [('lung','Lung'),('breast','Breast format-test')]:
 rows=[r for r in budgets if r['dataset']==ds];dsbud=sorted({r['budget_um2'] for r in rows},key=F)
 t += [r'\begin{table}[ht]\centering',r'\begin{tabular}{rrrrrr}\toprule',r'Area ($\mu\mathrm m^2$)&Adaptive&0.5&1&2&3\\\midrule']
 for d in dsbud:
  rr=[next(r for r in rows if r['setting']==st and r['budget_um2']==d) for st in settings]
  t.append(f'{float(F(d)):g} & '+' & '.join(str(r['non_singleton']) for r in rr)+r'\\')
 t += [r'\bottomrule\end{tabular}',f'\\caption{{{name}: non-singleton intervals on the common 17-budget grid. Column labels after Adaptive are required physical support halfwidths in micrometres. The numbers of question entries for the five families are '+', '.join(str(next(r for r in summary if r['dataset']==ds and r['setting']==st)['count_questions']) for st in settings)+'.}',r'\end{table}']
 t += [r'\begin{table}[ht]\centering',r'\begin{tabular}{rrrrrr}\toprule',r'Area ($\mu\mathrm m^2$)&Adaptive&0.5&1&2&3\\\midrule']
 for d in dsbud:
  rr=[next(r for r in rows if r['setting']==st and r['budget_um2']==d) for st in settings]
  t.append(f'{float(F(d)):g} & '+' & '.join(f"{r['tightened_questions']}/{r['eligible_marginal_width_ge2']}" for r in rr)+r'\\')
 t += [r'\bottomrule\end{tabular}',f'\\caption{{{name}: tightened cell--gene questions / questions with marginal width at least two at each budget. Zero over zero denotes no eligible comparison, not a measured zero effect rate. Each question can recur at multiple budgets.}}',r'\end{table}']
(R/'manuscript/reach_budget_tables.tex').write_text('\n'.join(t)+'\n')
# Original-grid budgets and change quantities, using exact stored results.
import csv
orig=[]
for ds,base in [('lung',R),('breast',R/'external')]:
 rows=list(csv.DictReader((base/'results/cell_expression_intervals.csv').open()))
 orig.append((ds,rows))
t=[r'\begin{table}[ht]\centering',r'\begin{tabular}{rrrrrrr}\toprule',r'&\multicolumn{3}{c}{Lung}&\multicolumn{3}{c}{Breast format-test}\\',r'Area&Variable&Width sum&Endpoints changed&Variable&Width sum&Endpoints changed\\\midrule']
prev={}
for d in ['0','1/20','1/10','1/4','1/2','1','2','4']:
 vals=[]
 for ds,rows in orig:
  z=[r for r in rows if F(r['budget_um2'])==F(d)]
  now={(r['roi'],r['cell'],r['gene']):(int(r['lower']),int(r['upper'])) for r in z}
  prior=prev.get(ds,now)
  vals += [sum(a<b for a,b in now.values()),sum(b-a for a,b in now.values()),sum(v!=prior[k] for k,v in now.items())]
  prev[ds]=now
 t.append(f'{float(F(d)):g} & '+' & '.join(map(str,vals))+r'\\')
t += [r'\bottomrule\end{tabular}',r'\caption{Primary eight-budget occupancy, summed widths and number of cell--gene intervals whose endpoints change from the preceding displayed budget. A constant non-singleton total can coexist with further endpoint widening.}',r'\end{table}']
(R/'manuscript/original_saturation_table.tex').write_text('\n'.join(t)+'\n')

print("Regenerated complete supplementary budget tables.")
