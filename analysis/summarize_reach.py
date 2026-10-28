"""Aggregate every locked width and budget, keeping denominators explicit."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import csv,json,math,sys
R=Path(__file__).resolve().parents[1]
WIDTHS=['cap4','cap8','cap16','cap20']
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def key(r):return (r['roi'],int(r['cell']),r['gene'])
def main():
 table=[];saturation=[];nuclear={};protocol=json.loads((R/'docs/integration/PROTOCOL_LOCK.json').read_text())
 for w in WIDTHS:
  for ds in ['lung','breast']:
   d=R/'reach_study'/w/ds; rows=read(d/'results/cell_expression_intervals.csv');cov=json.loads((d/'results/coverage_summary.json').read_text()); specs=json.loads((d/'results/cell_spectra.json').read_text())
   by=defaultdict(list);rich=set();tight=set();affected=set();widths=Counter();npatch=nedges=nRNA=ntot=ncand=0;maxband=0
   for r in rows:
    for field in ['lower','upper','independent_lower','independent_upper','cell']:r[field]=int(r[field])
    b=F(r['budget_um2']);by[b].append(r)
    if r['independent_upper']-r['independent_lower']>=2:rich.add(key(r))
    if (r['lower'],r['upper'])!=(r['independent_lower'],r['independent_upper']):tight.add(key(r));affected.add((r['roi'],r['cell']))
   for p in (d/'data').glob('atlas_*.json'):
    a=json.loads(p.read_text());pum=float(F(a['pixel_size_um']));npatch+=len(a['patches']);nedges+=a['summary']['selected_source_edges'];nRNA+=a['summary']['unique_editable_RNA'];ntot+=len(a['records']);ncand+=a['summary']['candidates']
    for pat in a['patches']:
     band=float(F(pat['band_halfwidth_pixels']))*pum/math.sqrt(2);maxband=max(maxband,band);widths[pat['halfwidth_pixels']]+=1
   # Exact largest minimum cost of a count endpoint is the endpoint-saturation budget.
   endpoint_sat=F(0);onset_sat=F(0);nchanged_afterquarter=0
   for s in specs:
    cost={int(k):F(v) for k,v in s['minima'].items()};lo,hi=min(cost),max(cost);endpoint_sat=max(endpoint_sat,cost[lo],cost[hi])
    if lo!=hi:
     onset=min(v for k,v in cost.items() if k!=s['source_count']);onset_sat=max(onset_sat,onset)
     quarter=[k for k,v in cost.items() if v<=F(1,4)];nchanged_afterquarter+=(min(quarter),max(quarter))!=(lo,hi)
   row={'width':w,'dataset':ds,'chart_cap':int(w[3:]),'max_support_halfwidth_um':maxband/0.75,'max_allowed_displacement_um':maxband,'width_histogram':dict(widths),
    'candidates':ncand,'patches':npatch,'editable_edges':nedges,'incident_edges':sum(x['all_incident_unique_edges'] for x in cov),'editable_RNA':nRNA,'retained_RNA':ntot,
    'cells':sum(x['crop_contained_connected_reference_cells'] for x in cov),'questions':len(specs),'multi_reachable_questions':len(rich),'tightened_questions':len(tight),'affected_cells':len(affected),
    'interval_rows':len(rows),'questions_still_widening_after_quarter':nchanged_afterquarter,'exact_endpoint_saturation_um2':str(endpoint_sat),'exact_nonsingleton_onset_saturation_um2':str(onset_sat),'by_budget':[]}
   for b,rr in sorted(by.items()):
    row['by_budget'].append({'budget_um2':str(b),'rows':len(rr),'non_singleton':sum(x['lower']<x['upper'] for x in rr),'marginal_width_ge2':sum(x['independent_upper']-x['independent_lower']>=2 for x in rr),
      'tightened':sum((x['lower'],x['upper'])!=(x['independent_lower'],x['independent_upper']) for x in rr),'max_width':max(x['upper']-x['lower'] for x in rr)})
   table.append(row);nuclear[(ds,w)]={key(r):r for r in read(d/'results/nucleus_reference_counts.csv') if int(r['source_count'])!=int(r['reference_count'])}
 proseg={w:{} for w in WIDTHS}
 for r in read(R/'results/proseg_reach_counts.csv'):
  if r['changed']=='True':proseg[r['width']][key(r)]=r
 comparison=[]
 for ds,method,rr in [('lung','nucleus',{w:nuclear[('lung',w)] for w in WIDTHS}),('breast','nucleus',{w:nuclear[('breast',w)] for w in WIDTHS}),('lung','proseg',proseg)]:
  common=set.intersection(*(set(q) for q in rr.values()))
  for w,questions in rr.items():
   vals=list(questions.values());cv=[questions[k] for k in sorted(common)]
   rec={'width':w,'dataset':ds,'reference':method,'changed_available':len(vals),'changed_common':len(cv),
     'attainable_available':sum(r['minimum_area_um2']!='inf' for r in vals),'attainable_common':sum(r['minimum_area_um2']!='inf' for r in cv),'by_budget':[]}
   for b in sorted(F(q['budget_um2']) for q in table[0]['by_budget']):
    rec['by_budget'].append({'budget_um2':str(b),'attainable_common':sum(r['minimum_area_um2']!='inf' and F(r['minimum_area_um2'])<=b for r in cv),'changed_common':len(cv)})
   comparison.append(rec)
 result={'protocol_lock':protocol,'table':table,'reference_comparison':comparison,'scope':'Fixed geometry-only width schedules; common changed-count keys across all schedules are primary. Per-width available denominators are retained. The 17-budget grid differs from the original eight-budget denominator audit.'}
 (R/'results/reach_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 for q in table:print(q['width'],q['dataset'],'patch',q['patches'],'coverage',q['editable_edges']/q['incident_edges'],'questions',q['questions'],'tight/rich',q['tightened_questions'],q['multi_reachable_questions'],'sat',float(F(q['exact_endpoint_saturation_um2'])))
 for q in comparison:print(q['width'],q['dataset'],q['reference'],q['attainable_common'],q['changed_common'])
if __name__=='__main__':main()
