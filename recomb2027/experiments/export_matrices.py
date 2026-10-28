"""Sparse full count-matrix export, retaining the corridor conditioning key.

Absent entries are zero only within the recorded observed-gene universe and
retained molecule set. Fixed positive entries are algebraic copies, not additional
optimized experiments. Different corridors must not be silently merged.
"""
from pathlib import Path
import json,csv,collections
ROOT=Path(__file__).resolve().parents[1]
fields=['corridor','roi','setting','cell_id','gene','budget_um2','source_count','lower','upper','status','method']
rows=[];models=[];fixed_rows=0;query_rows=0
observed=json.loads((ROOT/'data/observed_crop_counts.json').read_text())
genes=sorted({r['gene'] for d in observed['rois'].values() for r in d['records']})
for p in sorted((ROOT/'results/native').glob('*.json')):
 d=json.loads(p.read_text());first=d['rows'][0];budgets=list(dict.fromkeys(r['budget_um2'] for r in d['rows']))
 cells=sorted({r['cell_id'] for r in d['rows']})
 models.append({'corridor':first['corridor'],'roi':first['roi'],'setting':first['setting'],'cells':cells,'budgets_um2':budgets})
 for r in d['rows']:
  rows.append({k:r[k] for k in fields if k not in ['status','method']}|{'status':'exact','method':'complete_target_projection'});query_rows+=1
 local_genes={r['gene'] for r in d['rows']}
 fixed=[{'cell_id':cell,'gene':gene,'count':count} for cell in cells for gene,count in observed['rois'][first['roi']]['counts'].get(str(cell),{}).items() if gene not in local_genes and count>0]
 for r in fixed:
  for budget in budgets:
   rows.append({'corridor':first['corridor'],'roi':first['roi'],'setting':first['setting'],
                'cell_id':r['cell_id'],'gene':r['gene'],'budget_um2':budget,'source_count':r['count'],
                'lower':r['count'],'upper':r['count'],'status':'exact','method':'fixed_outside_editable_region'});fixed_rows+=1
keys=[tuple(r[k] for k in ['corridor','setting','cell_id','gene','budget_um2']) for r in rows]
assert len(keys)==len(set(keys))
with (ROOT/'results/expression_matrix_intervals.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
meta={'observed_gene_universe':genes,'models':models,'optimized_rows':query_rows,'fixed_positive_rows':fixed_rows,'sparse_rows':len(rows),
      'absent_entry_rule':'For a listed model/cell/budget and listed observed gene, absent CSV entries have count exactly zero under that model. This says nothing about undetected RNA or unmeasured genes.',
      'conditioning':'Separate matrices per corridor and geometric setting; no averaging or union over overlapping corridors.',
      'count_scope':observed['scope'],'raster_label_semantics':{'0':'reference background/unassigned compartment','-1':'unknown exterior candidate'}}
(ROOT/'results/expression_matrix_index.json').write_text(json.dumps(meta,indent=2)+'\n')
print({k:meta[k] for k in ['optimized_rows','fixed_positive_rows','sparse_rows']});print('observed genes',len(genes))
