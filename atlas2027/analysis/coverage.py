"""Boundary and observation denominators, including frozen interfaces."""
from pathlib import Path
import sys,os,json,csv,statistics
from fractions import Fraction as F
from collections import Counter,defaultdict
import numpy as np
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
rows=[];summary=[];indexes=[]
for p in sorted((ROOT/'data').glob('atlas_*.json')):
 a=json.loads(p.read_text());cells=set(a['complete_cells']);m=np.load(ROOT/'data/source_crops'/a['roi']/'masks.npz')['cell_mask']
 all_incident=0;perimeter=Counter();editable=Counter();trans=defaultdict(set);genes=defaultdict(Counter)
 for left,right in [(m[:,:-1],m[:,1:]),(m[:-1,:],m[1:,:])]:
  ix=np.where(left!=right)
  for u,v in zip(left[ix].tolist(),right[ix].tolist()):
   if u in cells or v in cells:all_incident+=1
   if u in cells:perimeter[u]+=1
   if v in cells:perimeter[v]+=1
 for pat in a['patches']:
  for c in pat['cell_labels']:
   if c in cells:
    editable[c]+=pat['length_edges']
    for r in pat['points']:trans[c].add(r['transcript_id']);genes[c][r['gene']]+=1
 for c in sorted(cells):
  rows.append({'atlas_id':a['atlas_id'],'roi':a['roi'],'cell':c,'all_reference_boundary_edges':perimeter[c],
   'editable_reference_boundary_edges':editable[c],'editable_fraction':editable[c]/perimeter[c],
   'observed_source_count':sum(a['source_counts'].get(str(c),{}).values()),'incident_editable_RNA':len(trans[c])})
 hist=Counter(n for vals in genes.values() for n in vals.values())
 summary.append({'roi':a['roi'],'crop_contained_connected_reference_cells':len(cells),'all_incident_unique_edges':all_incident,
                 'atlas_editable_unique_edges':a['summary']['selected_source_edges'],
                 'all_incident_edge_fraction':a['summary']['selected_source_edges']/all_incident,
                 'cell_coverage_min':min(editable[c]/perimeter[c] for c in cells),
                 'cell_coverage_median':statistics.median(editable[c]/perimeter[c] for c in cells),
                 'cell_coverage_max':max(editable[c]/perimeter[c] for c in cells),
                 'editable_cell_gene_multiplicity_histogram':dict(sorted(hist.items()))})
 indexes.append({'roi':a['roi'],'atlas_id':a['atlas_id'],'cells':sorted(cells),'observed_genes':sorted({r['gene'] for r in a['records']}),
   'input_observations':len(a['records']),'editable_observations':a['summary']['unique_editable_RNA'],
   'implicit_zero_rule':'A missing CSV key is [0,0] only for a listed cell and observed gene in this atlas; outside observations remain fixed. Rows with source counts but no mutable incident gene are explicitly present.',
   'joint_scope':'All rows share one atlas and one physical-area budget. Their marginal endpoints need not be jointly attainable.',
   'population_scope':'Counts of the supplied retained measured molecules for crop-contained reference regions; no undetected RNA, calibrated anatomy, or unobserved cell parts are inferred.'})
(ROOT/'results/coverage_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(ROOT/'results/matrix_index.json').write_text(json.dumps(indexes,indent=2)+'\n')
with (ROOT/'results/cell_boundary_coverage.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(summary)
