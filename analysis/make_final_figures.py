"""Publication figures from original masks and frozen exact records."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
R=Path(__file__).resolve().parents[1];O=R/'manuscript/figures';O.mkdir(exist_ok=True)
def save(fig,name):
 fig.tight_layout();fig.savefig(O/(name+'.pdf'));fig.savefig(O/(name+'.png'),dpi=160);plt.close(fig)
# Original source image, all visible observations, fixed axis extent.
a=json.loads((R/'data/atlas_roi_01.json').read_text());m=np.load(R/'data/source_crops/roi_01/masks.npz')['cell_mask'];p=float(F(a['pixel_size_um']))
fig,ax=plt.subplots(figsize=(6.4,3.55));ids=np.unique(m)
ax.imshow(np.searchsorted(ids,m),origin='upper',extent=(0,m.shape[1]*p,m.shape[0]*p,0),interpolation='nearest',alpha=.25)
lines=[np.asarray(q['source_path_pixels'])*p for q in a['patches']]
ax.add_collection(LineCollection(lines,linewidths=.8,label='Editable reference arcs'))
meta=json.loads((R/'data/source_crops/roi_01/metadata.json').read_text());ox,oy=meta['bbox_pixels'][:2]
pts=[(float(r['x_um'])-ox*p,float(r['y_um'])-oy*p) for r in a['records']]
pts=[(x,y) for x,y in pts if 45<=x<=125 and 0<=y<=44]
ax.scatter([x for x,y in pts],[y for x,y in pts],s=5,marker='.',label='Measured transcripts')
ax.plot([111,121],[40,40],linewidth=3,label='10 µm scale bar')
ax.set_xlim(45,125);ax.set_ylim(44,0);ax.set_aspect('equal');ax.set_xlabel('Crop coordinate (µm)',fontsize=11.5);ax.set_ylabel('Crop coordinate (µm)',fontsize=11.5);ax.tick_params(labelsize=11.5)
ax.legend(fontsize=11.5,loc='upper center',bbox_to_anchor=(.5,1.20),ncol=2,frameon=False);save(fig,'atlas_reference')
if (R/'results/reach_summary.json').exists():
 data=json.loads((R/'results/reach_summary.json').read_text());rows={(q['dataset'],q['width']):q for q in data['table']}
 fig,ax=plt.subplots(figsize=(6.4,3.25))
 for ds,ref,label,marker in [('lung','nucleus','Lung: nuclear reference','o'),('breast','nucleus','Breast test: nuclear reference','s'),('lung','proseg','Lung: recorded Proseg','^')]:
  rr=sorted([q for q in data['reference_comparison'] if q['dataset']==ds and q['reference']==ref],key=lambda q:rows[(ds,q['width'])]['chart_cap'])
  x=[100*rows[(ds,q['width'])]['editable_edges']/rows[(ds,q['width'])]['incident_edges'] for q in rr]
  y=[100*q['attainable_common']/q['changed_common'] for q in rr]
  ax.plot(x,y,marker=marker,label=label)
  for j,(xx,yy) in enumerate(zip(x,y)):
   ax.annotate(str(rows[(ds,rr[j]['width'])]['chart_cap']),(xx,yy),xytext=(4,4 if j%2==0 else -12),textcoords='offset points',fontsize=11.5)
 ax.set_xlabel('Editable incident boundary length (%)',fontsize=11.5);ax.set_ylabel('Alternative changed counts attainable (%)',fontsize=11.5);ax.tick_params(labelsize=11.5)
 ax.legend(fontsize=11.5,loc='upper center',bbox_to_anchor=(.5,1.32),frameon=False,ncol=1)
 save(fig,'reach_coverage')
 # Saturation of widths is different from saturation of indicator counts.
 fig,ax=plt.subplots(figsize=(6.4,3.0))
 for width,marker in [('cap4','o'),('cap8','s'),('cap16','^'),('cap20','D')]:
  q=rows[('breast',width)];rr=q['by_budget'];ax.plot([float(F(x['budget_um2'])) for x in rr if F(x['budget_um2'])>0],
   [x['non_singleton']/q['questions'] for x in rr if F(x['budget_um2'])>0],marker=marker,markersize=3,label='B cap '+str(q['chart_cap']))
 ax.set_xscale('log');ax.set_xlabel('Shared deformation budget (µm²)',fontsize=11.5);ax.set_ylabel('Non-singleton fraction (breast test)',fontsize=11.5);ax.tick_params(labelsize=11.5);ax.legend(fontsize=11.5,ncol=2);save(fig,'budget_reach')
print('Figures saved from fixed source data and completed summaries.')
