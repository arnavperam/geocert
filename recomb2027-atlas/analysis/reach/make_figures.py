"""Plot the complete fixed-width sensitivity census and reference-count reach."""
from pathlib import Path
import json
from fractions import Fraction as F
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'manuscript/figures';O.mkdir(parents=True,exist_ok=True)
summary=json.loads((R/'results/reach/summary.json').read_text());d={(r['dataset'],r['setting']):r for r in summary};order=['baseline','0.5','1','2','3']
fig,ax=plt.subplots(figsize=(6.4,3.4))
for ds,ref,label,marker in [('lung','proseg','Lung: Proseg','o'),('lung','nuclear','Lung: nuclear expansion','s'),('breast','nuclear','Breast format test: nuclear expansion','^')]:
 rows=[d[ds,k] for k in order];xx=[100*r['edge_coverage'] for r in rows];yy=[100*r[f'{ref}_asymptotic_reach'] for r in rows]
 ax.plot(xx,yy,marker=marker,linewidth=1.6,label=label)
 if ref=='proseg' or ds=='breast':
  for x,y,k in zip(xx,yy,order):
   ax.annotate('A' if k=='baseline' else k,(x,y),xytext=(4,5 if ds=='breast' else -13),textcoords='offset points',fontsize=12)
ax.set_xlim(0,78);ax.set_ylim(0,22)
ax.set_xlabel('Editable incident boundary (%)',fontsize=13)
ax.set_ylabel('Attainable changed counts (%)',fontsize=13)
ax.tick_params(labelsize=12)
fig.legend(*ax.get_legend_handles_labels(),fontsize=12,loc='upper center',bbox_to_anchor=(.53,1.01),frameon=False)
fig.subplots_adjust(left=.14,right=.98,bottom=.20,top=.69)
fig.savefig(O/'reach_coverage.pdf');fig.savefig(O/'reach_coverage.png',dpi=180);plt.close(fig)
bs=json.loads((R/'results/reach/budget_summary.json').read_text())
fig,ax=plt.subplots(figsize=(6.4,3.4))
for k in order:
 rows=sorted([q for q in bs if q['dataset']=='breast' and q['setting']==k and F(q['budget_um2'])>0],key=lambda q:F(q['budget_um2']))
 ax.plot([float(F(q['budget_um2'])) for q in rows],[100*q['non_singleton']/q['questions'] for q in rows],marker='o',markersize=3,label='Adaptive' if k=='baseline' else k+' µm support')
ax.set_xscale('log');ax.set_xlabel('Shared physical-area budget (µm²)',fontsize=13);ax.set_ylabel('Non-singleton count entries (%)',fontsize=13)
ax.tick_params(labelsize=12);ax.legend(fontsize=12,loc='upper left');fig.tight_layout();fig.savefig(O/'width_budget.pdf');fig.savefig(O/'width_budget.png',dpi=180);plt.close(fig)
print('Saved reach/coverage and all-budget width figures from complete recorded results.')
