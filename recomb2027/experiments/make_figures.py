"""Recreate manuscript figures from saved scientific records.

Requires Matplotlib only for rendering. No simulated image is presented as tissue.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'manuscript/figures';OUT.mkdir(parents=True,exist_ok=True)

def save(fig,name):
    fig.savefig(OUT/(name+'.pdf'))
    fig.savefig(OUT/(name+'.png'),dpi=180)
    plt.close(fig)

def main():
    data=json.loads((ROOT/'results/native_expression_summary.json').read_text())
    fig,ax=plt.subplots(figsize=(6.4,3.1))
    for setting,label,marker in [('baseline','Baseline geometry','o'),('maximal','Enlarged geometry','s')]:
        rows=[r for r in data['table'] if r['setting']==setting]
        ax.plot(range(len(rows)),[r['nontrivial_intervals'] for r in rows],marker=marker,label=label)
    ax.set_xticks(range(8),['0','0.05','0.1','0.25','0.5','1','2','4'])
    ax.set_xlabel('Declared deformation budget (µm²)',fontsize=11)
    ax.set_ylabel('Non-singleton count intervals / 899',fontsize=11)
    ax.tick_params(labelsize=10);ax.legend(fontsize=10);fig.tight_layout();save(fig,'native_intervals')
    data=json.loads((ROOT/'results/dual_gap.json').read_text())
    points=sorted((float(F(area)),int(count)) for count,area in data['spectrum'].items())
    fig,ax=plt.subplots(figsize=(6.4,3.1))
    ax.scatter([x for x,y in points],[y for x,y in points],s=35,label='Attainable count and minimum area')
    ax.plot([4,22],[2,8],'--',label='Supporting penalty bound')
    ax.axvline(10,linestyle=':',label='Hard budget = 10')
    ax.scatter([10],[4],marker='x',s=60,label='Best penalty upper bound')
    ax.annotate('Exact maximum: 3',(9,3),xytext=(2,5),arrowprops={'arrowstyle':'->'},fontsize=10.5)
    ax.set_xlabel('Minimum deformation (designed area units)',fontsize=11)
    ax.set_ylabel('Count in target cell',fontsize=11);ax.tick_params(labelsize=10)
    ax.legend(fontsize=10,loc='lower right');fig.tight_layout();save(fig,'hard_budget_gap')
    chart=next(c for c in json.loads((ROOT/'data/charts.json').read_text()) if c['id']=='4037abe9fea38e')
    record=json.loads((ROOT/'results/native/4037abe9fea38e_baseline.json').read_text())
    query=next(q for q in record['queries'] if q['gene']=='PECAM1' and q['cell_index']==0)
    witness=query['spectrum']['witnesses']['1']
    fig,ax=plt.subplots(figsize=(6.4,3.4));base=chart['baselines'][0]
    ax.plot([float(F(t)) for t,_ in base],[float(F(x)) for _,x in base],'--',label='Reference boundary')
    ts=[];xs=[]
    for a,b,u,v in witness['boundaries'][0]:ts.extend([float(F(a)),float(F(b))]);xs.extend([float(F(u)),float(F(v))])
    ax.plot(ts,xs,linewidth=1.8,label='Attaining count-one boundary')
    for is_target,marker,size,label in [(False,'.',22,'Other retained RNA'),(True,'x',65,'PECAM1 RNA')]:
        pts=[p for p in chart['points'] if (p['gene']=='PECAM1')==is_target]
        ax.scatter([float(F(p['t'])) for p in pts],[float(F(p['x'])) for p in pts],marker=marker,s=size,label=label)
    ax.set_xlabel('Longitudinal coordinate (µm)',fontsize=11);ax.set_ylabel('Transverse coordinate (µm)',fontsize=11)
    ax.tick_params(labelsize=10);ax.legend(fontsize=10,loc='upper center',bbox_to_anchor=(0.5,1.25),ncol=2)
    fig.subplots_adjust(top=.76,bottom=.18,left=.12,right=.98);save(fig,'native_count_witness')
    print('Three figures recreated from saved exact data and boundary witnesses.')
if __name__=='__main__':main()
