"""Figures from the fixed atlas and exact recorded count spectra."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R));O=R/'manuscript/figures';O.mkdir(exist_ok=True)
from atlas.model import compose
from atlas.transport import chart_to_pixel
from atlas.geometry import evaluate_path

def save(fig,name):
 fig.tight_layout();fig.savefig(O/(name+'.pdf'));fig.savefig(O/(name+'.png'),dpi=170);plt.close(fig)

# Source reference and selected patch support: actual cell mask, no generated imagery.
a=json.loads((R/'data/atlas_roi_01.json').read_text());z=np.load(R/'data/source_crops/roi_01/masks.npz');m=z['cell_mask'];p=float(F(a['pixel_size_um']))
fig,ax=plt.subplots(figsize=(6.3,3.5));ids=np.unique(m);show=np.searchsorted(ids,m)
ax.imshow(show,origin='upper',extent=(0,m.shape[1]*p,m.shape[0]*p,0),interpolation='nearest',alpha=.30)
paths=[np.asarray(q['source_path_pixels'])*p for q in a['patches']]
ax.add_collection(LineCollection(paths,linewidths=.9,label='Editable source arcs'))
x0,y0=a['bbox_pixels'][:2]
points=np.asarray([[float(F(r['x_um']))-x0*p,float(F(r['y_um']))-y0*p] for r in a['records']])
inside=(points[:,0]>=45)&(points[:,0]<=125)&(points[:,1]>=0)&(points[:,1]<=44)
ax.scatter(points[inside,0],points[inside,1],marker='.',s=6,alpha=.65,label='Observed RNA')
ax.annotate('',xy=(62,39),xytext=(52,39),arrowprops={'arrowstyle':'|-|','linewidth':1.4})
ax.text(57,36.5,'10 µm',ha='center',fontsize=12)
ax.set_xlim(45,125);ax.set_ylim(44,0);ax.set_aspect('equal');ax.set_xlabel('Crop coordinate (µm)',fontsize=12);ax.set_ylabel('Crop coordinate (µm)',fontsize=12);ax.tick_params(labelsize=12)
ax.set_title('Lung crop: fixed-seam boundary atlas',fontsize=12)
ax.legend(fontsize=11,loc='upper center',bbox_to_anchor=(.5,1.02),ncol=2,markerscale=2)
save(fig,'atlas_reference')

# Exact cell-level example selected after complete census; original source cells.
D=R/'external';a=json.loads((D/'data/atlas_breast_01.json').read_text());z=np.load(D/'data/source_crops/breast_01/masks.npz');m=z['cell_mask'];p=float(F(a['pixel_size_um']));cell=81;gene='LDHB'
s=next(q for q in json.loads((D/'results/cell_spectra.json').read_text()) if q['roi']=='breast_01' and q['cell']==cell and q['gene']==gene)
pm={q['id']:q for q in a['patches']};ls={(q['patch'],q['gene']):q for q in json.loads((D/'results/patch_spectra.json').read_text()) if q['roi']=='breast_01'}
fig,ax=plt.subplots(figsize=(6.3,3.6));yy,xx=np.nonzero(m==cell);x0,x1=max(0,xx.min()-15),min(m.shape[1],xx.max()+16);y0,y1=max(0,yy.min()-15),min(m.shape[0],yy.max()+16)
sub=m[y0:y1,x0:x1];ids=np.unique(sub);ax.imshow(np.searchsorted(ids,sub),origin='upper',extent=(x0*p,x1*p,y1*p,y0*p),interpolation='nearest',alpha=.22)
firstsource=firstalt=firstRNA=firstcollateral=True;record=[]; patch_number=0
for pid,delta in s['witness_traces']['17'].items():
 q=pm[pid];sp=ls[(pid,gene)];left=sp['source_left']+delta*(1 if q['cell_labels'][0]==cell else -1);w=sp['witnesses'][str(left)]
 if not delta:continue
 patch_number+=1
 path=np.asarray(q['source_path_pixels'])*p;ax.plot(path[:,0],path[:,1],'--',linewidth=1.2,label=f'Reference {patch_number}');firstsource=False
 seg=[list(map(F,v)) for v in w['boundaries'][0]];ts={u for v in seg for u in v[:2]};sx,sy=q['sx'],q['sy'];org=sx*q['source_path_pixels'][0][0]+sy*q['source_path_pixels'][0][1]
 ts|={F(sx*x+sy*y-org) for x,y in q['source_path_pixels']}
 def at(t):
  for u,v,x,y in seg:
   if u<=t<=v:return x+(y-x)*(t-u)/(v-u)
  raise ValueError(t)
 points=np.asarray([[float(v)*p for v in chart_to_pixel(q,t,at(t))] for t in sorted(ts)])
 ax.plot(points[:,0],points[:,1],linewidth=1.8,label=f'Attaining {patch_number}');firstalt=False
 moved=[]
 for obs,o in zip(q['points'],w['owners']):
  if o==obs['source']:continue
  xy=[float(v)*p for v in chart_to_pixel(q,F(obs['t']),F(obs['x']))];moved.append((obs['gene'],obs['transcript_id'],xy))
  ax.scatter(*xy,marker='x' if obs['gene']==gene else '+',s=40,label=('Moved LDHB' if firstRNA else None) if obs['gene']==gene else ('KRT8 collateral' if firstcollateral else None))
  if obs['gene']==gene:firstRNA=False
  else:firstcollateral=False
 record.append({'patch':pid,'delta':delta,'area_um2':w['physical_area_um2'],'moves':moved})
ax.set_xlabel('Crop coordinate (µm)',fontsize=12);ax.set_ylabel('Crop coordinate (µm)',fontsize=12);ax.tick_params(labelsize=12);ax.set_aspect('equal');ax.legend(fontsize=12,loc='upper left',bbox_to_anchor=(1,1));ax.set_title('Count 22 → 17 across two compatible patches',fontsize=12)
save(fig,'cell_witness');(R/'results/figure_example.json').write_text(json.dumps({'dataset':'external format-test crop','cell':cell,'gene':gene,'count':17,'min_area_um2':s['minima']['17'],'patch_witnesses':record},indent=2)+'\n')

add={(q['patch'],q['gene']):q for q in json.loads((D/'results/additive_local_spectra.json').read_text()) if q['roi']=='breast_01'}
parts=[]
for pid in s['patches']:
 q=pm[pid];sp=add[(pid,gene)];sign=1 if q['cell_labels'][0]==cell else -1
 parts.append((pid,{str((int(k)-sp['source_left'])*sign):v for k,v in sp['additive_minima'].items()}))
ad=compose(parts,s['source_count']);exact={int(k):F(v) for k,v in s['minima'].items()}
thresholds=sorted({F(0),F(13,20)}|{v for v in exact.values() if v<=F(13,20)}|{v for v,tr in ad.values() if v<=F(13,20)});xx=[float(v) for v in thresholds];lo=[min(k for k,v in exact.items() if v<=x) for x in thresholds];hi=[max(k for k,v in exact.items() if v<=x) for x in thresholds];al=[min(k for k,(v,tr) in ad.items() if v<=x) for x in thresholds]
fig,ax=plt.subplots(figsize=(6.3,3.0));ax.fill_between(xx,lo,hi,alpha=.18,step='post',label='Exact count interval');ax.step(xx,lo,where='post',label='Exact lower endpoint');ax.step(xx,al,where='post',linestyle='--',label='Additive-cost inner endpoint');ax.axvline(.5,linestyle=':',label='Budget 0.5 µm²')
ax.set_ylim(16.5,23.5);ax.set_yticks(range(17,24));ax.set_xlabel('One global deformation budget (µm²)',fontsize=12);ax.set_ylabel('Observed LDHB count',fontsize=12);ax.tick_params(labelsize=12);ax.legend(fontsize=12,loc='lower left');save(fig,'count_curve')

b=json.loads((R/'results/chain_benchmark_summary.json').read_text())['table'];fig,ax=plt.subplots(figsize=(6.3,2.7));ax.plot([q['n'] for q in b],[q['median_total_seconds'] for q in b],marker='o',label='Compile, all counts, and witness checks');ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks([q['n'] for q in b],[str(q['n']) for q in b]);ax.set_xlabel('Explicit point RNAs',fontsize=12);ax.set_ylabel('Median seconds',fontsize=12);ax.tick_params(labelsize=12);save(fig,'chain_scaling')
print('Saved four figures from original masks and exact recorded results.')
