from pathlib import Path
from fractions import Fraction
import json,numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=next(p for p in Path(__file__).resolve().parents if (p/'atlas_cli.py').is_file());O=Path(__file__).parent/'figures'
a=json.loads((R/'data/atlas_roi_01.json').read_text());m=np.load(R/'data/source_crops/roi_01/masks.npz')['cell_mask'];p=float(Fraction(a['pixel_size_um']));bbox=a['bbox_pixels']
fig,ax=plt.subplots(figsize=(6.3,3.35));ids=np.unique(m);ax.imshow(np.searchsorted(ids,m),origin='upper',extent=(0,m.shape[1]*p,m.shape[0]*p,0),interpolation='nearest',alpha=.25,cmap='tab20')
for q in a['patches']:
 xy=np.asarray(q['source_path_pixels'])*p;ax.plot(xy[:,0],xy[:,1],color='#1966ac',linewidth=.9,zorder=3)
xy=np.array([[float(Fraction(r['x_um']))-bbox[0]*p,float(Fraction(r['y_um']))-bbox[1]*p] for r in a['records']]);inside=(xy[:,0]>=45)&(xy[:,0]<=125)&(xy[:,1]>=0)&(xy[:,1]<=44);xy=xy[inside]
ax.scatter(xy[:,0],xy[:,1],s=4.5,color='#292929',alpha=.65,zorder=4,linewidths=0)
ax.plot([49,59],[40,40],color='black',lw=2.6,zorder=5);ax.text(54,38.7,'10 µm',ha='center',va='bottom',fontsize=11,bbox=dict(fc='white',alpha=.8,ec='none',pad=.5))
ax.set_xlim(45,125);ax.set_ylim(44,0);ax.set_aspect('equal');ax.set_xlabel('Crop x (µm)',fontsize=12);ax.set_ylabel('Crop y (µm)',fontsize=12);ax.tick_params(labelsize=11)
ax.set_title('Lung reference mask, transcripts and editable arcs',fontsize=12)
handles=[Line2D([0],[0],color='#1966ac',lw=1.5,label='Editable source arc'),Line2D([0],[0],color='#292929',marker='o',linestyle='None',markersize=3,label='Retained transcript')]
ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=11)
fig.tight_layout();fig.savefig(O/'atlas_reference.pdf');fig.savefig(O/'atlas_reference.png',dpi=170);plt.close(fig)
print('plotted',len(xy),'retained transcripts at original measured coordinates')
