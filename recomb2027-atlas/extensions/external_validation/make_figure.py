from pathlib import Path
import argparse,json
import numpy as np
import pyarrow.parquet as pq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def main(raw,analysis,out):
 plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
 fig,axs=plt.subplots(1,2,figsize=(8.2,4.8),constrained_layout=True)
 for ax,donor,title in zip(axs,['IDC','ILC'],['Ductal carcinoma donor','Lobular carcinoma donor']):
  t=pq.read_table(raw/donor/'cells.parquet',columns=['x_centroid','y_centroid']);x=t['x_centroid'].to_numpy();y=t['y_centroid'].to_numpy();ax.scatter(x[::20],y[::20],s=.35,c='#b8b8b8',rasterized=True)
  select=json.loads((analysis/donor/'selection.json').read_text());from fractions import Fraction
  p=float(Fraction(select['pixel_size_um']))
  for i,r in enumerate(select['regions'],1):
   for key,color,style,label in [('context_bbox_pixels','#266b9b','--','800 µm context'),('bbox_pixels','#c74a20','-','200 µm core')]:
    x0,y0,x1,y1=[p*v for v in r[key]];ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,edgecolor=color,linewidth=1.3,linestyle=style,label=label if i==1 else None))
   ax.annotate(str(i),(p*r['bbox_pixels'][0],p*r['bbox_pixels'][1]),xytext=(3,-12),textcoords='offset points',weight='bold',color='#83331a')
  ax.set(title=title,xlabel='x (µm)',ylabel='y (µm)');ax.set_aspect('equal');ax.invert_yaxis();ax.ticklabel_format(axis='both',style='sci',scilimits=(3,3));ax.spines[['top','right']].set_visible(False)
 axs[1].legend(loc='lower right',fontsize=8,frameon=True);fig.savefig(out);plt.close(fig)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.raw,a.analysis,a.output)
