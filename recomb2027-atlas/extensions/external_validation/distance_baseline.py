"""Conservative distance-count outer bound on the same raster reference.

Every boundary point lies at most p/2 from a midpoint of its unit raster edge.
A movement capped at r can cross only molecules within r+p/2 of a reference
edge midpoint. Counting these independently relaxes the ribbons, nuclear pins,
slope and shared area budget; it is an outer bound, not a competing exact solver.
"""
from pathlib import Path
from collections import defaultdict,Counter
from fractions import Fraction as F
import argparse,json,csv,gzip,math
import numpy as np
from scipy.spatial import cKDTree

def main(root):
 a=json.loads(next((root/'data').glob('atlas_*.json')).read_text());p=float(F(a['pixel_size_um']));m=np.load(root/'data/source_crops'/a['roi']/'masks.npz')['cell_mask'];edges=defaultdict(list)
 for axis in [0,1]:
  first,second=(m[:-1,:],m[1:,:]) if axis==0 else (m[:,:-1],m[:,1:]);rr,cc=np.nonzero(first!=second)
  for y,x in zip(rr,cc):
   point=(x+.5,y+1) if axis==0 else (x+1,y+.5)
   for k in [int(first[y,x]),int(second[y,x])]:
    if k>0:edges[k].append(point)
 x0,y0=a['bbox_pixels'][:2];coords=np.array([[float(F(float(r['x_um']))/F(a['pixel_size_um']))-x0,float(F(float(r['y_um']))/F(a['pixel_size_um']))-y0] for r in a['records']]);tree=cKDTree(coords);rmax=3*p*4/(4*math.sqrt(2));outer={}
 for cell in a['complete_cells']:
  found=set(i for ids in tree.query_ball_point(edges[cell],rmax/p+.5+1e-10) for i in ids);outgoing=Counter();incoming=Counter()
  for i in found:
   rec=a['records'][i];(outgoing if rec['cell']==cell else incoming)[rec['gene']]+=1
  outer[cell]=(outgoing,incoming)
 rows=[];tight=0
 for row in csv.DictReader((root/'results/cell_expression_intervals.csv').open()):
  k,g,n=int(row['cell']),row['gene'],int(row['source_count']);outgoing,incoming=outer[k];low,high=(n,n) if F(row['budget_um2'])==0 else (n-outgoing[g],n+incoming[g]);assert low<=int(row['lower'])<=int(row['upper'])<=high
  row.update(distance_lower=low,distance_upper=high);tight+=(low,high)!=(int(row['lower']),int(row['upper']));rows.append(row)
 with gzip.open(root/'results/distance_comparison.csv.gz','wt',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (root/'results/distance_summary.json').write_text(json.dumps({'roi':a['roi'],'displacement_cap_um':rmax,'edge_midpoint_allowance_um':p/2,'rows':len(rows),'exact_strictly_tighter_rows':tight,'failures':0,'interpretation':'Distance-only outer relaxation. It ignores the shared area budget for positive budgets and relaxes frozen geometry and nuclear constraints. Tightening relative to this baseline alone does not demonstrate shared-cost effects.'},indent=2)+'\n')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('roots',type=Path,nargs='+');args=ap.parse_args()
 for r in args.roots:main(r)
