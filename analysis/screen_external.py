"""Locked transfer screen: unchanged corridor rules, standard CSV ingestion.
The extraction logic follows the existing development adapter.
"""
from pathlib import Path
import sys,json,hashlib
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'analysis/external_inputs'
sys.path.insert(0,str(ROOT/'core'))
import pandas as pd
from layered import value

def runs(a):
 st=np.r_[0,np.where(a[1:]!=a[:-1])[0]+1]; en=np.r_[st[1:],len(a)]
 return [(int(a[s]),int(s),int(e)) for s,e in zip(st,en)]

def squash(knots):
 out=[]
 for t,x in knots:
  out.append((t,x))
  while len(out)>2:
   a,u=out[-3];b,v=out[-2];c,w=out[-1]
   if (v-u)*(c-b)==(w-v)*(b-a):out.pop(-2)
   else:break
 return out

def main(height=24,ncells=3,keep_nuc=True):
 (ROOT/'results/external').mkdir(parents=True,exist_ok=True)
 records=[];all_candidates=[];audit={}
 for roi in ['breast_01','breast_02']:
  p=INPUT/roi;meta=json.loads((p/'metadata.json').read_text());z=np.load(p/'masks.npz');mask=z['cell_mask'];nuc=np.where(z['nucleus_mask']>0,mask,0);df=pd.read_csv(p/'transcripts.csv.gz',float_precision='round_trip')
  # Preserve original float coordinates exactly as binary rationals, and stored transform.
  scale=F(float(meta['transform'][0][0]));pix=1/scale;x0,y0,*_=meta['bbox_pixels'];
  xs=np.array([float(F(float(x))*scale-x0) for x in df.x_location]);ys=np.array([float(F(float(y))*scale-y0) for y in df.y_location]);
  checks=int(np.sum(mask[np.floor(ys).astype(int),np.floor(xs).astype(int)]!=df.mask_label.to_numpy()))
  if checks:raise AssertionError('native source labels mismatch')
  stats={'source_RNA':len(df),'source_owner_mismatches':checks,'raw_triples':0,'topology_reject':0,'rail_reject':0,'nuclear_band_reject':0,'nuclear_interpolant_reject':0,'RNA_count_reject':0,'interpolation_owner_reject':0,'clearance_reject':0,'eligible':0}
  candidates=[]
  for orientation in [0,1]:
   mm=mask if orientation==0 else mask.T;nn=nuc if orientation==0 else nuc.T;xx,tt=(xs,ys) if orientation==0 else (ys,xs)
   for row in range(0,mm.shape[0]-height,height):
    rr=[runs(mm[j]) for j in range(row,row+height)]
    first=rr[0]
    for k in range(len(first)-ncells+1):
     trip=tuple(x[0] for x in first[k:k+ncells])
     if 0 in trip or len(set(trip))!=ncells:continue
     stats['raw_triples']+=1
     seq=[]
     for r in rr:
      matches=[r[j:j+ncells] for j in range(len(r)-ncells+1) if tuple(x[0] for x in r[j:j+ncells])==trip]
      if len(matches)!=1:break
      seq.append(matches[0])
     if len(seq)!=height:stats['topology_reject']+=1;continue
     left=max(s[0][1] for s in seq);right=min(s[-1][2] for s in seq)
     arrays=[[s[k][1] for s in seq] for k in range(1,ncells)];b1=arrays[0];b2=arrays[-1]
     pad=F(1,2)/pix
     if min(b1)-pad<=left or max(b2)+pad>=right:stats['rail_reject']+=1;continue
     bs=[]
     for arr in arrays:bs.append(squash([(F(j,1)*pix,(F(v)-left)*pix) for j,v in enumerate(arr)]))
     bands=[(min(x for _,x in b)-F(1,2),max(x for _,x in b)+F(1,2)) for b in bs]
     # Conservatively require every pixel intersecting either variable band to be nonnuclear.
     bad=False
     for lo,hi in bands:
      xa=max(0,int(np.floor(float(lo/pix+left))));xb=min(mm.shape[1],int(np.ceil(float(hi/pix+left))))
      if np.any(nn[row:row+height,xa:xb]):bad=True;break
     if bad and not keep_nuc:stats['nuclear_band_reject']+=1;continue
     inds=np.where((xx>=left)&(xx<right)&(tt>=row+.5)&(tt<=row+height-.5))[0]
     if not 3<=len(inds)<=24:stats['RNA_count_reject']+=1;continue
     points=[];mismatch=clear=False
     for i in inds:
      # Recover exact physical coordinates from the original floating-point values.
      rawx=F(float(df.iloc[i].x_location));rawy=F(float(df.iloc[i].y_location))
      gx,gt=(rawx,rawy) if orientation==0 else (rawy,rawx)
      offx=F(x0 if orientation==0 else y0)+left;offt=F(y0 if orientation==0 else x0)+F(2*row+1,2)
      x=gx-offx/scale;t=gt-offt/scale;v=[value(b,t) for b in bs]
      owner=sum(x>u for u in v)
      label=int(df.iloc[i].mask_label)
      if label not in trip or label!=trip[owner]:mismatch=True;break
      if min(abs(x-u) for u in v)<F(1,1000):clear=True;break
      points.append({'x':str(x),'t':str(t),'source':owner,'transcript_id':str(df.iloc[i].transcript_id),'gene':str(df.iloc[i].feature_name),'mask_label':label,'nuclear':bool(df.iloc[i].overlaps_nucleus)})
     if mismatch:stats['interpolation_owner_reject']+=1;continue
     if clear:stats['clearance_reject']+=1;continue
     gaps=[[(F(b)-F(a))*pix for a,b in zip(aa,bb)] for aa,bb in zip(arrays,arrays[1:])]
     slopes=[max([F(1)]+[abs(v-u)/(b-a) for (a,u),(b,v) in zip(bs0,bs0[1:])]) for bs0 in bs]
     nuclear_boxes=[[] for _ in trip];nuc_bad=False
     if keep_nuc:
      # Nuclear pixels are square physical regions, clipped at chart endpoints.
      for jj in range(height):
       tlo=max(F(0),(F(jj)-F(1,2))*pix);thi=min(F(height-1)*pix,(F(jj)+F(1,2))*pix)
       if tlo>=thi:continue
       for label,aa,bb in runs(nn[row+jj,left:right]):
        if label==0:continue
        if label not in trip:nuc_bad=True;break
        owner=trip.index(label);xx0=F(aa)*pix;xx1=F(bb)*pix
        for kk,b in enumerate(bs):
         vals=[value(b,t) for t in [tlo,thi]+[u for u,_ in b if tlo<u<thi]]
         if (owner<=kk and min(vals)<xx1+F(1,1000)) or (owner>kk and max(vals)>xx0-F(1,1000)):nuc_bad=True;break
        if nuc_bad:break
        nuclear_boxes[owner].append([str(xx0),str(xx1),str(tlo),str(thi)])
       if nuc_bad:break
     if nuc_bad:stats['nuclear_interpolant_reject']+=1;continue
     key=f'{roi}_{orientation}_{row}_{left}_{right}_{trip}'
     c={'id':hashlib.sha256(key.encode()).hexdigest()[:14],'roi':roi,'orientation':orientation,'row':row,'left':left,'right':right,'cell_labels':trip,'pixel_size_um':str(pix),'baselines':[[[str(t),str(x)] for t,x in b] for b in bs],'bands':[[str(a),str(b)] for a,b in bands],'slopes':list(map(str,slopes)),'minimum_gaps':[str(min(gg)/2) for gg in gaps],'maximum_gaps':[str(max(gg)*2) for gg in gaps],'points':points,'source_rows':list(range(row,row+height))}
     c['nuclear_boxes']=nuclear_boxes;c['height']=height;c['cells']=ncells;candidates.append(c);stats['eligible']+=1
  # Complete eligibility ledger retained; disjoint rectangular support selected by seeded hash ordering.
  chosen=[];rects=[]
  for c in sorted(candidates,key=lambda z:z['id']):
   if c['orientation']==0:rect=(c['left'],c['row'],c['right'],c['row']+height)
   else:rect=(c['row'],c['left'],c['row']+height,c['right'])
   if any(max(rect[0],r[0])<min(rect[2],r[2]) and max(rect[1],r[1])<min(rect[3],r[3]) for r in rects):continue
   chosen.append(c);rects.append(rect)
   if len(chosen)==12:break
  stats['selected']=len(chosen);stats['selection_rule']='lexicographic SHA256, nonoverlapping, at most 12 per ROI'
  records.extend(chosen);all_candidates.extend(candidates);audit[roi]=stats
 (ROOT/f'results/external/charts_h{height}_c{ncells}.json').write_text(json.dumps(records,indent=2))
 (ROOT/f'results/external/eligibility_h{height}_c{ncells}.json').write_text(json.dumps(audit,indent=2))
 (ROOT/f'results/external/all_eligible_h{height}_c{ncells}.json').write_text(json.dumps(all_candidates,indent=2))
 print(json.dumps(audit,indent=2));print('selected',len(records), 'sites',sum(len(r['points']) for r in records))
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--height',type=int,default=24);a.add_argument('--cells',type=int,default=3);a.add_argument('--keep-nuc',action='store_true',default=True);v=a.parse_args();main(v.height,v.cells,v.keep_nuc)
