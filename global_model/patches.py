"""Disjoint, trace-matched two-region patches from a labelled raster.

A single global reference is fitted once. Each accepted patch is a rectangle
containing exactly one ordered interface. Its continuous reference has the same
boundary trace as the raster cross-cut, retains every observed owner and fixed
nuclear region, and is changed only inside that rectangle. Uncovered interfaces
are frozen, not silently treated as uncertain.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict,Counter
import hashlib,json,math,sys
import numpy as np
from scipy.ndimage import label

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'core'))
from mixed import make_problem,verify_mixed

VERSION='global-patches-1'
PARAMETERS={'lengths_in_pixels':[96,48,24,12,6], 'stride_fraction':'1/2',
            'band_padding_um':'1/2','side_guard_pixels':'1/2',
            'minimum_slope':'1','clearance_um':'1/1000',
            'selection':'descending longitudinal length, then ascending box area, then orientation/position/labels',
            'support_overlap':'closed pixel cover disjoint; no RNA counted twice',
            'include_background':True,'fixed_complement':True}

def encode(o):
    if isinstance(o,F):return str(o)
    if isinstance(o,np.integer):return int(o)
    raise TypeError(type(o))

def raster_runs(mask):
    """All rowwise ordered adjacent region pairs, including explicit background."""
    bypair=defaultdict(lambda:defaultdict(list))
    for r,row in enumerate(mask):
        edges=np.flatnonzero(row[1:]!=row[:-1])+1
        limits=np.r_[0,edges,len(row)]
        labs=row[limits[:-1]]
        for j,c in enumerate(edges):
            a,b=int(labs[j]),int(labs[j+1])
            bypair[a,b][r].append((int(c),int(limits[j]),int(limits[j+2])))
    return bypair

def tracks(mask):
    for pair,rows in sorted(raster_runs(mask).items()):
        track=[];prev=-2
        for r,items in sorted(rows.items()):
            if len(items)!=1 or r!=prev+1:
                if len(track)>1:yield pair,track
                track=[]
            if len(items)==1:
                track.append((r,*items[0]));prev=r
            else:prev=-2
        if len(track)>1:yield pair,track

def candidates(mask,h):
    pad=F(1,2)/h
    for orient,arr in [('xy',mask),('yx',mask.T)]:
        for pair,tr in tracks(arr):
            for height in PARAMETERS['lengths_in_pixels']:
                if len(tr)<height:continue
                starts=set(range(0,len(tr)-height+1,max(1,height//2)))
                starts.add(len(tr)-height)
                for start in sorted(starts):
                    sub=tr[start:start+height]
                    rr=[z[0] for z in sub];xx=[z[1] for z in sub]
                    lo=max(F(max(z[2] for z in sub))+F(1,4),F(min(xx))-pad-F(1,2))
                    hi=min(F(min(z[3] for z in sub))-F(1,4),F(max(xx))+pad+F(1,2))
                    if lo>=min(xx) or hi<=max(xx):continue
                    t0=F(rr[0])+F(1,2);t1=F(rr[-1])+F(1,2)
                    box=(lo,hi,t0,t1) if orient=='xy' else (t0,t1,lo,hi)
                    ix0=math.floor(box[0]);ix1=math.ceil(box[1])
                    iy0=math.floor(box[2]);iy1=math.ceil(box[3])
                    yield {'orientation':orient,'pair':pair,'sub':sub,'rails':[lo,hi],
                           'box_pixels':list(box),'cover':[ix0,ix1,iy0,iy1],
                           'priority':(-height,float((hi-lo)*(t1-t0)),orient,rr[0],xx[0],pair)}

def source_val(base,t):
    import bisect
    j=min(len(base)-2,max(0,bisect.bisect_right([q[0] for q in base],t)-1))
    a,u=base[j];b,v=base[j+1]
    return u+(v-u)*(t-a)/(b-a)

def build_patch(cand,roi,mask,nuc,records,h,origin,rec_rows,rec_cols):
    orient=cand['orientation'];arr=mask if orient=='xy' else mask.T
    nm=nuc if orient=='xy' else nuc.T
    pair=cand['pair'];low,high=cand['rails'];sub=cand['sub']
    tbase=F(sub[0][0])+F(1,2);end=F(sub[-1][0])+F(1,2)
    xbase=low
    base=[((F(r)+F(1,2)-tbase)*h,(F(x)-xbase)*h) for r,x,_,_ in sub]
    band=[max(F(0)+h/F(4),min(v for _,v in base)-F(1,2)),
          min((high-low)*h-h/F(4),max(v for _,v in base)+F(1,2))]
    if not band[0]<min(v for _,v in base) or not band[1]>max(v for _,v in base):
        if band[0]>min(v for _,v in base) or band[1]<max(v for _,v in base):
            return None,'band'
    slope=max([F(1)]+[abs(v-u)/(b-a) for (a,u),(b,v) in zip(base,base[1:])])
    box=list(map(F,cand['box_pixels']))
    points=[]
    ix0,ix1,iy0,iy1=cand['cover']
    possible=np.flatnonzero((rec_rows>=iy0)&(rec_rows<iy1)&(rec_cols>=ix0)&(rec_cols<ix1))
    for ri in possible:
        z=records[int(ri)]
        x=(F(z['x_um'])-origin[0])/h;y=(F(z['y_um'])-origin[1])/h
        if not(box[0]<=x<=box[1] and box[2]<=y<=box[3]):continue
        if z['cell'] not in pair:return None,'third_owner'
        tx,tt=(x,y) if orient=='xy' else (y,x)
        tx=(tx-xbase)*h;tt=(tt-tbase)*h
        side=pair.index(z['cell']);fv=source_val(base,tt)
        if (side==0 and tx+F(1,1000)>fv) or (side==1 and tx-F(1,1000)<fv):
            return None,'source_rna'
        points.append({'transcript_id':z['id'],'gene':z['gene'],'source':side,
                       'x':str(tx),'t':str(tt),'pixel_col':z['pixel_col'],'pixel_row':z['pixel_row']})
    boxes=[[],[]]
    for r in range(math.floor(tbase),math.ceil(end)):
        for c in range(math.floor(low),math.ceil(high)):
            if not nm[r,c]:continue
            owner=int(arr[r,c])
            if owner not in pair:return None,'third_nucleus'
            x0=max(low,F(c));x1=min(high,F(c+1))
            y0=max(tbase,F(r));y1=min(end,F(r+1))
            if x1<=x0 or y1<=y0:continue
            x0=(x0-xbase)*h;x1=(x1-xbase)*h;y0=(y0-tbase)*h;y1=(y1-tbase)*h
            vals=[source_val(base,t) for t in [y0,y1]+[t for t,_ in base if y0<t<y1]]
            side=pair.index(owner)
            if (side==0 and x1+F(1,1000)>min(vals)) or (side==1 and x0-F(1,1000)<max(vals)):
                return None,'source_nucleus'
            boxes[side].append(list(map(str,[x0,x1,y0,y1])))
    raw={'roi':roi,'orientation':orient,'cell_labels':list(pair),'pixel_size_um':str(h),
         'baselines':[[list(map(str,z)) for z in base]],'slopes':[str(slope)],'bands':[list(map(str,band))],
         'minimum_gaps':[],'maximum_gaps':[],'points':points,'nuclear_boxes':boxes,
         'global_box_pixels':list(map(str,box)),'global_box_um':list(map(str,[origin[0]+box[0]*h,origin[0]+box[1]*h,origin[1]+box[2]*h,origin[1]+box[3]*h])),
         'local_frame_origin_pixels':list(map(str,[xbase,tbase])),
         'cover':cand['cover'],'source_rows':[[r,x] for r,x,_,_ in sub],
         'longitudinal_length_um':str((end-tbase)*h),
         'reference_arclength_um_float':sum(float(((b-a)**2+(v-u)**2))**.5 for (a,u),(b,v) in zip(base,base[1:])),
         'source_fit_area_um2':str(sum(abs(F(b[1])-F(a[1]))*h*h/F(4) for a,b in zip(sub,sub[1:])))}
    raw['id']=hashlib.sha256(json.dumps(raw,sort_keys=True).encode()).hexdigest()[:16]
    return raw,None

def construct(roi_dir,observed):
    roi=roi_dir.name;meta=json.loads((roi_dir/'metadata.json').read_text());npz=np.load(roi_dir/'masks.npz')
    mask=npz['cell_mask'];nuc=npz['nucleus_mask']
    h=1/F(str(meta['transform'][0][0]));bbox=meta['bbox_pixels']
    origin=(F(bbox[0])*h,F(bbox[1])*h)
    records=observed['records']
    for z in records:
        assert int(mask[z['pixel_row'],z['pixel_col']])==z['cell']
    rec_rows=np.array([z['pixel_row'] for z in records]);rec_cols=np.array([z['pixel_col'] for z in records])
    cand=list(candidates(mask,h));cand.sort(key=lambda z:z['priority'])
    occupied=np.zeros(mask.shape,bool);accepted=[];rejections=Counter();used=set()
    for c in cand:
        ix0,ix1,iy0,iy1=c['cover']
        if occupied[iy0:iy1,ix0:ix1].any():rejections['overlap']+=1;continue
        p,why=build_patch(c,roi,mask,nuc,records,h,origin,rec_rows,rec_cols)
        if why:rejections[why]+=1;continue
        ids={z['transcript_id'] for z in p['points']}
        if used&ids:raise AssertionError('RNA support repeated')
        used|=ids;accepted.append(p);occupied[iy0:iy1,ix0:ix1]=True
    border=set(map(int,np.unique(np.r_[mask[0],mask[-1],mask[:,0],mask[:,-1]])))
    cells=[]
    for k in sorted(set(map(int,np.unique(mask)))-{0}):
        comp,n=label(mask==k)
        cells.append({'cell':k,'pixels':int((mask==k).sum()),'border_clipped':k in border,
                      'components_4':n,'complete_observed_cell':k not in border and n==1})
    boundary_edges=int(np.count_nonzero(mask[1:]!=mask[:-1])+np.count_nonzero(mask[:,1:]!=mask[:,:-1]))
    covered=sum(len(p['source_rows'])-1+sum(abs(b[1]-a[1]) for a,b in zip(p['source_rows'],p['source_rows'][1:])) for p in accepted)
    stats={'roi':roi,'candidate_rectangles':len(cand),'accepted_patches':len(accepted),'rejected':dict(rejections),
           'transcripts':len(records),'patch_transcripts':len(used),'all_boundary_grid_edges':boundary_edges,
           'covered_grid_edge_length':covered,'covered_fraction_grid_length':covered/boundary_edges,
           'total_reference_fit_area_um2':str(sum(F(p['source_fit_area_um2']) for p in accepted)),
           'connected_complete_cells':sum(x['complete_observed_cell'] for x in cells),
           'cells':cells,'pixel_size_um':str(h),'source_origin_um':list(map(str,origin)),
           'zero_rna_patches':sum(not p['points'] for p in accepted)}
    return accepted,stats

if __name__=='__main__':
    observed=json.loads((ROOT/'data/observed_crop_counts.json').read_text())['rois']
    out=[];stats=[]
    for roi in sorted(observed):
        pp,ss=construct(ROOT/'data/source_crops'/roi,observed[roi]);out+=pp;stats.append(ss)
        print({k:v for k,v in ss.items() if k!='cells'},flush=True)
    payload={'schema':VERSION,'parameters':PARAMETERS,'patches':out,'census':stats,
             'scope':'One globally glued patch family per ROI; all complement fixed; background is an explicit compartment.'}
    text=json.dumps(payload,indent=2,default=encode)+'\n'
    (ROOT/'data/global_patches.json').write_text(text)
    print('global family sha256',hashlib.sha256(text.encode()).hexdigest(),flush=True)
