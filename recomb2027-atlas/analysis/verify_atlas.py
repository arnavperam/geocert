"""Recheck atlas acceptance using rational polygon/pixel operations.

Uses the supplied original masks and transcript coordinates, not the proposal
adapter's GEOS decisions. Checks membership of *every* retained transcript.
"""
from pathlib import Path
import sys,json,hashlib,math,time
from fractions import Fraction as F
from collections import Counter
import numpy as np
import os
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
from atlas.geometry import quad,check_pixels,positive_overlap,intersection,area2,evaluate_path

def bounds(points):
 return min(x for x,y in points),min(y for x,y in points),max(x for x,y in points),max(y for x,y in points)

def main():
 t0=time.time();total={'patches':0,'rational_half_ribbon_pixel_checks':0,'rational_patch_pair_checks':0,'raw_RNA_memberships':0,
                    'editable_RNA':0,'nuclear_fragments_rechecked':0,'failures':0};per=[]
 for af in sorted((ROOT/'data').glob('atlas_*.json')):
  a=json.loads(af.read_text());aid=a.pop('atlas_id');assert hashlib.sha256(json.dumps(a,sort_keys=True,separators=(',',':')).encode()).hexdigest()==aid;a['atlas_id']=aid
  d=ROOT/'data/source_crops'/a['roi'];masks=np.load(d/'masks.npz');mask=masks['cell_mask'];nuc=masks['nucleus_mask'];meta=json.loads((d/'metadata.json').read_text())
  for fname,sha in a['source_hashes'].items():assert hashlib.sha256((d/fname).read_bytes()).hexdigest()==sha
  p=1/F(float(meta['transform'][0][0]));assert F(a['pixel_size_um'])==p
  seen=[];stored={};row_candidates={}
  for pat in a['patches']:
   pts=pat['source_path_pixels'];sx,sy=pat['sx'],pat['sy'];width=F(pat['halfwidth_pixels']);sh=(sx*width/2,-sy*width/2)
   assert abs(sx)==abs(sy)==1 and F(pat['area_scale_um2'])==p*p/2
   assert 0<F(pat['band_halfwidth_pixels'])<width
   origin=sx*pts[0][0]+sy*pts[0][1];qseq=[];nuc_check=[]
   for v,w in zip(pts,pts[1:]):
    assert abs(w[0]-v[0])+abs(w[1]-v[1])==1
    ta=sx*v[0]+sy*v[1];tb=sx*w[0]+sy*w[1];assert tb>ta
    ua=sx*v[0]-sy*v[1];ub=sx*w[0]-sy*w[1];s=F(ub-ua,tb-ta);off=ua-s*ta
    for lo,hi,side in [(-1,0,0),(0,1,1)]:
     q=quad((v,w),sh,lo,hi)
     assert check_pixels(q,mask,{pat['cell_labels'][side]})
     total['rational_half_ribbon_pixel_checks']+=1
     x0,y0,x1,y1=bounds(q)
     for rr in range(math.floor(y0),math.ceil(y1)):
      for cc in range(math.floor(x0),math.ceil(x1)):
       if not (0<=rr<nuc.shape[0] and 0<=cc<nuc.shape[1]) or not nuc[rr,cc]:continue
       clip=intersection(q,[(cc,rr),(cc+1,rr),(cc+1,rr+1),(cc,rr+1)])
       if len(clip)<3 or area2(clip)==0:continue
       tv=[sx*x+sy*y for x,y in clip];hv=[sx*x-sy*y-(s*t+off) for (x,y),t in zip(clip,tv)]
       value=max(hv) if side==0 else min(hv);left,right=min(tv)-origin,max(tv)-origin
       # A stored constant obligation covering this fragment must be as strong.
       covering=[(F(l),F(r)) for k,b,l,r in pat['nuclear_constraints'] if k==side and (F(b)>=value if side==0 else F(b)<=value)]
       pos=left
       for l,r in sorted(covering):
        if l<=pos:pos=max(pos,r)
       assert pos>=right, ('nucleus constraint omitted',pat['id'],side,value,left,right)
       total['nuclear_fragments_rechecked']+=1
    qseq.append(quad((v,w),sh,-1,1))
   bbox=bounds([v for q in qseq for v in q])
   for oldbbox,oldqs in seen:
    if bbox[2]<=oldbbox[0] or oldbbox[2]<=bbox[0] or bbox[3]<=oldbbox[1] or oldbbox[3]<=bbox[1]:continue
    total['rational_patch_pair_checks']+=1
    for q in qseq:
     for qq in oldqs:assert not positive_overlap(q,qq)
   seen.append((bbox,qseq));total['patches']+=1
   for row in range(math.floor(bbox[1]),math.floor(bbox[3])+1):row_candidates.setdefault(row,[]).append(pat)
   for r in pat['points']:
    assert r['transcript_id'] not in stored;stored[r['transcript_id']]=(pat,r)
  editable=0;x0,y0=a['bbox_pixels'][:2]
  for r in a['records']:
   x=F(float(r['x_um']))/p-x0;y=F(float(r['y_um']))/p-y0
   assert int(mask[math.floor(y),math.floor(x)])==r['cell']
   found=[]
   for pat in row_candidates.get(math.floor(y),[]):
    sx,sy=pat['sx'],pat['sy'];pts=pat['source_path_pixels'];t=sx*x+sy*y;origin=sx*pts[0][0]+sy*pts[0][1]
    if not origin<t<origin+pat['domain_pixels'][1]:continue
    u=sx*x-sy*y-evaluate_path(pts,sx,sy,t)
    if abs(u)<F(pat['halfwidth_pixels']):found.append((pat['id'],t-origin,u,0 if u<0 else 1))
   assert len(found)<=1
   if found:
    pp,rr=stored[r['id']];pid,t,u,k=found[0]
    assert (pid,t,u,k)==(pp['id'],F(rr['t']),F(rr['x']),rr['source'])
    assert r['cell']==pp['cell_labels'][k];editable+=1
   else:assert r['id'] not in stored
   total['raw_RNA_memberships']+=1
  assert editable==len(stored);total['editable_RNA']+=editable
  per.append({'roi':a['roi'],'atlas_id':aid,'eligible_reference_edges':a['summary']['eligible_shared_edges'],
     'editable_edges':a['summary']['selected_source_edges'],'fraction':a['summary']['selected_source_edges']/a['summary']['eligible_shared_edges'],
     'editable_length_um':float(p*a['summary']['selected_source_edges']),
     'max_patch_arclength_um':float(p*max(z['length_edges'] for z in a['patches']))})
 total['seconds']=time.time()-t0;total['coverage']=per
 (ROOT/'results/atlas_geometric_validation.json').write_text(json.dumps(total,indent=2)+'\n');print(total)
if __name__=='__main__':main()
