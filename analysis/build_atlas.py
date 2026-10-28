"""Geometry-only atlas extraction from all interfaces of the supplied crops.

Proposals use rasterio/GEOS; exact rational checks independently certify the
accepted pixels and disjoint interiors. No molecular outcomes select patches.
"""
from pathlib import Path
import sys, json, hashlib, math, time
from fractions import Fraction as F
from collections import defaultdict, Counter
import numpy as np
from rasterio.features import shapes
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union
import os
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
from atlas.geometry import quad, ribbon_pieces, check_pixels, positive_overlap, evaluate_path, intersection, area2


def paths_from_edges(edges):
    adj=defaultdict(list)
    for i,(a,b) in enumerate(edges):adj[a].append((b,i));adj[b].append((a,i))
    seen=set();out=[]
    for start in sorted(adj, key=lambda p:(len(adj[p])==2,p)):
        for other,i in sorted(adj[start]):
            if i in seen:continue
            seen.add(i);path=[start,other];previous=start;current=other
            while len(adj[current])==2:
                nxt=[(p,j) for p,j in adj[current] if j not in seen]
                if not nxt:break
                p,j=nxt[0];seen.add(j);path.append(p);previous,current=current,p
                if current==start:break
            out.append(path)
    return out


def monotone_runs(path):
    choices={(1,1),(1,-1),(-1,1),(-1,-1)}
    run=[path[0]];ok=set(choices)
    for point in path[1:]:
        prev=run[-1];dx=point[0]-prev[0];dy=point[1]-prev[1]
        new={s for s in ok if s[0]*dx+s[1]*dy>0}
        if not new or len(run)>=129:
            if len(run)>=7:yield run, sorted(ok)
            run=[prev];ok=set(choices)
            new={s for s in ok if s[0]*dx+s[1]*dy>0}
        run.append(point);ok=new
    if len(run)>=7:yield run,sorted(ok)


def polygon_for(points, sx,sy,b,side=0):
    dx=sx*float(b)/2;dy=-sy*float(b)/2
    minus=[(x-dx,y-dy) for x,y in points]
    plus=[(x+dx,y+dy) for x,y in points]
    if side<0:return Polygon(minus+list(reversed(points)))
    if side>0:return Polygon(list(points)+list(reversed(plus)))
    return Polygon(minus+list(reversed(plus)))


def main():
    data=json.loads((ROOT/'data/observed_crop_counts.json').read_text())
    all_rois=[];start=time.time()
    for d in sorted((ROOT/'data/source_crops').iterdir()):
        if os.environ.get('GEOCERT_ROI') and d.name!=os.environ['GEOCERT_ROI']:continue
        print('start',d.name,flush=True)
        masks=np.load(d/'masks.npz');mask=masks['cell_mask'];nuc=masks['nucleus_mask']
        meta=json.loads((d/'metadata.json').read_text());h,w=mask.shape
        p=1/F(float(meta['transform'][0][0]))
        border=set(map(int,np.unique(np.r_[mask[0,:],mask[-1,:],mask[:,0],mask[:,-1]])))
        cells=defaultdict(list)
        for geo,val in shapes(mask.astype(np.int32)):
            cells[int(val)].append(shape(geo))
        eligible={k for k,v in cells.items() if k>0 and len(v)==1 and v[0].is_valid and k not in border}
        cellpoly={k:v[0] for k,v in cells.items() if k in eligible}
        cellpoly[0]=unary_union(cells.get(0,[]))
        nps=[shape(g) for g,v in shapes((nuc>0).astype(np.uint8)) if v]
        nuclear=unary_union(nps)
        edges=defaultdict(list)
        all_shared=0
        rr,cc=np.nonzero(mask[:,:-1]!=mask[:,1:])
        for r,c in zip(rr.tolist(),cc.tolist()):
            a,b=int(mask[r,c]),int(mask[r,c+1])
            if a!=b:all_shared+=1
            if (a in eligible and (b in eligible or b==0)) or (b in eligible and a==0):
                edges[tuple(sorted((a,b)))].append(((c+1,r),(c+1,r+1)))
        rr,cc=np.nonzero(mask[:-1,:]!=mask[1:,:])
        for r,c in zip(rr.tolist(),cc.tolist()):
            a,b=int(mask[r,c]),int(mask[r+1,c])
            if a!=b:all_shared+=1
            if (a in eligible and (b in eligible or b==0)) or (b in eligible and a==0):
                edges[tuple(sorted((a,b)))].append(((c,r+1),(c+1,r+1)))
        candidates=[];chain_count=0
        for pair,ee in sorted(edges.items()):
            for path in paths_from_edges(ee):
                chain_count+=1
                for run,orientations in monotone_runs(path):
                    # A frozen edge at both ends separates patches and junctions.
                    points=run[1:-1]
                    if len(points)>=5:candidates.append((pair,points,orientations))
        candidates.sort(key=lambda z:(-len(z[1]),z[0],z[1],z[2]))
        accepted=[];regions=[];geoms=[];reject=Counter();trials=Counter();exact_checks=0
        print('geometry_candidates',len(candidates),flush=True)
        for ci,(pair,points,orientations) in enumerate(candidates):
            if ci%100==0:print('candidate',d.name,ci,'accepted',len(accepted),'seconds',time.time()-start,flush=True)
            chosen=None;reason='no_orientation'
            for sx,sy in orientations:
                a,b=points[0],points[1];mid=((a[0]+b[0])/2,(a[1]+b[1])/2)
                lower=int(mask[math.floor(mid[1]+sy*.1),math.floor(mid[0]-sx*.1)])
                upper=int(mask[math.floor(mid[1]-sy*.1),math.floor(mid[0]+sx*.1)])
                if {lower,upper}!=set(pair):raise AssertionError('edge-side labels')
                for width in tuple(F(x) for x in os.environ.get('GEOCERT_WIDTHS','4,2,1,1/2').split(',')):
                    trials['widths']+=1
                    poly=polygon_for(points,sx,sy,width)
                    low=polygon_for(points,sx,sy,width,-1);high=polygon_for(points,sx,sy,width,1)
                    if not poly.is_valid or poly.area<=0:reason='invalid_ribbon';continue
                    if low.difference(cellpoly[lower]).area>1e-9 or high.difference(cellpoly[upper]).area>1e-9:
                        reason='third_label_or_non_graph_source';continue
                    # Nuclear ownership is retained by explicit zero-clearance structural constraints.
                    overlaps=[j for j,g in enumerate(geoms) if poly.intersects(g)]
                    if any(poly.intersection(geoms[j]).area>1e-9 for j in overlaps):reason='patch_overlap';continue
                    # Exact acceptance checks. Widths and mask vertices are dyadic.
                    shift=(F(sx)*width/2,-F(sy)*width/2)
                    quads=[]
                    for edge in zip(points,points[1:]):
                        for lo,hi,owner in [(-1,0,lower),(0,1,upper)]:
                            q=quad(edge,shift,lo,hi)
                            if not check_pixels(q,mask,{owner}):raise AssertionError('GEOS/rational pixel discrepancy')
                            exact_checks+=1
                        quads.append(quad(edge,shift,-1,1))
                    for j in overlaps:
                        for q in quads:
                            for other in regions[j]:
                                if positive_overlap(q,other):raise AssertionError('GEOS/rational disjointness discrepancy')
                    chosen=(sx,sy,width,lower,upper,poly,quads);break
                if chosen:break
            if chosen is None:reject[reason]+=1;continue
            sx,sy,width,lower,upper,poly,quads=chosen
            pid=hashlib.sha256(json.dumps([d.name,pair,points,sx,sy,str(width)],separators=(',',':')).encode()).hexdigest()[:16]
            accepted.append({'id':pid,'source_path_pixels':points,'sx':sx,'sy':sy,'halfwidth_pixels':str(width),
                'band_halfwidth_pixels':str(F(3,4)*width),'slope':'1','area_scale_um2':str(p*p/2),
                'cell_labels':[lower,upper],'length_edges':len(points)-1,'points':[], 'nuclear_constraints':[]})
            geoms.append(poly);regions.append(quads)
            constraints=[]
            for pa,pb in zip(points,points[1:]):
                sh=(F(sx)*width/2,-F(sy)*width/2)
                ta=sx*pa[0]+sy*pa[1];tb=sx*pb[0]+sy*pb[1]
                ua=sx*pa[0]-sy*pa[1];ub=sx*pb[0]-sy*pb[1]
                slope=F(ub-ua,tb-ta);intercept=F(ua)-slope*ta
                origin=sx*points[0][0]+sy*points[0][1]
                for l,h,side in [(-1,0,0),(0,1,1)]:
                    q=quad((pa,pb),sh,l,h)
                    for rr in range(math.floor(min(y for x,y in q)),math.ceil(max(y for x,y in q))):
                        for cc in range(math.floor(min(x for x,y in q)),math.ceil(max(x for x,y in q))):
                            if not (0<=rr<nuc.shape[0] and 0<=cc<nuc.shape[1]) or not nuc[rr,cc]:continue
                            clip=intersection(q,[(cc,rr),(cc+1,rr),(cc+1,rr+1),(cc,rr+1)])
                            if len(clip)<3 or area2(clip)==0:continue
                            tv=[sx*x+sy*y for x,y in clip]
                            hv=[sx*x-sy*y-(slope*t+intercept) for (x,y),t in zip(clip,tv)]
                            bound=max(hv) if side==0 else min(hv)
                            assert (bound<=0 if side==0 else bound>=0)
                            constraints.append((side,bound,min(tv)-origin,max(tv)-origin))
            # Exact interval-union and dominance compression for constant bounds.
            grouped=defaultdict(list)
            for side,bound,a,b in constraints:grouped[(side,bound)].append((a,b))
            merged=[]
            for (side,bound),intervals in sorted(grouped.items()):
                tmp=[]
                for a,b in sorted(intervals):
                    if tmp and a<=tmp[-1][1]:tmp[-1]=(tmp[-1][0],max(b,tmp[-1][1]))
                    else:tmp.append((a,b))
                merged.extend((side,bound,a,b) for a,b in tmp)
            retained=[]
            for side,bound,a,b in merged:
                if any(ss==side and aa<=a and bb>=b and ((vv>bound) if side==0 else (vv<bound)) for ss,vv,aa,bb in merged):continue
                retained.append([side,str(bound),str(a),str(b)])
            accepted[-1]['nuclear_constraints']=retained
            accepted[-1]['unmerged_nuclear_fragments']=len(constraints)
        # Molecular input is consulted only after every candidate has been selected.
        records=data['rois'][d.name]['records'];matched=set();near_boundary=[]
        buckets=defaultdict(list)
        for j,poly in enumerate(geoms):
            bx0,by0,bx1,by1=poly.bounds
            for by in range(math.floor(by0/16),math.floor(by1/16)+1):
                for bx in range(math.floor(bx0/16),math.floor(bx1/16)+1):buckets[(bx,by)].append(j)
        x0,y0=meta['bbox_pixels'][:2]
        for ri,r in enumerate(records):
            if ri%10000==0:print('observations',d.name,ri,'seconds',time.time()-start,flush=True)
            x=F(float(r['x_um']))/p-x0;y=F(float(r['y_um']))/p-y0
            found=[]
            for j in buckets[(math.floor(x/16),math.floor(y/16))]:
                a=accepted[j]
                sx,sy=a['sx'],a['sy'];t=sx*x+sy*y
                pts=a['source_path_pixels'];ta=sx*pts[0][0]+sy*pts[0][1];tb=sx*pts[-1][0]+sy*pts[-1][1]
                if not ta<t<tb:continue
                u=sx*x-sy*y-evaluate_path(pts,sx,sy,t)
                if not -F(a['halfwidth_pixels'])<u<F(a['halfwidth_pixels']):continue
                if u==0:raise AssertionError('RNA exactly on reference boundary; explicit tie policy needed')
                source=0 if u<0 else 1
                if a['cell_labels'][source]!=r['cell']:raise AssertionError('source ownership differs from raster')
                rho=F(1,1000000)/p
                if abs(u)<rho:near_boundary.append([r['id'],a['id']])
                a['points'].append({'transcript_id':r['id'],'gene':r['gene'],'source':source,'t':str(t-ta),'x':str(u),
                                    'global_x_um':r['x_um'],'global_y_um':r['y_um']})
                found.append(j)
            if len(found)>1:raise AssertionError('duplicate molecule in patch atlas')
            if found:matched.add(r['id'])
        if near_boundary:raise AssertionError(f'clearance not source-feasible: {near_boundary}')
        for a in accepted:
            pts=a['source_path_pixels'];sx,sy=a['sx'],a['sy'];a['domain_pixels']=[0,sx*(pts[-1][0]-pts[0][0])+sy*(pts[-1][1]-pts[0][1])]
            a['clearance_pixels']=str(F(1,1000000)/p)
            a['points'].sort(key=lambda q:int(q['transcript_id']))
        summary={'roi':d.name,'positive_labels':len(cells)-(0 in cells),'complete_connected_labels':len(eligible),'border_labels':len(border-{0}),
                 'all_shared_edges':all_shared,'eligible_shared_edges':sum(map(len,edges.values())),
                 'shared_pairs':len(edges),'chains':chain_count,'candidates':len(candidates),'selected_patches':len(accepted),
                 'selected_source_edges':sum(a['length_edges'] for a in accepted),'rejections':dict(reject),'width_trials':trials['widths'],
                 'exact_half_ribbon_pixel_checks':exact_checks,'retained_RNA':len(records),'unique_editable_RNA':len(matched),
                 'RNA_bearing_patches':sum(bool(a['points']) for a in accepted),
                 'widest_patch_RNA':max([len(a['points']) for a in accepted],default=0),
                 'cells_incident_to_patches':len({k for a in accepted for k in a['cell_labels'] if k>0}),'background_interface_patches':sum(0 in a['cell_labels'] for a in accepted),
                 'width_histogram':dict(Counter(a['halfwidth_pixels'] for a in accepted)), 'structural_constraints':sum(len(a['nuclear_constraints']) for a in accepted), 'patches_with_nuclei':sum(bool(a['nuclear_constraints']) for a in accepted)}
        atlas={'roi':d.name,'source_archive_sha256':meta['source_archive_sha256'],'pixel_size_um':str(p),
               'bbox_pixels':meta['bbox_pixels'],'complete_cells':sorted(eligible),'source_counts':data['rois'][d.name]['counts'],
               'records':records,'patches':accepted,'summary':summary,
               'source_hashes':data['rois'][d.name]['source_hashes']}
        atlas['atlas_id']=hashlib.sha256(json.dumps(atlas,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        (ROOT/f'data/atlas_{d.name}.json').write_text(json.dumps(atlas,indent=2)+'\n')
        all_rois.append(summary);print(summary,flush=True)
    (ROOT/('results/atlas_extraction_'+os.environ['GEOCERT_ROI']+'.json' if os.environ.get('GEOCERT_ROI') else 'results/atlas_extraction.json')).write_text(json.dumps({'protocol_sha256':hashlib.sha256((ROOT/'docs/ATLAS_PROTOCOL.json').read_bytes()).hexdigest(),
           'rois':all_rois,'seconds':time.time()-start},indent=2)+'\n')
if __name__=='__main__':main()
