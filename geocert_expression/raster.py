"""Dimension-independent candidate sets from a labeled image in physical units.

Requires NumPy only when called. Points and spacing use ARRAY AXIS order, e.g.
(y,x) in 2D or (z,y,x) in 3D, after any external coordinate transformations.
This is a displacement-based outer relaxation, not a continuous-area certificate.
"""
from fractions import Fraction as F
from math import floor

def label_candidates(labels, points, genes, *, spacing, radius, origin=None,
                     localization_radius=0, source_labels=None, outside_label=-1):
    import numpy as np
    labels=np.asarray(labels)
    if labels.ndim not in (2,3) or not np.issubdtype(labels.dtype,np.integer):raise ValueError('2D or 3D integer label image required')
    d=labels.ndim;spacing=[F(str(x)) for x in spacing]
    origin=[F(0)]*d if origin is None else [F(str(x)) for x in origin]
    if len(spacing)!=d or len(origin)!=d or any(x<=0 for x in spacing):raise ValueError('one positive physical spacing per image axis required')
    r=F(str(radius));e=F(str(localization_radius))
    if min(r,e)<0 or len(points)!=len(genes):raise ValueError('invalid radius or point/gene lengths')
    all_labels=set(map(int,np.unique(labels)))|{outside_label};records=[];outside=0
    for i,(point,gene) in enumerate(zip(points,genes)):
        if len(point)!=d:raise ValueError('point dimensionality mismatch')
        idx=[(F(str(x))-o)/h for x,o,h in zip(point,origin,spacing)]
        central=tuple(floor(x) for x in idx)
        if any(j<0 or j>=n for j,n in zip(central,labels.shape)):raise ValueError('observed point lies outside label image')
        source=int(labels[central])
        if source_labels is not None and int(source_labels[i])!=source:raise ValueError('provided owner disagrees with reference label image')
        lo=[floor(x-(r+e)/h) for x,h in zip(idx,spacing)];hi=[floor(x+(r+e)/h) for x,h in zip(idx,spacing)]
        if any(a<0 or b>=n for a,b,n in zip(lo,hi,labels.shape)):
            candidates=all_labels;outside+=1
        else:
            candidates=set(map(int,np.unique(labels[tuple(slice(a,b+1) for a,b in zip(lo,hi))])))
        assert source in candidates
        records.append({'source':source,'gene':gene,'candidates':sorted(candidates)})
    return {'records':records,'outside_context_points':outside,
            'model':'L-infinity label-transport displacement plus localization radius in the supplied physical coordinate system',
            'axis_order':'array axes, not automatically spatial-world axes',
            'scope':'candidate relaxation; no connected-boundary or topology witness; unknown exterior treated conservatively'}
