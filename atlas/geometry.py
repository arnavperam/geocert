"""Rational planar geometry for independently checking proposed ribbons.

The pixel adapter may use GEOS to propose candidates. Acceptance is additionally
checked by exact convex clipping against reference pixels and earlier ribbons.
"""
from fractions import Fraction as F
from math import floor, ceil


def cross(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def area2(poly):
    return sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(poly,poly[1:]+poly[:1]))


def ccw(poly):
    p=[(F(x),F(y)) for x,y in poly]
    return p if area2(p)>=0 else list(reversed(p))


def intersection(subject, clip):
    """Convex Sutherland--Hodgman intersection, all comparisons rational."""
    out=ccw(subject)
    c=ccw(clip)
    for a,b in zip(c,c[1:]+c[:1]):
        inp=out;out=[]
        if not inp:break
        prev=inp[-1];dp=cross(a,b,prev)
        for point in inp:
            dq=cross(a,b,point)
            if (dq>=0)!=(dp>=0):
                ratio=dp/(dp-dq)
                out.append((prev[0]+ratio*(point[0]-prev[0]), prev[1]+ratio*(point[1]-prev[1])))
            if dq>=0:out.append(point)
            prev,dp=point,dq
    return out


def positive_overlap(a,b):
    ax=[x for x,y in a];ay=[y for x,y in a]
    bx=[x for x,y in b];by=[y for x,y in b]
    if max(ax)<=min(bx) or max(bx)<=min(ax) or max(ay)<=min(by) or max(by)<=min(ay):return False
    p=intersection(a,b)
    return len(p)>=3 and area2(p)!=0


def quad(edge, shift, lower, upper):
    a,b=edge
    return ccw([(a[0]+lower*shift[0],a[1]+lower*shift[1]),
                (b[0]+lower*shift[0],b[1]+lower*shift[1]),
                (b[0]+upper*shift[0],b[1]+upper*shift[1]),
                (a[0]+upper*shift[0],a[1]+upper*shift[1])])


def check_pixels(poly, mask, allowed, nucleus=None):
    """Every positive-area pixel intersection must have the requested label."""
    h,w=mask.shape
    x0=floor(min(p[0] for p in poly));x1=ceil(max(p[0] for p in poly))
    y0=floor(min(p[1] for p in poly));y1=ceil(max(p[1] for p in poly))
    for r in range(y0,y1):
        for c in range(x0,x1):
            valid=0<=r<h and 0<=c<w and int(mask[r,c]) in allowed
            if valid and (nucleus is None or not nucleus[r,c]):continue
            square=[(c,r),(c+1,r),(c+1,r+1),(c,r+1)]
            if positive_overlap(poly,square):return False
    return True


def evaluate_path(points, sx, sy, t):
    # t is guaranteed in the monotone path's domain.
    from bisect import bisect_right
    ts=[sx*x+sy*y for x,y in points]
    j=max(0,min(len(points)-2,bisect_right(ts,t)-1))
    a,b=points[j:j+2];ta,tb=ts[j:j+2]
    ua=sx*a[0]-sy*a[1];ub=sx*b[0]-sy*b[1]
    return F(ua)+(F(t)-ta)*F(ub-ua,tb-ta)


def ribbon_pieces(points, sx, sy, width):
    sh=(F(sx)*width/2,-F(sy)*width/2)
    return [quad((a,b),sh,-1,1) for a,b in zip(points,points[1:])]
