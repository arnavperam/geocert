"""Exact chart homeomorphisms used by the global gluing theorem.

These maps are constructive witnesses of topology preservation. They are not
estimated tissue motion or a membrane probability model.
"""
from fractions import Fraction as F
from .geometry import evaluate_path


def transverse_forward(v, boundary, rail):
    v,b,rail=map(F,(v,boundary,rail))
    if rail<=0 or not -rail<b<rail or not -rail<=v<=rail:
        raise ValueError('strict interior boundary and in-ribbon coordinate required')
    return -rail+(v+rail)*(b+rail)/rail if v<=0 else b+v*(rail-b)/rail


def transverse_inverse(w,boundary,rail):
    w,b,rail=map(F,(w,boundary,rail))
    if rail<=0 or not -rail<b<rail or not -rail<=w<=rail:
        raise ValueError('strict interior boundary and in-ribbon coordinate required')
    return rail*(w+rail)/(b+rail)-rail if w<=b else rail*(w-b)/(rail-b)


def chart_to_pixel(patch,t,h):
    t,h=F(t),F(h);sx,sy=patch['sx'],patch['sy'];path=patch['source_path_pixels']
    origin=sx*path[0][0]+sy*path[0][1];tt=origin+t
    u=evaluate_path(path,sx,sy,tt)+h
    return (F(sx)*(tt+u)/2,F(sy)*(tt-u)/2)


def pixel_to_chart(patch,x,y):
    x,y=F(x),F(y);sx,sy=patch['sx'],patch['sy'];path=patch['source_path_pixels']
    tt=sx*x+sy*y;origin=sx*path[0][0]+sy*path[0][1]
    return (tt-origin,sx*x-sy*y-evaluate_path(path,sx,sy,tt))
