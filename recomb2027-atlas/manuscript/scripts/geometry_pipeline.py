"""Draw the illustrative geometry-to-spectrum diagram used in the manuscript.

This is an analytic schematic, not a tissue-data plot. Coordinates obey the
staircase shear chart in the manuscript; the lower tents and canonical clamp in the
middle panel use s=0.65, point requirements (2,0.6), (4,0.8), fixed endpoints,
and an upper fixed-owner interval. The final panel illustrates DP states.
Run from the project directory: python scripts/geometry_pipeline.py
"""
from pathlib import Path
import os
import tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.labelsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'lines.solid_capstyle': 'round',
    'savefig.facecolor': 'white',
})
NAVY, BLUE, ORANGE, GREEN = '#203449', '#2375a7', '#bf6d1a', '#60784c'
GRAY, LIGHT = '#9da8af', '#e9f1f5'
fig = plt.figure(figsize=(6.5, 2.3))
gs = fig.add_gridspec(1, 3, width_ratios=(1.65, 2.25, 1.85),
                      left=.025, right=.985, bottom=.19, top=.82, wspace=.36)
a,b,c = [fig.add_subplot(gs[0,i]) for i in range(3)]
for ax, title in zip((a,b,c), ('(a) Raster ribbon', '(b) Exact envelopes', '(c) Count states')):
    fig.text(ax.get_position().x0, .87, title, fontsize=10, fontweight='bold', color=NAVY)

# A: an exact coordinate-monotone pixel staircase, s_x=s_y=1, p=1.
t = np.linspace(0,6,1801)
r = np.interp(t, np.arange(7), [0,1,0,1,0,1,0])
def shear(t0,h):
    r0=np.interp(t0,np.arange(7),[0,1,0,1,0,1,0])
    return .5*(t0+r0+h), .5*(t0-r0-h)
B=1.4
xlo,ylo=shear(t,-B); xhi,yhi=shear(t,B)
a.fill(np.r_[xlo,xhi[::-1]],np.r_[ylo,yhi[::-1]],color=LIGHT,zorder=0)
a.plot(xlo,ylo,color=BLUE,lw=.8);a.plot(xhi,yhi,color=BLUE,lw=.8)
a.plot(*shear(t,0),color=NAVY,lw=1.8)
# Nuclear interval under the same chart; its owner stays above the boundary.
nt=np.r_[np.linspace(2.75,3.25,40), np.linspace(3.25,2.75,40)]
nh=np.r_[np.full(40,.65),np.full(40,1.15)]
nx,ny=shear(nt,nh)
a.fill(nx,ny,color=GREEN,alpha=.62,ec=GREEN,lw=.6)
for ti,hi in ((2,.55),(4,.75)):
    a.plot(*shear(ti,hi),'o',color=ORANGE,ms=3.7,mec='white',mew=.5,zorder=5)
a.annotate('$v=-B$', xy=shear(5.5,-B),xytext=(-.6,3.65),fontsize=10,
           arrowprops={'arrowstyle':'-','color':BLUE,'lw':.6},color=BLUE)
a.text(1.10,-1.10,r'$v=+B$',color=BLUE,fontsize=10)
a.annotate('nucleus',xy=shear(3,1),xytext=(2.70,.05),fontsize=10,color=GREEN,
           arrowprops={'arrowstyle':'-','color':GREEN,'lw':.6})
a.annotate('RNA',xy=shear(4,.75),xytext=(3.35,2.15),fontsize=10,color=ORANGE,
           arrowprops={'arrowstyle':'-','color':ORANGE,'lw':.6})
a.set_xlim(-.9,3.95); a.set_ylim(-1.3,3.95);a.set_aspect('equal');a.axis('off')

# B: Lipschitz lower/upper envelopes and their nearest-to-source realization.
s=.65
lower1=.6-s*np.abs(t-2);lower2=.8-s*np.abs(t-4)
L=np.maximum.reduce([np.full_like(t,-.9),-s*t,-s*(6-t),lower1,lower2])
U=np.minimum.reduce([np.full_like(t,.9),s*t,s*(6-t),
                     .6+s*np.maximum.reduce([2.75-t,np.zeros_like(t),t-3.25])])
f=np.minimum(U,np.maximum(L,0))
b.fill_between(t,0,f,color=ORANGE,alpha=.15,zorder=0)
b.plot(t,lower1,color=BLUE,lw=.75,ls=':');b.plot(t,lower2,color=BLUE,lw=.75,ls=':')
b.plot(t,L,color=BLUE,lw=1.2,label='$L$');b.plot(t,U,color=GREEN,lw=1.1,label='$U$')
b.axhline(0,color=GRAY,lw=.75)
b.plot(t,f,color=ORANGE,lw=2.1,zorder=4)
b.scatter([2,4],[.55,.75],s=14,c=ORANGE,edgecolors='white',linewidth=.45,zorder=6)
b.fill_between([2.75,3.25],[.65,.65],[.88,.88],color=GREEN,alpha=.35,zorder=0)
# Canonical events include tent zeros, intersection, point sites and endpoints.
events=[0,2-.6/s,2,3-(.2/(2*s)),4,4+.8/s,6]
for e in events:
    b.axvline(e,color=GRAY,lw=.45,ls=(0,(2,3)),zorder=-1)
    b.plot(e,np.interp(e,t,f),'o',ms=2.9,color=ORANGE,zorder=5)
b.annotate('nucleus', xy=(3,.78), xytext=(2.4,1.23), fontsize=10,
           ha='center', color=GREEN,
           arrowprops={'arrowstyle':'-','color':GREEN,'lw':.6})
b.text(4.75,1.03,'$U$',color=GREEN,fontsize=10)
b.annotate('$L$', xy=(.55,-.35), xytext=(.4,-.72), fontsize=10,
           color=BLUE, arrowprops={'arrowstyle':'-','color':BLUE,'lw':.6})
b.text(3.65,.24,'$f^*$',color=ORANGE,fontsize=10)
b.annotate('$f^0=0$', xy=(4.6,0), xytext=(4.45,-.72), fontsize=10,
           color='#626c75', arrowprops={'arrowstyle':'-','color':GRAY,'lw':.6})
b.set(xlim=(-.05,6.1),ylim=(-.88,1.45),xticks=[0,2,4,6],yticks=[])
b.text(1.035, -.015, '$t$', transform=b.transAxes, ha='left', va='center')
b.set_ylabel('$h$',rotation=0,labelpad=5)
b.spines[['top','right','left']].set_visible(False);b.spines['bottom'].set_color(GRAY)
b.tick_params(axis='x',length=2,color=GRAY,pad=1)

# C: a small abstract state diagram. Vertical nodes are retained critical heights;
# count q is an extra discrete state at each node, not a second geometric axis.
heights=[[0,.3,.6],[0,.3,.6,.9],[.3,.6,.9]]
for j in (0,1):
    for h in heights[j]:
        for h2 in heights[j+1]:
            if abs(h2-h)<=.3000001:
                c.plot([j,j+1],[h,h2],color='#c8d0d5',lw=.6,zorder=0)
for j,hs in enumerate(heights):
    c.scatter([j]*len(hs),hs,s=16,c=LIGHT,edgecolors=BLUE,linewidths=.8,zorder=2)
path=[0,.3,.6]
c.plot([0,1,2],path,color=ORANGE,lw=1.7,zorder=3)
c.scatter([0,1,2],path,s=18,c=ORANGE,edgecolors='white',linewidths=.45,zorder=4)
c.text(.06,1.31,r'Area $F_j(h,q)$',fontsize=10,color=NAVY)
c.text(.05,1.00,'critical heights',fontsize=10,color=BLUE)
c.text(1,-.54,'add RNA count',ha='center',fontsize=10,color=ORANGE)
c.annotate('',xy=(1,.3),xytext=(1.65,-.41),
           arrowprops={'arrowstyle':'-','color':ORANGE,'lw':.65})
c.text(0,-.26,'$t_{j-1}$',ha='center',fontsize=10)
c.text(1,-.26,'$t_j$',ha='center',fontsize=10)
c.text(2,-.26,'$t_{j+1}$',ha='center',fontsize=10)
c.set(xlim=(-.3,2.4),ylim=(-.64,1.57));c.axis('off')
# Shared understated footer establishes what the schematic teaches.
fig.text(.5,.035,r'RNA constraints $\longrightarrow$ exact envelopes $\longrightarrow$ minimum area by count',
         fontsize=10,ha='center',color=NAVY)
for extension in ('pdf', 'png'):
    with tempfile.NamedTemporaryFile(suffix='.'+extension, delete=False) as temporary:
        temporary_path = Path(temporary.name)
    try:
        fig.savefig(temporary_path, dpi=300)
        os.replace(temporary_path, OUT / ('geometry_pipeline.'+extension))
    finally:
        temporary_path.unlink(missing_ok=True)
plt.close(fig)

