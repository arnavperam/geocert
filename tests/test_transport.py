"""Rational inverse, monotonicity, fixed-trace and physical-area tests."""
import sys,json,random
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from atlas.transport import transverse_forward,transverse_inverse,chart_to_pixel,pixel_to_chart
from atlas.geometry import area2
rng=random.Random(401)
inverse=monotone=rails=coordinates=areas=0
for b in [F(1,2),F(1),F(2),F(4)]:
 for k in range(-7,8):
  h=b*F(k,8);values=[b*F(i,20) for i in range(-20,21)]
  mapped=[transverse_forward(v,h,b) for v in values]
  assert all(u<v for u,v in zip(mapped,mapped[1:]));monotone+=1
  for v,w in zip(values,mapped):assert transverse_inverse(w,h,b)==v;inverse+=1
  assert mapped[0]==-b and mapped[-1]==b and transverse_forward(0,h,b)==h;rails+=3
  assert all(transverse_forward(v,0,b)==v for v in values)
for path in sorted((ROOT/'data').glob('atlas_roi_*.json')):
 a=json.loads(path.read_text());p=F(a['pixel_size_um'])
 for patch in a['patches']:
  b=F(patch['halfwidth_pixels']);T=F(patch['domain_pixels'][1])
  for _ in range(8):
   t=T*F(rng.randrange(101),100);h=b*F(rng.randrange(-100,101),100)
   q=chart_to_pixel(patch,t,h);assert pixel_to_chart(patch,*q)==(t,h);coordinates+=1
  # Every source raster edge corresponds to a unit longitudinal interval. A
  # one-sided trapezoid in h has physical area p^2/2 times its chart integral.
  for j in range(len(patch['source_path_pixels'])-1):
   u=b*F(rng.randrange(1,8),8);v=b*F(rng.randrange(1,8),8)
   poly=[chart_to_pixel(patch,j,0),chart_to_pixel(patch,j+1,0),
         chart_to_pixel(patch,j+1,v),chart_to_pixel(patch,j,u)]
   assert abs(area2(poly))*p*p/2==(u+v)/2*F(patch['area_scale_um2']);areas+=1
for args in [(0,1,1),(0,-1,1),(2,0,1),(0,0,0)]:
 try:transverse_forward(*args)
 except ValueError:pass
 else:raise AssertionError('invalid homeomorphism accepted')
out=dict(inverse_points=inverse,strict_monotone_sections=monotone,rail_and_interface_equalities=rails,
         native_coordinate_roundtrips=coordinates,independent_physical_polygon_areas=areas,
         rejected_degenerate_maps=4,failures=0)
(ROOT/'results/transport_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
