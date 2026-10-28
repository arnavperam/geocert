"""Optional arbitrary-precision C++ flow, always checked by Python integers."""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess,time
from fractions import Fraction as F
from pinned_compiler import graph_for
from check_flow import check_certificate
from layered import realize
from mixed import verify_mixed

def flow(graph,binary,timeout=60):
    with TemporaryDirectory(prefix='geocert-flow-') as tmp:
        a=Path(tmp)/'input';b=Path(tmp)/'output'
        with a.open('w') as f:
            f.write(f"{graph['n']} {len(graph['edges'])} {graph['s']} {graph['t']}\n")
            for u,v,c in graph['edges']:f.write(f'{u} {v} {c}\n')
        subprocess.run([str(Path(binary).resolve()),str(a),str(b)],check=True,timeout=timeout,capture_output=True,text=True)
        lines=b.read_text().splitlines();value=int(lines[0]);side=list(map(int,lines[1].split()))
        if side[0]!=len(side)-1 or len(lines)!=2+len(graph['edges']):raise ValueError('malformed flow output')
        cert={'n':graph['n'],'s':graph['s'],'t':graph['t'],'value':value,'source_side':side[1:],
              'scale':graph.get('scale',1),'constant':str(graph.get('constant',0)),
              'arcs':[(u,v,c,int(f)) for (u,v,c),f in zip(graph['edges'],lines[2:])]}
        check_certificate(cert)
        return cert

def solve(compiled,costs,multiplier,binary,*,reconstruct=True,timeout=60):
    t0=time.perf_counter();graph=graph_for(compiled,costs,multiplier);proof=flow(graph,binary,timeout)
    out={'feasible':proof['value']<graph['big'],'certificate':proof,'backend':'cpp-arbitrary-precision-dinic'}
    if out['feasible']:
        side=set(proof['source_side']);p=compiled.problem;owners=p.decode([int(i in side) for i in range(p.n*p.m)])
        area=compiled.area(owners);objective=F(multiplier)*area+sum(F(r[o]) for r,o in zip(costs,owners))
        if objective!=F(proof['value'],graph['scale'])+graph['constant']:raise AssertionError('C++ cut objective mismatch')
        out.update(owners=owners,area=area,objective=objective)
        if reconstruct:
            w=realize(p,owners);assert w['feasible'] and w['area']==area;verify_mixed(p,owners,w['boundaries'],area);out['boundaries']=w['boundaries']
    out['seconds']=time.perf_counter()-t0;return out
