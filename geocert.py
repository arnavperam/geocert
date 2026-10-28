#!/usr/bin/env python3
"""Reviewer front door. Exact reproduction uses Python's standard library only."""
from pathlib import Path
import sys,json,argparse,time,hashlib,subprocess,os
if not __debug__:raise SystemExit('Run without -O: exact scientific assertions must remain enabled.')
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'core'))
from fractions import Fraction as F
from data_api import load_cases,reassignment_costs,encode
from scc_solver import prepare,solve_compiled

def emit(obj,out=None):
    text=json.dumps(obj,indent=2,default=encode)
    if out:
        q=Path(out);q.parent.mkdir(parents=True,exist_ok=True);q.write_text(text+'\n')
    print(text)

def showcase(out=None):
    case=next(z for z in load_cases() if z[1]['id']=='4037abe9fea38e' and z[2]=='baseline')
    _,ch,tag,p,pins,n,expected=case;results={}
    for rep in ['direct','height']:
        c,plan,stages=prepare(p,pins,rep,True)
        a=solve_compiled(c,rep,reassignment_costs(p,n),F(1),plan,certificate=True)
        if a['objective']!=F(expected['objectives']['1']):raise AssertionError('showcase mismatch')
        results[rep]={'objective':a['objective'],'deformation':a['area'],'changed_rna':sum(a['owners'][g]!=p.source[g] for g in range(n)),
                     'nodes':a['nodes'],'arcs':a['arcs'],'owners':a['owners'],'boundaries':a['boundaries'],
                     'certificate':a['certificate'],'old_to_contracted':a['mapping'],'seconds':stages['prepare_seconds']+a['seconds']}
    audit=json.loads((ROOT/'data/showcase_scalar_audit.json').read_text());q=audit['queries'][1]
    for kind,col,winner in [('joint','D',q['exact_owners']),('additive','D_add',q['additive_owners'])]:
        w=next(r for r in audit['rows'] if r['owners']==winner)
        for lam in [F(67,100),F(32,25)]:
            assert min(lam*(F(r[col])-F(w[col]))-(r['R']-w['R']) for r in audit['rows'] if r['owners']!=winner)>0
    result={'id':ch['id'],'scope':'native-derived fixed-geometry sensitivity, not corrected anatomy','exact_costs_agree':True,'penalty_interval':['67/100','32/25'],'results':results}
    if out:Path(out).parent.mkdir(parents=True,exist_ok=True);Path(out).write_text(json.dumps(result,indent=2,default=encode)+'\n')
    emit({'id':ch['id'],'exact_costs_agree':True,'verified_penalty_interval':['67/100','32/25'],'results':{k:{a:b for a,b in r.items() if a not in ['owners','boundaries','certificate','old_to_contracted']} for k,r in results.items()},'full_certificate_file':str(out) if out else None})

def compare(rep,use_scc,out):
    rows=[];outcomes=[]
    for i,ch,tag,p,pins,n,expected in load_cases():
        t=time.perf_counter();c,plan,stages=prepare(p,pins,rep,use_scc);rr=[]
        for lam in [F(1,4),F(1),F(4)]:
            a=solve_compiled(c,rep,reassignment_costs(p,n),lam,plan)
            if not a['feasible'] or a['objective']!=F(expected['objectives'][str(lam)]):raise AssertionError('native objective differs')
            outcomes.append([ch['id'],tag,str(lam),str(a['objective'])]);rr.append({k:a[k] for k in ['objective','area','nodes','arcs','seconds']})
        rows.append({'id':ch['id'],'setting':tag,**stages,'queries':rr,'seconds':time.perf_counter()-t})
    result={'representation':rep,'scc':use_scc,'models':len(rows),'queries':len(outcomes),'objective_hash':hashlib.sha256(json.dumps(outcomes).encode()).hexdigest(),'total_seconds':sum(r['seconds'] for r in rows),'rows':rows}
    if out:Path(out).parent.mkdir(parents=True,exist_ok=True);Path(out).write_text(json.dumps(result,indent=2,default=encode)+'\n')
    emit({k:v for k,v in result.items() if k!='rows'})

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    sh=sub.add_parser('showcase');sh.add_argument('--output',default='run/showcase.json')
    c=sub.add_parser('compare');c.add_argument('--representation',choices=['direct','height'],required=True);c.add_argument('--scc',action='store_true');c.add_argument('--output',required=True)
    b=sub.add_parser('benchmark');b.add_argument('--repetitions',type=int,default=1);b.add_argument('--scc',action='store_true');b.add_argument('--output-dir',default='run/comparison')
    sub.add_parser('verify-manifest')
    a=p.parse_args()
    if a.command=='showcase':showcase(a.output)
    elif a.command=='compare':compare(a.representation,a.scc,a.output)
    elif a.command=='benchmark':
        if not 1<=a.repetitions<=20:raise SystemExit('Choose 1 to 20 repetitions.')
        folder=Path(a.output_dir);folder.mkdir(exist_ok=True,parents=True)
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONHASHSEED='0')
        for rep in range(a.repetitions):
            for method in (['direct','height'] if rep%2==0 else ['height','direct']):
                cmd=[sys.executable,str(ROOT/'geocert.py'),'compare','--representation',method,'--output',str(folder/f'{method}_{rep}.json')]
                if a.scc:cmd.append('--scc')
                subprocess.run(cmd,check=True,env=env)
        hashes={json.loads(p.read_text())['objective_hash'] for p in folder.glob('*.json')}
        if len(hashes)!=1:raise AssertionError('benchmark result mismatch')
        emit({'completed_processes':2*a.repetitions,'identical_exact_objectives':True,'output_directory':str(folder)})
    else:
        path=ROOT/'SHA256SUMS.txt';count=0
        for line in path.read_text().splitlines():
            h,name=line.split('  ',1);q=ROOT/name
            if hashlib.sha256(q.read_bytes()).hexdigest()!=h:raise SystemExit('checksum mismatch: '+name)
            count+=1
        emit({'verified_files':count})
if __name__=='__main__':main()
