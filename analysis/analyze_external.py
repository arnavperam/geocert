"""Apply frozen joint/additive and count queries to the fixed external screen.

The tissue subset is vendor format-test data. This is a technical transfer check,
not a biological validation study. Every selected chart and resource limit is retained.
"""
from collections import Counter
from contextlib import contextmanager
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import argparse, json, signal, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'core'))
from data_api import variant, reassignment_costs, encode
from mixed import make_problem
from pin_closure import PinClosure
from scc_solver import prepare, solve_compiled

class Limit(Exception):
    pass

@contextmanager
def limited(seconds=15):
    def timeout(signum, frame):
        raise Limit('15-second per-stage/optimization cap')
    previous = signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)

def run(output):
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    charts = []
    for h, c in [(24,3),(12,3),(24,2),(12,2)]:
        charts.extend(json.loads((ROOT/f'results/external/charts_h{h}_c{c}.json').read_text()))
    (ROOT/'results/external/charts.json').write_text(json.dumps(charts,indent=2)+'\n')
    count_file=ROOT/'analysis/external_source_counts.json'
    if not count_file.exists():
        import pandas as pd
        counts={'counts':{},'complete_cell_labels':{},'source':{}}
        for roi in ['breast_01','breast_02']:
            folder=ROOT/'analysis/external_inputs'/roi
            meta=json.loads((folder/'metadata.json').read_text())
            df=pd.read_csv(folder/'transcripts.csv.gz',float_precision='round_trip')
            totals=Counter(zip(df.mask_label,df.feature_name))
            counts['complete_cell_labels'][roi]=meta['complete_source_labels']
            counts['counts'][roi]={}
            for chart in charts:
                if chart['roi']!=roi:continue
                for label in chart['cell_labels']:
                    if label not in meta['complete_source_labels']:continue
                    target=counts['counts'][roi].setdefault(str(label),{})
                    for gene in {z['gene'] for z in chart['points']}:
                        target[gene]=totals[label,gene]
            counts['source'][roi]={'RNA':len(df),'assignment_disagreements':meta['source_assignment_disagreements']}
        count_file.write_text(json.dumps(counts,separators=(',',':'))+'\n')
    counts = json.loads(count_file.read_text())
    rows, readouts, failures = [], [], []
    start = time.perf_counter()
    for ch in charts:
        p, pins, n = make_problem(ch)
        for setting, q in [('baseline',p),('maximal',variant(p))]:
            try:
                with limited():
                    compiled, plan, _ = prepare(q,pins,'direct',True)
                    closure = PinClosure(q,pins)
                single = [[F(0)]*(q.m+1) for _ in range(q.n)]
                for g in range(n):
                    for owner in range(q.m+1):
                        if owner == q.source[g]:continue
                        with limited():answer=closure.query({g:owner})
                        if answer['feasible']:single[g][owner]=answer['area']
            except Limit as e:
                failures.append({'chart':ch['id'],'setting':setting,'phase':'preparation','reason':str(e)})
                continue
            rewards=reassignment_costs(q,n)
            B=(q.hi-q.lo)*sum(hi-lo for lo,hi in q.bands)
            for penalty in [F(1,4),F(1),F(4)]:
                try:
                    with limited():j=solve_compiled(compiled,'direct',rewards,penalty,plan)
                    ac=[[penalty*d+r for d,r in zip(ds,rs)] for ds,rs in zip(single,rewards)]
                    with limited():a=solve_compiled(compiled,'direct',ac,F(0),plan)
                    scale=lcm(*(v.denominator for row in ac for v in row));M=int(penalty*B+n)+2
                    tc=[[M*scale*d+r for d,r in zip(ds,rs)] for ds,rs in zip(ac,rewards)]
                    with limited():tie=solve_compiled(compiled,'direct',tc,penalty,plan)
                    assert j['feasible'] and a['feasible'] and tie['feasible']
                    assert sum(ac[g][v] for g,v in enumerate(tie['owners']))==a['objective']
                    gain=sum(tie['owners'][g]!=q.source[g] for g in range(n))
                    regret=penalty*tie['area']-gain-j['objective'];assert regret>=0
                    rows.append({'chart':ch['id'],'setting':setting,'penalty':str(penalty),
                        'joint_objective':str(j['objective']),'joint_area':str(j['area']),
                        'joint_owners':j['owners'][:n],
                        'joint_changed':sum(j['owners'][g]!=q.source[g] for g in range(n)),
                        'additive_primary':str(a['objective']), 'additive_area':str(tie['area']),
                        'additive_changed':gain,'additive_owners':tie['owners'][:n],
                        'regret':str(regret)})
                except Limit as e:failures.append({'chart':ch['id'],'setting':setting,'penalty':str(penalty),'phase':'optional','reason':str(e)})
            local=Counter((z['source'],z['gene']) for z in ch['points'])
            for k,label in enumerate(ch['cell_labels']):
                if label not in counts['complete_cell_labels'][ch['roi']]:continue
                for gene in sorted({z['gene'] for z in ch['points']}):
                    total=counts['counts'][ch['roi']][str(label)].get(gene,0)
                    fixed=total-local[k,gene];assert fixed>=0
                    targets=[g for g,z in enumerate(ch['points']) if z['gene']==gene]
                    endpoints=[]
                    try:
                        for sign in [1,-1]:
                            costs=[[0]*(q.m+1) for _ in range(q.n)]
                            for g in targets:costs[g][k]=sign*(int(B)+2)
                            with limited():ans=solve_compiled(compiled,'direct',costs,F(1),plan)
                            assert ans['feasible']
                            v=fixed+sum(ans['owners'][g]==k for g in targets)
                            endpoints.append({'count':v,'area':str(ans['area']),'owners':ans['owners'][:n]})
                        assert endpoints[0]['count']<=total<=endpoints[1]['count']
                        readouts.append({'chart':ch['id'],'setting':setting,'roi':ch['roi'],
                            'cell_label':label,'cell_index':k,'gene':gene,'source_count':total,
                            'fixed_count':fixed,'minimum':endpoints[0],'maximum':endpoints[1]})
                    except Limit as e:failures.append({'chart':ch['id'],'setting':setting,'cell':label,'gene':gene,'phase':'count','reason':str(e)})
            print(ch['id'],setting,len(rows),len(readouts),len(failures),round(time.perf_counter()-start,2),flush=True)
            output.parent.mkdir(parents=True,exist_ok=True)
            output.write_text(json.dumps({'complete':False,'optional':rows,'readouts':readouts,'failures':failures},indent=2)+'\n')
    summary={}
    for setting in ['baseline','maximal']:
        a=[r for r in rows if r['setting']==setting];b=[r for r in readouts if r['setting']==setting]
        summary[setting]={'optional_questions':len(a),'positive_regret':sum(F(r['regret'])>0 for r in a),
           'readouts':len(b),'variable_readouts':sum(r['minimum']['count']<r['maximum']['count'] for r in b),
           'zero_positive_ambiguity':sum(r['minimum']['count']==0<r['maximum']['count'] for r in b)}
    output.write_text(json.dumps({'complete':True,'charts':len(charts),'summary':summary,
         'optional':rows,'readouts':readouts,'failures':failures,'elapsed_seconds':time.perf_counter()-start},indent=2)+'\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,default=ROOT/'results/external/analysis.json');run(a.parse_args().output)
