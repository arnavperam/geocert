"""Recompute the complete same-feasibility ablation and its gene-count readouts.

Uses all 42 existing corridors, both settings, and the three original penalties.
No geometry, penalty, gene panel, or case is selected using a new result.
The exact and tie-favorable additive answers are deterministic returned optima;
readout changes outside the unique showcase are not asserted for every tie.
"""
from pathlib import Path
from fractions import Fraction as F
from math import lcm,ceil
from collections import Counter
import argparse,csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'core'))
from data_api import load_cases,reassignment_costs,encode
from pin_closure import PinClosure
from scc_solver import prepare,solve_compiled

def readout(ch,owners,inputs):
    cells=inputs['rois'][ch['roi']]['cells'];labels=ch['cell_labels']
    before={(c,g):v for c in labels for g,v in cells[str(c)]['counts'].items()}
    after=Counter(before);moves=[]
    for i,pt in enumerate(ch['points']):
        old=labels[pt['source']];new=labels[owners[i]];gene=pt['gene']
        if old!=new:
            after[(old,gene)]-=1;after[(new,gene)]+=1
            moves.append({'transcript_id':pt['transcript_id'],'gene':gene,'from':old,'to':new})
    assert all(v>=0 for v in after.values())
    genes={g for _,g in set(before)|set(after)}
    assert all(sum(before.get((c,g),0) for c in labels)==sum(after[(c,g)] for c in labels) for g in genes)
    changed=[{'cell':c,'gene':g,'before':before.get((c,g),0),'after':after[(c,g)],
              'touches_crop_edge':cells[str(c)]['touches_crop_edge']}
             for c,g in sorted(set(before)|set(after)) if before.get((c,g),0)!=after[(c,g)]]
    toggles=[r for r in changed if (r['before']>0)!=(r['after']>0)]
    return {'moves':moves,'changed_counts':changed,'detection_toggles':toggles,
            'complete_cell_detection_toggles':[r for r in toggles if not r['touches_crop_edge']]}

def run(output):
    inputs=json.loads((ROOT/'data/readout_inputs.json').read_text())
    expected={(r['id'],r['setting'],r['lambda']):r for r in json.loads((ROOT/'data/expected_additive.json').read_text())}
    rows=[];singleton_queries=0
    for index,ch,setting,p,pins,n,prior in load_cases():
        pc=PinClosure(p,pins);dist=[[F(0)]*(p.m+1) for _ in range(p.n)]
        for g in range(n):
            for o in range(p.m+1):
                if o!=p.source[g]:
                    a=pc.query({g:o});singleton_queries+=1
                    dist[g][o]=a['area'] if a['feasible'] else None
        compiled,plan,_=prepare(p,pins,'direct',True);costs=reassignment_costs(p,n)
        for lam in [F(1,4),F(1),F(4)]:
            ex=solve_compiled(compiled,'direct',costs,lam,plan)
            additive=[[(F(0) if dist[g][o] is None else lam*dist[g][o])+costs[g][o] for o in range(p.m+1)] for g in range(p.n)]
            primary=solve_compiled(compiled,'direct',additive,F(0),plan)
            L=lcm(*(z.denominator for row in additive for z in row))
            B=(p.hi-p.lo)*sum(hi-lo for lo,hi in p.bands);M=ceil(lam*B+n)+1
            lex=[[M*L*additive[g][o]+costs[g][o] for o in range(p.m+1)] for g in range(p.n)]
            ad=solve_compiled(compiled,'direct',lex,lam,plan)
            ea=sum(additive[g][o] for g,o in enumerate(ad['owners']))
            assert ea==primary['objective']
            assert all(dist[g][o] is not None for ans in [ex,ad] for g,o in enumerate(ans['owners']))
            def E(ans):return lam*ans['area']+sum(costs[g][o] for g,o in enumerate(ans['owners']))
            regret=E(ad)-E(ex);assert regret>=0
            old=expected[(ch['id'],setting,str(lam))]
            assert (E(ex),E(ad),ea,regret)==tuple(F(old[k]) for k in ['joint_objective','additive_joint_objective','additive_primary_objective','regret'])
            rr={method:readout(ch,ans['owners'],inputs) for method,ans in [('joint',ex),('additive',ad)]}
            diff={(r['cell'],r['gene']) for r in rr['joint']['changed_counts']+rr['additive']['changed_counts']}
            def count(method,key):
                base=inputs['rois'][ch['roi']]['cells'][str(key[0])]['counts'].get(key[1],0)
                return next((r['after'] for r in rr[method]['changed_counts'] if (r['cell'],r['gene'])==key),base)
            disagreements=[{'cell':c,'gene':g,'joint':count('joint',(c,g)),'additive':count('additive',(c,g))}
                           for c,g in sorted(diff) if (count('joint',(c,g))>0)!=(count('additive',(c,g))>0)]
            rows.append({'id':ch['id'],'roi':ch['roi'],'setting':setting,'lambda':str(lam),
                         'joint_objective':str(E(ex)),'additive_joint_objective':str(E(ad)),'regret':str(regret),
                         'joint_area':str(ex['area']),'additive_area':str(ad['area']),
                         'joint_owners':ex['owners'],'additive_owners':ad['owners'],
                         'joint':rr['joint'],'additive':rr['additive'],'detection_disagreements':disagreements})
        print(index,setting,'done',flush=True)
    summary={}
    for setting in ['baseline','maximal']:
        ss=[r for r in rows if r['setting']==setting];strict=[r for r in ss if F(r['regret'])>0]
        summary[setting]={'questions':len(ss),'strict_regret':len(strict),
                          'joint_answers_changing_detection':sum(bool(r['joint']['detection_toggles']) for r in ss),
                          'joint_answers_changing_complete_cell_detection':sum(bool(r['joint']['complete_cell_detection_toggles']) for r in ss),
                          'strict_regret_with_detection_disagreement':sum(bool(r['detection_disagreements']) for r in strict),
                          'all_answers_with_detection_disagreement':sum(bool(r['detection_disagreements']) for r in ss)}
    result={'protocol':'Post-hoc readout analysis of all unchanged 252 queries; every model and gene retained. Detection means observed Q20 count > 0. Alternatives are evaluated separately; no stitching or cell-type claim.',
            'source_rna_records':inputs['source_rna_records'],'source_chart_records_checked':inputs['chart_records_checked'],
            'singleton_queries':singleton_queries,'fresh_optimization_queries':3*len(rows),'archived_objective_comparisons':4*len(rows),
            'summary':summary,'rows':rows}
    output.mkdir(parents=True,exist_ok=True);(output/'readout_audit.json').write_text(json.dumps(result,indent=2,default=encode)+'\n')
    with (output/'readout_census.csv').open('w',newline='') as f:
        keys=['id','roi','setting','lambda','joint_objective','additive_joint_objective','regret','joint_area','additive_area','joint_detection_toggles','additive_detection_toggles','detection_disagreements']
        w=csv.DictWriter(f,fieldnames=keys,lineterminator="\n");w.writeheader()
        for r in rows:w.writerow({**{k:r[k] for k in keys[:9]},'joint_detection_toggles':len(r['joint']['detection_toggles']),'additive_detection_toggles':len(r['additive']['detection_toggles']),'detection_disagreements':len(r['detection_disagreements'])})
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,default=ROOT/'run/readout');a=ap.parse_args();run(a.output_dir)
