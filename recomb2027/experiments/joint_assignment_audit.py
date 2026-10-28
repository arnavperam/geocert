"""Exhaustive, post-review native pair audit; does not select favorable pairs.

Run from any directory. Every individually feasible non-source owner move is
included, then every unordered pair on distinct molecular groups. The original
84 model settings and eight-budget grid are fixed. No group multiplicity weights
are used. These are conditional local assignment questions, not biological truth.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import json,sys,hashlib,csv,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from mixed import make_problem
from data_api import variant

def main():
    protocol={
        'origin':'Post-review exhaustive audit, specified before pair-area outcomes; not a preregistered biological study.',
        'universe':'All 42 retained corridors, baseline and inherited maximal geometry; all pairs of distinct individually movable transcript groups.',
        'move':'One observed transcript pinned to any individually feasible owner different from its source; all other unpinned transcripts may move as collateral.',
        'no_weights':'One explicit transcript per move; no binary-encoded bundles.',
        'comparators':{'marginal':'max of the two exact singleton costs; a sound lower bound on joint cost, never a realizability guarantee.',
                       'additive':'sum of exact singleton costs plus explicit common feasibility; an inner budget approximation for these source-feasible models.'},
        'budgets':['0','1/20','1/10','1/4','1/2','1','2','4'],
        'selection':'No pair, gene, model or outcome filtering; nonmovable groups excluded by the stated move definition.',
        'scope':'Development data from one donor, overlapping local corridor families; no error-prediction, whole-cell or downstream biological claim.'}
    out=ROOT/'results/joint_assignment_audit';out.mkdir(exist_ok=True)
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    t0=time.perf_counter();rows=[];models=[];witnesses=[];moves_all=[]
    for ch in json.loads((ROOT/'data/charts.json').read_text()):
        p0,pins,n=make_problem(ch)
        for setting,p in [('baseline',p0),('maximal',variant(p0))]:
            en=ExpressionEngine(p,pins);moves=[]
            for g in range(n):
                for o in range(p.m+1):
                    if o==p.source[g]:continue
                    a=en.clipped({g:o})
                    if a['feasible']:
                        moves.append((g,o,F(a['area'])))
                        moves_all.append({'corridor':ch['id'],'setting':setting,'group':g,'source':p.source[g],'owner':o,'gene':ch['points'][g]['gene'],'area':str(a['area'])})
            pair_count=0
            for a,b in combinations(moves,2):
                g,o,ca=a;h,v,cb=b
                if g==h:continue
                pair_count+=1
                ans=en.clipped({g:o,h:v});feasible=ans['feasible'];joint=F(ans['area']) if feasible else None
                low=max(ca,cb);add=ca+cb
                if feasible:assert low<=joint<=add
                row={'corridor':ch['id'],'roi':ch['roi'],'setting':setting,
                     'group_1':g,'source_1':p.source[g],'owner_1':o,'gene_1':ch['points'][g]['gene'],
                     'group_2':h,'source_2':p.source[h],'owner_2':v,'gene_2':ch['points'][h]['gene'],
                     'singleton_1':str(ca),'singleton_2':str(cb),'marginal_lower':str(low),'additive_sum':str(add),
                     'joint_feasible':feasible,'joint_area':str(joint) if feasible else '',
                     'joint_above_marginal':feasible and joint>low,'shared_cost_saving':feasible and joint<add,
                     'same_gene':ch['points'][g]['gene']==ch['points'][h]['gene']}
                rows.append(row)
                witnesses.append({'corridor':ch['id'],'setting':setting,'pins':{g:o,h:v},'answer':ans})
            models.append({'corridor':ch['id'],'setting':setting,'moves':len(moves),'pairs':pair_count})
            print(ch['id'],setting,len(moves),pair_count,flush=True)
    tables=[];budgets=[]
    for setting in ['baseline','maximal']:
        rr=[r for r in rows if r['setting']==setting];ff=[r for r in rr if r['joint_feasible']]
        tables.append({'setting':setting,'models':42,'models_with_pairs':sum(m['pairs']>0 for m in models if m['setting']==setting),
                       'moves':sum(m['moves'] for m in models if m['setting']==setting),'pairs':len(rr),'jointly_feasible':len(ff),
                       'incompatible':len(rr)-len(ff),'joint_above_marginal':sum(r['joint_above_marginal'] for r in ff),
                       'joint_equals_marginal':sum(not r['joint_above_marginal'] for r in ff),
                       'shared_cost_saving':sum(r['shared_cost_saving'] for r in ff),
                       'same_gene_pairs':sum(r['same_gene'] for r in rr),
                       'models_with_incompatibility':len({r['corridor'] for r in rr if not r['joint_feasible']}),
                       'models_with_cost_saving':len({r['corridor'] for r in rr if r['shared_cost_saving']})})
        for bs in protocol['budgets']:
            d=F(bs);marg=[r for r in rr if F(r['marginal_lower'])<=d]
            exact=[r for r in rr if r['joint_feasible'] and F(r['joint_area'])<=d]
            inner=[r for r in rr if r['joint_feasible'] and F(r['additive_sum'])<=d]
            budgets.append({'setting':setting,'budget':bs,'pairs':len(rr),'individually_admitted':len(marg),'jointly_admitted':len(exact),
                            'marginal_false_admissions':len(marg)-len(exact),'additive_admitted':len(inner),'additive_misses':len(exact)-len(inner),
                            'marginal_incompatible':sum(not r['joint_feasible'] for r in marg)})
    native=list(csv.DictReader((ROOT/'results/expression_intervals.csv').open()))
    questions={tuple(r[k] for k in ['corridor','setting','cell_id','gene']):r for r in native}
    stats={'rows':len(native),'statuses':dict(Counter(r['status'] for r in native)),
           'questions':len(questions),'local_gene_multiplicity':dict(Counter(r['local_gene_transcripts'] for r in questions.values())),
           'row_interval_widths':dict(Counter(int(r['upper'])-int(r['lower']) for r in native))}
    result={'protocol_sha256':hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest(),'tables':tables,'budget_table':budgets,
            'native_count_verification':stats,'model_counts':models,'seconds':time.perf_counter()-t0,
            'interpretation':'Joint local assignments exhibit compatibility and shared-deformation interactions; this does not imply useful marginal gene-count improvements or calibrated biological uncertainty.'}
    with (out/'pairs.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'moves.json').write_text(json.dumps(moves_all,indent=2)+'\n')
    (out/'witnesses.json').write_text(json.dumps(witnesses,default=str,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='model_counts'},indent=2))
if __name__=='__main__':main()
