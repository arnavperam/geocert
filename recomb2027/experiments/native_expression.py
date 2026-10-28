"""Predeclared native cell-by-gene expression-interval census. Resumable by model."""
from pathlib import Path
import sys,json,csv,hashlib,collections,time,argparse
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.budget import ExpressionEngine
from mixed import make_problem
from data_api import variant

def enc(v):
    if isinstance(v,F):return str(v)
    raise TypeError(type(v))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-new-models',type=int,default=84);args=ap.parse_args()
    protocol=json.loads((ROOT/'docs/EXPERIMENT_PROTOCOL.json').read_text())
    budgets=list(map(F,protocol['budgets_square_micrometres']));mus=list(map(F,protocol['lagrange_multipliers']))
    charts=json.loads((ROOT/'data/charts.json').read_text());obs=json.loads((ROOT/'data/observed_crop_counts.json').read_text())
    rows=[];models=[];new=0;t0=time.perf_counter();directory=ROOT/'results/native';directory.mkdir(exist_ok=True)
    for ci,ch in enumerate(charts):
        p0,pins,n=make_problem(ch)
        for setting,p in [('baseline',p0),('maximal',variant(p0))]:
            path=directory/f'{ch["id"]}_{setting}.json'
            if path.exists():
                rec=json.loads(path.read_text());rows.extend(rec['rows']);models.append(rec['summary']);continue
            if new>=args.max_new_models:
                print(json.dumps({'batch_completed':new,'models_present':len(list(directory.glob('*.json'))),'complete_census':False}));return
            t1=time.perf_counter();engine=ExpressionEngine(p,pins);records=[];metadata=[];cuts=patterns=0
            genes=sorted({z['gene'] for z in ch['points']});counts=obs['rois'][ch['roi']]['counts']
            for gene in genes:
                target=[i for i,z in enumerate(ch['points']) if z['gene']==gene]
                for k,cell in enumerate(ch['cell_labels']):
                    local=sum(ch['points'][i]['source']==k for i in target)
                    total=counts.get(str(cell),{}).get(gene,0);outside=total-local
                    if outside<0:raise AssertionError('negative exterior count')
                    scores=[[int(g in target and owner==k) for owner in range(p.m+1)] for g in range(p.n)]
                    sp=engine.spectrum(scores,outside)
                    if not sp.complete:raise RuntimeError('native target limit')
                    patterns+=sp.examined
                    dual=engine.dual_queries(scores,mus,certificate=(len(metadata)==0));cuts+=len(dual)
                    sd=[{key:v for key,v in a.items() if key in ['multiplier','sign','feasible','score','area','objective','owners','certificate']} for a in dual]
                    metadata.append({'gene':gene,'cell_index':k,'cell_id':cell,'local_indices':target,'fixed_outside_count':outside,
                                     'source_count':total,'spectrum':sp.serializable(),'dual_queries':sd})
                    for budget in budgets:
                        ex=sp.interval(budget);bd=engine.enclosure(scores,budget,outside,queries=dual)
                        assert bd['lower']<=ex['lower']<=ex['upper']<=bd['upper']
                        if budget==0:assert ex['lower']==ex['upper']==total
                        records.append({'corridor':ch['id'],'roi':ch['roi'],'setting':setting,'cell_id':cell,'cell_index':k,
                            'gene':gene,'budget_um2':str(budget),'source_count':total,'fixed_outside_count':outside,
                            'local_gene_transcripts':len(target),'lower':ex['lower'],'upper':ex['upper'],'status':'exact',
                            'dual_lower':bd['lower'],'dual_upper':bd['upper'],
                            'dual_closed':(bd['lower'],bd['upper'])==(ex['lower'],ex['upper']),
                            'dual_minimum_witness_value':bd['minimum_witness_value'],'dual_maximum_witness_value':bd['maximum_witness_value'],
                            'dual_self_closed':bd['status']=='exact',
                            'minimum_cost_for_lower':str(sp.minima[ex['lower']]),'minimum_cost_for_upper':str(sp.minima[ex['upper']]),
                            'cell_touches_crop_border':cell in obs['rois'][ch['roi']]['border_labels']})
            summary={'corridor':ch['id'],'setting':setting,'count_queries':len(metadata),'intervals':len(records),
                     'target_patterns_with_reuse':patterns,'unique_geometric_completions':len(engine._geometric_cache),
                     'exact_penalized_cuts':cuts,'fixed_nonzero_entries_outside_query_genes':sum(len(set(counts.get(str(c),{}))-set(genes)) for c in ch['cell_labels']),
                     'local_gene_max_multiplicity':max(collections.Counter(z['gene'] for z in ch['points']).values()),
                     'seconds':time.perf_counter()-t1}
            path.write_text(json.dumps({'input_chart':ch['id'],'summary':summary,'rows':records,'queries':metadata},default=enc,indent=2)+'\n')
            rows.extend(records);models.append(summary);new+=1
            print(f'{ci+1}/42 {setting}: {len(metadata)} count queries, {time.perf_counter()-t1:.2f}s',flush=True)
    with (ROOT/'results/expression_intervals.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    table=[]
    for setting in ['baseline','maximal']:
        for budget in budgets:
            rr=[r for r in rows if r['setting']==setting and r['budget_um2']==str(budget)]
            table.append({'setting':setting,'budget_um2':str(budget),'intervals':len(rr),
                'nontrivial_intervals':sum(r['lower']<r['upper'] for r in rr),'dual_closed':sum(r['dual_closed'] for r in rr),
                'dual_self_closed':sum(r['dual_self_closed'] for r in rr),
                'upper_excess':sum(r['dual_upper']-r['upper'] for r in rr),'lower_excess':sum(r['lower']-r['dual_lower'] for r in rr),
                'source_detection_can_disappear':sum(r['source_count']>0 and r['lower']==0 for r in rr),
                'source_zero_can_become_positive':sum(r['source_count']==0 and r['upper']>0 for r in rr),
                'maximum_interval_width':max(r['upper']-r['lower'] for r in rr)})
    lookup={(r['corridor'],r['setting'],r['cell_id'],r['gene'],r['budget_um2']):r for r in rows};nesting=0
    for r in rows:
        if r['setting']=='baseline':
            q=lookup[(r['corridor'],'maximal',r['cell_id'],r['gene'],r['budget_um2'])]
            assert q['lower']<=r['lower']<=r['upper']<=q['upper'];nesting+=1
    out={'scope':'42 selected corridors from one donor; local exact counts with observed crop exterior fixed.',
         'protocol_sha256':hashlib.sha256((ROOT/'docs/EXPERIMENT_PROTOCOL.json').read_bytes()).hexdigest(),
         'models':len(models),'count_queries':sum(r['count_queries'] for r in models),'intervals':len(rows),
         'exact_penalized_cuts':sum(r['exact_penalized_cuts'] for r in models),'target_patterns':sum(r['target_patterns_with_reuse'] for r in models),
         'geometric_completions':sum(r['unique_geometric_completions'] for r in models),'family_nesting_checks':nesting,
         'table':table,'model_records':models,'sum_model_seconds':sum(r['seconds'] for r in models)}
    (ROOT/'results/native_expression_summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='model_records'},indent=2))
if __name__=='__main__':main()
