"""Recheck count accounting and continuous witnesses without running optimization.

This checker verifies feasibility and accounting, not the upper-optimality claim.
The saved height-reference checks and full recomputation test those separately.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'core'))
from data_api import load_cases
from mixed import verify_mixed
from layered import realize

def check(path):
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    data=json.loads(path.read_text())
    models={(ch['id'],setting):(ch,p,pins,n) for _,ch,setting,p,pins,n,_ in load_cases()}
    cache={};records=0
    for record in data['rows']:
        ch,p,pins,n=models[record['chart'],record['setting']]
        k=record['cell_index'];gene=record['gene']
        assert ch['cell_labels'][k]==record['cell_label']
        original=sum(z['source']==k and z['gene']==gene for z in ch['points'])
        assert original==record['source_local_count']
        assert record['source_count']==record['fixed_count']+original
        for name in ['minimum','maximum']:
            answer=record[name];owners=list(p.source)
            owners[:n]=answer['owners']
            for group,label in pins.items():assert owners[group]==label
            observed=record['fixed_count']+sum(owners[g]==k and z['gene']==gene for g,z in enumerate(ch['points']))
            assert observed==answer['count']
            key=(record['chart'],record['setting'],tuple(owners))
            if key not in cache:
                realized=realize(p,owners)
                assert realized['feasible']
                verify_mixed(p,owners,realized['boundaries'],realized['area'])
                cache[key]=realized['area']
            assert cache[key]==F(answer['area'])
            records+=1
        assert record['minimum']['count']<=record['source_count']<=record['maximum']['count']
    expected={(r['chart'],r['setting'],r['cell_index'],r['gene']):r for r in data['rows']}
    for record in data['penalized']:
        ch,p,pins,n=models[record['chart'],'baseline'];k=ch['cell_labels'].index(record['cell'])
        old=expected[ch['id'],'baseline',k,record['gene']]
        owners=list(p.source);owners[:n]=record['owners']
        count=old['fixed_count']+sum(owners[g]==k and z['gene']==record['gene'] for g,z in enumerate(ch['points']))
        assert count==record['count']
        key=(ch['id'],'baseline',tuple(owners))
        if key not in cache:
            a=realize(p,owners);assert a['feasible'];verify_mixed(p,owners,a['boundaries'],a['area']);cache[key]=a['area']
        assert cache[key]==F(record['area'])
        sign=1 if record['direction']=='decrease' else -1
        assert F(record['lambda'])*cache[key]+sign*count==F(record['objective_with_fixed'])
        records+=1
    print(json.dumps({'extrema_pairs':len(data['rows']),'penalized':len(data['penalized']),
                      'witness_records':records,'distinct_model_witnesses':len(cache),'failures':0}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);check(p.parse_args().input)
