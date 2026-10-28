"""Integration, metric and failure-path checks for supplied paired truth."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,copy,gzip
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'core'),str(ROOT/'experiments')]
from paired_validation import run_validation,ranking_metrics,validate_dataset
from make_paired_fixture import generate
from data_api import encode


def main():
    d=generate(); outputs={}
    for name,data in d.items():
        assert json.loads(json.dumps(data,default=encode))==json.loads(gzip.decompress((ROOT/'data/paired_validation'/f'{name}.json.gz').read_bytes()))
        out=run_validation(data);outputs[name]=out
        dest=ROOT/'run/paired_validation'/f'{name}.json.gz';dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:
                g.write((json.dumps(out,indent=2,default=encode)+'\n').encode())
    normal,shift,failure=[outputs[k] for k in ['nominal','evaluation_shift','calibration_failure']]
    assert normal['calibration']==shift['calibration']
    assert failure['calibration']['status']=='abstain' and failure['calibration']['budget'] is None
    assert sum(r['family_error'] is not None for r in normal['boundary_audit'])==0
    assert sum(r['family_error'] is not None for r in shift['boundary_audit'])==12
    assert sum(r['family_error'] is not None for r in failure['boundary_audit'])==5
    for out in outputs.values():
        assert all(r['covered'] for r in out['intervals'] if r['family_contained'])
        assert len(out['evaluation_units'])==60
        assert all(s['query_rows']==480 for s in out['summaries'])
    assert ranking_metrics([0,1],[0,1])['auroc']==1
    assert ranking_metrics([0,1],[1,0])['auroc']==0
    tied=ranking_metrics([0,1,0,1],[1,1,1,1])
    assert tied['auroc']==F(1,2) and tied['average_precision']==F(1,2)
    assert ranking_metrics([1,1],[1,2])['status']=='single_class'
    bad=copy.deepcopy(d['nominal']);bad['records'][-1]['unit_id']=bad['records'][0]['unit_id']
    bad['records'][-1]['sample_id']='different-sample'
    try: validate_dataset(bad);raise AssertionError('Cross-split unit leakage accepted')
    except ValueError as e: assert 'split' in str(e)
    bad=copy.deepcopy(d['nominal']);bad['training_units']=[bad['records'][0]['unit_id']]
    try: validate_dataset(bad);raise AssertionError('Training leakage accepted')
    except ValueError as e: assert 'Training' in str(e)
    bad=copy.deepcopy(d['nominal']);bad['records'][0]['model']['nuclear_boxes']=[]
    try: validate_dataset(bad);raise AssertionError('Ignored structural regions')
    except ValueError as e: assert 'Unsupported' in str(e)
    bad=copy.deepcopy(d['nominal']);bad['records'][0]['truth_boundaries']=[[]]
    try: validate_dataset(bad);raise AssertionError('Missing truth accepted')
    except ValueError as e: assert 'boundary' in str(e)
    bad=copy.deepcopy(d['nominal']);bad.pop('coordinate_unit')
    try: validate_dataset(bad);raise AssertionError('Unknown units accepted')
    except ValueError as e: assert 'unit' in str(e)
    # Independent rank formula check against explicit positive-negative pairs.
    y=[0,1,1,0,1];s=[1,2,2,3,1]
    expected=sum(F(int(s[i]>s[j]))+F(int(s[i]==s[j]),2) for i in range(5) if y[i]
                 for j in range(5) if not y[j])/6
    assert ranking_metrics(y,s)['auroc']==expected
    summary={'passed':True,'scope':'synthetic software fixtures only','scenarios':{}}
    for k,o in outputs.items():
        cal=o['calibration'];level=cal['budget']
        summary['scenarios'][k]={'calibration':cal,'evaluation_units':len(o['evaluation_units']),
             'at_calibrated_budget':next((s for s in o['summaries'] if s['budget']==level),None),
             'family_violations':sum(r['family_error'] is not None for r in o['boundary_audit']),
             'ranking':o['pooled_error_ranking']}
    (ROOT/'run/paired_validation/tests.json').write_text(json.dumps(summary,indent=2,default=encode)+'\n')
    print(json.dumps(summary,indent=2,default=encode))


if __name__=='__main__':main()
