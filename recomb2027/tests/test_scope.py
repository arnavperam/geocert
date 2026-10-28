"""Prevent accidental whole-cell claims or gluing of overlapping corridors."""
from pathlib import Path
import copy,itertools,json,random,sys,tempfile
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.scope import corridor_scope,union_intervals
from geocert_expression.budget import ScoreSpectrum,compose_spectra,ExpressionEngine
from geocert_expression.cli import main
from mixed import make_problem

def record(family, values):
    return {'scope':{'kind':'single_corridor_conditional','reference_id':'observations-1','family_id':family},
            'query':{'roi':'r','cell_id':1,'gene':'G','budget':'1','budget_metric':'area'},
            'source_count':3,'status':'exact','lower':min(values),'upper':max(values),
            'minimum_witness':{'value':min(values)},'maximum_witness':{'value':max(values)}}

rng=random.Random(713);checks=0
for n in range(1,31):
    sets=[set(rng.sample(range(11),rng.randrange(1,10))) for _ in range(n)]
    rows=[record(str(i),v) for i,v in enumerate(sets)]
    result=union_intervals(rows);truth=set().union(*sets)
    assert (result['lower'],result['upper'])==(min(truth),max(truth))
    assert result['minimum_witness']['value']==min(truth)
    assert result['maximum_witness']['value']==max(truth)
    assert result['scope']['simultaneous_corridor_deformation'] is False
    assert result['scope']['whole_cell_boundary_uncertainty'] is False
    assert union_intervals(list(reversed(rows)))['scope']['family_id']==result['scope']['family_id']
    checks+=1
rows=[record('a',[2]),record('b',[2,3])]
assert (union_intervals(rows)['lower'],union_intervals(rows)['upper'])==(2,3)
# Intersection would silently discard an admitted assignment.
assert max(r['lower'] for r in rows)==min(r['upper'] for r in rows)==2
bad=[]
for field,value in [('cell_id',2),('gene','H'),('budget','2'),('roi','other'),('budget_metric','pixels')]:
    r=copy.deepcopy(rows);r[1]['query'][field]=value;bad.append(r)
for field,value in [('reference_id','other'),('reference_id',None),('family_id','a'),('kind','whole_cell')]:
    r=copy.deepcopy(rows);r[1]['scope'][field]=value;bad.append(r)
r=copy.deepcopy(rows);r[1]['source_count']=4;bad.append(r)
r=copy.deepcopy(rows);r[1]['status']='feasibility_unresolved';bad.append(r)
r=copy.deepcopy(rows);r[1]['lower']=5;bad.append(r)
bad.extend([[],[dict(rows[0],scope=dict(rows[0]['scope'],reference_id=None))]])
for field in ('roi','cell_id','gene','budget_metric'):
    r=copy.deepcopy(rows)
    for x in r:del x['query'][field]
    bad.append(r)
r=copy.deepcopy(rows)
for x in r:x['query']['budget']='-1'
bad.append(r)
r=copy.deepcopy(rows);r[1]['lower']='2';bad.append(r)
r=copy.deepcopy(rows)
for x in r:x['source_count']=-1
bad.append(r)
rejected=0
for rows_bad in bad:
    try:union_intervals(rows_bad)
    except ValueError:rejected+=1
    else:raise AssertionError('invalid union accepted')
sp=ScoreSpectrum({0:F(0),1:F(1)},{},True,2,2,0)
try:compose_spectra([sp,sp])
except ValueError:rejected+=1
else:raise AssertionError('unverified geometric composition accepted')
fixture=json.loads((ROOT/'geocert_expression/examples/native_count.json').read_text())
fixture['reference_id']='native-test-reference'
with tempfile.TemporaryDirectory() as d:
    p=Path(d);(p/'input.json').write_text(json.dumps(fixture))
    main(['certify',str(p/'input.json'),'--output',str(p/'output.json')])
    out=json.loads((p/'output.json').read_text())
    assert out['scope']['whole_cell_boundary_uncertainty'] is False
    assert all(q['scope']==out['scope'] for q in out['intervals'])
    assert all(q['binding']['model_key'] and q['binding']['score_key'] for q in out['intervals'])
    changed=copy.deepcopy(fixture['chart']);changed['baselines'][0][0][1]='17/2'
    assert corridor_scope(changed)['family_id']!=corridor_scope(fixture['chart'])['family_id']
    assert corridor_scope(fixture['chart'],'maximal')['family_id']!=corridor_scope(fixture['chart'])['family_id']
    assert corridor_scope(fixture['chart'],rna_radius='0.1')['family_id']!=corridor_scope(fixture['chart'])['family_id']
    assert corridor_scope(fixture['chart'],rna_radius='0.1')['family_id']==corridor_scope(fixture['chart'],rna_radius='1/10')['family_id']
summary={'exact_union_set_checks':checks,'unsafe_operations_rejected':rejected,
         'native_cli_scope_and_binding_verified':True,'family_changes_detected':3,'failures':0}
(ROOT/'results/scope_tests.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
