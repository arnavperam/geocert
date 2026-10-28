"""Strict-width inputs retain the same original-source integrity contract."""
from pathlib import Path
import copy, json, os, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from atlas.binding import verify_prepared_atlas

valid=[]
for family in sorted((ROOT/'reach').glob('*/*')):
    audit=json.loads((family/'results/atlas_geometric_validation.json').read_text())
    for path in sorted((family/'data').glob('atlas_*.json')):
        atlas=json.loads(path.read_text())
        report=verify_prepared_atlas(atlas,family/'data/source_crops'/atlas['roi'],audit)
        valid.append({'family':str(family.relative_to(ROOT)),**report})
assert len(valid)==16

family=ROOT/'reach/lung/r1'
atlas=json.loads((family/'data/atlas_roi_01.json').read_text())
audit=json.loads((family/'results/atlas_geometric_validation.json').read_text())
try:
    verify_prepared_atlas(atlas,ROOT/'external/data/source_crops/breast_01',audit)
except ValueError as exc:
    assert 'source-file hash mismatch' in str(exc)
else:raise AssertionError('wrong original source collection accepted')

with tempfile.TemporaryDirectory() as td:
    temp=Path(td);(temp/'data').mkdir();(temp/'results').mkdir()
    bad=copy.deepcopy(atlas);bad['patches'][0]['slope']='2'
    (temp/'data/atlas_bad.json').write_text(json.dumps(bad))
    (temp/'results/atlas_geometric_validation.json').write_text(json.dumps(audit))
    env=dict(os.environ,GEOCERT_DATA_ROOT=str(temp),GEOCERT_RESULT_DIR=str(temp/'output'))
    proc=subprocess.run([sys.executable,str(ROOT/'analysis/reach/certify_width.py')],
                        cwd=ROOT,env=env,capture_output=True,text=True)
    assert proc.returncode and 'atlas content digest mismatch' in proc.stderr
    assert not (temp/'output/cell_expression_intervals.csv').exists()

proc=subprocess.run([sys.executable,'-O',str(ROOT/'analysis/reach/certify_width.py')],
                    cwd=ROOT,capture_output=True,text=True)
assert proc.returncode and 'do not run Python with -O' in proc.stderr
result={'valid_strict_width_atlases':len(valid),'validations':valid,
        'invalid_cases_rejected':['wrong_original_source_collection',
                                  'width_entry_refuses_tampered_atlas_before_output',
                                  'width_entry_refuses_optimized_python'],
        'failures':0,
        'scope':'Existing saved geometry audits and original source bytes are bound to every strict-width input. This is an integrity check, not a new geometry or count experiment.'}
(ROOT/'results/width_binding_tests.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='validations'},indent=2))
