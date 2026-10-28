"""Prepared atlas identity and composition-precondition regression tests."""
from pathlib import Path
import copy, hashlib, json, os, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from atlas.binding import atlas_digest, verify_prepared_atlas

base=json.loads((ROOT/'data/atlas_roi_01.json').read_text())
audit=json.loads((ROOT/'results/atlas_geometric_validation.json').read_text())
sources=ROOT/'data/source_crops/roi_01'
report=verify_prepared_atlas(base,sources,audit)
assert report['patches']==276 and report['editable_transcripts']==627
rejected=[]

def reject(name, a, rec=None, source=sources):
    try:verify_prepared_atlas(a,source,audit if rec is None else rec)
    except (ValueError, KeyError):rejected.append(name)
    else:raise AssertionError('accepted invalid case '+name)

def rebind_for_structural_test(a):
    # This deliberately fabricated audit record is ONLY a unit-test fixture.
    # It exercises structural guards after the independent identity guards.
    a['atlas_id']=atlas_digest(a);rec=copy.deepcopy(audit)
    next(r for r in rec['coverage'] if r['roi']=='roi_01')['atlas_id']=a['atlas_id']
    return rec

x=copy.deepcopy(base);x['patches'][0]['slope']='2';reject('stale_content_digest',x)
x=copy.deepcopy(base);x['patches'][0]['slope']='2';x['atlas_id']=atlas_digest(x);reject('changed_digest_absent_from_geometry_audit',x)
rec=copy.deepcopy(audit);rec['failures']=1;reject('failed_geometry_audit',base,rec)
rec=copy.deepcopy(audit);rec['coverage'].append(copy.deepcopy(rec['coverage'][0]));reject('duplicated_geometry_audit_identity',base,rec)
x=copy.deepcopy(base);x['patches'].append(copy.deepcopy(x['patches'][0]));reject('duplicate_patch',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);x['records'].append(copy.deepcopy(x['records'][0]));reject('duplicate_observation',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);cell=next(iter(x['source_counts']));gene=next(iter(x['source_counts'][cell]));x['source_counts'][cell][gene]+=1;reject('incorrect_fixed_count',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);p=next(p for p in x['patches'] if p['points']);p['points'].append(copy.deepcopy(p['points'][0]));reject('duplicate_editable_observation',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);p=next(p for p in x['patches'] if p['points']);p['points'][0]['gene']='MUTATED';reject('gene_binding',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);p=next(p for p in x['patches'] if p['points']);p['points'][0]['x']='0';reject('coordinate_binding',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);x['patches'][0]['area_scale_um2']='1';reject('physical_area_scale',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);x['patches'][0]['band_halfwidth_pixels']=x['patches'][0]['halfwidth_pixels'];reject('non_strict_rails',x,rebind_for_structural_test(x))
x=copy.deepcopy(base);x['complete_cells'].append(x['complete_cells'][0]);reject('duplicate_complete_cell',x,rebind_for_structural_test(x))
with tempfile.TemporaryDirectory() as td:
    d=Path(td)
    for name in base['source_hashes']:(d/name).write_bytes((sources/name).read_bytes())
    (d/'metadata.json').write_text('{}')
    reject('original_source_hash',base,source=d)
    x=copy.deepcopy(base);metadata=json.loads((sources/'metadata.json').read_text());metadata['transform'][1][1]*=2
    (d/'metadata.json').write_text(json.dumps(metadata));x['source_hashes']['metadata.json']=hashlib.sha256((d/'metadata.json').read_bytes()).hexdigest()
    reject('anisotropic_source_transform',x,rebind_for_structural_test(x),d)
with tempfile.TemporaryDirectory() as td:
    d=Path(td);(d/'data').mkdir();(d/'results').mkdir()
    x=copy.deepcopy(base);x['patches'][0]['slope']='2'
    (d/'data/atlas_bad.json').write_text(json.dumps(x))
    (d/'results/atlas_geometric_validation.json').write_text(json.dumps(audit))
    env=dict(os.environ,GEOCERT_DATA_ROOT=str(d),GEOCERT_RESULT_DIR=str(d/'output'))
    proc=subprocess.run([sys.executable,str(ROOT/'analysis/certify_cells.py')],env=env,cwd=ROOT,capture_output=True,text=True)
    assert proc.returncode and 'atlas content digest mismatch' in proc.stderr
    assert not (d/'output/cell_expression_intervals.csv').exists()
    rejected.append('certification_entry_refuses_tampered_atlas_before_output')
for sub in [ROOT,ROOT/'external']:
    rec=json.loads((sub/'results/atlas_geometric_validation.json').read_text())
    for path in sorted((sub/'data').glob('atlas_*.json')):
        a=json.loads(path.read_text());verify_prepared_atlas(a,sub/'data/source_crops'/a['roi'],rec)
proc=subprocess.run([sys.executable,'-O',str(ROOT/'analysis/certify_cells.py')],capture_output=True,text=True)
assert proc.returncode and 'do not run Python with -O' in proc.stderr
rejected.append('optimized_python_mode')
result={'valid_prepared_atlases':4,'rejected_invalid_inputs':len(rejected),'cases':rejected,
        'certification_entry_guard':True,'failures':0,
        'scope':'Integrity checks bind to the saved exact geometry audit; they do not replace the rational geometry validator or defend against forged audit records.'}
(ROOT/'results/atlas_binding_tests.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
