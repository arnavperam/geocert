"""Replay the frozen scientific core in a disposable, checksum-verified copy.

Recorded measurements in the source release are never overwritten.
Model/query binding metadata was added in 0.2.0a2; separate scope tests verify it,
while historical numerical replay omits this new metadata field. Exact scientific
JSON values must match; only explicitly identified time/memory fields are omitted.
Optional source ingestion/raster checks require the analysis extras and are a
separate command, not silently treated as part of the dependency-free replay.
"""
from __future__ import annotations
from pathlib import Path
import argparse,copy,hashlib,json,os,shutil,subprocess,sys,tempfile,time,csv
ROOT=Path(__file__).resolve().parents[1]
TIMING={'seconds','elapsed','wall_seconds','sum_model_seconds','seconds_excluding_compilation',
        'median_process_seconds','max_peak_rss_mib','peak_rss_mib','peak_rss_kib','max_rss_kib',
        'cpp_peak_rss_kib','native_seconds','total_seconds','process_seconds','binding'}
def normalize(x):
    if isinstance(x,dict):return {k:normalize(v) for k,v in x.items() if k not in TIMING and not k.endswith('_seconds')}
    if isinstance(x,list):return [normalize(v) for v in x]
    return x

def verify_manifest(root):
    count=0
    for line in (root/'SHA256SUMS.txt').read_text().splitlines():
        h,p=line.split('  ',1);path=root/p
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=h:raise AssertionError('checksum mismatch: '+p)
        count+=1
    return count

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');ap.add_argument('--cpp',action='store_true');ap.add_argument('--output',type=Path,default=ROOT/'run/replay.json');a=ap.parse_args()
    started=time.perf_counter();n=verify_manifest(ROOT)
    commands=[
      ('tests/test_budget.py',['budget_tests.json','dual_gap.json']),
      ('tests/test_joint.py',['joint_tests.json','joint_presence_example.json']),
      ('tests/test_binding.py',['binding_tests.json']),
      ('tests/test_resources.py',['resource_tests.json']),
      ('tests/test_composition.py',['composition_tests.json']),
      ('tests/test_parametric_scope.py',['parametric_scope_tests.json']),
      ('tests/test_relaxation.py',['relaxation_tests.json']),
      ('tests/test_evaluation.py',['evaluation_tests.json']),
      ('tests/test_scope.py',['scope_tests.json']),
    ]
    if a.full:
        commands += [
         ('experiments/synthetic_calibration.py',['synthetic_calibration.json']),
         ('experiments/probability_diagnostic.py',['probability_diagnostic.json']),
         ('experiments/native_expression.py',['native_expression_summary.json']+['native/'+p.name for p in sorted((ROOT/'results/native').glob('*.json'))]),
         ('validation/check_native_intervals.py',['native_interval_replay.json']),
         ('experiments/compare_marginal_bounds.py',['marginal_comparison.json']),
         ('experiments/native_joint_counts.py',['native_joint_counts.json']),
         ('experiments/export_matrices.py',['expression_matrix_index.json']),
        ]
    if a.cpp:
        if not shutil.which('g++'):raise SystemExit('g++ requested but unavailable')
        commands += [('tests/test_cpp_flow.py',['cpp_flow_tests.json'])]
    reference={rel:normalize(json.loads((ROOT/'results'/rel).read_text())) for _,files in commands for rel in files}
    output=a.output.resolve();output.parent.mkdir(parents=True,exist_ok=True)
    logdir=output.parent/(output.stem+'_logs');logdir.mkdir(exist_ok=True)
    stages=[];compared=[]
    with tempfile.TemporaryDirectory(prefix='geocert-expression-replay-') as td:
        work=Path(td)/'source'
        shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('__pycache__','run','build','*.egg-info','.venv'))
        assert verify_manifest(work)==n
        (work/'run').mkdir(exist_ok=True)
        if a.full:
            shutil.rmtree(work/'results/native')
        for index,(script,files) in enumerate(commands):
            t0=time.perf_counter();proc=subprocess.run([sys.executable,script],cwd=work,text=True,capture_output=True,timeout=300)
            (logdir/f'{index:02d}.log').write_text(proc.stdout+'\n'+proc.stderr)
            rec={'script':script,'returncode':proc.returncode,'seconds':time.perf_counter()-t0,'compared':[]}
            if proc.returncode:
                output.write_text(json.dumps({'manifest_entries':n,'stages':stages+[rec],'complete':False},indent=2)+'\n')
                raise RuntimeError('replay failed: '+script+'\n'+proc.stdout[-1500:]+'\n'+proc.stderr[-1500:])
            for rel in files:
                got=normalize(json.loads((work/'results'/rel).read_text()))
                if got!=reference[rel]:
                    (logdir/'comparison_failure.json').write_text(json.dumps({'file':rel,'expected':reference[rel],'actual':got},indent=2))
                    raise AssertionError('scientific mismatch: '+rel)
                rec['compared'].append(rel);compared.append(rel)
            stages.append(rec);print(json.dumps({'passed':script,'scientific_files':len(files),'seconds':rec['seconds']}),flush=True)
        if a.full:
            for rel in ['expression_intervals.csv','expression_matrix_intervals.csv']:
                with (ROOT/'results'/rel).open() as f:left=list(csv.DictReader(f))
                with (work/'results'/rel).open() as f:right=list(csv.DictReader(f))
                assert left==right,'CSV mismatch: '+rel
                compared.append(rel)
    assert verify_manifest(ROOT)==n
    result={'manifest_entries':n,'python':sys.version,'mode':'full' if a.full else 'core','cpp':a.cpp,
            'stages':stages,'scientific_files_matched':len(compared),'matched_paths':compared,
            'source_manifest_unchanged_after_replay':True,'complete':True,'seconds':time.perf_counter()-started,
            'scope':'Internal same-environment reproduction from a checksum-verified disposable copy. Full mode regenerates all84 native models. Analysis-dependency source/raster studies and the3-run streaming timing census are separate.'}
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['stages','matched_paths']},indent=2))
if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O')
    main()
