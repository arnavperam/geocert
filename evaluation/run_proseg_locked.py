"""Run the predeclared published-method comparison; never tune on its outputs."""
from pathlib import Path
import hashlib, json, os, platform, subprocess, time, shutil
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'proseg_locked'
OUT.mkdir(parents=True, exist_ok=True)
binary = shutil.which('proseg')
if not binary:
    raise RuntimeError('Proseg binary missing')
version = subprocess.check_output([binary, '--version'], text=True).strip()
if '3.2.0' not in version:
    raise RuntimeError(version)
(OUT/'help.txt').write_text(subprocess.check_output([binary,'--help'],text=True))
record = {'version':version,'binary_sha256':hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
          'platform':platform.platform(),'protocol':'evaluation/global_patch_protocol.md',
          'scope':'published-method comparison on existing development crops; not truth or calibration', 'runs':[]}
for roi in ('roi_01','roi_02'):
    src = ROOT/'data'/'native'/roi/'transcripts.parquet'
    df = pq.read_table(src).to_pandas()
    inp = OUT/(roi+'_input.csv.gz')
    df.to_csv(inp,index=False,compression='gzip')
    for mode in ('fixed_xy','default'):
        target=OUT/(roi+'_'+mode)
        target.mkdir(exist_ok=True)
        cmd=[binary,'--xenium','--nthreads','2','--output-spatialdata','spatial.zarr',
             '--output-transcript-metadata','transcripts.csv.gz',
             '--output-cell-metadata','cells.csv.gz','--output-gene-metadata','genes.csv.gz',
             '--output-counts','counts.mtx.gz','--output-cell-polygons','polygons.geojson.gz']
        if mode=='fixed_xy':
            cmd += ['--ignore-z-coord','--voxel-layers','1','--no-diffusion']
        cmd += [str(inp)]
        entry={'roi':roi,'mode':mode,'command':cmd,'input_rows':len(df),
               'input_parquet_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
               'sampling_defaults':{'burnin':200,'samples':200,'hillclimb':50},
               'stochastic':'one realization; upstream does not expose a seed option'}
        t=time.monotonic()
        with (target/'stdout.log').open('w') as so,(target/'stderr.log').open('w') as se:
            try:
                p=subprocess.run(cmd,cwd=target,stdout=so,stderr=se,timeout=900)
                entry['returncode']=p.returncode
            except subprocess.TimeoutExpired:
                entry['returncode']='timeout_900s'
        entry['elapsed_seconds']=time.monotonic()-t
        entry['outputs']={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                          for p in target.iterdir() if p.is_file()}
        record['runs'].append(entry)
        (OUT/'run_record.json').write_text(json.dumps(record,indent=2)+'\n')
        # Preserve compact original tables and logs, not the larger duplicated zarr tree.
        shutil.rmtree(target/'spatial.zarr',ignore_errors=True)
(OUT/'run_record.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
if not any(r['returncode']==0 for r in record['runs']):
    raise SystemExit('No Proseg run completed; inspect preserved logs.')
