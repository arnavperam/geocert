"""Documented small-crop amendment: --no-factorization on BOTH supplied crops."""
from pathlib import Path
import hashlib,json,subprocess,shutil,time
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/proseg_unfactorized';OUT.mkdir(parents=True,exist_ok=True)
binary=shutil.which('proseg');assert binary
version=subprocess.check_output([binary,'--version'],text=True).strip();assert '3.2.0' in version
record={'version':version,'binary_sha256':hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
        'amendment':'Both crops use upstream --no-factorization after ROI1 feature selection returned EmptyInput. No source modification.', 'runs':[]}
for roi in ('roi_01','roi_02'):
    inp=ROOT/'results/proseg_locked'/(roi+'_input.csv.gz')
    for mode in ('fixed_xy','default'):
        folder=OUT/(roi+'_'+mode);folder.mkdir(exist_ok=True)
        cmd=[binary,'--xenium','--nthreads','2','--no-factorization',
             '--output-spatialdata','spatial.zarr','--output-transcript-metadata','transcripts.csv.gz',
             '--output-cell-metadata','cells.csv.gz','--output-gene-metadata','genes.csv.gz',
             '--output-counts','counts.mtx.gz','--output-cell-polygons','polygons.geojson.gz']
        if mode=='fixed_xy':cmd+=['--ignore-z-coord','--voxel-layers','1','--no-diffusion']
        cmd+=[str(inp)]
        t=time.monotonic()
        try:
            r=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=900)
            entry={'roi':roi,'mode':mode,'command':cmd,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
        except subprocess.TimeoutExpired as e:
            entry={'roi':roi,'mode':mode,'command':cmd,'returncode':'timeout_900s','stdout':str(e.stdout),'stderr':str(e.stderr)}
        entry['seconds']=time.monotonic()-t
        shutil.rmtree(folder/'spatial.zarr',ignore_errors=True)
        entry['outputs']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir() if p.is_file()}
        record['runs'].append(entry);(OUT/'run_record.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
