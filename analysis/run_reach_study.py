from pathlib import Path
import os,sys,json,subprocess,time
R=Path(__file__).resolve().parents[1]
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,default=R/'run/reach');args=ap.parse_args()
S=args.output_dir.resolve();S.mkdir(parents=True,exist_ok=True)
BUDGETS='0,.01,.025,.05,.075,.1,.15,.2,.25,.35,.5,.75,1,2,4,8,16'
ledger=[]
for name,widths in [('cap4','4,2,1,1/2'),('cap8','8,4,2,1,1/2'),('cap16','16,8,4,2,1,1/2'),('cap20','20,16,8,4,2,1,1/2')]:
 for dataset,src in [('lung',R),('breast',R/'external')]:
  d=S/name/dataset; (d/'data').mkdir(parents=True,exist_ok=True);(d/'results').mkdir(exist_ok=True);(d/'logs').mkdir(exist_ok=True)
  for x,y in [(src/'data/source_crops',d/'data/source_crops'),(src/'data/observed_crop_counts.json',d/'data/observed_crop_counts.json'),(src/'docs',d/'docs')]:
   if not y.exists():y.symlink_to(os.path.relpath(x,y.parent),target_is_directory=x.is_dir())
  env=os.environ.copy();env.update(GEOCERT_DATA_ROOT=str(d),GEOCERT_WIDTHS=widths,GEOCERT_BUDGETS=BUDGETS,GEOCERT_COUNT_METHOD='enumeration',PYTHONHASHSEED='0')
  for cmd in ['build_atlas.py','verify_atlas.py','certify_cells.py','coverage.py','nucleus_reference.py']:
   t=time.time(); log=d/'logs'/f'{cmd}.log';print(name,dataset,cmd,'start',flush=True)
   with log.open('w') as f:
    try:proc=subprocess.run([sys.executable,str(R/'analysis'/cmd)],env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800); code=proc.returncode
    except subprocess.TimeoutExpired:code=124
   rec={'width':name,'dataset':dataset,'stage':cmd,'returncode':code,'seconds':time.time()-t,'log':str(log.relative_to(S))};ledger.append(rec)
   (S/'execution.json').write_text(json.dumps(ledger,indent=2)+'\n');print(rec,flush=True)
   if code:print(log.read_text()[-4000:],flush=True);raise SystemExit(code)
print('ALL COMPLETE: compare outputs with reach_study using verify_reach_replay.py.',flush=True)
