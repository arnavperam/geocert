"""Run frozen analysis stages for each predefined region, retaining all failures."""
from pathlib import Path
import argparse,subprocess,os,json,time,concurrent.futures
ROOT=Path(__file__).resolve().parents[2]
def run(roi):
 env={**os.environ,'GEOCERT_DATA_ROOT':str(roi.resolve()),'GEOCERT_COUNT_METHOD':'enumeration','OPENBLAS_NUM_THREADS':'1'};logs=roi/'logs';logs.mkdir(exist_ok=True);stages=[]
 for stage in ['build_atlas','verify_atlas','certify_cells','additive_baseline','check_global_counts','nucleus_reference']:
  start=time.time()
  try:
   with (logs/(stage+'.log')).open('w') as log:
    p=subprocess.run(['python3','-u',str(ROOT/'analysis'/f'{stage}.py')],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=3*3600)
   row={'stage':stage,'returncode':p.returncode,'seconds':time.time()-start}
  except subprocess.TimeoutExpired:row={'stage':stage,'returncode':None,'status':'timed_out','seconds':time.time()-start}
  stages.append(row);(roi/'run_status.json').write_text(json.dumps(stages,indent=2)+'\n');print(roi.name,row,flush=True)
  if row['returncode']!=0:break
 return roi.name
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('roots',type=Path,nargs='+');ap.add_argument('--workers',type=int,default=3);a=ap.parse_args()
 with concurrent.futures.ThreadPoolExecutor(a.workers) as ex:list(ex.map(run,a.roots))
