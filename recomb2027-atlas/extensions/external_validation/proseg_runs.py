"""Run matched core-only and buffered Proseg comparisons, retaining default failures."""
from pathlib import Path
import argparse,subprocess,os,json,time,hashlib,concurrent.futures,shutil

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def execute(job):
 roi,side,variant,binary=job;dest=roi/'proseg'/f'{side}_{variant}';dest.mkdir(parents=True,exist_ok=True)
 inp=roi/'context_transcripts.csv.gz' if side==800 else roi/'data/source_crops'/roi.name/'transcripts.csv.gz'
 source=inp;inp=dest/'input.csv.gz';shutil.copyfile(source,inp)
 cmd=[str(binary),'--xenium',str(inp.resolve()),'--nthreads','2','--output-path',str(dest.resolve()),'--output-transcript-metadata','transcript-metadata.csv.gz','--output-cell-metadata','cell-metadata.csv.gz','--output-gene-metadata','gene-metadata.csv.gz','--output-counts','counts.mtx.gz']
 if variant=='no_factorization':cmd.append('--no-factorization')
 rec={'command':cmd,'source_input_path':str(source.resolve()),'input_sha256':sha(inp),'binary_sha256':sha(binary),'input_side_um':side,'variant':variant,'version':subprocess.check_output([str(binary),'--version'],text=True).strip(),'settings':'200 burnin, 200 sampling, 50 hillclimb; 2 threads; stochastic upstream run, no seed option'}
 start=time.time()
 try:
  with (dest/'run.log').open('w') as log:p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=5400)
  rec.update(returncode=p.returncode,seconds=time.time()-start)
 except subprocess.TimeoutExpired:rec.update(returncode=None,status='timed_out',seconds=time.time()-start)
 (dest/'run.json').write_text(json.dumps(rec,indent=2)+'\n');print(roi.name,side,variant,rec['returncode'],round(rec['seconds'],2),flush=True)
 return rec
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path,required=True);ap.add_argument('roots',type=Path,nargs='+');ap.add_argument('--workers',type=int,default=2);a=ap.parse_args();jobs=[(r,side,v,a.binary.resolve()) for r in a.roots for side in [200,800] for v in ['default','no_factorization']]
 with concurrent.futures.ThreadPoolExecutor(a.workers) as ex:list(ex.map(execute,jobs))
