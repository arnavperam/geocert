"""Read paired boundary records; calibrate and evaluate without fitting on test labels."""
from pathlib import Path
import argparse,json,sys,hashlib,gzip
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'core'))
from paired_validation import run_validation
from data_api import encode
if __name__=='__main__':
    if not __debug__: raise SystemExit('Keep assertions enabled; do not use -O')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path);ap.add_argument('--alpha',default='1/10')
    ap.add_argument('--budgets',nargs='+',default=['0','1/100','1/20','1/4'])
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    data=json.loads(gzip.decompress(a.input.read_bytes()) if a.input.suffix=='.gz' else a.input.read_bytes())
    result=run_validation(data,a.alpha,a.budgets)
    result['input_sha256']=hashlib.sha256(a.input.read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,default=encode)+'\n')
    print(json.dumps({k:result[k] for k in ['evidence_kind','calibration','summaries','pooled_error_ranking']},indent=2,default=encode))
