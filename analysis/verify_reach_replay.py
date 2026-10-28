"""Check all frozen interval and minimum-area outputs after a new width-study run."""
from pathlib import Path
import argparse, hashlib, json, csv
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input-dir',type=Path,required=True);a=ap.parse_args()
 expected=json.loads((R/'results/reach_scientific_hashes.json').read_text())
 for rel,h in expected['byte_hashes'].items():
  p=a.input_dir/rel
  if not p.is_file() or sha(p)!=h:raise ValueError('scientific replay mismatch: '+rel)
 for rel,h in expected['minimum_spectrum_hashes'].items():
  p=a.input_dir/rel;data=json.loads(p.read_text())
  canonical=[{'roi':x['roi'],'cell':x['cell'],'gene':x['gene'],'source_count':x['source_count'],'minima':x['minima']} for x in data]
  actual=hashlib.sha256(json.dumps(canonical,sort_keys=True,separators=(',',':')).encode()).hexdigest()
  if actual!=h:raise ValueError('minimum-area spectrum mismatch: '+rel)
 print(json.dumps({'status':'passed','byte_identical_files':len(expected['byte_hashes']),'exact_minimum_spectra_files':len(expected['minimum_spectrum_hashes']),'scope':'Timings and equally optimal witness choices are not compared.'}))
if __name__=='__main__':main()
