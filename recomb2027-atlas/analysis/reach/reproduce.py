"""Recompute one frozen width family without overwriting stored results."""
from pathlib import Path
import argparse, csv, json, os, subprocess, sys
R=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',choices=['lung','breast'],required=True)
    p.add_argument('--width',choices=['0.5','1','2','3'],required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args(); base=R/'reach'/a.dataset/('r'+a.width.replace('.','p'))
    out=a.output_dir.resolve()
    if out==(base/'results').resolve():raise ValueError('Choose a fresh output directory')
    env=os.environ.copy();env.update(GEOCERT_DATA_ROOT=str(base),GEOCERT_RESULT_DIR=str(out),GEOCERT_COUNT_METHOD='chain')
    subprocess.run([sys.executable,str(R/'analysis/reach/certify_width.py')],cwd=R,env=env,check=True)
    def spectra(path):return {(r['roi'],r['cell'],r['gene']):r['minima'] for r in json.loads(path.read_text())}
    if spectra(base/'results/cell_spectra.json')!=spectra(out/'cell_spectra.json'):raise ValueError('Exact spectra differ')
    with (base/'results/cell_expression_intervals.csv').open() as expected, (out/'cell_expression_intervals.csv').open() as actual:
        from itertools import zip_longest
        count=0
        for e,v in zip_longest(csv.DictReader(expected),csv.DictReader(actual)):
            if e!=v:raise ValueError(f'Interval row {count} differs')
            count+=1
    print(json.dumps({'status':'exact match','dataset':a.dataset,'width':a.width,'rows':count}))
if __name__=='__main__':main()
