"""Recompute each local minimum independently; compose by Cartesian truth.

The source geometry and prepared observation coordinates are shared inputs.
Rational geometry acceptance and raw-data regeneration are separate stages.
"""
from pathlib import Path
import sys,json,csv,time
from fractions import Fraction as F
from itertools import product
from collections import defaultdict,Counter
import os
CODE_ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE_ROOT))
ROOT=Path(os.environ.get('GEOCERT_DATA_ROOT',str(CODE_ROOT))).resolve()
from atlas.scalar_reference import spectrum
from atlas.model import verify_witness

def main():
 t0=time.time();aa={a['roi']:a for a in [json.loads(p.read_text()) for p in sorted((ROOT/'data').glob('atlas_*.json'))]}
 pp={(roi,p['id']):p for roi,a in aa.items() for p in a['patches']}
 raw_lookup={roi:{q['id']:q for q in a['records']} for roi,a in aa.items()}
 raw_counts={roi:Counter((q['cell'],q['gene']) for q in a['records']) for roi,a in aa.items()}
 local=json.loads((ROOT/'results/patch_spectra.json').read_text());bykey={}
 comparisons=0;owners=0;witnesses=0
 for rec in local:
  p=pp[(rec['roi'],rec['patch'])];truth,n=spectrum(p,rec['gene']);owners+=n
  expected={int(k):F(v) for k,v in rec['minima'].items()}
  if truth!=expected:raise AssertionError(('local spectrum',rec['roi'],rec['patch'],rec['gene'],truth,expected))
  comparisons+=1;bykey[(rec['roi'],rec['patch'],rec['gene'])]=rec
  for key,w in rec['witnesses'].items():
   verify_witness(p,w['owners'],w['boundaries'],w['physical_area_um2'])
   count=sum(r['gene']==rec['gene'] and o==0 for r,o in zip(p['points'],w['owners']))
   if count!=int(key):raise AssertionError('local score witness')
   witnesses+=1
 spectra=json.loads((ROOT/'results/cell_spectra.json').read_text());global_checked=0;cartesian=0;global_witnesses=0
 for r in spectra:
  options=[]
  for pid in r['patches']:
   lp=bykey[(r['roi'],pid,r['gene'])];p=pp[(r['roi'],pid)];sign=1 if p['cell_labels'][0]==r['cell'] else -1
   options.append([(sign*(int(k)-lp['source_left']),F(v)) for k,v in lp['minima'].items()])
  truth={}
  for outcomes in product(*options):
   value=r['source_count']+sum(v for v,d in outcomes);cost=sum((d for v,d in outcomes),F(0));cartesian+=1
   truth[value]=min(truth.get(value,cost),cost)
  if truth!={int(v):F(d) for v,d in r['minima'].items()}:raise AssertionError('global spectrum differs from Cartesian product')
  for value,trace in r['witness_traces'].items():
   changed={};area=F(0)
   for pid,delta in trace.items():
    lp=bykey[(r['roi'],pid,r['gene'])];p=pp[(r['roi'],pid)];sign=1 if p['cell_labels'][0]==r['cell'] else -1
    localvalue=lp['source_left']+sign*delta;w=lp['witnesses'][str(localvalue)];area+=F(w['physical_area_um2'])
    for pt,owner in zip(p['points'],w['owners']):
     if pt['transcript_id'] in changed:raise AssertionError('duplicate observation across patches')
     changed[pt['transcript_id']]=p['cell_labels'][owner]
   actual=raw_counts[r['roi']][(r['cell'],r['gene'])]
   assert actual==r['source_count']
   for ident,newowner in changed.items():
    q=raw_lookup[r['roi']][ident]
    if q['gene']==r['gene']:actual+=int(newowner==r['cell'])-int(q['cell']==r['cell'])
   if actual!=int(value) or area!=F(r['minima'][value]):raise AssertionError('global physical witness/count')
   global_witnesses+=1
  global_checked+=1
 lookup={(r['roi'],r['cell'],r['gene']):r for r in spectra};rows=0;seen=set()
 with (ROOT/'results/cell_expression_intervals.csv').open() as f:
  for r in csv.DictReader(f):
   key=(r['roi'],int(r['cell']),r['gene']);budget=F(r['budget_um2']);ident=(*key,budget)
   if ident in seen:raise AssertionError('duplicated cell/gene/budget')
   seen.add(ident);s=lookup[key];v=[int(n) for n,d in s['minima'].items() if F(d)<=budget]
   if (min(v),max(v))!=(int(r['lower']),int(r['upper'])):raise AssertionError('invalid interval')
   rows+=1
 result={'independent_local_spectra':comparisons,'full_target_binary_patterns':owners,'local_witness_records':witnesses,
         'global_spectra':global_checked,'Cartesian_combinations':cartesian,'global_attaining_count_witnesses':global_witnesses,
         'unique_interval_rows':rows,'failures':0,'seconds':time.time()-t0,
         'independence':'Scalar reference imports no production implication, envelope, compiler or optimizer. Prepared extraction inputs remain shared; no independent donor or theorem review.'}
 (ROOT/'results/global_count_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
