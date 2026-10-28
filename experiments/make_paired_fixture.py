"""Generate frozen-protocol fixtures, not realistic segmentation benchmarks."""
from pathlib import Path
import json,random,hashlib,sys,copy,gzip
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'core'))
from data_api import encode


def generate():
    raw=(ROOT/'docs/validation/SYNTHETIC_PROTOCOL.json').read_bytes(); plan=json.loads(raw)
    rng=random.Random(plan['seed']); records=[]
    sites=[{'x':str(x),'t':str(t),'gene':gene} for x,t,gene in
           [(F(-5,64),1,'SYN_A'),(F(5,64),3,'SYN_A'),
            (F(-13,64),1,'SYN_B'),(F(13,64),3,'SYN_B')]]
    for split,n in [('calibration',plan['calibration_units']),('evaluation',plan['evaluation_units'])]:
        for u in range(n):
            for j in range(plan['records_per_unit']):
                while True:
                    h=F(rng.randint(-900,900),2000)
                    if all(abs(F(s['x'])-h/2)>=F(1,100) for s in sites): break
                records.append({'unit_id':f'synthetic-{split}-{u:03d}','sample_id':f'pair-{j}',
                    'split':split,'sites':copy.deepcopy(sites),
                    'model':{'baselines':[[[0,0],[4,0]]],'source':[0,1,0,1],'slopes':[1],
                             'minimum_gaps':[],'maximum_gaps':[],'clearance':'1/100',
                             'bands':[['-1/2','1/2']],'mode':'nominal'},
                    'truth_owners':[int(F(s['x'])>h/2) for s in sites],
                    'truth_boundaries':[[[0,2,0,h],[2,4,h,0]]],
                    'targets':[{'gene':g,'owner':o,'outside':rng.randrange(4)}
                               for g in ['SYN_A','SYN_B'] for o in [0,1]]})
    base={'schema':'geocert.paired_boundaries.v1','protocol_id':plan['protocol_id'],
          'protocol_sha256':hashlib.sha256(raw).hexdigest(),'evidence_kind':'synthetic_fixture',
          'coordinate_unit':'synthetic_length_unit','outside_assignments_fixed':True,'training_units':[],
          'provenance':'Generated two-cell triangular-boundary software fixture. No donors, microscopy, real segmenter or tissue shape model.',
          'records':records}
    output={}
    for scenario in plan['scenarios']:
        d=copy.deepcopy(base); d['scenario']=scenario
        for r in d['records']:
            idx=int(r['unit_id'].rsplit('-',1)[1])
            bad=(scenario=='evaluation_shift' and r['split']=='evaluation' and idx%5==0 or
                 scenario=='calibration_failure' and r['split']=='calibration' and idx<5)
            if bad and r['sample_id']=='pair-0':
                r['truth_boundaries']=[[[0,2,0,F(3,4)],[2,4,F(3,4),0]]]
                r['truth_owners']=[int(F(s['x'])>F(3,8)) for s in r['sites']]
        output[scenario]=d
    return output


if __name__=='__main__':
    dest=ROOT/'data/paired_validation';dest.mkdir(exist_ok=True)
    for name,d in generate().items():
        with (dest/(name+'.json.gz')).open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:
                g.write((json.dumps(d,indent=2,default=encode)+'\n').encode())
    print('Three labeled synthetic fixtures written')
