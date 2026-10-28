"""Descriptive C++ accumulation benchmark; candidate construction is excluded."""
from pathlib import Path
import sys,json,csv,subprocess,time,hashlib,statistics,random,platform
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.relaxation import candidate_count_intervals
run=ROOT/'run/streaming';run.mkdir(parents=True,exist_ok=True)
binary=run/'count_bounds'
subprocess.run(['g++','-std=c++17','-O2',str(ROOT/'cpp/count_bounds.cpp'),'-o',str(binary)],check=True)
# Exact cross-language comparison, not a timing claim.
r=random.Random(3911);small=[]
for i in range(5000):
    source=r.randrange(200);gene=r.randrange(16);cand=sorted({source,r.randrange(200),r.randrange(200)})
    small.append({'source':source,'gene':gene,'candidates':cand})
def write(path,records):
    with path.open('w') as f:
        for x in records:f.write(f"{x['gene']} {x['source']} {len(x['candidates'])} "+' '.join(map(str,x['candidates']))+'\n')
write(run/'small.tsv',small)
checks=0
for budget in [0,3,100,'all']:
    out=run/'small_output.tsv';subprocess.run([str(binary),str(run/'small.tsv'),str(out),str(budget)],check=True,capture_output=True)
    py=candidate_count_intervals(small,None if budget=='all' else budget)
    with out.open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==len(py)
    for rr in rows:
        exp=py[(int(rr['cell']),int(rr['gene']))]
        assert (int(rr['source']),int(rr['lower']),int(rr['upper']))==(exp['source_count'],exp['lower'],exp['upper']);checks+=1
records=[]
for n in [100000,1000000]:
    inp=run/f'generated_{n}.tsv'
    with inp.open('w') as f:
        for i in range(n):
            s=i%2000;g=(i//2000)%32;c=[s] if i%4==0 else [s,(s+1)%2000,(s+1999)%2000]
            f.write(f'{g} {s} {len(c)} '+' '.join(map(str,c))+'\n')
    expected=None
    for rep in range(3):
        out=run/'large_output.tsv';timing=run/'time.txt';t0=time.perf_counter()
        p=subprocess.run(['/usr/bin/time','-f','%e %M','-o',str(timing),str(binary),str(inp),str(out),'all'],check=True,capture_output=True,text=True)
        elapsed=time.perf_counter()-t0;external,rss=timing.read_text().split();digest=hashlib.sha256(out.read_bytes()).hexdigest()
        if expected is not None:assert digest==expected
        expected=digest
        records.append({'RNA_records':n,'repetition':rep,'wall_seconds':elapsed,'maximum_RSS_KiB':int(rss),
                        'output_sha256':digest,'input_bytes':inp.stat().st_size,'output_bytes':out.stat().st_size,'stdout':p.stderr.strip()})
    inp.unlink()
summary={'scope':'C++ candidate-bound accumulation including text input and output; excludes candidate extraction, geometry, inference, and network I/O.',
         'cross_language_entries_checked':checks,'runs':records,
         'table':[{'RNA_records':n,'median_wall_seconds':statistics.median(r['wall_seconds'] for r in records if r['RNA_records']==n),
                   'max_RSS_MiB':max(r['maximum_RSS_KiB'] for r in records if r['RNA_records']==n)/1024} for n in [100000,1000000]],
         'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
         'platform':platform.platform(),'no_whole_slide_runtime_claim':True}
(ROOT/'results/streaming_benchmark.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='runs'},indent=2))
