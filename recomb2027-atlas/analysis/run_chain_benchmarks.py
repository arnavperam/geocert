from pathlib import Path
import subprocess,sys,json,random,os,time
R=Path(__file__).resolve().parents[1];D=R/'results/chain_benchmarks';D.mkdir(exist_ok=True)
schedule=[(n,r) for n in [16,32,64,128] for r in range(3)];random.Random(214).shuffle(schedule)
(D/'schedule.json').write_text(json.dumps(schedule))
for n,r in schedule:
 p=D/f'n{n}_r{r}.json';log=D/f'n{n}_r{r}.log'
 with log.open('w') as f:
  try:
   x=subprocess.run([sys.executable,str(R/'analysis/benchmark_chain.py'),'--n',str(n),'--output',str(p)],stdout=f,stderr=subprocess.STDOUT,timeout=160,cwd=R)
   if x.returncode:print('ERROR',n,r,x.returncode,flush=True)
  except subprocess.TimeoutExpired:
   p.write_text(json.dumps({'n':n,'status':'wall_timeout','seconds':160}));print('TIMEOUT',n,r,flush=True)
 print(n,r,'complete',flush=True)
(D/'status').write_text('0')
