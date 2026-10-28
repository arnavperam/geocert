"""Exact cell-total count spectra for ONE validated disjoint-patch family."""
from __future__ import annotations
from fractions import Fraction as F
from collections import Counter,defaultdict
from pathlib import Path
import sys,json,hashlib,csv,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from global_model.engine import LocalEngine as ExpressionEngine
from mixed import make_problem,verify_mixed

BUDGETS=list(map(F,['0','0.05','0.1','0.25','0.5','1','2','4']))

def _encode(x):
    if isinstance(x,F):return str(x)
    raise TypeError(type(x))

def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':'),default=_encode)
def digest(x):return hashlib.sha256(canonical(x).encode()).hexdigest()

def patch_spectra(patch):
    p,pins,n=make_problem(patch)
    eng=ExpressionEngine(p,pins)
    if eng.source_engine is None:raise ValueError('source not admissible')
    src=eng.clipped({})
    if not src['feasible'] or F(src['area'])!=0:raise ValueError('invalid source')
    genes=sorted(set(z['gene'] for z in patch['points']))
    spectra={};witnesses={};singleton=[]
    for gene in genes:
        scores=[[int(g<n and patch['points'][g]['gene']==gene and o==0) for o in range(2)] for g in range(p.n)]
        s=eng.spectrum(scores,pattern_limit=65536)
        if not s.complete:raise RuntimeError('target enumeration limit reached: '+patch['id']+'/'+gene)
        ds={};ws={}
        for count,cost in s.minima.items():
            w=s.witnesses[count];key=digest(w['owners'])[:20]
            old=witnesses.get(key)
            if old is not None and old['area']!=str(cost):raise AssertionError('owner cost inconsistency')
            witnesses[key]={'owners':w['owners'],'area':str(cost),'boundaries':w['boundaries']}
            ds[str(count)]=str(cost);ws[str(count)]=key
        spectra[gene]={'local_total':sum(z['gene']==gene for z in patch['points']),
                       'source_left':sum(z['gene']==gene and z['source']==0 for z in patch['points']),
                       'minima_left':ds,'witness_ids':ws,'examined':s.examined}
    for g,z in enumerate(patch['points']):
        a=eng.clipped({g:1-z['source']})
        singleton.append({'transcript_id':z['transcript_id'],'gene':z['gene'],'source':z['source'],
                          'alternative_area':str(a['area']) if a['feasible'] else None})
    return {'id':patch['id'],'roi':patch['roi'],'cell_labels':patch['cell_labels'],
            'spectra':spectra,'witnesses':witnesses,'singleton':singleton,
            'source_zero_verified':True,'patch_sha256':digest(patch)}

def convolve(parts,fixed=0):
    """Min-plus product with one shared budget; retain complete backpointers."""
    state={int(fixed):(F(0),[])}
    for pid,spec in parts:
        nxt={}
        for a,(cost,wit) in state.items():
            for b,(extra,reference) in spec.items():
                value=a+b;newcost=cost+F(extra)
                if value not in nxt or newcost<nxt[value][0]:
                    nxt[value]=(newcost,wit+[[pid,reference]])
        state=nxt
    return state

def run():
    t0=time.perf_counter();family=json.loads((ROOT/'data/global_patches.json').read_text())
    observed=json.loads((ROOT/'data/observed_crop_counts.json').read_text())['rois']
    family_hash=digest(family);outdir=ROOT/'results/patch_spectra';outdir.mkdir(exist_ok=True)
    patchres={};new=0
    for j,patch in enumerate(family['patches']):
        path=outdir/(patch['id']+'.json')
        if path.exists():
            r=json.loads(path.read_text())
            if r['patch_sha256']!=digest(patch):raise ValueError('stale patch cache')
        else:
            r=patch_spectra(patch);path.write_text(json.dumps(r,indent=2,default=_encode)+'\n');new+=1
        patchres[patch['id']]=r
        if j%25==0:print('patch',j,'/',len(family['patches']),'elapsed',round(time.perf_counter()-t0,2),flush=True)
    memberships=defaultdict(list)
    complete={}
    for ss in family['census']:
        for cell in ss['cells']:complete[ss['roi'],cell['cell']]=cell['complete_observed_cell']
    for patch in family['patches']:
        for side,cell in enumerate(patch['cell_labels']):
            if cell:memberships[patch['roi'],cell].append((patch,side,patchres[patch['id']]))
    records=[];rows=[]
    for roi,obs in sorted(observed.items()):
        allcells=sorted(map(int,obs['counts']))
        for cell in allcells:
            if cell==0:continue
            raw=obs['counts'][str(cell)];incident=memberships[roi,cell]
            genes=set(raw)
            for patch,side,pr in incident:genes.update(pr['spectra'])
            for gene in sorted(genes):
                source=int(raw.get(gene,0));inside=0;parts=[];loss=[];gain=[]
                for patch,side,pr in incident:
                    if gene not in pr['spectra']:continue
                    s=pr['spectra'][gene];total=s['local_total'];left=s['source_left']
                    inside+=left if side==0 else total-left
                    d={}
                    for k,v in s['minima_left'].items():
                        count=int(k) if side==0 else total-int(k)
                        d[count]=(F(v),s['witness_ids'][k])
                    parts.append((patch['id'],d))
                    for mol in pr['singleton']:
                        if mol['gene']!=gene or mol['alternative_area'] is None:continue
                        (loss if mol['source']==side else gain).append(F(mol['alternative_area']))
                fixed=source-inside
                if fixed<0:raise AssertionError('transcript double counting')
                sp=convolve(parts,fixed)
                if source not in sp or sp[source][0]!=0:raise AssertionError('source missing from global spectrum')
                key=f'{roi}:{cell}:{gene}'
                rec={'key':key,'roi':roi,'cell':cell,'gene':gene,'complete_observed_cell':bool(complete.get((roi,cell),False)),
                     'source_count':source,'fixed_count':fixed,'patch_count':len(parts),
                     'model_id':family_hash,'minima':{str(k):str(v[0]) for k,v in sorted(sp.items())},
                     'witnesses':{str(k):v[1] for k,v in sorted(sp.items())}}
                records.append(rec)
                for delta in BUDGETS:
                    feasible=sorted(k for k,(cost,_) in sp.items() if cost<=delta)
                    lo,hi=feasible[0],feasible[-1]
                    marginal_lo=source-sum(c<=delta for c in loss)
                    marginal_hi=source+sum(c<=delta for c in gain)
                    if not 0<=marginal_lo<=lo<=source<=hi<=marginal_hi:raise AssertionError('invalid bound')
                    rows.append({'model_id':family_hash,'roi':roi,'cell':cell,'gene':gene,
                                 'complete_observed_cell':rec['complete_observed_cell'],
                                 'budget_um2':str(delta),'source_count':source,'lower':lo,'upper':hi,
                                 'status':'exact','incident_gene_patches':len(parts),
                                 'marginal_lower':marginal_lo,'marginal_upper':marginal_hi,
                                 'minimum_witness_area_um2':str(sp[lo][0]),'maximum_witness_area_um2':str(sp[hi][0])})
    keys=[(r['roi'],r['cell'],r['gene'],r['budget_um2']) for r in rows]
    if len(keys)!=len(set(keys)):raise AssertionError('duplicate cell/gene/budget')
    (ROOT/'results/cell_spectra.json').write_text(json.dumps(records,indent=2)+'\n')
    with (ROOT/'results/cell_intervals.csv').open('w') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    bybudget=[]
    for delta in BUDGETS:
        rr=[r for r in rows if r['budget_um2']==str(delta) and r['complete_observed_cell']]
        bybudget.append({'budget':str(delta),'complete_cell_entries':len(rr),
                         'nontrivial':sum(r['lower']!=r['upper'] for r in rr),
                         'width_sum':sum(r['upper']-r['lower'] for r in rr),
                         'max_width':max(r['upper']-r['lower'] for r in rr),
                         'tighter_than_marginals':sum(r['lower']>r['marginal_lower'] or r['upper']<r['marginal_upper'] for r in rr),
                         'detection_variable':sum(r['lower']==0<r['upper'] for r in rr)})
    summary={'model_id':family_hash,'new_patch_computations':new,'patches':len(patchres),
             'local_gene_spectra':sum(len(r['spectra']) for r in patchres.values()),
             'cell_gene_spectra':len(records),'interval_rows':len(rows),
             'complete_cells':sum(complete.values()),
             'complete_cell_gene_spectra':sum(r['complete_observed_cell'] for r in records),
             'multi_patch_gene_queries':sum(r['patch_count']>1 for r in records),
             'maximum_gene_patches':max(r['patch_count'] for r in records),
             'budgets':bybudget,'seconds':time.perf_counter()-t0,
             'scope':'One interval per observed cell/gene/budget in one global family. Only disjoint accepted patches move; complement fixed.'}
    (ROOT/'results/global_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':run()
