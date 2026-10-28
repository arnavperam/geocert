"""Command-line access to explicit, model-conditioned count certificates."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,sys
from .budget import ExpressionEngine
from .calibration import threshold
from .relaxation import candidate_count_intervals
from mixed import make_problem
from data_api import variant

def encode(v):
    if isinstance(v,F):return str(v)
    raise TypeError(type(v))

def main(argv=None):
    parser=argparse.ArgumentParser(prog='geocert-expression',description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    q=sub.add_parser('certify',help='gene-count intervals for a supplied corridor')
    q.add_argument('input',type=Path);q.add_argument('--output',type=Path,required=True)
    q.add_argument('--method',choices=['auto','bounds'],default='auto')
    q.add_argument('--max-patterns',type=int,default=100000);q.add_argument('--cpp-binary',type=Path)
    q=sub.add_parser('demo',help='run the packaged native PECAM1 example');q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('calibrate',help='split-conformal threshold from held-out region scores')
    q.add_argument('input',type=Path);q.add_argument('--alpha',default='1/10');q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('candidate-bounds',help='outer intervals from supplied possible cell labels')
    q.add_argument('input',type=Path);q.add_argument('--rna-budget',type=int);q.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if not __debug__:parser.error('run without -O; the inherited exact engine uses scientific assertions')
    if args.command in ('certify','demo'):
        path=Path(__file__).parent/'examples/native_count.json' if args.command=='demo' else args.input
        raw=path.read_bytes();r=json.loads(raw);ch=r['chart'];p,pins,n=make_problem(ch,r.get('mode','nominal'),F(r.get('rna_radius','0')))
        if r.get('setting','baseline')=='maximal':p=variant(p)
        elif r.get('setting','baseline')!='baseline':raise ValueError('unknown geometry setting')
        cell=r['cell_index'];fixed=r['fixed_outside_count']
        if type(cell) is not int or not 0<=cell<=p.m or type(fixed) is not int or fixed<0:raise ValueError('invalid cell or exterior count')
        gene=str(r['gene']);target=[g for g,z in enumerate(ch['points']) if z['gene']==gene]
        scores=[[int(g in target and o==cell) for o in range(p.m+1)] for g in range(p.n)]
        engine=ExpressionEngine(p,pins,cpp_binary=getattr(args,'cpp_binary',None));method=getattr(args,'method','auto')
        spectrum=None
        if method=='auto' and engine.source_engine is not None:
            spectrum=engine.spectrum(scores,fixed,pattern_limit=getattr(args,'max_patterns',100000))
        if spectrum is not None and spectrum.complete:
            intervals=[spectrum.interval(F(d)) for d in r['budgets']];certificate={'kind':'complete_target_projection','spectrum':spectrum.serializable()}
        else:
            queries=engine.dual_queries(scores,[0,F(1,16),F(1,4),1,4,16,64])
            intervals=[engine.enclosure(scores,F(d),fixed,queries=queries) for d in r['budgets']]
            certificate={'kind':'rational_lagrange_cuts','queries':queries,'exact_refinement_complete':False}
        out={'schema':'geocert-expression-0.2','input_sha256':hashlib.sha256(raw).hexdigest(),'corridor':ch['id'],
             'gene':gene,'cell_index':cell,'cell_id':ch['cell_labels'][cell],
             'budget_metric':'summed_interface_deformation_um2','fixed_outside_count':fixed,
             'source_count':fixed+sum(ch['points'][g]['source']==cell for g in target),
             'intervals':intervals,'certificate':certificate,
             'interpretation':'deterministic count bounds conditional on this geometric family; not a calibrated probability or corrected-expression estimate'}
    elif args.command=='calibrate':
        out=threshold(json.loads(args.input.read_text())['scores'],F(args.alpha))
    else:
        records=json.loads(args.input.read_text())['records'];rows=candidate_count_intervals(records,args.rna_budget)
        out={'scope':'independent candidate-owner relaxation','rows':[{'cell':c,'gene':g,**v} for (c,g),v in rows.items()]}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,default=encode,indent=2)+'\n')
    print(str(args.output))
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,KeyError,TypeError,TimeoutError) as e:
        print(f'geocert-expression: {e}',file=sys.stderr);raise SystemExit(2)
