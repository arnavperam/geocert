"""Reindex archived interval results by their actual family; no study rerun.

The union output quantifies over a choice of ONE corridor family. It cannot be
used as a whole-cell deformation guarantee or as simultaneous corridor gluing.
"""
from collections import defaultdict,Counter
from pathlib import Path
from fractions import Fraction
import csv,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geocert_expression.scope import corridor_scope,digest,union_intervals

def read_csv(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))

def write_csv(path, rows):
    rows=list(rows)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

charts={c['id']:c for c in json.loads((ROOT/'data/charts.json').read_text())}
source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted((ROOT/'data/source_crops').glob('*/*')) if p.is_file()}
ref=digest(source_hashes)
index=json.loads((ROOT/'results/expression_matrix_index.json').read_text())
scopes={(m['corridor'],m['setting']):corridor_scope(charts[m['corridor']],m['setting'],reference_id=ref)
        for m in index['models']}
rows=read_csv(ROOT/'results/expression_intervals.csv')
by_question=defaultdict(list);cells=defaultdict(set)
for row in rows:
    by_question[tuple(row[k] for k in ('roi','setting','cell_id','gene','budget_um2'))].append(row)
    cells[(row['roi'],row['cell_id'])].add(row['corridor'])
differences=[]
for key,items in sorted(by_question.items()):
    endpoints={(r['lower'],r['upper']) for r in items}
    if len(endpoints)>1:
        differences.append(dict(zip(('roi','setting','cell_id','gene','budget_um2'),key),
                                family_count=len(items),
                                corridor_endpoints=json.dumps({r['corridor']:[int(r['lower']),int(r['upper'])] for r in items},sort_keys=True),
                                union_lower=min(int(r['lower']) for r in items),
                                union_upper=max(int(r['upper']) for r in items)))
write_csv(ROOT/'results/corridor_scope_disagreements.csv',differences)
# Preserve old CSVs as source records; scoped replacements retain all endpoints.
for name in ('expression_intervals','expression_matrix_intervals'):
    source=read_csv(ROOT/f'results/{name}.csv')
    for row in source:
        scope=scopes[row['corridor'],row['setting']]
        row.update(family_id=scope['family_id'],reference_id=ref,
                   scope='single_corridor_conditional',exterior_policy='all_observations_outside_corridor_fixed',
                   whole_cell_boundary_uncertainty='False',simultaneous_corridor_deformation='False')
    write_csv(ROOT/f'results/{name}_scoped.csv',source)
# A union index explicitly includes every model for each cell/setting/budget.
# The sparse source's documented absent-entry rule provides exact zeros only for
# the listed observed genes and declared model/cell/budget. Never infer new genes.
matrix=read_csv(ROOT/'results/expression_matrix_intervals.csv')
lookup={tuple(r[k] for k in ('corridor','setting','cell_id','gene','budget_um2')):r for r in matrix}
models=defaultdict(list)
for m in index['models']:
    for cell in m['cells']:
        for budget in m['budgets_um2']:
            models[(m['roi'],m['setting'],str(cell),budget)].append(m['corridor'])
observed=json.loads((ROOT/'data/observed_crop_counts.json').read_text())
unions=[]
for (roi,setting,cell,budget),corridors in sorted(models.items()):
    for gene in index['observed_gene_universe']:
        source_count=observed['rois'][roi]['counts'].get(cell,{}).get(gene,0)
        records=[]
        for c in sorted(corridors):
            r=lookup.get((c,setting,cell,gene,budget))
            if r is None:
                # The input declares a missing row to be exact zero, not fixed source.
                if source_count!=0:raise AssertionError('positive source missing from sparse matrix')
                lower=upper=0;status='exact'
            else:
                assert int(r['source_count'])==source_count
                lower,upper,status=int(r['lower']),int(r['upper']),r['status']
            records.append({'scope':scopes[c,setting],
                            'query':{'roi':roi,'cell_id':cell,'gene':gene,'budget':budget,'budget_metric':'summed_interface_deformation_um2'},
                            'source_count':source_count,'status':status,'lower':lower,'upper':upper})
        union=union_intervals(records)
        # Retain a sparse union matrix: the full domain and zero rule are explicit below.
        if source_count or union['lower'] or union['upper']:
            reverse={scopes[c,setting]['family_id']:c for c in corridors}
            unions.append({'roi':roi,'setting':setting,'cell_id':cell,'gene':gene,'budget_um2':budget,
                           'source_count':source_count,'lower':union['lower'],'upper':union['upper'],'status':union['status'],
                           'family_id':union['scope']['family_id'],'member_corridors':';'.join(sorted(corridors)),
                           'lower_attaining_corridor':reverse[union['lower_attaining_family_id']],
                           'upper_attaining_corridor':reverse[union['upper_attaining_family_id']],
                           'scope':'union_of_single_corridor_families','whole_cell_boundary_uncertainty':'False',
                           'simultaneous_corridor_deformation':'False'})
write_csv(ROOT/'results/corridor_union_intervals.csv',unions)
summary={'input_rows':len(rows),'distinct_cell_gene_setting_budget_keys':len(by_question),
         'keys_with_different_conditional_intervals':len(differences),
         'differences_by_setting':dict(Counter(r['setting'] for r in differences)),
         'roi_qualified_cells':len(cells),'cells_in_multiple_corridors':sum(len(cs)>1 for cs in cells.values()),
         'maximum_corridors_per_cell':max(map(len,cells.values())),
         'source_file_hashes':source_hashes,'reference_id':ref,
         'original_endpoint_values_unchanged':True,'union_sparse_rows':len(unions),
         'union_full_domain_rows':len(models)*len(index['observed_gene_universe']),
         'interpretation':'Different conditional families can have different exact intervals. Union chooses one corridor family at a time; it does not permit simultaneous corridor edits or unrestricted whole-cell boundary changes.',
         'study_rerun':False}
(ROOT/'results/corridor_scope_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
(ROOT/'results/corridor_union_index.json').write_text(json.dumps({
    'reference_id':ref,'scope':'union_of_single_corridor_families',
    'whole_cell_boundary_uncertainty':False,'simultaneous_corridor_deformation':False,
    'observed_gene_universe':index['observed_gene_universe'],
    'models':[{'roi':k[0],'setting':k[1],'cell_id':k[2],'budget_um2':k[3],'member_corridors':sorted(v)} for k,v in sorted(models.items())],
    'absent_entry_rule':'Only for listed models and observed genes, absent union CSV entries have exact count zero. Unlisted genes and unmeasured transcripts are unresolved.',
    'endpoint_witness_lookup':'The attaining corridor identifies an original local count frontier/witness. If its gene is absent from that corridor, use the unchanged source geometry.',
    'exterior_policy':'choose one corridor family; fix its exterior'},indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='source_file_hashes'},indent=2))
