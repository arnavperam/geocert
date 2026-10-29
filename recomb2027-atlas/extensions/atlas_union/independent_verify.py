"""Independent verification of a finite union of recorded exact atlas spectra.

This module imports no production atlas/union optimizer. It rebuilds question
censuses and min-over-family costs, checks common observation bindings, and can
compare a production output with exact rational arithmetic. Component geometry
certificates remain explicit shared inputs; selected boundary witnesses receive
an additional direct piecewise-linear check below.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time

SETTINGS = ('baseline', '0.5', '1', '2', '3')
BUDGETS = tuple(map(F, ('0','.005','.01','.025','.05','.1','.15','.2','.25','.375','.5','.75','1','2','4','8','16')))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def locate(path):
    if path.is_file():
        return path
    zipped = Path(str(path) + '.gz')
    require(zipped.is_file(), 'missing input: ' + str(path))
    return zipped


def read_json(path):
    path = locate(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
        return json.load(stream)


def csv_rows(path):
    path = locate(path)
    with (gzip.open(path, 'rt', newline='') if path.suffix == '.gz' else path.open(newline='')) as stream:
        yield from csv.DictReader(stream)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def root_for(data, dataset, setting):
    if setting == 'baseline':
        return data if dataset == 'lung' else data / 'external'
    return data / 'reach' / dataset / ('r' + setting.replace('.', 'p'))


def exact_key(row):
    return row['roi'], int(row['cell']), row['gene']


def source_count(atlases, key):
    roi, cell, gene = key
    return atlases[roi]['source_counts'].get(str(cell), {}).get(gene, 0)


def load_components(data, dataset):
    loaded = {}
    hashes_checked = 0
    geometry_ids = []
    common = None
    for setting in SETTINGS:
        root = root_for(data, dataset, setting)
        files = sorted({Path(str(p)[:-3]) if p.suffix == '.gz' else p for p in (root/'data').glob('atlas_*.json*')})
        atlases = {}
        audit = read_json(root/'results/atlas_geometric_validation.json')
        require(audit.get('failures') == 0, 'component geometry audit has failures')
        audit_ids = [(r['roi'], r['atlas_id']) for r in audit.get('coverage', [])]
        for file in files:
            a = read_json(file)
            roi = a['roi']
            require(roi not in atlases, 'duplicate atlas ROI')
            actual = hashlib.sha256(json.dumps({k:v for k,v in a.items() if k != 'atlas_id'}, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            require(actual == a['atlas_id'], 'atlas content identity mismatch')
            require(audit_ids.count((roi, actual)) == 1, 'component audit identity missing/duplicated')
            geometry_ids.append({'dataset':dataset, 'setting':setting, 'roi':roi, 'atlas_id':actual})
            require(len(set(a['complete_cells'])) == len(a['complete_cells']), 'duplicate eligible cells')
            records = {r['id']:r for r in a['records']}
            require(len(records) == len(a['records']), 'duplicate observed transcript')
            count = defaultdict(Counter)
            for rec in records.values():
                count[str(rec['cell'])][rec['gene']] += 1
            expected_counts = {cell:dict(genes) for cell,genes in count.items()}
            supplied_counts = {cell:{g:n for g,n in genes.items() if n} for cell,genes in a['source_counts'].items() if any(genes.values())}
            require(expected_counts == supplied_counts, 'source counts do not reconstruct from observations')
            if common is None:
                source_dir = root_for(data, dataset, 'baseline')/'data/source_crops'/roi
                for name, value in a['source_hashes'].items():
                    require(digest(source_dir/name) == value, 'source file digest mismatch: '+name)
                    hashes_checked += 1
            else:
                c = common[roi]
                for field in ('source_archive_sha256','pixel_size_um','bbox_pixels','complete_cells','source_counts','source_hashes'):
                    require(a[field] == c[field], 'family source mismatch: '+field)
                require(records == {r['id']:r for r in c['records']}, 'family retained observation universe differs')
            seen_points = set()
            for patch in a['patches']:
                require(F(patch['area_scale_um2']) == F(a['pixel_size_um'])**2/2, 'patch physical area units differ')
                for point in patch['points']:
                    ident = point['transcript_id']
                    require(ident in records and ident not in seen_points, 'invalid or repeated editable transcript')
                    seen_points.add(ident)
                    rec = records[ident]
                    require(point['gene'] == rec['gene'] and patch['cell_labels'][point['source']] == rec['cell'], 'editable point source differs')
            atlases[roi] = a
        if common is None:
            common = atlases
        else:
            require(set(atlases) == set(common), 'family ROI set differs')
        rows = read_json(root/'results/cell_spectra.json')
        index = {exact_key(r):r for r in rows}
        require(len(index) == len(rows), 'duplicate component spectrum key')
        census = set()
        for roi,a in atlases.items():
            cells = set(a['complete_cells'])
            for cell,genes in a['source_counts'].items():
                if int(cell) in cells:
                    census.update((roi,int(cell),g) for g,n in genes.items() if n)
            for patch in a['patches']:
                for cell in patch['cell_labels']:
                    if cell in cells:
                        census.update((roi,cell,p['gene']) for p in patch['points'])
        require(census == set(index), 'component missing or extraneous positive/editable question')
        for key,row in index.items():
            a = atlases[key[0]]
            require(row['atlas_id'] == a['atlas_id'], 'spectrum atlas binding mismatch')
            require(row['source_count'] == source_count(atlases, key), 'spectrum fixed source count mismatch')
            costs = {int(q):F(v) for q,v in row['minima'].items()}
            require(costs and all(q >= 0 and v >= 0 for q,v in costs.items()), 'invalid count spectrum')
            require(costs.get(row['source_count']) == 0, 'source state absent or charged')
            require(set(row['minima']) == set(row['witness_traces']), 'missing component trace')
        loaded[setting] = {'root':root, 'atlases':atlases, 'spectra':index}
    return loaded, hashes_checked, geometry_ids


def union_truth(loaded):
    keys = set().union(*(set(v['spectra']) for v in loaded.values()))
    costs = {}
    winners = {}
    for key in sorted(keys):
        n = source_count(loaded['baseline']['atlases'], key)
        values, choices = {}, {}
        for setting in SETTINGS:
            row = loaded[setting]['spectra'].get(key)
            spec = row['minima'] if row else {str(n):'0'}
            for q,value in spec.items():
                q, value = int(q), F(value)
                if q not in values or value < values[q]:
                    values[q], choices[q] = value, {setting}
                elif value == values[q]:
                    choices[q].add(setting)
        costs[key], winners[key] = values, choices
    return costs, winners


def verify_boundary(patch, witness):
    require(len(witness['boundaries']) == 1, 'expected scalar boundary')
    pieces = [tuple(map(F,p)) for p in witness['boundaries'][0]]
    require(pieces, 'empty boundary witness')
    lo,hi = map(F,patch['domain_pixels'])
    require(pieces[0][0] == lo and pieces[-1][1] == hi, 'boundary domain mismatch')
    require(pieces[0][2] == pieces[-1][3] == 0, 'endpoint seam moved')
    band, slope = F(patch['band_halfwidth_pixels']), F(patch['slope'])
    require(0 < band < F(patch['halfwidth_pixels']), 'displacement does not stay inside fixed rails')
    integral = F(0)
    for i,(a,b,u,v) in enumerate(pieces):
        require(a < b and abs(v-u) <= slope*(b-a), 'invalid segment/slope')
        require(max(abs(u),abs(v)) <= band, 'boundary outside band')
        if i:
            require(pieces[i-1][1] == a and pieces[i-1][3] == u, 'boundary discontinuity')
        integral += (b-a)*(abs(u)+abs(v))/2 if u*v >= 0 else (b-a)*(u*u+v*v)/(2*(abs(u)+abs(v)))
    def value(t):
        for a,b,u,v in pieces:
            if a <= t <= b:
                return u+(v-u)*(t-a)/(b-a)
        raise ValueError('point outside witness domain')
    nuclear_sides = list(dict.fromkeys(int(row[0]) for row in patch['nuclear_constraints']))
    require(len(witness['owners']) == len(patch['points'])+len(nuclear_sides), 'owner vector length mismatch')
    require(witness['owners'][len(patch['points']):] == nuclear_sides, 'structural group owner mismatch')
    for point,owner in zip(patch['points'],witness['owners']):
        require(owner in (0,1), 'nonbinary point owner')
        h,x,rho = value(F(point['t'])),F(point['x']),F(patch['clearance_pixels'])
        require(h >= x+rho if owner == 0 else h <= x-rho, 'owner clearance violated')
    for side,bound,start,end in patch['nuclear_constraints']:
        start,end,bound = F(start),F(end),F(bound)
        heights = [value(start),value(end)] + [u for a,b,u,v in pieces if start <= a <= end]
        require(min(heights) >= bound if side == 0 else max(heights) <= bound, 'nuclear ownership violated')
    physical = integral*F(patch['area_scale_um2'])
    require(physical == F(witness['physical_area_um2']), 'physical witness integral mismatch')
    return physical


def summaries(dataset, loaded, costs):
    baseline = set(loaded['baseline']['spectra'])
    out = []
    for population,keys in [('baseline_common',baseline),('union_all',set(costs))]:
        previous = {}
        for budget in BUDGETS:
            intervals = {key:(min(q for q,v in costs[key].items() if v <= budget),max(q for q,v in costs[key].items() if v <= budget)) for key in keys}
            constant = sum(l == u for l,u in intervals.values())
            out.append({'dataset':dataset,'population':population,'budget_um2':str(budget),'questions':len(keys),'singleton':constant,'non_singleton':len(keys)-constant,'summed_interval_width':sum(u-l for l,u in intervals.values()),'endpoint_changes_from_previous_budget':sum(intervals[k] != previous[k] for k in keys) if previous else 0})
            previous = intervals
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--union-output',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    start = time.perf_counter()
    output = {'status':'independent_truth_computed','geometry_exactness_scope':'Component geometry audits are shared inputs; atlas/source hashes and selected continuous witnesses are independently checked.','datasets':[], 'budget_summary':[],'reference_summary':[]}
    all_costs = {}
    for dataset in ('lung','breast'):
        loaded, hashes, ids = load_components(args.data_root, dataset)
        costs,winners = union_truth(loaded)
        all_costs[dataset] = costs
        output['budget_summary'].extend(summaries(dataset, loaded, costs))
        output['datasets'].append({'dataset':dataset,'union_questions':len(costs),'baseline_questions':len(loaded['baseline']['spectra']),'minimum_area_states':sum(map(len,costs.values())),'source_files_hashed':hashes,'member_atlas_ids':ids})
        print(dataset,output['datasets'][-1]['union_questions'],'questions',flush=True)
        # Production comparisons are deliberately independent of the producing code.
        if args.union_output:
            verify_production(args.union_output,dataset,loaded,costs,winners,output)
            verify_intervals(args.data_root,args.union_output,dataset,loaded,costs,output)
    grouped = defaultdict(lambda:defaultdict(dict))
    counts = {}
    for row in csv_rows(args.data_root/'results/reach/reference_counts.csv.gz'):
        if row['changed'] != 'True':
            continue
        group = row['dataset'],row['reference']
        key = exact_key(row)
        metadata = int(row['source_count']),int(row['reference_count'])
        require((group,key) not in counts or counts[group,key] == metadata,'reference key count universe differs')
        counts[group,key] = metadata
        require(row['setting'] not in grouped[group][key],'duplicate member reference query')
        grouped[group][key][row['setting']] = None if row['minimum_area_um2'] == 'inf' else F(row['minimum_area_um2'])
    for (dataset,reference),keys in sorted(grouped.items()):
        baseline = reached = 0
        for key,members in keys.items():
            require(set(members) == set(SETTINGS),'reference comparison member keys differ')
            candidate = min((v for v in members.values() if v is not None),default=None)
            source,target = counts[(dataset,reference),key]
            expected = all_costs[dataset].get(key,{source:F(0)}).get(target)
            require(candidate == expected,'reference query is not the exact discrete union count state')
            baseline += members['baseline'] is not None
            reached += candidate is not None
        output['reference_summary'].append({'dataset':dataset,'reference':reference,'changed_count_keys':len(keys),'baseline_attainable':baseline,'union_attainable':reached})
    if args.union_output:
        verify_summaries(args.union_output,output)
    output['seconds'] = time.perf_counter()-start
    if args.union_output:
        output['status'] = 'passed'
    (args.output_dir/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':output['status'],'reference_summary':output['reference_summary'],'seconds':output['seconds']},indent=2))


def verify_production(folder,dataset,loaded,costs,winners,report):
    # Filled against the producer's published file schema; no producer import.
    identity_fields = ('source_archive_sha256','pixel_size_um','bbox_pixels','complete_cells','source_counts','records','source_hashes')
    common = {roi:{field:a[field] for field in identity_fields} for roi,a in loaded['baseline']['atlases'].items()}
    canonical = lambda obj: json.dumps(obj,sort_keys=True,separators=(',', ':')).encode()
    common_digest = hashlib.sha256(canonical(common)).hexdigest()
    member_ids = {setting:{roi:a['atlas_id'] for roi,a in member['atlases'].items()} for setting,member in loaded.items()}
    model_id = hashlib.sha256(canonical({'schema':'geocert-atlas-union-1','dataset':dataset,'reference_digest':common_digest,'members':member_ids})).hexdigest()
    rows = read_json(folder/'union_spectra.json')
    rows = [r for r in rows if r['dataset'] == dataset]
    actual = {exact_key(r):r for r in rows}
    require(len(actual) == len(rows) and set(actual) == set(costs),'union question census mismatch')
    local = {}
    patch_index = {}
    checked = set()
    traces = 0
    implicit = 0
    for key,row in actual.items():
        require(row['model_id'] == model_id,'union model identity mismatch')
        require(row['baseline_question'] == (key in loaded['baseline']['spectra']),'baseline question flag differs')
        require({int(q):F(v) for q,v in row['minima'].items()} == costs[key], 'union minimum cost differs')
        require(row['source_count'] == source_count(loaded['baseline']['atlases'],key),'union source count differs')
        require(set(row['witnesses']) == set(row['minima']),'union count missing witness')
        for q,w in row['witnesses'].items():
            count = int(q)
            setting = w['setting']
            require(setting in winners[key][count],'union witness member is not attaining minimum')
            member = loaded[setting]
            atlas = member['atlases'][key[0]]
            require(w['atlas_id'] == atlas['atlas_id'],'union witness atlas identity differs')
            expected_path = ('results/cell_spectra.json' if dataset == 'lung' else 'external/results/cell_spectra.json') if setting == 'baseline' else f"reach/{dataset}/r{setting.replace('.', 'p')}/results/cell_spectra.json"
            require(w['source_spectrum_path'] == expected_path,'union witness source path differs')
            original = member['spectra'].get(key)
            if original is None:
                require(w['implicit_fixed'] and w['trace'] == {} and count == row['source_count'] and costs[key][count] == 0,'invalid implicit fixed state')
                implicit += 1
                continue
            require(not w['implicit_fixed'] and w['trace'] == original['witness_traces'][q],'union witness trace differs from named member')
            if setting not in local:
                lr = read_json(member['root']/'results/patch_spectra.json')
                local[setting] = {(r['roi'],r['patch'],r['gene']):r for r in lr}
                require(len(local[setting]) == len(lr),'duplicate local spectrum')
                patch_index[setting] = {(roi,p['id']):p for roi,a in member['atlases'].items() for p in a['patches']}
            require(count == row['source_count'] + sum(w['trace'].values()),'union count trace sum differs')
            area = F(0)
            for pid,delta in w['trace'].items():
                patch = patch_index[setting][key[0],pid]
                require(key[1] in patch['cell_labels'],'witness patch not incident to requested cell')
                loc = local[setting][key[0],pid,key[2]]
                sign = 1 if patch['cell_labels'][0] == key[1] else -1
                target = loc['source_left'] + sign*delta
                state = loc['witnesses'][str(target)]
                require(sum(p['gene'] == key[2] and owner == 0 for p,owner in zip(patch['points'],state['owners'])) == target,'local owner count differs')
                value = F(state['physical_area_um2'])
                require(value == F(loc['minima'][str(target)]),'local trace minimum differs')
                check_id = setting,key[0],pid,key[2],target
                if check_id not in checked:
                    verify_boundary(patch,state)
                    checked.add(check_id)
                area += value
            require(area == costs[key][count],'union trace area differs')
            traces += 1
    report.setdefault('production_checks',[]).append({'dataset':dataset,'model_id':model_id,'spectra_checked':len(actual),'traces_checked':traces,'implicit_fixed_states_checked':implicit,'distinct_continuous_local_witnesses_checked':len(checked)})


def verify_intervals(data,folder,dataset,loaded,costs,report):
    """Check every budget row and the outer bound for the same union family."""
    outer = {}
    for key in costs:
        n = source_count(loaded['baseline']['atlases'],key)
        outer[key] = [[n,n] for _ in BUDGETS]
    budget_index = {d:i for i,d in enumerate(BUDGETS)}
    component_rows = 0
    for setting in SETTINGS:
        member = loaded[setting]
        path = data/'results/reach'/f'baseline_{dataset}_fine_intervals.csv' if setting == 'baseline' else member['root']/'results/cell_expression_intervals.csv'
        seen = set()
        for row in csv_rows(path):
            key,delta = exact_key(row),F(row['budget_um2'])
            require(key in member['spectra'] and delta in budget_index,'unexpected component interval query')
            ident = key,delta
            require(ident not in seen,'duplicate component interval row')
            seen.add(ident)
            require(row['atlas_id'] == member['atlases'][key[0]]['atlas_id'],'component interval atlas mismatch')
            spec = member['spectra'][key]['minima']
            feasible = [int(q) for q,v in spec.items() if F(v) <= delta]
            exact = min(feasible),max(feasible)
            require(row['status'] == 'exact' and (int(row['lower']),int(row['upper'])) == exact,'component interval endpoint differs from spectrum')
            lower,upper = int(row['independent_lower']),int(row['independent_upper'])
            require(lower <= exact[0] <= exact[1] <= upper,'component marginal is not an outer bound')
            bounds = outer[key][budget_index[delta]]
            bounds[0],bounds[1] = min(bounds[0],lower),max(bounds[1],upper)
            component_rows += 1
        require(len(seen) == len(member['spectra'])*len(BUDGETS),'incomplete component interval grid')
    expected_model_id = next(r['model_id'] for r in report['production_checks'] if r['dataset'] == dataset)
    seen = set()
    tightened = defaultdict(set)
    eligible = defaultdict(set)
    for row in csv_rows(folder/'intervals.csv'):
        if row['dataset'] != dataset:
            continue
        key,delta = exact_key(row),F(row['budget_um2'])
        require(row['model_id'] == expected_model_id,'interval model identity differs')
        require(row['baseline_question'] == str(key in loaded['baseline']['spectra']),'interval baseline question flag differs')
        require(key in costs and delta in budget_index,'unexpected union interval key')
        require((key,delta) not in seen,'duplicate union interval row')
        seen.add((key,delta))
        feasible = [q for q,v in costs[key].items() if v <= delta]
        exact = min(feasible),max(feasible)
        require(row['status'] == 'exact' and (int(row['lower']),int(row['upper'])) == exact,'union interval endpoints differ')
        require(int(row['source_count']) == source_count(loaded['baseline']['atlases'],key),'union interval source count differs')
        require(F(row['lower_cost_um2']) == costs[key][exact[0]] and F(row['upper_cost_um2']) == costs[key][exact[1]],'union endpoint costs differ')
        for side,q in [('lower',exact[0]),('upper',exact[1])]:
            setting = row[side+'_setting']
            require(setting in SETTINGS,'interval witness setting is unknown')
            member = loaded[setting]
            require(row[side+'_atlas_id'] == member['atlases'][key[0]]['atlas_id'],'interval witness atlas differs')
            recorded = member['spectra'].get(key)
            candidate = {int(n):F(v) for n,v in recorded['minima'].items()} if recorded else {source_count(member['atlases'],key):F(0)}
            require(candidate.get(q) == costs[key][q],'interval witness member does not attain endpoint cost')
        expected_outer = outer[key][budget_index[delta]]
        observed_outer = [int(row['independent_lower']),int(row['independent_upper'])]
        require(observed_outer == expected_outer,'union marginal bound differs from union of member bounds')
        if expected_outer[1]-expected_outer[0] >= 2:
            eligible['union_all'].add(key)
            if key in loaded['baseline']['spectra']:
                eligible['baseline_common'].add(key)
        if list(exact) != expected_outer:
            tightened['union_all'].add(key)
            if key in loaded['baseline']['spectra']:
                tightened['baseline_common'].add(key)
        # The baseline family is included even when its key is an implicit zero.
        row0 = loaded['baseline']['spectra'].get(key)
        spec0 = {int(q):F(v) for q,v in row0['minima'].items()} if row0 else {source_count(loaded['baseline']['atlases'],key):F(0)}
        f0 = [q for q,v in spec0.items() if v <= delta]
        require(exact[0] <= min(f0) <= max(f0) <= exact[1],'union violates baseline interval inclusion')
        require((int(row['baseline_lower']),int(row['baseline_upper'])) == (min(f0),max(f0)),'baseline endpoint columns differ')
    require(len(seen) == len(costs)*len(BUDGETS),'union interval grid incomplete')
    report.setdefault('interval_checks',[]).append({'dataset':dataset,'component_rows_recomputed':component_rows,'union_rows_recomputed':len(seen),'eligible_marginal_width_ge2':{k:len(v) for k,v in eligible.items()},'tightened_questions':{k:len(v) for k,v in tightened.items()}})


def verify_summaries(folder,report):
    rows = list(csv_rows(folder/'budget_summary.csv'))
    observed = {(r['dataset'],r['population'],F(r['budget_um2'])):r for r in rows}
    require(len(observed) == len(rows) == 68,'budget summary census mismatch')
    population_names = {'baseline_common':'baseline_questions','union_all':'union_questions'}
    for expected in report['budget_summary']:
        key = expected['dataset'],population_names[expected['population']],F(expected['budget_um2'])
        row = observed[key]
        for source,target in [('questions','questions'),('singleton','invariant'),('non_singleton','non_singleton'),('summed_interval_width','summed_width'),('endpoint_changes_from_previous_budget','endpoint_changes_from_previous_budget')]:
            require(int(row[target]) == expected[source],'budget summary value mismatch: '+target)
        require(abs(float(row['invariant_fraction'])-expected['singleton']/expected['questions']) < 1e-14,'budget summary fraction mismatch')
    reference = read_json(folder/'reference_summary.json')
    unbounded = {(r['dataset'],r['reference']):r for r in reference if r['budget_um2'] == 'unbounded'}
    require(len(unbounded) == 3,'reference summary comparison census differs')
    for expected in report['reference_summary']:
        row = unbounded[expected['dataset'],expected['reference']]
        for field in ('changed_count_keys','baseline_attainable','union_attainable'):
            require(row[field] == expected[field],'reference summary count mismatch')
    report['summary_checks'] = {'budget_rows':len(rows),'unbounded_reference_comparisons':len(unbounded)}


if __name__ == '__main__':
    main()
