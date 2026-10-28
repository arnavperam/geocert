"""Paired-boundary validation on declared independent units.

This accepts supplied truth; it cannot establish that labels are anatomical truth,
that units are independent, or that a recorded protocol was actually preregistered.
"""
from fractions import Fraction as F
from collections import defaultdict
from math import sqrt
from layered import Layers, Box, value
from verify import verify
from calibration import split_budget
from pin_closure import PinClosure
from budget_counts import CountTarget, enumerate_gene_states, count_frontier, exact_interval


def make_problem(record):
    m = record['model']
    allowed = {'baselines','source','slopes','minimum_gaps','maximum_gaps','clearance','bands','mode'}
    if set(m) - allowed:
        raise ValueError('Unsupported model fields; structural regions are not implemented by this importer')
    if m.get('mode', 'nominal') != 'nominal':
        raise ValueError('Paired importer currently requires nominal point sites')
    boxes = [Box(z['x'], z['x'], z['t'], z['t'], i) for i, z in enumerate(record['sites'])]
    p = Layers(m['baselines'], boxes, m['source'], m['slopes'], m['minimum_gaps'],
               m.get('maximum_gaps'), clearance=m['clearance'], bands=m['bands'], mode='nominal')
    if p.n != len(record['sites']) or any(type(o) is not int or not 0 <= o <= p.m for o in p.source):
        raise ValueError('Each site needs one source owner')
    return p


def ranking_metrics(labels, scores):
    """Exact tie-aware AUROC and threshold-grouped average precision."""
    if len(labels) != len(scores) or any(type(y) is not int or y not in (0, 1) for y in labels):
        raise ValueError('Matched binary labels and scores required')
    groups = defaultdict(lambda: [0, 0])
    for y, s in zip(labels, scores):
        groups[F(s)][y] += 1
    positives = sum(labels); negatives = len(labels) - positives
    if not positives or not negatives:
        return {'status': 'single_class', 'positives': positives, 'negatives': negatives,
                'auroc': None, 'average_precision': None}
    below = 0; concordance = F(0)
    for s, (n, p) in sorted(groups.items()):
        concordance += p * (below + F(n, 2)); below += n
    tp = fp = 0; ap = F(0)
    for s, (n, p) in sorted(groups.items(), reverse=True):
        tp += p; fp += n
        ap += F(p, positives) * F(tp, tp + fp)
    return {'status': 'computed', 'positives': positives, 'negatives': negatives,
            'auroc': concordance / (positives * negatives), 'average_precision': ap}


def wilson(successes, total):
    if not total:
        return None
    z = 1.959963984540054; phat = successes / total; den = 1 + z*z/total
    center = (phat + z*z/(2*total)) / den
    half = z * sqrt(phat*(1-phat)/total + z*z/(4*total*total)) / den
    return [0.0 if successes == 0 else max(0.0, center-half),
            1.0 if successes == total else min(1.0, center+half)]


def validate_dataset(data):
    if data.get('schema') != 'geocert.paired_boundaries.v1':
        raise ValueError('Unsupported paired data schema')
    if data.get('outside_assignments_fixed') is not True:
        raise ValueError('Current intervals require fixed outside assignments')
    if data.get('evidence_kind') not in ('synthetic_fixture', 'externally_supplied_paired_boundaries'):
        raise ValueError('Explicit evidence kind required')
    if not data.get('protocol_id') or not data.get('provenance'):
        raise ValueError('Protocol identifier and provenance required')
    units = {}; seen = set(); training = set(data.get('training_units', []))
    for r in data['records']:
        if not r['unit_id'] or not r['sample_id'] or r['split'] not in ('calibration', 'evaluation'):
            raise ValueError('Named unit and calibration/evaluation split required')
        key = (r['unit_id'], r['sample_id'])
        if key in seen: raise ValueError('Duplicate paired sample')
        seen.add(key)
        if r['unit_id'] in training: raise ValueError('Training unit leaked into calibration/evaluation')
        if r['unit_id'] in units and units[r['unit_id']] != r['split']:
            raise ValueError('Unit crosses calibration/evaluation split')
        units[r['unit_id']] = r['split']
        p = make_problem(r)
        if len(r['truth_owners']) != p.n or any(type(o) is not int or not 0 <= o <= p.m for o in r['truth_owners']):
            raise ValueError('Complete valid supplied truth owners required')
        boundaries = r['truth_boundaries']
        if len(boundaries) != p.m or any(not b for b in boundaries):
            raise ValueError('Complete supplied boundary arrays required')
        for b in boundaries:
            for seg in b:
                if len(seg) != 4: raise ValueError('Boundary segment needs four rational values')
                list(map(F, seg))
        keys = [(g['gene'], g['owner']) for g in r['targets']]
        if len(keys) != len(set(keys)): raise ValueError('Duplicate count target')
        for g in r['targets']:
            CountTarget(tuple(i for i, s in enumerate(r['sites']) if s['gene'] == g['gene']),
                        g['owner'], g['outside']).checked_weights(p)
            if not any(s['gene'] == g['gene'] for s in r['sites']):
                raise ValueError('Targets must name a represented gene')
        if not r['targets']: raise ValueError('At least one target required')
        for name, scores in r.get('external_error_scores', {}).items():
            if name in ('singleton_fragility', 'distance_to_boundary') or name not in data.get('external_method_provenance', {}):
                raise ValueError('External scores need distinct method names and provenance')
            if len(scores) != p.n: raise ValueError('External scores must cover every site')
            list(map(F, scores))
    if set(units.values()) != {'calibration', 'evaluation'}:
        raise ValueError('Nonempty disjoint calibration and evaluation units required')
    return units


def truth_score(p, record):
    """Actual supplied boundary area, not minimum area for its owner vector."""
    try:
        boundaries = [[tuple(map(F, seg)) for seg in boundary] for boundary in record['truth_boundaries']]
        return verify(p, record['truth_owners'], boundaries), None
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError) as e:
        return None, str(e)


def run_validation(data, alpha=F(1, 10), budgets=(F(0), F(1, 100), F(1, 20), F(1, 4)), max_assignments=100000):
    units = validate_dataset(data)
    budgets = tuple(sorted(set(map(F, budgets))))
    if any(b < 0 for b in budgets): raise ValueError('Nonnegative probe budgets required')
    prepared = []; scores = defaultdict(list)
    for r in data['records']:
        p = make_problem(r)
        engine = PinClosure(p, {})  # Includes the source-admissibility check.
        if engine.inconsistent: raise ValueError('Reference model is infeasible')
        area, error = truth_score(p, r)
        scores[r['unit_id']].append(area)
        prepared.append((r, p, engine, area, error))
    unit_scores = {u: None if any(x is None for x in s) else max(s) for u, s in scores.items()}
    cal = sorted(u for u in units if units[u] == 'calibration')
    test = sorted(u for u in units if units[u] == 'evaluation')
    fit = split_budget([unit_scores[u] for u in cal], alpha, cal, test)
    levels = sorted(set(budgets + (() if fit['budget'] is None else (fit['budget'],))))
    rows = []; boundary_audit = []; labels = []; rank_scores = defaultdict(list)
    ext_methods = None
    for r, p, engine, area, error in prepared:
        boundary_audit.append({'unit': r['unit_id'], 'sample': r['sample_id'], 'split': r['split'],
                               'area_score': area, 'family_error': error})
        if r['split'] != 'evaluation': continue
        methods = set(r.get('external_error_scores', {}))
        if ext_methods is None: ext_methods = methods
        elif ext_methods != methods: raise ValueError('External methods need identical evaluation coverage')
        minimum = []
        for i in range(p.n):
            costs = {}
            for o in range(p.m+1):
                a = engine.query({i:o})
                if a['feasible']: costs[o] = a['area']
            minimum.append(costs)
            other = [v for o,v in costs.items() if o != p.source[i]]
            # A bounded monotone transform also handles impossible alternatives.
            rank_scores['singleton_fragility'].append(F(0) if not other else F(1)/(1+min(other)))
            z = r['sites'][i]
            distance = min(abs(F(z['x'])-value(b,F(z['t']))) for b in p.baselines)
            rank_scores['distance_to_boundary'].append(-distance)
            labels.append(int(r['truth_owners'][i] != p.source[i]))
        for name in methods: rank_scores[name].extend(map(F, r['external_error_scores'][name]))
        cache = {}
        for t in r['targets']:
            groups = tuple(i for i,s in enumerate(r['sites']) if s['gene'] == t['gene'])
            if groups not in cache: cache[groups] = enumerate_gene_states(engine, groups, max_assignments)
            target = CountTarget(groups, t['owner'], t['outside'])
            frontier = count_frontier(p,target,cache[groups]); truth = target.count(r['truth_owners'])
            for delta in levels:
                interval = exact_interval(frontier,delta)
                lo,hi = interval['lower'],interval['upper']
                allowed = [[o for o,c in minimum[i].items() if c<=delta] for i in groups]
                slo = target.outside + sum(all(o==target.owner for o in a) for a in allowed)
                shi = target.outside + sum(target.owner in a for a in allowed)
                assert slo <= lo <= hi <= shi
                contained = area is not None and area<=delta
                covered = lo<=truth<=hi
                if contained and not covered: raise AssertionError('Valid truth boundary escaped exact count interval')
                rows.append({'unit':r['unit_id'],'sample':r['sample_id'],'gene':t['gene'],
                             'owner':t['owner'],'budget':delta,'truth_count':truth,'lower':lo,'upper':hi,
                             'covered':covered,'family_contained':contained,'singleton_lower':slo,
                             'singleton_upper':shi,'singleton_covered':slo<=truth<=shi})
    summaries = []
    for delta in levels:
        rr = [r for r in rows if r['budget']==delta]
        per_unit = {u:[r for r in rr if r['unit']==u] for u in test}
        all_covered = sum(all(r['covered'] for r in rs) for rs in per_unit.values())
        all_contained = sum(unit_scores[u] is not None and unit_scores[u]<=delta for u in test)
        macro_coverage = sum(F(sum(r['covered'] for r in rs),len(rs)) for rs in per_unit.values())/len(test)
        macro_width = sum(F(sum(r['upper']-r['lower'] for r in rs),len(rs)) for rs in per_unit.values())/len(test)
        summaries.append({'budget':delta,'evaluation_units':len(test),'query_rows':len(rr),
            'simultaneous_count_covered_units':all_covered,'simultaneous_count_coverage':F(all_covered,len(test)),
            'unit_coverage_wilson_95':wilson(all_covered,len(test)),
            'family_contained_units':all_contained,'family_coverage':F(all_contained,len(test)),
            'family_coverage_wilson_95':wilson(all_contained,len(test)),
            'macro_count_coverage':macro_coverage,'macro_interval_width':macro_width,
            'joint_strictly_tighter_rows':sum(r['lower']>r['singleton_lower'] or r['upper']<r['singleton_upper'] for r in rr)})
    return {'schema':'geocert.paired_validation.v1','protocol_id':data['protocol_id'],
            'evidence_kind':data['evidence_kind'],'calibration':fit,'calibration_units':cal,'evaluation_units':test,
            'unit_scores':unit_scores,'summaries':summaries,'intervals':rows,'boundary_audit':boundary_audit,
            'pooled_error_ranking':{name:ranking_metrics(labels,ss) for name,ss in rank_scores.items()},
            'ranking_scope':'Pooled descriptive transcript ranking; no transcript-independent confidence intervals.',
            'interval_scope':'Wilson intervals use unit-level indicators and require independent comparable units; identifiers alone do not verify this.',
            'guarantee_scope':'Split-conformal coverage requires exchangeable calibration/new units, a frozen family, complete valid paired truth and the same sampling protocol. Out-of-family truth has infinite score. Counts outside modeled sites stay fixed.'}
