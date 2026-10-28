"""Identity and safe aggregation of explicitly conditional uncertainty families.

A union permits one member family at a time. It is not simultaneous deformation
of several corridors, nor a whole-cell geometric model.
"""
from fractions import Fraction
import hashlib
import json


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def corridor_scope(chart, setting='baseline', mode='nominal', rna_radius='0', reference_id=None):
    descriptor = {'chart': chart, 'setting': setting, 'mode': mode,
                  'rna_radius': str(Fraction(rna_radius)),
                  'exterior_policy': 'all_observations_outside_corridor_fixed'}
    return {'kind': 'single_corridor_conditional', 'family_id': digest(descriptor),
            'reference_id': reference_id, 'roi': chart['roi'],
            'corridor': chart['id'], 'setting': setting,
            'exterior_policy': descriptor['exterior_policy'],
            'whole_cell_boundary_uncertainty': False,
            'simultaneous_corridor_deformation': False}


def union_intervals(records):
    """Exact endpoints over a declared union when all constituent endpoints are exact.

    Records must have scope, query (cell_id, gene, budget, budget_metric), and
    source_count. A shared nonempty reference_id binds their observations. Missing
    families are never inferred. Each member is allowed independently; overlap is
    harmless for set union but forbids interpreting this as joint deformation.
    """
    records = list(records)
    if not records:
        raise ValueError('at least one fully scoped interval is required')
    keys = ('roi', 'cell_id', 'gene', 'budget_metric')
    first = records[0]
    ref = first['scope'].get('reference_id')
    if not isinstance(ref, str) or not ref.strip():
        raise ValueError('union requires a shared explicit reference_id')
    query = first['query']
    if any(k not in query or query[k] is None or query[k] == '' for k in keys):
        raise ValueError('complete region, cell, gene and budget-metric identity is required')
    budget = Fraction(query['budget'])
    if budget < 0:
        raise ValueError('negative budget')
    families = []
    for row in records:
        scope = row['scope']
        if scope.get('kind') != 'single_corridor_conditional':
            raise ValueError('union accepts explicitly scoped single-corridor intervals only')
        if scope.get('reference_id') != ref:
            raise ValueError('reference observations differ between families')
        if any(row['query'].get(k) != query.get(k) for k in keys) or Fraction(row['query']['budget']) != budget:
            raise ValueError('union requires the same region, cell, gene, metric and budget')
        if type(row['source_count']) is not int or row['source_count'] < 0:
            raise ValueError('nonnegative integer source counts required')
        if type(row['lower']) is not int or type(row['upper']) is not int:
            raise ValueError('integer interval endpoints required')
        if row['source_count'] != first['source_count']:
            raise ValueError('source counts disagree between families')
        if row['status'] not in ('exact', 'bounded'):
            raise ValueError('union requires resolved nonempty constituent intervals')
        family = scope.get('family_id')
        if not isinstance(family, str) or not family or family in families:
            raise ValueError('unique nonempty family identifiers are required')
        if row['lower'] > row['upper']:
            raise ValueError('reversed constituent interval')
        families.append(family)
    lower_row = min(records, key=lambda r: r['lower'])
    upper_row = max(records, key=lambda r: r['upper'])
    members = sorted(families)
    result = {'scope': {'kind': 'union_of_single_corridor_families',
              'family_id': digest({'reference_id': ref, 'members': members}),
              'reference_id': ref, 'member_family_ids': members,
              'whole_cell_boundary_uncertainty': False,
              'simultaneous_corridor_deformation': False,
              'exterior_policy': 'choose_one_corridor_family; fix_its_exterior'},
              'query': dict(query), 'source_count': first['source_count'],
              'status': 'exact' if all(r['status'] == 'exact' for r in records) else 'bounded',
              'lower': lower_row['lower'], 'upper': upper_row['upper'],
              'lower_attaining_family_id': lower_row['scope']['family_id'] if lower_row['status'] == 'exact' else None,
              'upper_attaining_family_id': upper_row['scope']['family_id'] if upper_row['status'] == 'exact' else None}
    for name, row in [('minimum_witness', lower_row), ('maximum_witness', upper_row)]:
        if row['status'] == 'exact' and name in row:
            result[name] = row[name]
    return result
