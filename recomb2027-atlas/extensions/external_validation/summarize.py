"""Recompute donor summaries from interval tables, with failures kept explicit."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import argparse
import csv
import json


def summarize(analysis, output):
    output.mkdir(parents=True, exist_ok=True)
    donor_rows = []
    region_rows = []
    budget_rows = []
    for donor in ['IDC', 'ILC']:
        total = Counter()
        per_budget = {}
        unique_cells = set()
        for region in sorted((analysis / donor).glob(donor.lower() + '_*')):
            results = region / 'results'
            atlas = json.loads(next((region / 'data').glob('atlas_*.json')).read_text())
            audits = {name: json.loads((results / name).read_text()) for name in [
                'global_count_validation.json', 'source_input_audit.json', 'distance_summary.json']}
            assert all(a['failures'] == 0 for a in audits.values())
            assert all(s['returncode'] == 0 for s in json.loads((region / 'run_status.json').read_text()))
            assert not unique_cells.intersection(atlas['complete_cells']), 'Overlapping complete-cell sets'
            unique_cells.update(atlas['complete_cells'])
            rows = list(csv.DictReader((results / 'cell_expression_intervals.csv').open()))
            assert all(r['status'] == 'exact' for r in rows)
            fields = set(rows[0])
            # The exporter labels the matched independent interval explicitly.
            lower = 'independent_lower'
            upper = 'independent_upper'
            assert {lower, upper} <= fields, fields
            eligible = set()
            tightened = set()
            for row in rows:
                key = (int(row['cell']), row['gene'])
                budget = str(Fraction(row['budget_um2']))
                group = per_budget.setdefault(budget, Counter())
                group['questions'] += 1
                group['non_singleton'] += int(row['lower']) != int(row['upper'])
                different = (int(row['lower']), int(row['upper'])) != (int(row[lower]), int(row[upper]))
                group['tightened'] += different
                if int(row[upper]) - int(row[lower]) >= 2:
                    eligible.add(key)
                if different:
                    tightened.add(key)
            assert tightened <= eligible
            at4 = [r for r in rows if Fraction(r['budget_um2']) == 4]
            extraction = json.loads((results / 'atlas_extraction.json').read_text())
            coverage = extraction['rois'][0]
            additive = json.loads((results / 'additive_summary.json').read_text())
            record = dict(donor=donor, region=region.name,
                cells=len(atlas['complete_cells']), RNAs=len(atlas['records']), patches=len(atlas['patches']),
                queries=len(at4), interval_rows=len(rows), invariant_at4=sum(r['lower'] == r['upper'] for r in at4),
                eligible=len(eligible), tightened=len(tightened), additive_missed=additive['affected_cell_gene_questions'],
                covered_edges=coverage['selected_source_edges'], eligible_edges=coverage['eligible_shared_edges'])
            region_rows.append(record)
            total.update({k: v for k, v in record.items() if isinstance(v, int)})
        donor_rows.append({'donor': donor, **total, 'edge_coverage': total['covered_edges'] / total['eligible_edges']})
        budget_rows.extend({'donor': donor, 'budget_um2': b, **per_budget[b]} for b in sorted(per_budget, key=Fraction))
    for name, records in [('donor_summary', donor_rows), ('region_summary', region_rows), ('budget_summary', budget_rows)]:
        (output / (name + '.json')).write_text(json.dumps(records, indent=2) + '\n')
        with (output / (name + '.csv')).open('w', newline='') as target:
            writer = csv.DictWriter(target, fieldnames=records[0].keys())
            writer.writeheader()
            writer.writerows(records)
    print(json.dumps(donor_rows, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summarize(args.analysis, args.output)
