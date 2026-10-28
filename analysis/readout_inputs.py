"""Build the outcome-independent cell/gene query universe from published inputs."""
from collections import Counter
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def build_inputs() -> dict:
    charts = json.loads((ROOT / 'data/charts.json').read_text())
    source = json.loads((ROOT / 'analysis/source_counts.json').read_text())
    cells, readouts = [], []
    for chart in charts:
        roi = chart['roi']
        complete = set(source['complete_cell_labels'][roi])
        genes = sorted({point['gene'] for point in chart['points']})
        local = Counter((point['source'], point['gene']) for point in chart['points'])
        for index, label in enumerate(chart['cell_labels']):
            cells.append({'chart': chart['id'], 'roi': roi, 'cell_index': index,
                          'mask_label': label, 'complete_cell': label in complete})
            if label not in complete:
                continue
            for gene in genes:
                total = source['counts'][roi][str(label)][gene]
                inside = local[index, gene]
                if total < inside:
                    raise ValueError(f'Invalid exterior accounting: {chart["id"]}, {label}, {gene}')
                readouts.append({'chart': chart['id'], 'roi': roi, 'cell_index': index,
                                 'cell_label': label, 'gene': gene, 'source_count': total,
                                 'fixed_count': total - inside, 'source_local_count': inside})
    return {'protocol_sha256': hashlib.sha256((ROOT / 'analysis/PROTOCOL.md').read_bytes()).hexdigest(),
            'source': source['source'], 'cells': cells, 'readouts': readouts}


if __name__ == '__main__':
    result = build_inputs()
    print(json.dumps({'cell_records': len(result['cells']), 'eligible_cell_records': sum(c['complete_cell'] for c in result['cells']),
                      'readouts_per_setting': len(result['readouts'])}, indent=2))
