"""Recheck the published evidence and reproduce the reported core comparisons."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def replay(analysis):
    for donor in ['IDC', 'ILC', 'LUNG']:
        for region in sorted((analysis / donor).iterdir()):
            if not region.is_dir():
                continue
            env = {**os.environ, 'GEOCERT_DATA_ROOT': str(region.resolve())}
            if donor != 'LUNG':
                for script in ['verify_atlas', 'check_global_counts']:
                    subprocess.run([sys.executable, str(ROOT / 'analysis' / (script + '.py'))], env=env, check=True)
            comparison = region / 'results/proseg_comparison.json'
            expected = json.loads(comparison.read_text())
            subprocess.run([sys.executable, str(HERE / 'evaluate_proseg.py'), '--saved-core', str(region)], check=True)
            assert json.loads(comparison.read_text()) == expected, region.name
    destination = analysis / 'replayed_summary'
    subprocess.run([sys.executable, str(HERE / 'summarize.py'), '--analysis', str(analysis), '--output', str(destination)], check=True)
    for name in ['donor_summary', 'region_summary', 'budget_summary']:
        assert json.loads((destination / (name + '.json')).read_text()) == json.loads((HERE / 'results' / (name + '.json')).read_text()), name
    print('Prepared geometry, exact spectra, witnesses, donor summaries and all saved core comparisons reproduced.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', type=Path, required=True)
    replay(parser.parse_args().analysis)
