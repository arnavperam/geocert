"""Diagnostic only: local Gumbel perturbations need not sample the owner Gibbs law.

Two coincident RNA sites must share an owner. Their non-source energy is 1.
This is an uncalibrated synthetic distribution, not a biological posterior.
"""
from pathlib import Path
import argparse, json, math, random, sys
from fractions import Fraction as F

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'core'))
from layered import Layers, Box
from scc_solver import prepare, solve_compiled


def main(samples, output):
    p = Layers([[(0, 0), (3, 0)]], [Box(F(1, 5), F(1, 5), 1, 1, g) for g in range(2)],
               [1, 1], F(1), [], [], clearance=F(1, 20), bands=[(F(-1, 2), F(1, 2))])
    compiled, plan, _ = prepare(p, {}, 'direct', True)
    rng = random.Random(20261028)
    changed = 0
    for _ in range(samples):
        rows = []
        for g in range(2):
            row = []
            for o in range(2):
                u = rng.randrange(1, 2 ** 53) / 2 ** 53
                gumbel = -math.log(-math.log(u))
                row.append(-F(round(gumbel * 10 ** 6), 10 ** 6))
            rows.append(row)
        a = solve_compiled(compiled, 'direct', rows, F(16), plan)
        assert tuple(a['owners']) in [(0, 0), (1, 1)]
        changed += a['owners'][0] == 0
    freq = changed / samples
    z = 1.959963984540054
    den = 1 + z * z / samples
    center = (freq + z * z / (2 * samples)) / den
    half = z * math.sqrt(freq * (1 - freq) / samples + z * z / (4 * samples ** 2)) / den
    exact_gibbs = 1 / (1 + math.exp(1))
    result = {'scope': 'Synthetic low-order perturb-and-MAP diagnostic; no posterior calibration claim.',
              'samples': samples, 'seed': 20261028, 'quantization': '1/1000000',
              'feasible_owner_vectors': [[0, 0], [1, 1]], 'unperturbed_energy': [1, 0],
              'gibbs_non_source_probability': exact_gibbs, 'perturb_map_non_source_frequency': freq,
              'wilson_95_interval_for_sampler_frequency': [center - half, center + half],
              'empirical_total_variation_from_gibbs': abs(freq - exact_gibbs),
              'all_sampled_configurations_geometrically_verified': True}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--samples', type=int, default=2000)
    ap.add_argument('--output', type=Path, default=ROOT / 'run/perturb_map.json')
    args = ap.parse_args()
    if args.samples <= 0: raise SystemExit('Positive sample count required')
    main(args.samples, args.output)
