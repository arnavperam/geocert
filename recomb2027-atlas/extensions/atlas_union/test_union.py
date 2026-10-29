"""Focused tests for union semantics and admission guards; no external packages."""
import json
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

import union


class UnionTests(unittest.TestCase):
    def test_union_keeps_unsupported_count_and_chooses_whole_member(self):
        # A reference count inside an endpoint interval need not be attainable.
        a = {0: F(0), 2: F(4), 3: F(9)}
        b = {0: F(0), 2: F(6), 5: F(8)}
        costs, selected = union.union_minima([('a', a), ('b', b)])
        self.assertEqual(costs, {0: F(0), 2: F(4), 3: F(9), 5: F(8)})
        self.assertNotIn(1, costs)
        self.assertNotIn(4, costs)
        self.assertEqual(selected, {0: 'a', 2: 'a', 3: 'a', 5: 'b'})
        self.assertEqual(sorted(q for q, d in costs.items() if d <= 8), [0, 2, 5])

    def test_ties_are_deterministic_and_exact_rationals_are_preserved(self):
        small = F(1, 10**24)
        costs, selected = union.union_minima([('baseline', {1: F(1, 3)}),
                                             ('wide', {1: F(1, 3) - small})])
        self.assertEqual(selected[1], 'wide')
        self.assertEqual(costs[1], F(1, 3) - small)
        _, tied = union.union_minima([('baseline', {1: F(1, 3)}), ('wide', {1: F(1, 3)})])
        self.assertEqual(tied[1], 'baseline')

    def test_expected_census_includes_zero_source_editable_destinations(self):
        atlas = {'complete_cells': [1, 2], 'source_counts': {'1': {'G': 2}, '9': {'X': 1}},
                 'patches': [{'id': 'p', 'cell_labels': [1, 2], 'points': [{'gene': 'G'}]}]}
        questions, parts = union.expected_questions({'roi': atlas})
        self.assertEqual(questions, {('roi', 1, 'G'), ('roi', 2, 'G')})
        self.assertEqual(parts[('roi', 2, 'G')], {'p'})

    def test_complete_reference_identity_checks_observations_not_just_metadata(self):
        a = {'roi': {'source_archive_sha256': 'x', 'pixel_size_um': '1', 'bbox_pixels': [0, 0, 2, 2],
                     'complete_cells': [1], 'source_counts': {'1': {'G': 1}},
                     'records': [{'id': 'm', 'cell': 1, 'gene': 'G', 'x_um': 0.1}],
                     'source_hashes': {'metadata.json': 'm'}}}
        b = json.loads(json.dumps(a))
        b['roi']['records'][0]['x_um'] = 0.2
        with self.assertRaisesRegex(union.ContractError, 'identical reference/count observations'):
            union.check_common_reference([{'atlases': a}, {'atlases': b}])

    def test_checksum_pin_rejects_edited_component_before_missing_other_records(self):
        pins = union.read_json(Path(union.__file__).with_name('frozen_inputs.json'))
        first = next(iter(pins['files']))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / first
            target.parent.mkdir(parents=True)
            target.write_text('changed component')
            with self.assertRaisesRegex(union.ContractError, 'source checksum mismatch'):
                union.check_pinned_inputs(root)

    def test_negative_or_noncanonical_count_and_negative_area_are_rejected(self):
        for value in ({'-1': '0'}, {'01': '0'}, {'1': '-1'}, {}):
            with self.assertRaises(ValueError):
                union.parse_minima(value, 'fixture')

    def test_output_cannot_overwrite_frozen_records(self):
        root = Path('/input/frozen')
        for destination in (root, root / 'results/reach', root.parent):
            with self.assertRaisesRegex(union.ContractError, 'separate from'):
                union.check_output_location(root, destination)
        union.check_output_location(root, Path('/output/union'))

    def test_fixed_fallback_requires_an_eligible_source_cell(self):
        member = {'spectra': {}, 'atlases': {'roi': {'complete_cells': [1], 'source_counts': {}}}}
        self.assertEqual(union.spectrum_or_fixed(member, ('roi', 1, 'G')), {0: F(0)})
        with self.assertRaisesRegex(union.ContractError, 'ineligible'):
            union.spectrum_or_fixed(member, ('roi', 2, 'G'))


if __name__ == '__main__':
    unittest.main()
