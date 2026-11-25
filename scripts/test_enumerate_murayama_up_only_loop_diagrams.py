"""Unit tests for enumerate_murayama_up_only_loop_diagrams.py (topology only, no vertex is evaluated).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_enumerate_murayama_up_only_loop_diagrams
"""
import math
import unittest

import enumerate_murayama_up_only_loop_diagrams as enumeration

# Group of diagrams -> diagram labels in the code list (all allowed by all rules, each once).
EXPECTED_LABELS = {
    'class 1': ['1.1', '1.2', '1.3'],
    'class 2': ['2.1', '2.2', '2.3', '2.4', '2.5'],
    'class 3': ['3.1', '3.2', '3.3'],
    'class 4 bubble': [f'4.{i}' for i in range(1, 7)],
    'class 4 triangle': [f'4.{i}' for i in range(7, 15)],
    'class 5 bubble': [f'5.{i}' for i in range(1, 5)],
    'class 5 triangle': [f'5.{i}' for i in range(5, 13)],
    'class 6 bubble': ['6.1'],
    'class 6 triangle': ['6.2', '6.3', '6.4'],
    'class 6 box': ['6.5', '6.6'],
    'class 6 sigma self-energy': ['6.7'],
}
# Cycles outside the estimate -> number of diagrams allowed by R0-R3 (non-zero ones dropped by the
# three-or-more-line truncation R4).
TRUNCATION = '(three-or-more-line truncation)'
EXPECTED_OTHER_COUNTS = {
    'A triangle': 0, 'A box': 0, 'two A lines, triangle': 0, 'two A lines, box': 0,
    'class 4 box': 0, 'class 5 box': 0, 'closed S line (rule R2)': 0, f'three S lines {TRUNCATION}': 0,
    f'two S lines and an A line, triangle {TRUNCATION}': 1, f'two S lines and an A line, box {TRUNCATION}': 2,
    f'two S lines and two A lines, box {TRUNCATION}': 0, f'two A lines and an S line, triangle {TRUNCATION}': 2,
    f'two A lines and an S line, box {TRUNCATION}': 0,
}
# The diagrams with a sigma^3 vertex: label -> (diagram, symmetry factor).
SIGMA_CUBED_DIAGRAMS = {
    '6.1': (((('beta',), ('0', 'A')), ('S', 'S')), 0.5),
    '6.3': (((('beta',), ('A',), ('0',)), ('S', 'M', 'S')), 1.0),
    '6.4': (((('gamma',), ('+',), ('-',)), ('S', 'M', 'S')), 1.0),
    '6.7': (((('beta',), ('gamma',)), ('S', 'S')), 0.5),
}


class TestDihedralCounting(unittest.TestCase):
    def test_box_with_four_single_legs_and_alternating_lines(self):
        # 4!/4: the alternating S, M pattern is kept by 4 of the 8 dihedral symmetries.
        found = enumeration.enumerate_diagrams(('S', 'S', 'M', 'M'), rules=('parity',))
        self.assertEqual(len(found), math.factorial(4) // 4)

    def test_symmetry_factors(self):
        self.assertEqual(enumeration.symmetry_factor(((('A', '+', '-', '0'),), ('A',))), 0.5)
        self.assertEqual(enumeration.symmetry_factor(((('A', '0'), ('+', '-')), ('A', 'A'))), 0.5)
        self.assertEqual(enumeration.symmetry_factor(((('A', '0'), ('+', '-')), ('A', 'M'))), 1.0)
        self.assertEqual(enumeration.symmetry_factor(((('A', '0'), ('+',), ('-',)), ('S', 'M', 'S'))), 1.0)

    def test_canonical_form_identifies_reflected_triangle(self):
        first = ((('A', '+'), ('-',), ('0',)), ('A', 'S', 'M'))
        reflected = ((('A', '+'), ('0',), ('-',)), ('M', 'S', 'A'))
        self.assertEqual(enumeration.canonical(first), enumeration.canonical(reflected))


class TestSelectionRules(unittest.TestCase):
    def test_charge_forbids_charged_vertex_between_neutral_lines(self):
        self.assertFalse(enumeration.charge_allowed((('A', '+'), ('-', '0')), ('A', 'A')))
        self.assertTrue(enumeration.charge_allowed((('A', '+'), ('-', '0')), ('A', 'M')))

    def test_radial_rule_needs_two_alp_carrying_fields(self):
        self.assertFalse(enumeration.radial_allowed((('+',), ('gamma', '-')), ('M', 'S')))
        self.assertTrue(enumeration.radial_allowed((('0',), ('beta', 'A')), ('M', 'S')))

    def test_sigma_cubed_vertex_has_no_further_field(self):
        self.assertTrue(enumeration.radial_allowed((('beta',), ('A', '0')), ('S', 'S')))
        self.assertTrue(enumeration.radial_allowed((('beta',), ('gamma',)), ('S', 'S')))
        self.assertFalse(enumeration.radial_allowed((('beta', 'A'), ('0',)), ('S', 'S')))
        self.assertFalse(enumeration.radial_allowed((('gamma',), ('+', '-')), ('S', 'S')))
        self.assertFalse(enumeration.radial_allowed((('beta', 'A', '0'),), ('S',)))

    def test_other_cycles(self):
        self.assertEqual(set(EXPECTED_OTHER_COUNTS), set(enumeration.OTHER_CYCLES))
        for name, line_multiset in enumeration.OTHER_CYCLES.items():
            with self.subTest(cycle=name):
                self.assertEqual(len(enumeration.enumerate_diagrams(line_multiset)), EXPECTED_OTHER_COUNTS[name])

    def test_excluded_sigma_self_energies_are_not_in_code_list(self):
        listed = {d for entries in enumeration.code_diagrams().values() for d, _, _ in entries}
        for line_multiset in enumeration.EXCLUDED_SIGMA_SELF_ENERGY_CYCLES.values():
            found = enumeration.enumerate_diagrams(line_multiset, object_sets=enumeration.SIGMA_SELF_ENERGY_OBJECTS)
            self.assertEqual(len(found), 1)
            self.assertFalse(found & listed)


class TestCodeList(unittest.TestCase):
    def setUp(self):
        self.code = enumeration.code_diagrams()

    def test_labels_and_counts(self):
        self.assertEqual(set(EXPECTED_LABELS), set(enumeration.CYCLE_GROUPS))
        self.assertEqual(set(self.code), set(enumeration.CYCLE_GROUPS))
        self.assertEqual(sum(len(labels) for labels in EXPECTED_LABELS.values()), 44)
        for group, labels in EXPECTED_LABELS.items():
            _, line_multiset, object_sets = enumeration.CYCLE_GROUPS[group]
            with self.subTest(group=group):
                self.assertEqual([label for _, _, label in self.code[group]], labels)
                self.assertEqual(len(enumeration.enumerate_diagrams(line_multiset, object_sets=object_sets)),
                                 len(labels))

    def test_code_list_is_exactly_the_allowed_diagrams_once(self):
        for group, (_, line_multiset, object_sets) in enumeration.CYCLE_GROUPS.items():
            with self.subTest(group=group):
                listed = [d for d, _, _ in self.code[group]]
                self.assertEqual(len(listed), len(set(listed)))
                self.assertEqual(enumeration.enumerate_diagrams(line_multiset, object_sets=object_sets), set(listed))

    def test_sigma_cubed_diagrams_listed_with_symmetry_factors(self):
        listed = {label: (d, w) for entries in self.code.values() for d, w, label in entries}
        for label, (diagram, weight) in SIGMA_CUBED_DIAGRAMS.items():
            with self.subTest(diagram=label):
                self.assertEqual(listed[label][0], enumeration.canonical(diagram))
                self.assertAlmostEqual(listed[label][1], weight)
                self.assertAlmostEqual(enumeration.symmetry_factor(diagram), weight)

    def test_code_symmetry_factors(self):
        for group, entries in self.code.items():
            for diagram, weight, label in entries:
                with self.subTest(group=group, diagram=label):
                    self.assertAlmostEqual(weight, enumeration.symmetry_factor(diagram))


if __name__ == '__main__':
    unittest.main()
