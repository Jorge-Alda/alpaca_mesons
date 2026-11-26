"""Unit tests for derive_sigma_pipi_loop_frame_invariance.py (exact symbolic identities only).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_frame_invariance
"""
import unittest

import sympy as sp

import derive_sigma_pipi_loop_frame_invariance as frame
import derive_sigma_pipi_loop_traces as traces


class FrameInvarianceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = frame.frames(radial_order=3, alp_order=3)

    def test_dressed_sources_unchanged(self):
        differences = frame.dressed_differences(self.data)
        self.assertEqual(len(differences), 32)
        self.assertTrue(all(not value for value in differences.values()))

    def test_raw_sources_change(self):
        original, raw = self.data['original'], self.data['raw new']
        changed = frame.poly_add(raw[1][(0, 1)], frame.poly_scale(original[1][(0, 1)], -1))
        self.assertTrue(changed)
        unchanged = frame.poly_add(raw[0][(1, 0)], frame.poly_scale(original[0][(1, 0)], -1))
        self.assertFalse(unchanged)

    def test_traces(self):
        entries = traces.kept_graph_traces()
        self.assertEqual(len(entries), 26 + 6 + 5 + 5)
        self.assertTrue(frame.sources_covered(entries))
        results = frame.trace_checks(entries, self.data)
        self.assertEqual(len(results), len(entries))
        for label, (unchanged, changed) in results.items():
            self.assertTrue(unchanged, label)
        changing = {label for label, (_, changed) in results.items() if changed}
        self.assertIn('tadpole {pi+,pi-,sigma}-A-', changing)
        self.assertIn('bubble {pi+,pi-}-M-{sigma}-M-', changing)
        self.assertNotIn('tadpole {pi+,pi-,sigma}-S-', changing)
        only_radial_sources = {entry['label'] for entry in entries
                               if all(name[2:] == '0' for name in entry['words'].letters() if name[0] in 'SP'
                                      and name[1:].isdigit())}
        self.assertEqual(changing, {entry['label'] for entry in entries} - only_radial_sources)

    def test_theta(self):
        theta = self.data['theta']
        self.assertEqual(sp.simplify(theta['new dressed'] - theta['original dressed']), 0)
        self.assertNotEqual(sp.simplify(theta['new raw'] - theta['original dressed']), 0)


if __name__ == '__main__':
    unittest.main()
