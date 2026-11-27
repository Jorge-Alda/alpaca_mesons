"""Unit tests for derive_sigma_pipi_cartesian_loop_frame_invariance.py (noncommuting symbols; no matrices).

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_frame_invariance
"""
import pathlib
import unittest

import sympy as sp

import derive_sigma_pipi_cartesian_loop_frame_invariance as frame
from test_derive_sigma_pipi_cartesian_loop_traces import FORBIDDEN_COMPONENTS, FORBIDDEN_LOOP_VALUES

RECORD = (pathlib.Path(__file__).resolve().parent.parent / 'memory' / 'alp' / 'radial-mode' / 'finite-mixing'
          / 'two-pion-loop' / 'cartesian-estimate' / 'cartesian-sigma-pipi-loop-frame-invariance.md')
Polynomial = frame.Polynomial


class AlgebraTest(unittest.TestCase):
    def test_transport_of_constant_connection_is_exponential(self):
        kappa = frame.redefinition_letter()
        w = frame.transport(frame.Series({(0, 0): kappa}))
        exponential = frame.exponential_along_alp(kappa, 1)
        for q in range(frame.ORDER + 1):
            self.assertEqual(w.coefficient(0, q), exponential.coefficient(0, q))

    def test_exponentials_are_inverse(self):
        kappa = frame.redefinition_letter()
        product = frame.exponential_along_alp(kappa, 1) * frame.exponential_along_alp(kappa, -1)
        self.assertEqual(product.coefficient(0, 0), Polynomial.unit())
        for q in range(1, frame.ORDER + 1):
            self.assertTrue(product.coefficient(0, q).is_zero())

    def test_dagger(self):
        word = Polynomial({('SC10', 'kappa'): sp.I})
        self.assertEqual(word.dagger(), Polynomial({('kappa', 'SC10'): -sp.I}))


class CheckOneTest(unittest.TestCase):
    def test_dressed_coefficients_agree(self):
        count, agree = frame.dressed_coefficients_agree()
        self.assertEqual(count, 32)
        self.assertTrue(agree)

    def test_raw_coefficients_differ(self):
        """The raw coefficients with q >= 1 change (first order in kappa), those with q = 0 do not."""
        shifts = frame.raw_first_order_shift()
        for key, shift in shifts.items():
            self.assertEqual(shift, frame.expected_first_order_shift(*key))
            self.assertEqual(shift.is_zero(), key[2] == 0)


class TraceTest(unittest.TestCase):
    def test_entering_letters_in_range(self):
        letters = frame.entering_letters()
        self.assertTrue(frame.within_covered_range(letters))
        self.assertTrue({'SC00', 'PC01', 'SC10', 'SC02', 'SC20'} <= frame.pole_data_letters())

    def test_coefficients_are_scalars(self):
        self.assertTrue(frame.coefficients_are_scalars())

    def test_traces_unchanged(self):
        checks = frame.trace_checks()
        self.assertEqual(len(checks), 37)
        self.assertTrue(all(value[0] for value in checks.values()))

    def test_negative_control(self):
        changed, with_alp, with_shifted = frame.negative_control()
        self.assertEqual(changed, with_shifted)
        self.assertEqual(frame.letters_with_vanishing_shift(), {'SC01', 'SC11'})
        self.assertEqual(with_alp - changed, {'C7', 'C9', 'C19', 'C21', 'C39', 'Cpi7', 'Cpi9'})
        self.assertEqual(len(changed), 20)
        self.assertEqual(changed, {'C2', 'C3', 'C11', 'C12', 'C13', 'C14', 'C15', 'C41', 'C42', 'C43', 'C44', 'C45',
                                   'C46', 'C47', 'C48', 'Csigma5', 'Csigma6', 'Csigma7', 'Csigma8', 'Csigma9'})
        self.assertEqual(len(with_alp), 27)

    def test_theta(self):
        raw, dressed, others = frame.theta_check()
        self.assertEqual(raw, -frame.TRACE_OF_REDEFINITION)
        self.assertEqual(dressed, 0)
        self.assertTrue(others)


class RecordTest(unittest.TestCase):
    def test_record(self):
        text = RECORD.read_text()
        for token in FORBIDDEN_COMPONENTS + FORBIDDEN_LOOP_VALUES:
            self.assertNotIn(token, text, token)
        self.assertIn('-21/(64\\pi^2)', text)
        self.assertIn('authoritative', text)
        self.assertNotIn('action', text.lower().replace('interaction', ''))
        changed, with_alp, _ = frame.negative_control()
        self.assertIn(f'{len(changed)} of the {len(frame.trace_checks())}', text)
        self.assertIn(f'{len(with_alp)} kept', text)


if __name__ == '__main__':
    unittest.main()
