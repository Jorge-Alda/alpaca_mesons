"""Tests for the decoupling of the charged mesons at LO (sources commuting with Q)."""

import unittest

import sympy as sp

from derive_alp_charged_meson_decoupling import QUARK_CHARGE, derive, symbolic_hermitian
from symbolic_utils import commutator


class ChargedMesonDecouplingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for label, passed in self.result["checks"].items():
            self.assertTrue(passed, label)

    def test_symbolic_sources_commute_with_charge(self):
        self.assertEqual(commutator(QUARK_CHARGE, symbolic_hermitian("X")), sp.zeros(3))
        self.assertNotEqual(commutator(QUARK_CHARGE, symbolic_hermitian("X", charge_conserving=False)), sp.zeros(3))

    def test_general_sources_couple_charged_mesons(self):
        self.assertGreater(self.result["general_summary"]["number of nonzero charged-neutral entries for general sources"], 0)

    def test_charged_kinetic_mixing_is_the_u_entry_of_axial_source(self):
        kinetic_entries = self.result["general"]["charged-neutral kinetic entries"]
        ratio = sp.Symbol("f_pi", positive=True) / sp.Symbol("f_a", positive=True)
        self.assertEqual(sp.expand(kinetic_entries[0, 5] + 2 * ratio * sp.Symbol("cA_ud_re", real=True)), 0)
        self.assertEqual(kinetic_entries[:, :5], sp.zeros(4, 5))


if __name__ == "__main__":
    unittest.main()
