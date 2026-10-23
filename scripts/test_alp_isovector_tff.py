"""Symbolic checks for the rho projection; run with .venv/bin/python -m unittest discover -s scripts -p test_alp_isovector_tff.py."""

import unittest

import sympy as sp

import derive_alp_doubly_virtual_tff as derivation


class IsovectorTffTests(unittest.TestCase):
    """Check Bose symmetry, flavour zeros, and ordered-sum multiplicity."""

    @classmethod
    def setUpClass(cls):
        cls.terms = derivation.isovector_contribution()

    def test_bose_symmetry(self):
        for name in ("vjp_scalar", "rho_pair_scalar", "mixed_scalar"):
            expression = self.terms[name]
            exchanged = expression.xreplace({derivation.tensor.s1: derivation.tensor.s2,
                                              derivation.tensor.s2: derivation.tensor.s1})
            self.assertEqual(sp.simplify(expression - exchanged), 0)

    def test_isospin_symmetric_alp_has_no_mixed_exchange(self):
        expression = self.terms["mixed_pair"]
        symbols = {str(symbol): symbol for symbol in expression.free_symbols}
        substitutions = {symbols["c_d^A"]: symbols["c_u^A"],
                         symbols["c_d^P"]: symbols["c_u^P"],
                         symbols["chi_d"]: symbols["chi_u"]}
        self.assertEqual(sp.simplify(expression.subs(substitutions)), 0)

    def test_projection_equals_full_sum_minus_isoscalar_sum(self):
        tensor = derivation.tensor
        rho_mass, isoscalar_mass, rho_residue, isoscalar_residue = sp.symbols("M_rho M_I R_rho R_I")
        terms = self.terms
        isoscalar_vjp, isoscalar_pair = sp.symbols("V_I P_II")
        full = sum(derivation.assemble_scalar(
            [rho_mass, isoscalar_mass], [rho_residue, isoscalar_residue],
            [terms["vjp"], isoscalar_vjp],
            [[terms["rho_pair"], terms["mixed_pair"]], [terms["mixed_pair"], isoscalar_pair]]))
        isoscalar = sum(derivation.assemble_scalar(
            [isoscalar_mass], [isoscalar_residue], [isoscalar_vjp], [[isoscalar_pair]]))
        projected = terms["vjp_scalar"] + terms["rho_pair_scalar"] + terms["mixed_scalar"]
        actual_rho_mass = next(symbol for symbol in projected.free_symbols if str(symbol) == "M_rho")
        source_constant = next(symbol for symbol in projected.free_symbols if str(symbol) == "F_V")
        solved_source_constant = sp.solve(terms["residue"] - rho_residue, source_constant)[0]
        projected = projected.subs(actual_rho_mass, rho_mass).subs(source_constant, solved_source_constant)
        projected = projected.subs({sp.Function("W_I")(tensor.s1): isoscalar_residue / (isoscalar_mass**2 - tensor.s1),
                                    sp.Function("W_I")(tensor.s2): isoscalar_residue / (isoscalar_mass**2 - tensor.s2)})
        self.assertEqual(sp.simplify(full - isoscalar - projected), 0)

    def test_no_strange_or_double_trace_source_in_rho_residue(self):
        self.assertEqual({str(symbol) for symbol in self.terms["residue"].free_symbols},
                         {"F_V", "z_V", "x_l"})


if __name__ == "__main__":
    unittest.main()
