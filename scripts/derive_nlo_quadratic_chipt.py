#!/usr/bin/env python3
"""Derive the displayed NLO ChiPT quadratic terms with SymPy.

Run: .venv/bin/python scripts/derive_nlo_quadratic_chipt.py
Tests: .venv/bin/python scripts/derive_nlo_quadratic_chipt.py --test

Inputs: memory/nlo/lagrangian.md and the indexed Alda:2026nrd Appendix A
and topological-invariant extracts. Natural units; f_pi, f_a and fields
have mass dimension one, chi_0=2 B_0 m_q has dimension two. Arbitrary
noncommuting flavour matrices are retained. No loops or p^6 terms are
included. Only Phi^2 and a Phi coefficients are returned; a^2 is omitted.
"""

import sys
import unittest

import sympy as sp


def trace(expression):
    """Convert expanded noncommuting products to standard matrix traces.

    SymPy 1.14 quantum Tr fails on cyclic words with repeated symbols.
    Its matrix Trace lacks cyclic collection, so rotate each finite word
    using default_sort_key before constructing Trace; no commutation is used.
    """
    dimension = sp.Symbol("N_f", integer=True, positive=True)
    terms = []
    for term in sp.Add.make_args(sp.expand(expression)):
        commuting, noncommuting = term.args_cnc()
        word = []
        for factor in noncommuting:
            base, power = factor.as_base_exp()
            word.extend([base] * int(power))
        if not word:
            terms.append(sp.Mul(*commuting) * dimension)
            continue
        rotations = [word[index:] + word[:index] for index in range(len(word))]
        ordered = min(rotations, key=lambda values: tuple(sp.default_sort_key(value) for value in values))
        matrices = [sp.MatrixSymbol(str(value), dimension, dimension) for value in ordered]
        terms.append(sp.Mul(*commuting) * sp.Trace(sp.MatMul(*matrices)))
    return sp.Add(*terms)


def derive():
    """Expand Appendix A building blocks and extract each operator coefficient."""
    meson_order, alp_order = sp.symbols("t_Phi t_a", real=True)
    f_pi, f_a, a, derivative_a, c_G = sp.symbols(
        "f_pi f_a a partial_mu_a c_G", real=True
    )
    Phi, derivative_Phi, chi_0, c_A, c_P, c_S, c_V = sp.symbols(
        "Phi partial_mu_Phi chi_0 c_q^A c_q^P c_q^S c_q^V", commutative=False
    )
    phi = meson_order * Phi / (2 * f_pi)
    chi_S = chi_0 + sp.I * alp_order * a / f_a * (c_S * chi_0 - chi_0 * c_S)
    chi_P = alp_order * a / f_a * (c_P * chi_0 + chi_0 * c_P)

    def anticommutator(left, right):
        return left * right + right * left

    # The Lorentz index is contracted between repeated derivative factors.
    u_mu = (
        -meson_order * derivative_Phi / f_pi
        + 2 * alp_order * derivative_a / f_a * c_A
        - 2 * sp.I * alp_order * derivative_a / f_a * (phi * c_V - c_V * phi)
    )
    chi_plus = (
        2 * chi_S - anticommutator(phi, anticommutator(phi, chi_S))
        + 2 * anticommutator(phi, chi_P)
    )
    chi_minus = 2 * sp.I * (chi_P - anticommutator(phi, chi_S))
    X_TC = -meson_order * trace(Phi) / f_pi - 2 * c_G * alp_order * a / f_a
    operators = {
        "L_5": trace(u_mu * u_mu * chi_plus),
        "L_8": trace(chi_plus**2 + chi_minus**2) / 2,
        "Lambda_1": f_pi**2 * trace(u_mu)**2 / 12,
        "Lambda_2": sp.I * f_pi**2 * X_TC * trace(chi_minus) / 12,
        "L_4": trace(u_mu * u_mu) * trace(chi_plus),
        "L_6": trace(chi_plus)**2,
        "L_7": trace(chi_minus)**2,
        "L_18": trace(u_mu) * trace(u_mu * chi_plus),
        "L_25": -sp.I * X_TC * trace(chi_plus * chi_minus),
    }
    results = {}
    for name, expression in operators.items():
        expanded = sp.expand(expression)
        results[name] = {
            "PhiPhi": sp.simplify(expanded.coeff(meson_order, 2).coeff(alp_order, 0)),
            "aPhi": sp.simplify(expanded.coeff(meson_order, 1).coeff(alp_order, 1)),
        }
    return {
        "building_blocks": {"u_mu": u_mu, "chi_plus": chi_plus,
                            "chi_minus": chi_minus, "X_TC": X_TC},
        "operators": operators, "results": results,
        "symbols": locals(),
    }


class QuadraticDerivationTests(unittest.TestCase):
    """Independent scalar exponential and absent-coupling checks."""

    @classmethod
    def setUpClass(cls):
        cls.derivation = derive()

    def test_vector_and_scalar_couplings_absent(self):
        for sectors in self.derivation["results"].values():
            for expression in sectors.values():
                matrix_names = {str(matrix) for matrix in expression.atoms(sp.MatrixSymbol)}
                self.assertTrue(matrix_names.isdisjoint({"c_q^S", "c_q^V"}))

    def test_scalar_l8_exponential(self):
        # For one flavour, chi=chi_0*(1+2 i c_P a/f_a) and U=exp(i Phi/f_pi).
        d = self.derivation["symbols"]
        field, mass, coupling, t, q = sp.symbols("P M c t q", real=True)
        chi = mass * (1 + 2 * sp.I * coupling * q * d["a"] / d["f_a"])
        chi_conjugate = sp.conjugate(chi)
        exact = chi**2 * sp.exp(-2 * sp.I * t * field / d["f_pi"]) + chi_conjugate**2 * sp.exp(2 * sp.I * t * field / d["f_pi"])
        scalar_expansion = sp.series(exact, t, 0, 3).removeO().expand()
        for sector, powers in {"PhiPhi": (2, 0), "aPhi": (1, 1)}.items():
            result = self.derivation["results"]["L_8"][sector]
            # Tr reduces to identity in one flavour, including all products.
            result = result.replace(
                lambda expr: isinstance(expr, sp.Trace),
                lambda expr: sp.Mul(*[
                    {"Phi": field, "chi_0": mass, "c_q^P": coupling}.get(str(factor), factor)
                    for factor in sp.Mul.make_args(expr.args[0])
                ]),
            )
            result = result.replace(
                lambda expr: isinstance(expr, sp.MatPow),
                lambda expr: {"Phi": field, "chi_0": mass, "c_q^P": coupling}[str(expr.base)]**expr.exp,
            )
            expected = scalar_expansion.coeff(t, powers[0]).coeff(q, powers[1])
            self.assertEqual(sp.simplify(result - expected), 0)


def main():
    calculation = derive()
    print("Natural units; arbitrary noncommuting flavour matrices; Lorentz derivatives contracted.")
    print("Only Phi^2 and a Phi terms; no a^2, loops, or p^6 terms.\n")
    for label, expression in calculation["building_blocks"].items():
        print(f"Appendix A {label}:\n{sp.latex(expression)}\n")
    for name, sectors in calculation["results"].items():
        for sector, expression in sectors.items():
            print(f"Coefficient of {name}, {sector}:\n{sp.latex(expression)}\n")


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
