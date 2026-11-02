#!/usr/bin/env python3
"""Derive NLO neutral pseudoscalar matrices in the (3,8,0) basis.

Natural units.  The quadratic convention is
L = (1/2) dP^T delta_K dP - (1/2) P^T delta_M_squared P.
Run with .venv/bin/python scripts/derive_nlo_neutral_matrices.py.
Optional unit tests: append --test.
Inputs: memory/nlo/quadratic-lagrangian.md.
"""

import sys
import unittest

import sympy as sp


def derive_matrices():
    """Project the recorded quadratic Lagrangian by taking field Hessians."""
    b_zero, f_pi = sp.symbols("B_0 f_pi", nonzero=True, real=True)
    mass_up, mass_down, mass_strange = sp.symbols("m_u m_d m_s", real=True)
    pi_three, eta_eight, eta_zero = sp.symbols("pi_3 eta_8 eta_0", real=True)
    derivative_pi_three, derivative_eta_eight, derivative_eta_zero = sp.symbols(
        "partial_pi_3 partial_eta_8 partial_eta_0", real=True
    )
    fields = sp.Matrix([pi_three, eta_eight, eta_zero])
    derivatives = sp.Matrix(
        [derivative_pi_three, derivative_eta_eight, derivative_eta_zero]
    )
    generators = (
        sp.diag(1, -1, 0),
        sp.diag(1, 1, -2) / sp.sqrt(3),
        sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
    )
    phi_matrix = sum(
        (field * generator for field, generator in zip(fields, generators)),
        sp.zeros(3),
    )
    derivative_phi = sum(
        (
            derivative * generator
            for derivative, generator in zip(derivatives, generators)
        ),
        sp.zeros(3),
    )
    chi_zero = 2 * b_zero * sp.diag(mass_up, mass_down, mass_strange)
    trace = sp.trace

    # These expressions are the recorded Phi^2 terms before projection.
    kinetic_operators = {
        "L_5": 2 * trace(chi_zero * derivative_phi**2) / f_pi**2,
        "Lambda_1": trace(derivative_phi) ** 2 / 12,
        "L_4": 2 * trace(derivative_phi**2) * trace(chi_zero) / f_pi**2,
        "L_18": (
            2 * trace(chi_zero * derivative_phi) * trace(derivative_phi) / f_pi**2
        ),
    }
    nonderivative_operators = {
        "L_8": (
            -2
            * (
                trace(phi_matrix**2 * chi_zero**2)
                + trace(phi_matrix * chi_zero * phi_matrix * chi_zero)
            )
            / f_pi**2
        ),
        "Lambda_2": -trace(phi_matrix * chi_zero) * trace(phi_matrix) / 6,
        "L_6": -4 * trace(phi_matrix**2 * chi_zero) * trace(chi_zero) / f_pi**2,
        "L_7": -4 * trace(phi_matrix * chi_zero) ** 2 / f_pi**2,
        "L_25": 4 * trace(phi_matrix * chi_zero**2) * trace(phi_matrix) / f_pi**2,
    }
    kinetic_matrices = {
        coupling: sp.hessian(expression, derivatives).applyfunc(sp.factor)
        for coupling, expression in kinetic_operators.items()
    }
    mass_matrices = {
        coupling: (-sp.hessian(expression, fields)).applyfunc(sp.factor)
        for coupling, expression in nonderivative_operators.items()
    }
    couplings = {
        name: sp.Symbol(name, real=True)
        for name in (*kinetic_matrices, *mass_matrices)
    }
    total_kinetic = sum(
        (couplings[name] * matrix for name, matrix in kinetic_matrices.items()),
        sp.zeros(3),
    ).applyfunc(sp.factor)
    total_mass_squared = sum(
        (couplings[name] * matrix for name, matrix in mass_matrices.items()),
        sp.zeros(3),
    ).applyfunc(sp.factor)
    return {
        "B_0": b_zero,
        "f_pi": f_pi,
        "m_u": mass_up,
        "m_d": mass_down,
        "m_s": mass_strange,
        "fields": fields,
        "generators": generators,
        "Phi": phi_matrix,
        "chi_0": chi_zero,
        "kinetic_operators": kinetic_operators,
        "nonderivative_operators": nonderivative_operators,
        "kinetic": kinetic_matrices,
        "mass_squared": mass_matrices,
        "total_kinetic": total_kinetic,
        "total_mass_squared": total_mass_squared,
    }


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    result = derive_matrices()
    print("Natural units; real diagonal quark masses; no ALP.")
    print("Basis (pi_3, eta_8, eta_0); tr(lambda_i lambda_j)=2 delta_ij.")
    print("L_quad,NLO = 1/2 dP^T delta_K dP - 1/2 P^T delta_M_squared P.")
    print("Displayed local L^(1), L^(2) terms only; no loops or p^6 terms.\n")
    show("Phi", result["Phi"])
    show("chi_0", result["chi_0"])
    for coupling, expression in result["kinetic_operators"].items():
        show(f"Quadratic derivative coefficient of {coupling}", sp.factor(expression))
        show(f"Coefficient of {coupling} in delta_K", result["kinetic"][coupling])
    for coupling, expression in result["nonderivative_operators"].items():
        show(f"Quadratic nonderivative coefficient of {coupling}", sp.factor(expression))
        show(
            f"Coefficient of {coupling} in delta_M_squared",
            result["mass_squared"][coupling],
        )
    for name in ("total_kinetic", "total_mass_squared"):
        for row in range(3):
            for column in range(row, 3):
                show(f"{name} entry ({row + 1},{column + 1})", result[name][row, column])


class NeutralMatrixTests(unittest.TestCase):
    """Check normalization and the flavor symmetry limits."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_matrices()

    def test_generator_normalization(self):
        for row, left in enumerate(self.result["generators"]):
            for column, right in enumerate(self.result["generators"]):
                self.assertEqual(
                    sp.simplify(sp.trace(left * right)), 2 if row == column else 0
                )

    def test_isospin_limit_decouples_pi_three(self):
        substitution = {self.result["m_d"]: self.result["m_u"]}
        for sector in ("kinetic", "mass_squared"):
            for matrix in self.result[sector].values():
                self.assertEqual(sp.simplify(matrix[0, 1].subs(substitution)), 0)
                self.assertEqual(sp.simplify(matrix[0, 2].subs(substitution)), 0)

    def test_chiral_limit(self):
        substitution = {
            self.result["m_u"]: 0,
            self.result["m_d"]: 0,
            self.result["m_s"]: 0,
        }
        for matrix in self.result["mass_squared"].values():
            self.assertEqual(matrix.subs(substitution), sp.zeros(3))
        for coupling, matrix in self.result["kinetic"].items():
            expected = sp.diag(0, 0, 1) if coupling == "Lambda_1" else sp.zeros(3)
            self.assertEqual(matrix.subs(substitution), expected)


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
