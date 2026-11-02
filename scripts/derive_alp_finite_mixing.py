"""Derive finite ALP mixing, pole normalization and external amplitudes.

Run: .venv/bin/python scripts/derive_alp_finite_mixing.py
Tests: .venv/bin/python scripts/derive_alp_finite_mixing.py --test
Inputs: memory/nlo/lsz-invariant-split.md and kinetic-block-alignment.md.
Real fields have dimension one; K is dimensionless; M^2 and s have dimension
two. The local quadratic form is exact in mixing, at a specified EFT order.
K and M^2 are symmetric; K is positive definite. A simple real pole and an
invertible meson inverse propagator are required for the Schur expressions.
The full generalized eigenproblem also applies where that block is singular.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    dimension = sp.Symbol("n", integer=True, positive=True)
    s = sp.Symbol("s", real=True)
    kinetic_meson = sp.MatrixSymbol("K_P", dimension, dimension)
    mass_meson = sp.MatrixSymbol("M_P^2", dimension, dimension)
    kinetic_mixing = sp.MatrixSymbol("K_{Pa}", dimension, 1)
    mass_mixing = sp.MatrixSymbol("M_{Pa}^2", dimension, 1)
    kinetic_alp, mass_alp = sp.symbols("K_aa M_aa^2", real=True)
    inverse_meson = s * kinetic_meson - mass_meson
    inverse_mixed = s * kinetic_mixing - mass_mixing
    inverse_alp = sp.Matrix([[s * kinetic_alp - mass_alp]])
    inverse_full = sp.BlockMatrix([
        [inverse_meson, inverse_mixed],
        [inverse_mixed.T, inverse_alp],
    ])
    unnormalized_meson_overlap = -inverse_meson.I * inverse_mixed
    unnormalized_overlap = sp.BlockMatrix([
        [unnormalized_meson_overlap], [sp.Identity(1)]
    ])
    pole_equation = sp.block_collapse(inverse_full * unnormalized_overlap)
    schur_inverse = pole_equation.blocks[1, 0]
    # SymPy 1.14 differentiates a collapsed abstract MatPow via scalar Pow
    # and raises TypeError. Derive d(A^-1)/ds from A A^-1=1 first,
    # then apply the product rule, keeping standard MatrixExpr objects.
    derivative_inverse_meson = -inverse_meson.I * sp.diff(inverse_meson, s) * inverse_meson.I
    derivative_schur = (
        sp.diff(inverse_alp, s)
        - sp.diff(inverse_mixed.T, s) * inverse_meson.I * inverse_mixed
        - inverse_mixed.T * derivative_inverse_meson * inverse_mixed
        - inverse_mixed.T * inverse_meson.I * sp.diff(inverse_mixed, s)
    ).doit()
    kinetic_full = sp.BlockMatrix([
        [kinetic_meson, kinetic_mixing],
        [kinetic_mixing.T, sp.Matrix([[kinetic_alp]])],
    ])
    pole_norm = sp.block_collapse(
        unnormalized_overlap.T * kinetic_full * unnormalized_overlap
    )
    meson_vertex = sp.MatrixSymbol("V_P", 1, dimension)
    alp_vertex = sp.MatrixSymbol("V_a", 1, 1)
    amputated_vertex = sp.BlockMatrix([[meson_vertex, alp_vertex]])
    projected_amplitude = sp.block_collapse(amputated_vertex * unnormalized_overlap)
    meson_shift = sp.MatrixSymbol("rho", dimension, 1)
    field_conversion = sp.BlockMatrix([
        [sp.Identity(dimension), -meson_shift],
        [sp.ZeroMatrix(1, dimension), sp.Identity(1)],
    ])
    converted_inverse = sp.block_collapse(
        field_conversion.T * inverse_full * field_conversion
    )
    converted_alp = converted_inverse.blocks[1, 1].expand()

    # A two-state canonically normalized example derives both eigenvalues
    # directly from the characteristic determinant, without an expansion.
    pion_mass, alp_mass, mass_mixing_scalar = sp.symbols(
        "m_pi0_sq m_a0_sq delta_sq", real=True
    )
    two_state_mass = sp.Matrix([
        [pion_mass, mass_mixing_scalar],
        [mass_mixing_scalar, alp_mass],
    ])
    characteristic = (s * sp.eye(2) - two_state_mass).det()
    two_state_eigenvalues = sp.solve(characteristic, s)

    # Exact finite rotation can turn a multi-a operator into a vertex with
    # only one physical ALP. This polynomial is an explicitly stated example.
    physical_alp, physical_pion, mixing_angle, coupling = sp.symbols(
        "a_phys pi_phys theta g", real=True
    )
    original_alp = sp.cos(mixing_angle) * physical_alp - sp.sin(mixing_angle) * physical_pion
    original_pion = sp.sin(mixing_angle) * physical_alp + sp.cos(mixing_angle) * physical_pion
    original_interaction = coupling * original_alp**2 * original_pion**2
    one_alp_interaction = sp.expand(original_interaction).coeff(physical_alp, 1) * physical_alp
    return locals()


class FiniteMixingTests(unittest.TestCase):
    def test_exact_pole_and_residue(self):
        # A nonorthogonal field change produces kinetic and mass mixing while
        # preserving known, nondegenerate poles and their normalized residues.
        conversion = sp.Matrix([[1, sp.Rational(1, 3)], [sp.Rational(1, 2), 1]])
        kinetic = conversion.T * conversion
        mass = conversion.T * sp.diag(2, 7) * conversion
        s = sp.Symbol("s")
        inverse = s * kinetic - mass
        meson_inverse, mixing_inverse = inverse[0, 0], inverse[0, 1]
        schur = sp.factor(inverse[1, 1] - mixing_inverse**2 / meson_inverse)
        self.assertEqual(sp.factor(meson_inverse * schur - inverse.det()), 0)
        for pole in (2, 7):
            overlap = sp.Matrix([-mixing_inverse / meson_inverse, 1]).subs(s, pole)
            norm = (overlap.T * kinetic * overlap)[0]
            self.assertEqual(sp.simplify(sp.diff(schur, s).subs(s, pole) - norm), 0)
            self.assertEqual(inverse.subs(s, pole) * overlap, sp.zeros(2, 1))
            propagator_residue = inverse.inv().applyfunc(lambda entry: sp.limit((s - pole) * entry, s, pole))
            self.assertEqual(sp.simplify(propagator_residue - overlap * overlap.T / norm), sp.zeros(2))

    def test_amplitude_and_pole_norm_under_field_change(self):
        kinetic = sp.Matrix([[2, 1], [1, 3]])
        overlap = sp.Matrix([2, -1])
        vertex = sp.Matrix([[3, 5]])
        conversion = sp.Matrix([[1, 2], [0, 1]])
        converted_overlap = conversion.inv() * overlap
        self.assertEqual((vertex * overlap)[0], (vertex * conversion * converted_overlap)[0])
        self.assertEqual(
            (overlap.T * kinetic * overlap)[0],
            (converted_overlap.T * conversion.T * kinetic * conversion * converted_overlap)[0],
        )

    def test_exact_spectrum_removes_perturbative_degeneracy(self):
        result = derive()
        for eigenvalue in result["two_state_eigenvalues"]:
            self.assertEqual(sp.simplify(result["characteristic"].subs(result["s"], eigenvalue)), 0)
        equal_diagonals = result["two_state_mass"].subs(result["alp_mass"], result["pion_mass"])
        self.assertEqual(sp.factor(equal_diagonals.det()), result["pion_mass"]**2 - result["mass_mixing_scalar"]**2)


def main():
    result = derive()
    print("Exact local quadratic mixing; simple stable pole; no expansion in f_pi/f_a.")
    print("All masses and interactions must follow from the same specified covariant EFT.")
    for label, name in [
        ("Full inverse propagator", "inverse_full"),
        ("Meson components relative to original a component", "unnormalized_meson_overlap"),
        ("Schur inverse propagator; its zeros determine the mixed poles", "schur_inverse"),
        ("Derivative at a pole gives inverse original-a residue", "derivative_schur"),
        ("Same pole normalization from the full kinetic metric", "pole_norm"),
        ("Unnormalized one-pole amplitude; divide by square root of pole norm", "projected_amplitude"),
        ("ALP diagonal block after P=Pprime-rho a; keep every term", "converted_alp"),
        ("Two-state mass matrix", "two_state_mass"),
        ("Characteristic determinant", "characteristic"),
        ("Exact two-state eigenvalues", "two_state_eigenvalues"),
        ("Example g a_original^2 P_original^2 after finite rotation: one physical ALP", "one_alp_interaction"),
    ]:
        show(label, result[name])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
