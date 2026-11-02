"""Derive LO finite ALP mixing with a canonical bare ALP quadratic term.

Run: .venv/bin/python scripts/derive_alp_lo_finite_mixing.py
Validation: .venv/bin/python scripts/derive_alp_lo_finite_mixing.py --test
Inputs: memory/extracts/Alda:2026nrd/appendix-a-expansion.md and
lo-topological-action.md; memory/nlo/kinetic-block-alignment.md.
Assume a CP-conserving stationary vacuum and Hermitian flavor sources.
The quark spurion is defined to be linear in a in these coordinates; no
independent a^2 quark spurion is added. The LO kinetic and anomaly squares
are kept in full. Fields and decay constants have dimension one, chi_0 and
m_0^2 have dimension two. No expansion in f_pi/f_a is made.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def expanded_trace(expression):
    """Expand the matrix polynomial before SymPy applies trace linearity."""
    return sp.Trace(expression.expand()).doit()


def derive():
    flavor_count = sp.Symbol("N_f", integer=True, positive=True)
    f_pi, f_a = sp.symbols("f_pi f_a", positive=True)
    anomaly_mass, bare_alp_mass = sp.symbols("m_0_sq m_a0_sq", real=True)
    gluon_coupling = sp.Symbol("c_G", real=True)
    alp, meson_i, meson_j = sp.symbols("a P_i P_j", real=True)
    derivative_alp, derivative_i, derivative_j = sp.symbols("d_a d_P_i d_P_j", real=True)
    generator_i = sp.MatrixSymbol("lambda_i", flavor_count, flavor_count)
    generator_j = sp.MatrixSymbol("lambda_j", flavor_count, flavor_count)
    axial_source = sp.MatrixSymbol("c_q^A", flavor_count, flavor_count)
    pseudoscalar_source = sp.MatrixSymbol("c_q^P", flavor_count, flavor_count)
    quark_spurion = sp.MatrixSymbol("chi_0", flavor_count, flavor_count)
    meson_matrix = meson_i * generator_i + meson_j * generator_j
    derivative_meson = derivative_i * generator_i + derivative_j * generator_j
    vielbein = -derivative_meson / f_pi + 2 * derivative_alp * axial_source / f_a
    kinetic_lagrangian = derivative_alp**2 / 2 + f_pi**2 * expanded_trace(vielbein * vielbein) / 4
    topological_invariant = -expanded_trace(meson_matrix) / f_pi - 2 * gluon_coupling * alp / f_a
    chi_pseudoscalar = alp * (pseudoscalar_source * quark_spurion + quark_spurion * pseudoscalar_source) / f_a
    # This is the quadratic expansion of f_pi^2 <chi_+>/4, after cyclicity.
    mass_lagrangian = (
        -expanded_trace(meson_matrix * meson_matrix * quark_spurion) / 4
        + f_pi * expanded_trace(meson_matrix * chi_pseudoscalar) / 2
        - f_pi**2 * anomaly_mass * topological_invariant**2 / (4 * flavor_count)
        - bare_alp_mass * alp**2 / 2
    )
    # SymPy does not expand abstract matrix powers under a Trace, and its
    # Trace scalar derivative can fail there. Apply the quadratic product
    # rule before forming each Trace; matrix derivatives remain standard.
    kinetic_entries = [(
        sp.diff(derivative_alp**2 / 2, left, right)
        + f_pi**2 / 4 * sp.Trace(sp.diff(vielbein, left) * sp.diff(vielbein, right)
                                + sp.diff(vielbein, right) * sp.diff(vielbein, left))
    ).doit() for left, right in (
        (derivative_i, derivative_j), (derivative_i, derivative_alp), (derivative_alp, derivative_alp)
    )]
    topological_invariant = (-meson_i * sp.Trace(generator_i) - meson_j * sp.Trace(generator_j)) / f_pi - 2 * gluon_coupling * alp / f_a
    mass_entries = []
    for left, right in ((meson_i, meson_j), (meson_i, alp), (alp, alp)):
        mass_entries.append((
            sp.Trace((sp.diff(meson_matrix, left) * sp.diff(meson_matrix, right)
                      + sp.diff(meson_matrix, right) * sp.diff(meson_matrix, left)) * quark_spurion) / 4
            - f_pi / 2 * sp.Trace(sp.diff(meson_matrix, left) * sp.diff(chi_pseudoscalar, right)
                                   + sp.diff(meson_matrix, right) * sp.diff(chi_pseudoscalar, left))
            + f_pi**2 * anomaly_mass / (2 * flavor_count) * sp.diff(topological_invariant, left) * sp.diff(topological_invariant, right)
            + sp.diff(bare_alp_mass * alp**2 / 2, left, right)
        ).doit())

    # The normalized complete U(N_f) generator basis gives K_PP=I and K_Pa=-v.
    dimension = sp.Symbol("n", integer=True, positive=True)
    axial_overlap = sp.MatrixSymbol("v", dimension, 1)
    meson_mass = sp.MatrixSymbol("M_PP_sq", dimension, dimension)
    mass_mixing = sp.MatrixSymbol("M_Pa_sq", dimension, 1)
    alp_mass_entry = sp.MatrixSymbol("M_aa_sq", 1, 1)
    kinetic_full = sp.BlockMatrix([
        [sp.Identity(dimension), -axial_overlap],
        [-axial_overlap.T, sp.Identity(1) + axial_overlap.T * axial_overlap],
    ])
    mass_full = sp.BlockMatrix([[meson_mass, mass_mixing], [mass_mixing.T, alp_mass_entry]])
    field_conversion = sp.BlockMatrix([
        [sp.Identity(dimension), axial_overlap], [sp.ZeroMatrix(1, dimension), sp.Identity(1)]
    ])
    canonical_kinetic = sp.block_collapse(field_conversion.T * kinetic_full * field_conversion)
    canonical_mass = sp.block_collapse(field_conversion.T * mass_full * field_conversion)
    return locals()


class LoFiniteMixingTests(unittest.TestCase):
    def test_exact_canonicalization_and_poles(self):
        axial_overlap = sp.Rational(3, 2)
        kinetic = sp.Matrix([[1, -axial_overlap], [-axial_overlap, 1 + axial_overlap**2]])
        mass = sp.Matrix([[3, 1], [1, 5]])
        conversion = sp.Matrix([[1, axial_overlap], [0, 1]])
        self.assertEqual(conversion.T * kinetic * conversion, sp.eye(2))
        pole = sp.Symbol("s")
        self.assertEqual(sp.expand((pole * kinetic - mass).det() - (pole * sp.eye(2) - conversion.T * mass * conversion).det()), 0)

    def test_bare_term_does_not_remove_source_squares(self):
        result = derive()
        axial_source = sp.Matrix([[1, 2], [2, 3]])
        kinetic_alp = result["kinetic_entries"][2].subs(result["axial_source"], axial_source).doit()
        kinetic_alp = kinetic_alp.replace(lambda item: isinstance(item, sp.Trace), lambda item: sp.trace(item.arg.as_explicit()))
        self.assertEqual(sp.expand(kinetic_alp - 1 - 2 * result["f_pi"]**2 / result["f_a"]**2 * sp.trace(axial_source**2)), 0)
        self.assertEqual(sp.simplify(result["mass_entries"][2] - result["bare_alp_mass"] - 2 * result["f_pi"]**2 * result["anomaly_mass"] * result["gluon_coupling"]**2 / (result["flavor_count"] * result["f_a"]**2)), 0)


def main():
    result = derive()
    print("LO chiral theory; linear quark spurion; full kinetic/anomaly squares; canonical bare ALP term.")
    for label, name in [
        ("Linear vielbein", "vielbein"),
        ("Quadratic kinetic Lagrangian", "kinetic_lagrangian"),
        ("Quadratic nonderivative Lagrangian", "mass_lagrangian"),
        ("Hessian entries K_ij, K_ia, K_aa", "kinetic_entries"),
        ("Hessian entries M_ij^2, M_ia^2, M_aa^2", "mass_entries"),
        ("K with complete normalized generator basis", "kinetic_full"),
        ("Coordinate change P=Q+v a", "field_conversion"),
        ("Canonical kinetic matrix", "canonical_kinetic"),
        ("Canonical full mass matrix", "canonical_mass"),
    ]:
        show(label, result[name])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
