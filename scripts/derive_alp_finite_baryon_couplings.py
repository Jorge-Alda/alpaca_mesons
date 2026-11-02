"""Derive pointlike LO ALP-baryon couplings at an exact mixed pole.

Run: .venv/bin/python scripts/derive_alp_finite_baryon_couplings.py
Tests: .venv/bin/python scripts/derive_alp_finite_baryon_couplings.py --test

Inputs: memory/alp/shared/finite-alp-meson-mixing.md;
memory/nlo/lsz-invariant-split.md; Alda:2026nrd, Sec. V.
The relative pole eigenvector has original a component one. Use a complete
Hermitian meson basis tr(lambda_i lambda_j)=2 delta_ij, real symmetric
quadratic matrices, a simple stable pole, and the linear quark spurion of
the stated LO theory. No expansion in f_pi/f_a or the pole overlap is made.
Fields and f_pi,f_a have mass dimension one; chi_0,m_0^2,s have dimension
two; baryon coefficients and source matrices are dimensionless. Baryon
form factors are set to their LO pointlike values. Abstract flavor
matrices do not commute. Trace cyclicity is used only in the displayed
equivalent vector coefficient, and is checked with explicit matrices.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    flavor_count = sp.Symbol("N_f", integer=True, positive=True)
    f_pi, f_a = sp.symbols("f_pi f_a", positive=True)
    anomaly_mass, pole = sp.symbols("m_0^2 m_a^2", real=True)
    gluon_source = sp.Symbol("c_G", real=True)
    quark_spurion = sp.MatrixSymbol("chi_0", flavor_count, flavor_count)
    axial_source = sp.MatrixSymbol("c_q^A", flavor_count, flavor_count)
    pseudoscalar_source = sp.MatrixSymbol("c_q^P", flavor_count, flavor_count)
    vector_source = sp.MatrixSymbol("c_q^{VS}", flavor_count, flavor_count)
    meson_overlap = sp.MatrixSymbol("w_i^{(a)}lambda_i", flavor_count, flavor_count)
    effective_source = sp.MatrixSymbol("S_a", flavor_count, flavor_count)
    identity = sp.Identity(flavor_count)

    def mass_operator(flavor_matrix):
        return (
            (quark_spurion * flavor_matrix + flavor_matrix * quark_spurion) / 2
            + anomaly_mass / flavor_count * sp.Trace(flavor_matrix) * identity
            - pole * flavor_matrix
        )

    kinetic_mixing_flavor = -2 * f_pi / f_a * axial_source
    mass_mixing_flavor = (
        -f_pi / f_a * (quark_spurion * pseudoscalar_source
                        + pseudoscalar_source * quark_spurion)
        + 2 * f_pi * anomaly_mass * gluon_source / (flavor_count * f_a) * identity
    )
    pole_rhs = pole * kinetic_mixing_flavor - mass_mixing_flavor
    # Apply the linear operator to S_a=f_a W_a/(2 f_pi)-c_q^A.
    # The pole equation supplies M_a(W_a)=pole_rhs directly.
    source_from_pole = sp.expand(
        f_a / (2 * f_pi) * pole_rhs - mass_operator(axial_source)
    )
    invariant_source = (
        -(quark_spurion * (axial_source - pseudoscalar_source)
          + (axial_source - pseudoscalar_source) * quark_spurion) / 2
        - anomaly_mass / flavor_count * (gluon_source + sp.Trace(axial_source)) * identity
    )
    meson_overlap_substitution = 2 * f_pi / f_a * (axial_source + effective_source)
    vielbein_coefficient_before = -meson_overlap / f_pi + 2 * axial_source / f_a
    vielbein_coefficient = sp.expand(vielbein_coefficient_before.subs(
        meson_overlap, meson_overlap_substitution
    ))
    topological_coefficient_before = -sp.Trace(meson_overlap) / f_pi - 2 * gluon_source / f_a
    topological_coefficient = sp.expand(topological_coefficient_before.subs(
        meson_overlap, meson_overlap_substitution
    ).doit())
    # Completeness and tr(lambda_i lambda_j)=2 delta_ij imply that
    # 1+sum_i(w_i-v_i)^2 = 1+tr[(W_a-2 f_pi c_q^A/f_a)^2]/2.
    overlap_after_kinetic_shift = sp.expand(
        (meson_overlap - 2 * f_pi / f_a * axial_source).subs(
            meson_overlap, meson_overlap_substitution
        )
    )
    # Use MatMul's scalar/matrix decomposition because SymPy does not
    # distribute a scalar prefactor in an abstract matrix square.
    norm_scalar, norm_matrix = overlap_after_kinetic_shift.as_coeff_mmul()
    lo_pole_normalization = 1 + norm_scalar**2 / 2 * sp.Trace(norm_matrix**2)

    meson_count = sp.Symbol("n", integer=True, positive=True)
    kinetic_meson = sp.MatrixSymbol("K_PP", meson_count, meson_count)
    kinetic_mixing = sp.MatrixSymbol("K_Pa", meson_count, 1)
    kinetic_alp = sp.MatrixSymbol("K_aa", 1, 1)
    overlap_vector = sp.MatrixSymbol("w^{(a)}", meson_count, 1)
    kinetic_full = sp.BlockMatrix([
        [kinetic_meson, kinetic_mixing], [kinetic_mixing.T, kinetic_alp]
    ])
    coordinate_change = sp.BlockMatrix([
        [sp.Identity(meson_count), overlap_vector],
        [sp.ZeroMatrix(1, meson_count), sp.Identity(1)],
    ])
    kinetic_transformed = sp.block_collapse(coordinate_change.T * kinetic_full * coordinate_change)
    pole_normalization = kinetic_transformed.blocks[1, 1]

    baryon_first_adjoint = sp.Adjoint(sp.MatrixSymbol(r"\lambda_{B_1}", flavor_count, flavor_count))
    baryon_second = sp.MatrixSymbol(r"\lambda_{B_2}", flavor_count, flavor_count)
    d_baryon, f_baryon, s_baryon = sp.symbols("D_B F_B S_B", real=True)
    s_topological = sp.Symbol(r"\widetilde S_B", real=True)
    axial_coefficient = f_a * (
        d_baryon / 2 * sp.Trace(baryon_first_adjoint * (vielbein_coefficient * baryon_second
                                                      + baryon_second * vielbein_coefficient))
        + f_baryon / 2 * sp.Trace(baryon_first_adjoint * (vielbein_coefficient * baryon_second
                                                       - baryon_second * vielbein_coefficient))
        + s_baryon / 2 * sp.Trace(baryon_first_adjoint * baryon_second) * sp.Trace(vielbein_coefficient)
        + s_topological / 2 * sp.Trace(baryon_first_adjoint * baryon_second) * topological_coefficient
    )
    axial_coefficient = axial_coefficient.expand().doit().expand()
    vector_coefficient = sp.Trace(baryon_first_adjoint * (vector_source * baryon_second
                                                         - baryon_second * vector_source))
    vector_cyclic_form = -sp.Trace(vector_source * (baryon_first_adjoint * baryon_second
                                                   - baryon_second * baryon_first_adjoint))
    physical_axial_coefficient = axial_coefficient / sp.sqrt(lo_pole_normalization)
    physical_vector_coefficient = vector_cyclic_form / sp.sqrt(lo_pole_normalization)
    return locals()


class FiniteBaryonTests(unittest.TestCase):
    def test_source_equation_for_noncommuting_flavor_sources(self):
        result = derive()
        substitutions = {
            result["flavor_count"]: 3,
            result["quark_spurion"]: sp.diag(2, 3, 7),
            result["axial_source"]: sp.Matrix([[1, 2, 0], [2, 0, 1], [0, 1, 3]]),
            result["pseudoscalar_source"]: sp.Matrix([[0, 1, 2], [1, 3, 0], [2, 0, 1]]),
        }
        difference = (result["source_from_pole"] - result["invariant_source"]).subs(substitutions).doit()
        self.assertEqual(sp.simplify(difference.as_explicit()), sp.zeros(3))

    def test_ccwz_coefficients_and_pole_norm(self):
        result = derive()
        self.assertEqual(sp.expand(result["vielbein_coefficient"] + 2 * result["effective_source"] / result["f_a"]), sp.ZeroMatrix(result["flavor_count"], result["flavor_count"]))
        difference = result["topological_coefficient"] + 2 / result["f_a"] * (
            result["gluon_source"] + sp.Trace(result["axial_source"]) + sp.Trace(result["effective_source"])
        )
        self.assertEqual(sp.simplify(difference), 0)
        kinetic = sp.Matrix([[2, 1], [1, 3]])
        overlap = sp.Rational(4, 3)
        conversion = sp.Matrix([[1, overlap], [0, 1]])
        expected = (sp.Matrix([[overlap, 1]]) * kinetic * sp.Matrix([overlap, 1]))[0]
        self.assertEqual((conversion.T * kinetic * conversion)[1, 1], expected)

    def test_vector_trace_cyclicity(self):
        result = derive()
        substitutions = {
            result["baryon_first_adjoint"]: sp.Matrix([[1, 0, 1], [2, 0, 0], [0, 1, 1]]),
            result["baryon_second"]: sp.Matrix([[0, 1, 0], [1, 2, 1], [0, 0, 1]]),
            result["vector_source"]: sp.Matrix([[1, 2, 0], [2, 3, 1], [0, 1, 0]]),
        }
        self.assertEqual((result["vector_coefficient"] - result["vector_cyclic_form"]).subs(substitutions).doit(), 0)

    def test_lo_norm_in_invariant_source(self):
        result = derive()
        expected = 1 + 2 * result["f_pi"]**2 / result["f_a"]**2 * sp.Trace(result["effective_source"]**2)
        self.assertEqual(sp.simplify(result["lo_pole_normalization"] - expected), 0)


def main():
    result = derive()
    print("LO one-pseudoscalar vertices; exact simple stable pole; pointlike baryons; natural units.")
    for label, expression in [
        ("Pole equation M_a(W_a) equals", result["pole_rhs"]),
        ("Source obtained after S_a=f_a W_a/(2 f_pi)-c_q^A", result["source_from_pole"]),
        ("Equivalent invariant C", result["invariant_source"]),
        ("Thus M_a(S_a)=C", sp.Eq(result["mass_operator"](result["effective_source"]), result["source_from_pole"])),
        ("u_mu coefficient of partial_mu a", result["vielbein_coefficient"]),
        ("partial_mu X_TC coefficient of partial_mu a", result["topological_coefficient"]),
        ("Kprime_aa equals the exact pole normalization", result["pole_normalization"]),
        ("LO pole normalization in the invariant source S_a=M_a^-1(C)", result["lo_pole_normalization"]),
        ("Section V axial coefficient before physical field normalization", result["axial_coefficient"]),
        ("Section V vector coefficient by cyclic trace", result["vector_cyclic_form"]),
        ("Physical axial coefficient", result["physical_axial_coefficient"]),
        ("Physical vector coefficient", result["physical_vector_coefficient"]),
    ]:
        show(label, expression)


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
