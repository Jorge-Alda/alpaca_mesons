"""Project the local NLO ALP mixing and its invariant source.

Run .venv/bin/python scripts/derive_nlo_lsz_source_projection.py [--test].
Natural units; diagonal neutral ALP sources; basis (3,8,0).  The Lagrangian
is L_mix = partial(a) k^T partial(P) - a mu^T P.  Local LECs are retained
linearly.  Loops, O(1/f_a**2), and higher-order products of LECs are omitted.
"""

import sys
import unittest

import sympy as sp

from derive_nlo_neutral_matrices import derive_matrices


def derive():
    meson = derive_matrices()
    f_pi, b_zero = meson["f_pi"], meson["B_0"]
    f_a, m_zero_squared = sp.symbols("f_a m_0^2", nonzero=True, real=True)
    c_g, c_ga = sp.symbols("c_G c_GA", real=True)
    axial = sp.diag(*sp.symbols("c_u^A c_d^A c_s^A", real=True))
    pseudoscalar = sp.diag(*sp.symbols("c_u^P c_d^P c_s^P", real=True))
    invariant_ap = sp.diag(*sp.symbols("c_u^{AP} c_d^{AP} c_s^{AP}", real=True))
    phi, chi = meson["Phi"], meson["chi_0"]
    generators = meson["generators"]
    fields = meson["fields"]
    derivative_fields = sp.Matrix(sp.symbols("d_pi_3 d_eta_8 d_eta_0", real=True))
    derivative_phi = sum((value * matrix for value, matrix in zip(derivative_fields, generators)), sp.zeros(3))
    trace = sp.trace
    # Coefficients of partial(a) and a, respectively, copied from the
    # indexed quadratic-Lagrangian record before differentiating the fields.
    kinetic_operators = {
        "L_5": -4 * (trace(axial * chi * derivative_phi) + trace(axial * derivative_phi * chi)) / (f_a * f_pi),
        "Lambda_1": -f_pi * trace(axial) * trace(derivative_phi) / (3 * f_a),
        "L_4": -8 * trace(axial * derivative_phi) * trace(chi) / (f_a * f_pi),
        "L_18": -4 * (trace(axial * chi) * trace(derivative_phi) + trace(chi * derivative_phi) * trace(axial)) / (f_a * f_pi),
    }
    nonderivative_operators = {
        "L_8": 4 * (trace(phi * chi**2 * pseudoscalar) + trace(phi * pseudoscalar * chi**2) + 2 * trace(phi * chi * pseudoscalar * chi)) / (f_a * f_pi),
        "Lambda_2": f_pi * (-c_g * trace(phi * chi) + trace(pseudoscalar * chi) * trace(phi)) / (3 * f_a),
        "L_6": 8 * (trace(phi * pseudoscalar * chi) + trace(phi * chi * pseudoscalar)) * trace(chi) / (f_a * f_pi),
        "L_7": 16 * trace(phi * chi) * trace(pseudoscalar * chi) / (f_a * f_pi),
        "L_25": 8 * (c_g * trace(phi * chi**2) - trace(pseudoscalar * chi**2) * trace(phi)) / (f_a * f_pi),
    }
    kinetic_mixing = {name: sp.Matrix([sp.diff(expression, value) for value in derivative_fields]).applyfunc(sp.factor) for name, expression in kinetic_operators.items()}
    mass_mixing = {name: sp.Matrix([-sp.diff(expression, value) for value in fields]).applyfunc(sp.factor) for name, expression in nonderivative_operators.items()}
    axial_shift = sp.Matrix([f_pi * trace(generator * axial) / f_a for generator in generators])
    leading_mass = sp.Matrix(3, 3, lambda row, column: trace(generators[row] * generators[column] * chi) / 2 + m_zero_squared * trace(generators[row]) * trace(generators[column]) / 6)
    leading_kinetic_mixing = -axial_shift
    # Derive the LO mass mixing from f_pi^2 tr(chi_+)/4 and the recorded
    # topological term, using the linear a*Phi pieces in those building blocks.
    leading_nonderivative = f_pi * trace(phi * pseudoscalar * chi) / f_a - f_pi * m_zero_squared * c_g * trace(phi) / (3 * f_a)
    leading_mass_mixing = sp.Matrix([-sp.diff(leading_nonderivative, value) for value in fields])
    leading_invariant = (-leading_mass * axial_shift - leading_mass_mixing).applyfunc(sp.expand)
    source_substitution = {pseudoscalar[index, index]: axial[index, index] - invariant_ap[index, index] for index in range(3)}
    source_substitution[c_g] = c_ga - trace(axial)
    leading_invariant = leading_invariant.subs(source_substitution).applyfunc(sp.factor)
    invariant_corrections = {}
    for name, correction in mass_mixing.items():
        invariant_corrections[name] = (-meson["mass_squared"][name] * axial_shift - correction).subs(source_substitution).applyfunc(sp.factor)
    # Completeness reconstructs the diagonal flavor matrix whose generator
    # projections are (f_pi/f_a) times the invariant-source components.
    leading_flavor_source = sp.simplify(sum((generator * value for generator, value in zip(generators, leading_invariant)), sp.zeros(3)) * f_a / (2 * f_pi))
    flavor_corrections = {name: sp.simplify(sum((generator * value for generator, value in zip(generators, vector)), sp.zeros(3)) * f_a / (2 * f_pi)) for name, vector in invariant_corrections.items()}
    # Recover coefficients of the spurion structures present in each original
    # potential operator, using polynomial coefficient matching after projection.
    chi_abstract = sp.MatrixSymbol("chi_0", 3, 3)
    ap_abstract = sp.MatrixSymbol("c_q^{AP}", 3, 3)
    candidate_structures = {
        "L_8": ((chi**2 * invariant_ap, chi_abstract**2 * ap_abstract),),
        "Lambda_2": ((c_ga * chi, c_ga * chi_abstract), (trace(chi * invariant_ap) * sp.eye(3), sp.Trace(chi_abstract * ap_abstract) * sp.Identity(3))),
        "L_6": ((trace(chi) * chi * invariant_ap, sp.Trace(chi_abstract) * chi_abstract * ap_abstract),),
        "L_7": ((trace(chi * invariant_ap) * chi, sp.Trace(chi_abstract * ap_abstract) * chi_abstract),),
        "L_25": ((c_ga * chi**2, c_ga * chi_abstract**2), (trace(chi**2 * invariant_ap) * sp.eye(3), sp.Trace(chi_abstract**2 * ap_abstract) * sp.Identity(3))),
    }
    compact_flavor_corrections = {}
    polynomial_variables = (meson["m_u"], meson["m_d"], meson["m_s"], *invariant_ap.diagonal(), c_ga)
    for name, structures in candidate_structures.items():
        coefficients = sp.symbols(f"coefficient_0:{len(structures)}")
        difference = flavor_corrections[name] - sum((coefficient * structure[0] for coefficient, structure in zip(coefficients, structures)), sp.zeros(3))
        equations = []
        for expression in difference:
            equations.extend(sp.Poly(sp.expand(expression), *polynomial_variables).coeffs())
        solution = sp.solve(equations, coefficients, dict=True)[0]
        compact_flavor_corrections[name] = sum((solution[coefficient] * structure[1] for coefficient, structure in zip(coefficients, structures)), sp.ZeroMatrix(3, 3))
    return locals()


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


class SourceProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_kinetic_terms_are_axial_shift_of_meson_terms(self):
        result = self.result
        for name, vector in result["kinetic_mixing"].items():
            self.assertEqual((vector + result["meson"]["kinetic"][name] * result["axial_shift"]).applyfunc(sp.simplify), sp.zeros(3, 1))

    def test_invariant_source_contains_no_axial_coefficients(self):
        result = self.result
        for vector in (result["leading_invariant"], *result["invariant_corrections"].values()):
            self.assertFalse(vector.free_symbols.intersection(result["axial"].free_symbols))

    def test_lo_flavor_source_matches_recorded_C(self):
        result = self.result
        expected = -result["chi"] * result["invariant_ap"] - result["m_zero_squared"] * result["c_ga"] * sp.eye(3) / 3
        self.assertEqual((result["leading_flavor_source"] - expected).applyfunc(sp.simplify), sp.zeros(3))

    def test_flavor_reconstruction_recovers_all_projections(self):
        result = self.result
        for name, flavor_source in result["flavor_corrections"].items():
            reconstructed = sp.Matrix([result["f_pi"] * sp.trace(generator * flavor_source) / result["f_a"] for generator in result["generators"]])
            self.assertEqual((reconstructed - result["invariant_corrections"][name]).applyfunc(sp.simplify), sp.zeros(3, 1))


def main():
    result = derive()
    print("Natural units; real diagonal neutral sources; local LECs linearly; no loops.")
    print("L_mix = partial(a) k^T partial(P) - a mu^T P; basis (3,8,0).")
    show("Axial meson shift r_A", result["axial_shift"])
    show("LO meson mass matrix", result["leading_mass"])
    show("LO mass mixing mu", result["leading_mass_mixing"])
    for name, vector in result["kinetic_mixing"].items():
        show(f"Coefficient of {name} in delta k", vector)
    for name, vector in result["mass_mixing"].items():
        show(f"Coefficient of {name} in delta mu", vector)
    show("LO invariant flavor source", result["leading_flavor_source"])
    for name, matrix in result["flavor_corrections"].items():
        show(f"Coefficient of {name} in invariant flavor source", matrix)
        show(f"Coefficient of {name} after solving for spurion structures", result["compact_flavor_corrections"][name])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
