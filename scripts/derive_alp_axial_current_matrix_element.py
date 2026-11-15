"""Derive the vacuum-to-ALP matrix element of the light-quark axial currents.

Run .venv/bin/python scripts/derive_alp_axial_current_matrix_element.py [--test].

Natural units; metric (+---), epsilon^{0123}=+1. Conventions of Alda:2026nrd:
Phi = sum_P lambda_P P, u = exp(i Phi/(2 f_pi)), f_pi = 93 MeV,
tr(lambda_i lambda_j) = 2 delta_ij, lambda_0 = sqrt(2/3) 1. The quark-level
axial source is defined by L_QCD containing qbar gamma^mu gamma5 A_mu q, so
that qbar_f gamma^mu gamma5 q_f = dL/dA_mu^{ff} for a diagonal source. The ALP
derivative coupling (d_mu a/f_a) qbar gamma^mu c_q^A gamma5 q of Alda:2026nrd
has the same form, hence the source enters u_mu exactly as 2 c_q^A d_mu a/f_a
does: u_mu contains 2 A_mu.

Order and approximations (all explicit): leading-order chiral Lagrangian,
L_2 = (f_pi^2/4) tr(u_mu u^mu) + terms without explicit axial-source
dependence; linear order in one ALP, O(1/f_a); three light flavors; diagonal
sources. The ALP sources are those of the no-mixing basis of Alda:2026nrd:
U_a = diag(U_u, U_d, U_s) = -2 S/f_a is the coefficient of d_mu a in u_mu,
Y_a = 2 {c_q^{AP} + S, chi}/f_a, and the topological contact coefficient is
Z_a = -(2 g/f_a)(c_GA + tr S) with g = m_0^2/3. In that basis the meson
fields do not connect the vacuum to the ALP at O(1/f_a).

f_pi is the leading-order chiral coupling constant multiplying tr(u_mu u^mu).
It is not a meson decay constant used as an ALP input; no eta or eta' mass or
decay constant enters any result of this script.

The local short-distance coefficient is taken from
memory/alp/lepton-pair/alp-hard-current-ope.md:
F contains 2i/((q^2+i0) m_a^2) sum_f Q_f^2 p^sigma <0|qbar_f gamma_sigma
gamma5 q_f|a(p)>, with Q_f the quark charges and q = (q_1 - q_2)/2.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp


def show(label, expression):
    """Print expressions through the standard SymPy LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def lorentz_square(metric, lower_components):
    """Contract a list of lower-index flavor matrices with the inverse metric."""
    inverse_metric = metric.inv()
    return sum((inverse_metric[first, second] * lower_components[first] * lower_components[second]
                for first in range(4) for second in range(4)), sp.zeros(3))


def axial_currents(lagrangian, axial_source, metric):
    """Return dL/dA_mu^{ff} at vanishing source, as rows mu (upper), columns f."""
    return sp.Matrix(4, 3, lambda lorentz, flavor: sp.expand(
        sp.diff(lagrangian, axial_source[lorentz][flavor, flavor])
        .subs({axial_source[index][entry, entry]: 0 for index in range(4) for entry in range(3)})))


@lru_cache(maxsize=1)
def derive():
    """Carry out the five steps of the derivation and return all intermediate objects."""
    metric = sp.diag(1, -1, -1, -1)
    pion_coupling = sp.Symbol("f_pi", positive=True)
    alp_coupling_scale = sp.Symbol("f_a", positive=True)
    alp_mass = sp.Symbol("m_a", positive=True)
    flavor_names = ("u", "d", "s")

    # ---- Step 1: plane wave, axial current from L_2, vacuum-to-ALP element.
    position = sp.Matrix(sp.symbols("x^0 x^1 x^2 x^3", real=True))
    spatial_momentum = sp.symbols("p^1 p^2 p^3", real=True)
    energy = sp.sqrt(alp_mass**2 + sum(component**2 for component in spatial_momentum))
    momentum_upper = sp.Matrix([energy, *spatial_momentum])
    momentum_lower = metric * momentum_upper
    momentum_squared = sp.simplify((momentum_upper.T * metric * momentum_upper)[0])
    # <0|a(x)|a(p)> = exp(-i p.x) for the relativistic one-particle normalization.
    plane_wave = sp.exp(-sp.I * (momentum_lower.T * position)[0])
    plane_wave_gradient_lower = sp.Matrix([sp.diff(plane_wave, position[index]) for index in range(4)])
    plane_wave_gradient_upper_at_origin = (metric.inv() * plane_wave_gradient_lower).subs(
        {coordinate: 0 for coordinate in position})

    # Generic lower-index derivatives of the ALP field and of the diagonal
    # entries of Phi, and a generic diagonal axial source A_mu.
    alp_gradient = sp.symbols("d_0a d_1a d_2a d_3a")
    meson_gradient = [sp.diag(*sp.symbols(f"d_{index}Phi_uu d_{index}Phi_dd d_{index}Phi_ss"))
                      for index in range(4)]
    axial_source = [sp.diag(*sp.symbols(f"A_{index}^uu A_{index}^dd A_{index}^ss"))
                    for index in range(4)]
    derivative_sources = sp.symbols("U_u U_d U_s")
    derivative_source_matrix = sp.diag(*derivative_sources)

    # u_mu to linear order in the fields and sources. Only diagonal entries of
    # Phi multiply a diagonal source in the trace, so they are sufficient.
    vielbein_lower = [derivative_source_matrix * alp_gradient[index] + 2 * axial_source[index]
                      - meson_gradient[index] / pion_coupling for index in range(4)]
    # chi_+ and the topological term X_TC^2 depend on the axial source only
    # through the fields, not explicitly, so they do not contribute to
    # dL/dA_mu and are omitted here.
    kinetic_lagrangian = sp.expand(pion_coupling**2 / 4
                                   * sp.trace(lorentz_square(metric, vielbein_lower)))
    current_operator = axial_currents(kinetic_lagrangian, axial_source, metric)
    inverse_metric = metric.inv()
    expected_current_operator = sp.Matrix(4, 3, lambda lorentz, flavor: sp.expand(sum(
        inverse_metric[lorentz, lower] * (pion_coupling**2 * derivative_sources[flavor] * alp_gradient[lower]
                                          - pion_coupling * meson_gradient[lower][flavor, flavor])
        for lower in range(4))))

    # Vacuum-to-ALP element: d_mu a -> -i p_mu; <0|Phi|a(p)> = 0 at O(1/f_a)
    # in the no-mixing basis, so the meson-field term drops.
    alp_substitution = {alp_gradient[index]: plane_wave_gradient_lower[index].subs(
        {coordinate: 0 for coordinate in position}) for index in range(4)}
    alp_substitution.update({meson_gradient[index][flavor, flavor]: 0
                             for index in range(4) for flavor in range(3)})
    alp_matrix_element = current_operator.subs(alp_substitution).applyfunc(sp.simplify)
    alp_matrix_element_contracted = (momentum_lower.T * alp_matrix_element).applyfunc(sp.simplify)

    # ---- Step 2: flavor decomposition of the charge-squared weights.
    charges = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    charge_squared = charges * charges
    identity = sp.eye(3)
    lambda_isovector = sp.diag(1, -1, 0)
    lambda_octet = sp.diag(1, 1, -2) / sp.sqrt(3)
    lambda_singlet = sp.sqrt(sp.Rational(2, 3)) * identity
    generators = {"isovector": lambda_isovector, "octet": lambda_octet, "singlet": lambda_singlet}
    charge_coefficients = {name: sp.nsimplify(sp.trace(charge_squared * generator) / 2)
                           for name, generator in generators.items()}
    charge_identity_coefficient = sp.trace(charge_squared) / 3
    source_coefficients = {name: sp.simplify(sp.trace(derivative_source_matrix * generator) / 2)
                           for name, generator in generators.items()}
    source_identity_coefficient = sp.trace(derivative_source_matrix) / 3
    charge_reconstruction = sum((charge_coefficients[name] * generators[name] for name in generators),
                                sp.zeros(3))
    source_reconstruction = sum((source_coefficients[name] * generators[name] for name in generators),
                                sp.zeros(3)).applyfunc(sp.simplify)
    # sum_f Q_f^2 <0|A^mu_f|a(p)> = -i f_pi^2 p^mu tr(Q^2 U_a), split by generator.
    weighted_source_parts = {name: sp.expand(2 * charge_coefficients[name] * source_coefficients[name])
                             for name in generators}
    weighted_source_trace = sp.expand(sp.trace(charge_squared * derivative_source_matrix))
    weighted_matrix_element = (alp_matrix_element * sp.Matrix(
        [charge_squared[flavor, flavor] for flavor in range(3)])).applyfunc(sp.simplify)
    weighted_matrix_element_parts = {name: (-sp.I * pion_coupling**2 * part * momentum_upper)
                                     for name, part in weighted_source_parts.items()}

    # ---- Step 3: equation-of-motion form in terms of Y_a and Z_a.
    quark_mass_terms = sp.symbols("chi_u chi_d chi_s", positive=True)
    pseudoscalar_sources = sp.symbols("c_u^{AP} c_d^{AP} c_s^{AP}")
    gluon_coupling = sp.Symbol("c_GA")
    topological_mass_term = sp.Symbol("g", positive=True)  # g = m_0^2/3
    unknown_profile = sp.symbols("S_u S_d S_s")
    profile_equations = [
        (quark_mass_terms[flavor] - alp_mass**2) * unknown_profile[flavor]
        + topological_mass_term * sum(unknown_profile)
        + quark_mass_terms[flavor] * pseudoscalar_sources[flavor]
        + topological_mass_term * gluon_coupling for flavor in range(3)]
    (profile_solution,) = sp.linsolve(profile_equations, unknown_profile)
    profile = tuple(sp.simplify(component) for component in profile_solution)
    derivative_source_solved = sp.diag(*[-2 * component / alp_coupling_scale for component in profile])
    mass_matrix = sp.diag(*quark_mass_terms)
    shifted_source = sp.diag(*[pseudoscalar_sources[flavor] + profile[flavor] for flavor in range(3)])
    mass_source_solved = (2 * (shifted_source * mass_matrix + mass_matrix * shifted_source)
                          / alp_coupling_scale).applyfunc(sp.simplify)
    contact_solved = sp.simplify(-2 * topological_mass_term / alp_coupling_scale
                                 * (gluon_coupling + sum(profile)))
    equation_of_motion_residual = (mass_source_solved + 2 * alp_mass**2 * derivative_source_solved
                                   - 2 * contact_solved * identity).applyfunc(sp.simplify)

    mass_sources = sp.symbols("Y_u Y_d Y_s")
    contact_coefficient = sp.Symbol("Z_a")
    # m_a^2 U_ff = Z_a - Y_ff/2, obtained by solving the verified identity for U.
    derivative_sources_from_identity = tuple(sp.solve(
        mass_sources[flavor] + 2 * alp_mass**2 * derivative_sources[flavor] - 2 * contact_coefficient,
        derivative_sources[flavor])[0] for flavor in range(3))
    to_mass_and_contact = dict(zip(derivative_sources, derivative_sources_from_identity))
    alp_matrix_element_mass_form = alp_matrix_element.subs(to_mass_and_contact).applyfunc(sp.simplify)
    alp_matrix_element_contracted_mass_form = alp_matrix_element_contracted.subs(
        to_mass_and_contact).applyfunc(sp.simplify)

    # ---- Step 4: combine with the leading local coefficient of the current product.
    hard_momentum_squared = sp.Symbol("q^2", real=True, nonzero=True)
    euclidean_momentum_squared = sp.Symbol("Q^2", positive=True)

    def hard_form_factor(contracted_elements):
        """2i/(q^2 m_a^2) sum_f Q_f^2 p_sigma <0|A^sigma_f|a(p)>."""
        return sp.simplify(2 * sp.I / (hard_momentum_squared * alp_mass**2) * sum(
            charge_squared[flavor, flavor] * contracted_elements[flavor] for flavor in range(3)))

    hard_form_factor_derivative_form = hard_form_factor(alp_matrix_element_contracted)
    hard_form_factor_mass_form = hard_form_factor(alp_matrix_element_contracted_mass_form)
    hard_form_factor_lagrangian_sources = sp.simplify(hard_form_factor_derivative_form.subs(
        {derivative_sources[flavor]: derivative_source_solved[flavor, flavor] for flavor in range(3)}))
    hard_form_factor_euclidean = hard_form_factor_derivative_form.subs(
        hard_momentum_squared, -euclidean_momentum_squared)

    # ---- Step 5: the results may depend only on these symbols.
    allowed_symbols = {pion_coupling, alp_coupling_scale, alp_mass, hard_momentum_squared,
                       euclidean_momentum_squared, contact_coefficient, gluon_coupling,
                       topological_mass_term, *derivative_sources, *mass_sources,
                       *quark_mass_terms, *pseudoscalar_sources}
    final_form_factors = (hard_form_factor_derivative_form, hard_form_factor_mass_form,
                          hard_form_factor_lagrangian_sources, hard_form_factor_euclidean)
    final_symbols = set().union(*(expression.free_symbols for expression in final_form_factors))

    # ---- Method check only (never used for the ALP): the neutral pion.
    # With Phi = pi lambda_3 and <0|pi(x)|pi(p)> = exp(-i p.x), the same
    # functional derivative gives the textbook pion matrix element.
    pion_substitution = {alp_gradient[index]: 0 for index in range(4)}
    pion_substitution.update({meson_gradient[index][flavor, flavor]:
                              lambda_isovector[flavor, flavor] * (-sp.I * momentum_lower[index])
                              for index in range(4) for flavor in range(3)})
    pion_matrix_element = current_operator.subs(pion_substitution).applyfunc(sp.simplify)
    # The pion mass replaces m_a in the on-shell momentum of this check.
    pion_form_factor = hard_form_factor((momentum_lower.T * pion_matrix_element).applyfunc(sp.simplify))
    return locals()


class AxialCurrentMatrixElementTests(unittest.TestCase):
    """Checks of every step of the derivation."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_plane_wave_derivative_sign(self):
        result = self.result
        self.assertEqual((result["plane_wave_gradient_upper_at_origin"]
                          + sp.I * result["momentum_upper"]).applyfunc(sp.simplify), sp.zeros(4, 1))
        self.assertEqual(sp.simplify(result["momentum_squared"] - result["alp_mass"]**2), 0)

    def test_current_operator_and_alp_matrix_element(self):
        result = self.result
        self.assertEqual((result["current_operator"] - result["expected_current_operator"])
                         .applyfunc(sp.expand), sp.zeros(4, 3))
        for lorentz in range(4):
            for flavor in range(3):
                expected = (-sp.I * result["pion_coupling"]**2 * result["derivative_sources"][flavor]
                            * result["momentum_upper"][lorentz])
                self.assertEqual(sp.simplify(result["alp_matrix_element"][lorentz, flavor] - expected), 0)
        for flavor in range(3):
            expected = (-sp.I * result["pion_coupling"]**2 * result["alp_mass"]**2
                        * result["derivative_sources"][flavor])
            self.assertEqual(sp.simplify(result["alp_matrix_element_contracted"][flavor] - expected), 0)

    def test_flavor_weights_and_reconstruction(self):
        result = self.result
        self.assertEqual(result["charge_squared"],
                         sp.diag(sp.Rational(4, 9), sp.Rational(1, 9), sp.Rational(1, 9)))
        coefficients = result["charge_coefficients"]
        self.assertEqual(sp.simplify(coefficients["isovector"] - sp.Rational(1, 6)), 0)
        self.assertEqual(sp.simplify(coefficients["octet"] - 1 / (6 * sp.sqrt(3))), 0)
        self.assertEqual(sp.simplify(coefficients["singlet"] - 2 / (3 * sp.sqrt(6))), 0)
        self.assertEqual(result["charge_identity_coefficient"], sp.Rational(2, 9))
        for generator in result["generators"].values():
            self.assertEqual(sp.simplify(sp.trace(generator * generator)), 2)
        self.assertEqual((result["charge_reconstruction"] - result["charge_squared"])
                         .applyfunc(sp.simplify), sp.zeros(3))
        self.assertEqual((result["source_reconstruction"] - result["derivative_source_matrix"])
                         .applyfunc(sp.simplify), sp.zeros(3))
        up, down, strange = result["derivative_sources"]
        parts = result["weighted_source_parts"]
        self.assertEqual(sp.expand(parts["isovector"] - (up - down) / 6), 0)
        self.assertEqual(sp.expand(parts["octet"] - (up + down - 2 * strange) / 18), 0)
        self.assertEqual(sp.expand(parts["singlet"] - 2 * (up + down + strange) / 9), 0)
        self.assertEqual(sp.expand(sum(parts.values()) - result["weighted_source_trace"]), 0)
        summed_parts = sum(result["weighted_matrix_element_parts"].values(), sp.zeros(4, 1))
        self.assertEqual((summed_parts - result["weighted_matrix_element"]).applyfunc(sp.simplify),
                         sp.zeros(4, 1))

    def test_equation_of_motion_and_mass_form(self):
        result = self.result
        for equation in result["profile_equations"]:
            self.assertEqual(sp.simplify(equation.subs(
                dict(zip(result["unknown_profile"], result["profile"])))), 0)
        self.assertEqual(result["equation_of_motion_residual"], sp.zeros(3))
        mass, coupling = result["alp_mass"], result["pion_coupling"]
        for flavor in range(3):
            combination = result["contact_coefficient"] - result["mass_sources"][flavor] / 2
            self.assertEqual(sp.simplify(result["derivative_sources_from_identity"][flavor]
                                         - combination / mass**2), 0)
            self.assertEqual(sp.simplify(result["alp_matrix_element_contracted_mass_form"][flavor]
                                         + sp.I * coupling**2 * combination), 0)
            for lorentz in range(4):
                expected = -sp.I * coupling**2 * combination * result["momentum_upper"][lorentz] / mass**2
                self.assertEqual(sp.simplify(
                    result["alp_matrix_element_mass_form"][lorentz, flavor] - expected), 0)

    def test_hard_form_factor_and_symbol_content(self):
        result = self.result
        coupling, hard_squared = result["pion_coupling"], result["hard_momentum_squared"]
        expected = 2 * coupling**2 * result["weighted_source_trace"] / hard_squared
        self.assertEqual(sp.simplify(result["hard_form_factor_derivative_form"] - expected), 0)
        self.assertNotIn(result["alp_mass"], result["hard_form_factor_derivative_form"].free_symbols)
        self.assertEqual(sp.simplify(result["hard_form_factor_euclidean"] + 2 * coupling**2
                                     * result["weighted_source_trace"]
                                     / result["euclidean_momentum_squared"]), 0)
        weighted_mass_form = sum(result["charge_squared"][flavor, flavor]
                                 * (result["contact_coefficient"] - result["mass_sources"][flavor] / 2)
                                 for flavor in range(3))
        self.assertEqual(sp.simplify(result["hard_form_factor_mass_form"] - 2 * coupling**2
                                     * weighted_mass_form / (hard_squared * result["alp_mass"]**2)), 0)
        self.assertTrue(result["final_symbols"] <= result["allowed_symbols"])
        for symbol in result["final_symbols"]:
            self.assertNotIn("eta", symbol.name.lower())
        # The solved Lagrangian-source form stays finite at m_a = 0.
        massless = sp.simplify(result["hard_form_factor_lagrangian_sources"].subs(result["alp_mass"], 0))
        self.assertTrue(massless.is_finite is not False and not massless.has(sp.zoo, sp.nan))

    def test_neutral_pion_method_check(self):
        result = self.result
        coupling = result["pion_coupling"]
        for lorentz in range(4):
            for flavor in range(3):
                expected = (sp.I * coupling * result["lambda_isovector"][flavor, flavor]
                            * result["momentum_upper"][lorentz])
                self.assertEqual(sp.simplify(result["pion_matrix_element"][lorentz, flavor] - expected), 0)
        self.assertEqual(sp.simplify(result["pion_form_factor"]
                                     + 2 * coupling / (3 * result["hard_momentum_squared"])), 0)


def main():
    result = derive()
    print("LO chiral Lagrangian, O(1/f_a), three light flavors, diagonal sources, no-mixing basis;")
    print("metric (+---); <0|a(x)|a(p)> = exp(-i p.x); f_pi is the LO chiral coupling constant.")
    show("upper-index gradient of the plane wave at the origin", result["plane_wave_gradient_upper_at_origin"])
    show("axial current operators (rows mu upper, columns u,d,s)", result["current_operator"])
    show("<0|A^mu_f|a(p)> (rows mu upper, columns u,d,s)", result["alp_matrix_element"])
    show("p_mu <0|A^mu_f|a(p)>", result["alp_matrix_element_contracted"])
    show("charge-squared weights", result["charge_squared"])
    show("coefficients of Q^2 on lambda_3, lambda_8, lambda_0", sp.Matrix(list(result["charge_coefficients"].values())))
    show("coefficients of U_a on lambda_3, lambda_8, lambda_0", sp.Matrix(list(result["source_coefficients"].values())))
    for name, part in result["weighted_source_parts"].items():
        show(f"{name} part of tr(Q^2 U_a)", part)
    show("tr(Q^2 U_a)", result["weighted_source_trace"])
    show("solved profile S_u", result["profile"][0])
    show("Y_a + 2 m_a^2 U_a - 2 Z_a 1 on the solved profile", result["equation_of_motion_residual"])
    show("U_f from the identity", sp.Matrix(result["derivative_sources_from_identity"]))
    show("<0|A^mu_f|a(p)> in terms of Y_a, Z_a", result["alp_matrix_element_mass_form"])
    show("p_mu <0|A^mu_f|a(p)> in terms of Y_a, Z_a", result["alp_matrix_element_contracted_mass_form"])
    show("hard form factor, U_a form", result["hard_form_factor_derivative_form"])
    show("hard form factor, Euclidean q^2=-Q^2", result["hard_form_factor_euclidean"])
    show("hard form factor, Y_a and Z_a form", result["hard_form_factor_mass_form"])
    show("hard form factor, Lagrangian sources", result["hard_form_factor_lagrangian_sources"])
    print("The U_a form has no 1/m_a^2; the Y_a, Z_a form shows an explicit 1/m_a^2 that cancels")
    print("once Y_a and Z_a are evaluated on the solved profile.")
    show("symbols in the final results", sp.Matrix(sorted(result["final_symbols"], key=str)))
    print("Method check only, not an ALP input: neutral pion with Phi = pi lambda_3.")
    show("<0|A^mu_f|pi(p)>", result["pion_matrix_element"])
    show("pion form factor from the same local coefficient", result["pion_form_factor"])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
