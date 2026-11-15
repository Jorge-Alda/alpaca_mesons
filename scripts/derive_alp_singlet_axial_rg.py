"""Derive the scale dependence of the flavor-singlet axial-current matrix element.

Run .venv/bin/python scripts/derive_alp_singlet_axial_rg.py [--test].

Units and conventions: natural units, scales mu and mu_0 in GeV. MS-bar
scheme, massless QCD with n_f active flavors, coupling a = alpha_s/(4 pi) as
in Larin:1993tq. The singlet axial current sum_f qbar_f gamma_mu gamma5 q_f is
renormalized as in Larin:1993tq, including the finite constant Z_5^s, so that
the anomaly relation keeps its one-loop form. Inputs:

- memory/extracts/Larin:1993tq/singlet-axial-current-renormalization.md:
  anomalous dimension of the singlet current,
  mu^2 d ln J / d mu^2 = -6 C_F n_f a^2
      + (-142/3 C_F C_A n_f + 4/3 C_F n_f^2 + 18 C_F^2 n_f) a^3,
  and vanishing anomalous dimension of the non-singlet currents.
- memory/extracts/Agaev:2014wna/singlet-decay-constant-evolution.md:
  the linear solution 1 + (2 n_f/(pi beta_0)) [alpha_s(mu) - alpha_s(mu_0)]
  and mu d ln f/d mu = -4 n_f (alpha_s/(2 pi))^2.
- Beta function (stated assumption of the task, standard MS-bar values for
  C_A = 3): mu^2 da/dmu^2 = -beta_0 a^2 - beta_1 a^3, beta_0 = 11 - 2 n_f/3,
  beta_1 = 102 - 38 n_f/3.

Because the current is multiplicatively renormalized, every matrix element of
it has the same scale dependence; the script works with the ratio of the
singlet axial-current matrix element at scale mu to that at scale mu_0.

Approximations (all explicit):
- Leading running order: anomalous dimension truncated at a^2 and beta
  function truncated at a^2. Within this truncation the solution is exact.
- Next order: anomalous dimension to a^3 and beta function to a^3. The
  logarithm of the ratio is obtained as a power series in a(mu) and a(mu_0),
  dropping terms of order a^3. Here a(mu) is treated as a quantity obeying
  the two-loop equation; no explicit two-loop solution a(mu) is used, so this
  is a series check, not an exact two-loop solution.
- No quark-mass thresholds: n_f is fixed. Matching of the singlet current
  across thresholds is not treated here.
- The numerical table at the end only illustrates the size of the effect. It
  is not a benchmark and not an ALP input.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp


def show(label, expression, symbol_names=None):
    """Print expressions through the standard SymPy LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression, symbol_names=symbol_names or {})}\n")


def scale_derivative(expression, scale):
    """Return mu^2 d/dmu^2 of an expression, written as (mu/2) d/dmu."""
    return scale / 2 * sp.diff(expression, scale)


def integrate_between_couplings(integrand, integration_coupling, initial_coupling, final_coupling):
    """Integrate a polynomial in the coupling from its value at mu_0 to its value at mu."""
    return sp.integrate(integrand, (integration_coupling, initial_coupling, final_coupling))


@lru_cache(maxsize=1)
def derive():
    """Carry out the six steps of the derivation and return all intermediate objects."""
    scale = sp.Symbol("mu", positive=True)
    reference_scale = sp.Symbol("mu_0", positive=True)
    flavor_number = sp.Symbol("n_f", positive=True)
    casimir_fundamental = sp.Symbol("C_F", positive=True)
    casimir_adjoint = sp.Symbol("C_A", positive=True)
    beta_zero = sp.Symbol("beta_0", positive=True)
    beta_one = sp.Symbol("beta_1", real=True)
    # a = alpha_s/(4 pi) at the two scales, as independent symbols.
    coupling = sp.Symbol("coupling_at_scale", positive=True)
    coupling_reference = sp.Symbol("coupling_at_reference_scale", positive=True)
    strong_coupling = sp.Symbol("strong_coupling_at_scale", positive=True)
    strong_coupling_reference = sp.Symbol("strong_coupling_at_reference_scale", positive=True)
    symbol_names = {coupling: r"a(\mu)", coupling_reference: r"a(\mu_0)",
                    strong_coupling: r"\alpha_s(\mu)", strong_coupling_reference: r"\alpha_s(\mu_0)"}
    integration_coupling = sp.Dummy("integration_coupling", positive=True)
    quantum_chromodynamics = {casimir_fundamental: sp.Rational(4, 3), casimir_adjoint: 3}
    beta_values = {beta_zero: 11 - 2 * flavor_number / 3, beta_one: 102 - 38 * flavor_number / 3}
    to_strong_coupling = {coupling: strong_coupling / (4 * sp.pi),
                          coupling_reference: strong_coupling_reference / (4 * sp.pi)}

    # Inputs: coefficients of a^2 and a^3 in the singlet anomalous dimension.
    anomalous_dimension_second_order = -6 * casimir_fundamental * flavor_number
    anomalous_dimension_third_order = (-sp.Rational(142, 3) * casimir_fundamental * casimir_adjoint * flavor_number
                                       + sp.Rational(4, 3) * casimir_fundamental * flavor_number**2
                                       + 18 * casimir_fundamental**2 * flavor_number)

    def singlet_anomalous_dimension(value, highest_power):
        """Singlet anomalous dimension truncated at the given power of a."""
        result = anomalous_dimension_second_order * value**2
        if highest_power >= 3:
            result += anomalous_dimension_third_order * value**3
        return result

    def beta_function(value, highest_power):
        """mu^2 da/dmu^2 truncated at the given power of a."""
        result = -beta_zero * value**2
        if highest_power >= 3:
            result += -beta_one * value**3
        return result

    # ---- Step 1: leading running order, exact within the truncation.
    logarithm_of_scales = sp.log(scale**2 / reference_scale**2)
    one_loop_coupling = coupling_reference / (1 + beta_zero * coupling_reference * logarithm_of_scales)
    one_loop_coupling_residual = sp.simplify(scale_derivative(one_loop_coupling, scale)
                                             - beta_function(one_loop_coupling, 2))
    # d ln(ratio)/da = anomalous dimension / beta function; a^2 cancels exactly.
    leading_integrand = sp.cancel(singlet_anomalous_dimension(integration_coupling, 2)
                                  / beta_function(integration_coupling, 2))
    leading_log_ratio = sp.factor(integrate_between_couplings(
        leading_integrand, integration_coupling, coupling_reference, coupling))
    leading_ratio = sp.exp(leading_log_ratio)
    leading_ratio_of_scale = leading_ratio.subs(coupling, one_loop_coupling)
    leading_evolution_residual = sp.simplify(
        scale_derivative(sp.log(leading_ratio_of_scale), scale)
        - singlet_anomalous_dimension(one_loop_coupling, 2))
    leading_ratio_at_reference_scale = sp.simplify(leading_ratio_of_scale.subs(scale, reference_scale))

    # ---- Step 2: alpha_s form, first-order expansion, comparison with Agaev:2014wna.
    leading_log_ratio_strong_coupling = sp.factor(
        leading_log_ratio.subs(to_strong_coupling).subs(quantum_chromodynamics))
    expansion_counter = sp.Dummy("expansion_counter")
    leading_ratio_first_order = sp.series(
        sp.exp(expansion_counter * leading_log_ratio_strong_coupling), expansion_counter, 0, 2
    ).removeO().subs(expansion_counter, 1)
    linear_form_of_agaev = 1 + 2 * flavor_number / (sp.pi * beta_zero) * (
        strong_coupling - strong_coupling_reference)
    first_neglected_term = sp.factor(sp.series(
        sp.exp(expansion_counter * leading_log_ratio_strong_coupling), expansion_counter, 0, 3
    ).removeO().subs(expansion_counter, 1) - leading_ratio_first_order)
    # mu d/dmu = 2 mu^2 d/dmu^2, with the explicit one-loop alpha_s(mu).
    one_loop_strong_coupling = sp.simplify(4 * sp.pi * one_loop_coupling.subs(
        coupling_reference, strong_coupling_reference / (4 * sp.pi)))
    leading_log_derivative = sp.simplify(scale * sp.diff(
        leading_log_ratio_strong_coupling.subs(strong_coupling, one_loop_strong_coupling), scale))
    expected_log_derivative = -flavor_number * one_loop_strong_coupling**2 / sp.pi**2
    log_derivative_residual = sp.simplify(leading_log_derivative - expected_log_derivative)
    log_derivative_in_agaev_form = -4 * flavor_number * (one_loop_strong_coupling / (2 * sp.pi))**2
    # Derivative of the linear form itself, in units of the value at mu_0.
    linear_form_derivative_residual = sp.simplify(scale * sp.diff(
        linear_form_of_agaev.subs(strong_coupling, one_loop_strong_coupling), scale)
        - expected_log_derivative)

    # ---- Step 3: next order, as a series in the coupling.
    next_integrand = sp.series(sp.cancel(singlet_anomalous_dimension(integration_coupling, 3)
                                         / beta_function(integration_coupling, 3)),
                               integration_coupling, 0, 2).removeO()
    next_log_ratio = sp.expand(integrate_between_couplings(
        next_integrand, integration_coupling, coupling_reference, coupling))
    expected_next_log_ratio = (
        -anomalous_dimension_second_order / beta_zero * (coupling - coupling_reference)
        - sp.Rational(1, 2) * (anomalous_dimension_third_order / beta_zero
                               - anomalous_dimension_second_order * beta_one / beta_zero**2)
        * (coupling**2 - coupling_reference**2))
    # Chain rule: mu^2 d/dmu^2 acts on a(mu) through the two-loop beta function.
    next_evolution_residual = sp.expand(sp.diff(next_log_ratio, coupling) * beta_function(coupling, 3)
                                        - singlet_anomalous_dimension(coupling, 3))
    next_evolution_residual_coefficients = sp.Poly(next_evolution_residual, coupling).all_coeffs()[::-1]
    next_log_ratio_three_flavors = sp.factor(next_log_ratio.subs(quantum_chromodynamics)
                                             .subs(beta_values).subs(flavor_number, 3))

    # ---- Step 4: non-singlet currents have vanishing anomalous dimension.
    non_singlet_log_ratio = integrate_between_couplings(
        sp.cancel(sp.Integer(0) / beta_function(integration_coupling, 3)),
        integration_coupling, coupling_reference, coupling)
    non_singlet_ratio = sp.exp(non_singlet_log_ratio)

    # ---- Step 5: separation of the singlet part of the charge-squared matrix.
    charges = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    charge_squared = charges * charges
    identity = sp.eye(3)
    singlet_coefficient = sp.trace(charge_squared) / 3
    traceless_charge_squared = charge_squared - singlet_coefficient * identity
    generic_diagonal_matrix = sp.diag(*sp.symbols("M_u M_d M_s"))
    weighted_trace = sp.expand(sp.trace(charge_squared * generic_diagonal_matrix))
    traceless_weighted_trace = sp.expand(sp.trace(traceless_charge_squared * generic_diagonal_matrix))
    singlet_weighted_trace = sp.expand(singlet_coefficient * sp.trace(generic_diagonal_matrix))
    # Only the term with the trace of M multiplies the singlet current, so
    # only that term receives the ratio derived above.
    evolved_weighted_trace = traceless_weighted_trace + leading_ratio * singlet_weighted_trace

    # ---- Step 6: illustration of the size of the effect only.
    illustration_flavor_number = 3
    illustration_reference_scale = 1  # GeV
    illustration_reference_coupling = sp.Rational(1, 2)
    illustration_substitution = {flavor_number: illustration_flavor_number,
                                 reference_scale: illustration_reference_scale,
                                 strong_coupling_reference: illustration_reference_coupling}
    illustration_beta_zero = beta_values[beta_zero].subs(flavor_number, illustration_flavor_number)
    illustration_rows = []
    for illustration_scale in (2, 10, 100):  # GeV
        coupling_value = one_loop_strong_coupling.subs(beta_zero, illustration_beta_zero).subs(
            illustration_substitution).subs(scale, illustration_scale)
        point = {strong_coupling: coupling_value, beta_zero: illustration_beta_zero,
                 **illustration_substitution}
        illustration_rows.append([
            illustration_scale, sp.N(coupling_value, 6),
            sp.N(sp.exp(leading_log_ratio_strong_coupling).subs(point), 6),
            sp.N(linear_form_of_agaev.subs(point), 6)])
    illustration_table = sp.Matrix(illustration_rows)
    return locals()


class SingletAxialCurrentEvolutionTests(unittest.TestCase):
    """Checks of every step of the derivation."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_leading_order_differential_equations_hold_exactly(self):
        result = self.result
        self.assertEqual(result["one_loop_coupling_residual"], 0)
        self.assertEqual(result["leading_evolution_residual"], 0)
        self.assertEqual(result["leading_ratio_at_reference_scale"], 1)
        expected = (6 * result["casimir_fundamental"] * result["flavor_number"] / result["beta_zero"]
                    * (result["coupling"] - result["coupling_reference"]))
        self.assertEqual(sp.simplify(result["leading_log_ratio"] - expected), 0)

    def test_linear_form_and_derivative_of_agaev(self):
        result = self.result
        self.assertEqual(sp.simplify(result["leading_ratio_first_order"]
                                     - result["linear_form_of_agaev"]), 0)
        self.assertEqual(result["log_derivative_residual"], 0)
        self.assertEqual(sp.simplify(result["log_derivative_in_agaev_form"]
                                     - result["expected_log_derivative"]), 0)
        self.assertEqual(result["linear_form_derivative_residual"], 0)
        # Larin:1993tq leading coefficient for C_F = 4/3 is -8 n_f.
        self.assertEqual(result["anomalous_dimension_second_order"].subs(
            result["quantum_chromodynamics"]), -8 * result["flavor_number"])
        # The first neglected term is quadratic in the difference of couplings.
        difference = result["strong_coupling"] - result["strong_coupling_reference"]
        self.assertEqual(sp.simplify(result["first_neglected_term"] / difference**2).free_symbols,
                         {result["flavor_number"], result["beta_zero"]})

    def test_next_order_residual_starts_at_fourth_power(self):
        result = self.result
        self.assertEqual(sp.expand(result["next_log_ratio"] - result["expected_next_log_ratio"]), 0)
        coefficients = result["next_evolution_residual_coefficients"]
        self.assertEqual(len(coefficients), 5)
        for power in range(4):
            self.assertEqual(sp.simplify(coefficients[power]), 0)
        expected_fourth = (result["anomalous_dimension_third_order"] / result["beta_zero"]
                           - result["anomalous_dimension_second_order"] * result["beta_one"]
                           / result["beta_zero"]**2) * result["beta_one"]
        self.assertEqual(sp.simplify(coefficients[4] - expected_fourth), 0)
        self.assertNotEqual(sp.simplify(coefficients[4]), 0)
        # The leading-order logarithm alone leaves a residual already at a^3.
        coupling = result["coupling"]
        leading_only = sp.expand(sp.diff(result["leading_log_ratio"], coupling)
                                 * result["beta_function"](coupling, 3)
                                 - result["singlet_anomalous_dimension"](coupling, 3))
        self.assertNotEqual(sp.simplify(leading_only.coeff(coupling, 3)), 0)
        self.assertEqual(sp.simplify(leading_only.coeff(coupling, 2)), 0)

    def test_non_singlet_ratio_is_one(self):
        self.assertEqual(self.result["non_singlet_log_ratio"], 0)
        self.assertEqual(self.result["non_singlet_ratio"], 1)

    def test_flavor_projection(self):
        result = self.result
        self.assertEqual(result["singlet_coefficient"], sp.Rational(2, 9))
        self.assertEqual(sp.trace(result["traceless_charge_squared"]), 0)
        self.assertEqual(result["traceless_charge_squared"]
                         + result["singlet_coefficient"] * result["identity"], result["charge_squared"])
        self.assertEqual(sp.expand(result["weighted_trace"] - result["traceless_weighted_trace"]
                                   - result["singlet_weighted_trace"]), 0)
        # Without scale dependence the separated form reduces to the full trace.
        self.assertEqual(sp.expand(result["evolved_weighted_trace"].subs(
            result["coupling"], result["coupling_reference"]) - result["weighted_trace"]), 0)
        # A matrix proportional to the identity has no traceless part.
        proportional = {symbol: 1 for symbol in result["generic_diagonal_matrix"].free_symbols}
        self.assertEqual(result["traceless_weighted_trace"].subs(proportional), 0)

    def test_illustration_is_consistent(self):
        result = self.result
        self.assertEqual(result["illustration_beta_zero"], 9)
        previous_coupling = result["illustration_reference_coupling"]
        for _, coupling_value, exponential_ratio, linear_ratio in result["illustration_table"].tolist():
            self.assertLess(coupling_value, previous_coupling)
            previous_coupling = coupling_value
            self.assertLess(exponential_ratio, 1)
            self.assertLess(abs(exponential_ratio - linear_ratio), 0.01)
            self.assertAlmostEqual(float(linear_ratio),
                                   1 + 2 * 3 / (float(sp.pi) * 9) * (float(coupling_value) - 0.5), places=5)


def main():
    result = derive()

    def show_result(label, expression):
        show(label, expression, result["symbol_names"])

    print("MS-bar, massless QCD with n_f active flavors, a = alpha_s/(4 pi) (Larin:1993tq);")
    print("singlet axial current renormalized with the finite constant Z_5^s; no thresholds.")
    show_result("coefficient of a^2 in the singlet anomalous dimension", result["anomalous_dimension_second_order"])
    show_result("coefficient of a^3 in the singlet anomalous dimension", result["anomalous_dimension_third_order"])
    print("Step 1: leading running order (anomalous dimension and beta function truncated at a^2).")
    show_result("one-loop a(mu)", result["one_loop_coupling"])
    show_result("mu^2 da/dmu^2 + beta_0 a^2 on the one-loop a(mu)", result["one_loop_coupling_residual"])
    show_result("logarithm of the ratio of the singlet axial-current matrix element at mu to that at mu_0",
         result["leading_log_ratio"])
    show_result("ratio of the singlet axial-current matrix element at mu to that at mu_0", result["leading_ratio"])
    show_result("mu^2 d ln(ratio)/dmu^2 minus the anomalous dimension, on the one-loop a(mu)",
         result["leading_evolution_residual"])
    show_result("ratio at mu = mu_0", result["leading_ratio_at_reference_scale"])
    print("Step 2: alpha_s form for C_F = 4/3.")
    show_result("logarithm of the ratio in terms of alpha_s", result["leading_log_ratio_strong_coupling"])
    show_result("ratio truncated at first order in the exponent", result["leading_ratio_first_order"])
    show_result("linear form of Agaev:2014wna", result["linear_form_of_agaev"])
    show_result("first term neglected by the linear form", result["first_neglected_term"])
    show_result("one-loop alpha_s(mu)", result["one_loop_strong_coupling"])
    show_result("mu d ln(ratio)/dmu, with alpha_s(mu) the one-loop coupling",
                -result["flavor_number"] * result["strong_coupling"]**2 / sp.pi**2)
    show_result("mu d ln(ratio)/dmu + n_f alpha_s(mu)^2/pi^2", result["log_derivative_residual"])
    print("Step 3: next order (anomalous dimension and beta function to a^3). Series check in a(mu):")
    print("a(mu) is only assumed to obey the two-loop equation; this is not an exact two-loop solution.")
    show_result("logarithm of the ratio, terms of order a^3 dropped", result["next_log_ratio"])
    show_result("the same for C_F = 4/3, C_A = 3, n_f = 3", result["next_log_ratio_three_flavors"])
    show_result("mu^2 d ln(ratio)/dmu^2 minus the anomalous dimension (starts at a^4)",
         sp.factor(result["next_evolution_residual"]))
    print("Step 4: non-singlet currents.")
    show_result("ratio of the non-singlet axial-current matrix element at mu to that at mu_0",
         result["non_singlet_ratio"])
    print("Step 5: separation of the singlet part of the charge-squared matrix.")
    show_result("Q^2", result["charge_squared"])
    show_result("trace of Q^2 divided by the number of flavors", result["singlet_coefficient"])
    show_result("traceless part of Q^2", result["traceless_charge_squared"])
    show_result("trace of Q^2 M", result["weighted_trace"])
    show_result("trace of (traceless part of Q^2) M", result["traceless_weighted_trace"])
    show_result("(2/9) trace of M", result["singlet_weighted_trace"])
    show_result("trace of Q^2 M with the singlet term multiplied by the leading-order ratio",
         result["evolved_weighted_trace"])
    print("Step 6: ILLUSTRATION of the size of the effect only; not a benchmark, not an ALP input.")
    print("n_f = 3 fixed at all scales (no thresholds), mu_0 = 1 GeV, alpha_s(mu_0) = 0.5, one-loop running.")
    show_result("columns: mu in GeV, alpha_s(mu), leading-order ratio (exponential), first-order expansion",
         result["illustration_table"])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
