"""Derive the heavy-flavor treatment of the short-distance ALP two-photon coefficient.

Run .venv/bin/python scripts/derive_alp_heavy_flavor_threshold.py [--test].

Units and conventions: natural units, all masses and scales in GeV. Metric
(+---), epsilon^{0123}=+1, gamma5 = i gamma0 gamma1 gamma2 gamma3. Fourier
convention i integral d^4x exp(i q.x) T{J^mu(x/2) J^nu(-x/2)}, the
electromagnetic current excludes e, q=(q1-q2)/2. MS-bar scheme for the strong
coupling and the quark masses. Thresholds are symbols; no numerical quark
mass is used anywhere.

Inputs:
- memory/alp/lepton-pair/alp-hard-current-ope.md: massless free-quark
  derivation and its sign, reproduced here at m=0.
- memory/extracts/Chetyrkin:1997sg/alphas-threshold-matching.md: matching
  alpha_s^{(n_f-1)}/alpha_s^{(n_f)} with alpha_s/pi, mu_h=m_h(m_h), constant 11/72.
- memory/extracts/Chetyrkin:1997un/alphas-decoupling-constant.md: beta
  function mu^2 d(alpha_s/pi)/dmu^2 = -beta_0 (alpha_s/pi)^2 - beta_1 (alpha_s/pi)^3,
  beta_0=(11-2n_f/3)/4, beta_1=(102-38n_f/3)/16, and the matching written with
  the running mass m_h(mu) (coefficient -11/24).
- memory/extracts/Larin:1993tq/singlet-axial-current-renormalization.md:
  singlet anomalous dimension -6 C_F n_f (alpha_s/(4 pi))^2 in mu^2 d/dmu^2.
- scripts/derive_alp_singlet_axial_rg.py: leading-order ratio of the singlet
  axial-current matrix element at fixed n_f, extended here.

Part A. One free quark flavor of mass m, propagator
S_F(k)=i(kslash+m)/(k^2-m^2+i0), leading order in alpha_s, local term without
derivatives of the uncontracted quark fields, fixed soft external momentum.
Approximation: the same local expansion as in the massless record. For
|q^2| << m^2 the heavy quark is not a light field of the low-energy theory,
so the factor obtained there is the free-quark expression only.

Part B. Strong-coupling matching at a heavy-quark threshold. The scale
mu_h=m_h(m_h) does not depend on mu. The series in alpha_s/pi is truncated
at relative order (alpha_s/pi)^2.

Part C. Leading-order evolution of the singlet axial-current matrix element
with 3, 4, 5 active flavors. Approximations and assumptions (all explicit):
- one-loop running of alpha_s in each interval;
- alpha_s continuous at each threshold (the constant 11/72 of Part B is
  neglected; it is of higher order than the one-loop running);
- ASSUMPTION: the singlet axial-current matrix element is continuous across
  each threshold. No matching coefficient for the singlet current is
  available in the sources read, and Larin:1994va does not establish
  decoupling of the heavy quark from this current;
- the reference scale mu_0 lies below the charm threshold;
- the top-quark region is outside the supported domain, because electroweak
  effects are not treated. The symbol for the top threshold only marks the
  nominal upper end of the five-flavor interval.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp
from sympy.physics.matrices import mgamma


def show(label, expression, symbol_names=None):
    """Print expressions through the standard SymPy LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression, symbol_names=symbol_names or {})}\n")


def scale_derivative(expression, scale):
    """Return mu^2 d/dmu^2 of an expression, written as (mu/2) d/dmu."""
    return scale / 2 * sp.diff(expression, scale)


@lru_cache(maxsize=1)
def derive_massive_coefficient():
    """Part A: project both Wick contractions with a massive quark onto the Clifford basis."""
    metric = sp.diag(1, -1, -1, -1)
    gamma = tuple(mgamma(index) for index in range(4))
    gamma5 = mgamma(5)
    identity = sp.eye(4)
    quark_mass = sp.Symbol("m", positive=True)
    hard_momentum = sp.Matrix(sp.symbols("q_0 q_1 q_2 q_3", real=True))
    hard_lower = metric * hard_momentum
    hard_squared = (hard_momentum.T * metric * hard_momentum)[0]
    slash = sum((hard_lower[index] * gamma[index] for index in range(4)), sp.zeros(4))

    # The two single contractions are gamma^mu S_F(q) gamma^nu and
    # gamma^nu S_F(-q) gamma^mu. With S_F(k)=i(kslash+m)/(k^2-m^2+i0) and the
    # Fourier prefactor i, the overall factor is i*i/(q^2-m^2) = -1/(q^2-m^2)
    # times the matrix formed below.
    contraction_sum = tuple(tuple(sp.expand(
        gamma[first] * (slash + quark_mass * identity) * gamma[second]
        + gamma[second] * (-slash + quark_mass * identity) * gamma[first])
        for second in range(4)) for first in range(4))
    denominator = hard_squared - quark_mass**2

    # Full Clifford basis: 1, gamma5, gamma^sigma, gamma^sigma gamma5,
    # sigma^{rho sigma} = (i/2)[gamma^rho, gamma^sigma] with rho < sigma.
    clifford_basis = [("scalar", identity), ("pseudoscalar", gamma5)]
    clifford_basis += [(f"vector {index}", gamma[index]) for index in range(4)]
    clifford_basis += [(f"axial {index}", gamma[index] * gamma5) for index in range(4)]
    clifford_basis += [(f"tensor {first}{second}",
                        sp.I / 2 * (gamma[first] * gamma[second] - gamma[second] * gamma[first]))
                       for first in range(4) for second in range(first + 1, 4)]
    basis_labels = tuple(label for label, _ in clifford_basis)
    basis_matrices = tuple(matrix for _, matrix in clifford_basis)
    basis_norms = tuple(sp.trace(matrix * matrix) for matrix in basis_matrices)
    basis_overlaps = sp.Matrix(16, 16, lambda first, second:
                               sp.trace(basis_matrices[first] * basis_matrices[second]))

    def project(matrix):
        """Coefficients of a 4x4 matrix in the Clifford basis, from traces."""
        return tuple(sp.simplify(sp.trace(basis_matrix * matrix) / norm)
                     for basis_matrix, norm in zip(basis_matrices, basis_norms))

    # Coefficients of the full Fourier-transformed contraction, prefactor included.
    clifford_coefficients = tuple(tuple(tuple(
        sp.simplify(-coefficient / denominator) for coefficient in project(contraction_sum[first][second]))
        for second in range(4)) for first in range(4))
    reconstructed = tuple(tuple(sum(
        (clifford_coefficients[first][second][position] * basis_matrices[position]
         for position in range(16)), sp.zeros(4))
        for second in range(4)) for first in range(4))

    # Expected: (2i/(q^2-m^2)) epsilon^{mu nu rho sigma} q_rho gamma_sigma gamma5
    # and -(2m/(q^2-m^2)) g^{mu nu}; gamma_sigma = g_{sigma sigma} gamma^sigma.
    def expected_coefficient(first, second, position):
        label = basis_labels[position]
        if label == "scalar":
            return -2 * quark_mass * metric[first, second] / denominator
        if label.startswith("axial"):
            axial = int(label.split()[1])
            return 2 * sp.I * metric[axial, axial] * sum(
                sp.LeviCivita(first, second, momentum, axial) * hard_lower[momentum]
                for momentum in range(4)) / denominator
        return sp.Integer(0)

    expected_coefficients = tuple(tuple(tuple(
        sp.simplify(expected_coefficient(first, second, position)) for position in range(16))
        for second in range(4)) for first in range(4))

    # Terms linear in the mass, and their antisymmetric part in mu, nu.
    mass_linear = tuple(tuple(sp.diff(contraction_sum[first][second], quark_mass)
                              for second in range(4)) for first in range(4))
    mass_linear_antisymmetric = tuple(tuple(
        (mass_linear[first][second] - mass_linear[second][first]).applyfunc(sp.simplify)
        for second in range(4)) for first in range(4))

    # Scalar coefficient functions multiplying the two surviving structures.
    hard_squared_symbol = sp.Symbol("q^2", real=True)
    axial_coefficient = 2 * sp.I / (hard_squared_symbol - quark_mass**2)
    scalar_density_coefficient = -2 * quark_mass / (hard_squared_symbol - quark_mass**2)
    massless_axial_coefficient = 2 * sp.I / hard_squared_symbol
    heavy_flavor_factor = sp.simplify(axial_coefficient / massless_axial_coefficient)

    # Spacelike momentum: q^2 = -|q^2|.
    spacelike_magnitude = sp.Symbol("spacelike_magnitude", positive=True)
    symbol_names = {spacelike_magnitude: r"|q^2|"}
    heavy_flavor_factor_spacelike = sp.simplify(
        heavy_flavor_factor.subs(hard_squared_symbol, -spacelike_magnitude))
    # |q^2| >> m^2: expand in m^2/|q^2|.
    mass_ratio = sp.Symbol("mass_ratio", positive=True)
    large_momentum_series = sp.series(
        heavy_flavor_factor_spacelike.subs(quark_mass, sp.sqrt(mass_ratio * spacelike_magnitude)),
        mass_ratio, 0, 3).removeO().subs(mass_ratio, quark_mass**2 / spacelike_magnitude)
    # |q^2| << m^2: expand in |q^2|/m^2.
    momentum_ratio = sp.Symbol("momentum_ratio", positive=True)
    small_momentum_series = sp.series(
        heavy_flavor_factor_spacelike.subs(spacelike_magnitude, momentum_ratio * quark_mass**2),
        momentum_ratio, 0, 3).removeO().subs(momentum_ratio, spacelike_magnitude / quark_mass**2)
    large_momentum_limit = sp.limit(heavy_flavor_factor_spacelike, spacelike_magnitude, sp.oo)
    small_momentum_limit = sp.limit(heavy_flavor_factor_spacelike, spacelike_magnitude, 0)
    return locals()


@lru_cache(maxsize=1)
def derive_coupling_matching():
    """Part B: renormalization-group consistency of the alpha_s matching and the jump at mu_h."""
    flavor_number = sp.Symbol("n_f", positive=True)
    lighter_flavor_number = flavor_number - 1
    # alpha_s/pi of the theory with n_f flavors, at scale mu.
    coupling = sp.Symbol("coupling_full", positive=True)
    # ln(mu^2/mu_h^2) with mu_h = m_h(m_h) independent of mu.
    threshold_logarithm = sp.Symbol("threshold_logarithm", real=True)
    symbol_names = {coupling: r"\frac{\alpha_s^{(n_f)}(\mu)}{\pi}",
                    threshold_logarithm: r"\ln\frac{\mu^2}{\mu_h^2}"}

    def beta_zero(number):
        return (11 - sp.Rational(2, 3) * number) / 4

    def beta_one(number):
        return (102 - sp.Rational(38, 3) * number) / 16

    def beta_function(value, number):
        """mu^2 d(alpha_s/pi)/dmu^2 through (alpha_s/pi)^3."""
        return -beta_zero(number) * value**2 - beta_one(number) * value**3

    # General form of the matching through relative order (alpha_s/pi)^2.
    first_order_logarithm = sp.Symbol("first_order_logarithm_coefficient")
    first_order_constant = sp.Symbol("first_order_constant")
    second_order_squared_logarithm = sp.Symbol("second_order_squared_logarithm_coefficient")
    second_order_logarithm = sp.Symbol("second_order_logarithm_coefficient")
    second_order_constant = sp.Symbol("second_order_constant")
    general_ratio = (1 + coupling * (first_order_logarithm * threshold_logarithm + first_order_constant)
                     + coupling**2 * (second_order_squared_logarithm * threshold_logarithm**2
                                      + second_order_logarithm * threshold_logarithm
                                      + second_order_constant))
    general_lighter_coupling = sp.expand(general_ratio * coupling)

    # mu^2 d/dmu^2 acts on the explicit logarithm (derivative 1) and on the
    # n_f-flavor coupling through its beta function. The result must equal
    # the beta function of the (n_f - 1)-flavor theory.
    def consistency_residual(lighter_coupling):
        derivative = (sp.diff(lighter_coupling, threshold_logarithm)
                      + sp.diff(lighter_coupling, coupling) * beta_function(coupling, flavor_number))
        return sp.expand(derivative - beta_function(lighter_coupling, lighter_flavor_number))

    general_residual = consistency_residual(general_lighter_coupling)
    residual_second_order = sp.expand(general_residual.coeff(coupling, 2))
    residual_third_order = sp.expand(general_residual.coeff(coupling, 3))
    consistency_equations = (sp.Poly(residual_second_order, threshold_logarithm).all_coeffs()
                             + sp.Poly(residual_third_order, threshold_logarithm).all_coeffs())
    consistency_solution = sp.solve(
        consistency_equations,
        [first_order_logarithm, second_order_squared_logarithm, second_order_logarithm], dict=True)[0]
    derived_first_order_logarithm = sp.simplify(consistency_solution[first_order_logarithm])
    derived_second_order_squared_logarithm = sp.simplify(
        consistency_solution[second_order_squared_logarithm])
    derived_second_order_logarithm_general = sp.simplify(consistency_solution[second_order_logarithm])
    # Chetyrkin:1997sg has no constant at first order.
    derived_second_order_logarithm = derived_second_order_logarithm_general.subs(first_order_constant, 0)

    literature_first_order_logarithm = -sp.Rational(1, 6)
    literature_second_order_squared_logarithm = sp.Rational(1, 36)
    literature_second_order_logarithm = -sp.Rational(19, 24)
    literature_second_order_constant = sp.Rational(11, 72)
    coefficient_differences = (
        derived_first_order_logarithm - literature_first_order_logarithm,
        derived_second_order_squared_logarithm - literature_second_order_squared_logarithm,
        derived_second_order_logarithm - literature_second_order_logarithm)

    literature_values = {first_order_logarithm: literature_first_order_logarithm,
                         first_order_constant: 0,
                         second_order_squared_logarithm: literature_second_order_squared_logarithm,
                         second_order_logarithm: literature_second_order_logarithm}
    literature_ratio_general_constant = general_ratio.subs(literature_values)
    literature_ratio = literature_ratio_general_constant.subs(
        second_order_constant, literature_second_order_constant)
    literature_residual = consistency_residual(sp.expand(literature_ratio * coupling))
    literature_residual_low_orders = tuple(
        sp.expand(literature_residual.coeff(coupling, power)) for power in range(4))

    # Cross-check with the form of Chetyrkin:1997un written with m_h(mu):
    # ln(mu^2/m_h(mu)^2) = ln(mu^2/mu_h^2) (1 + 2 x alpha_s/pi) + higher orders
    # if mu^2 d m_h/dmu^2 = -x (alpha_s/pi) m_h. The value of x that
    # reconciles -11/24 with -19/24 is solved for; it is not an input.
    mass_anomalous_dimension = sp.Symbol("mass_anomalous_dimension_coefficient")
    running_mass_logarithm = threshold_logarithm * (1 + 2 * mass_anomalous_dimension * coupling)
    running_mass_form = (1 - coupling * running_mass_logarithm / 6
                         + coupling**2 * (sp.Rational(11, 72) - sp.Rational(11, 24) * running_mass_logarithm
                                          + running_mass_logarithm**2 / 36))
    running_mass_difference = sp.expand(running_mass_form - literature_ratio)
    running_mass_difference_second_order = sp.factor(running_mass_difference.coeff(coupling, 2))
    reconciling_mass_anomalous_dimension = sp.solve(
        running_mass_difference_second_order, mass_anomalous_dimension)

    # Jump at mu = mu_h.
    ratio_at_threshold = literature_ratio_general_constant.subs(threshold_logarithm, 0)
    first_order_term_at_threshold = ratio_at_threshold.coeff(coupling, 1)
    coupling_jump = sp.expand((ratio_at_threshold - 1) * coupling)
    # Coefficient 1 - alpha_s/pi of Ludtke:2020moa evaluated just below
    # (n_f - 1 flavors) and just above (n_f flavors) the threshold.
    coefficient_below = 1 - ratio_at_threshold * coupling
    coefficient_above = 1 - coupling
    coefficient_jump = sp.expand(coefficient_below - coefficient_above)
    coefficient_jump_numerical_constant = coefficient_jump.subs(
        second_order_constant, literature_second_order_constant)
    coefficient_jump_lowest_power = min(
        power for (power,), _ in sp.Poly(coefficient_jump_numerical_constant, coupling).terms())
    return locals()


@lru_cache(maxsize=1)
def derive_singlet_evolution():
    """Part C: leading-order singlet matrix-element ratio with 3, 4, 5 active flavors."""
    scale = sp.Symbol("mu", positive=True)
    reference_scale = sp.Symbol("mu_0", positive=True)
    charm_threshold = sp.Symbol("mu_c", positive=True)
    bottom_threshold = sp.Symbol("mu_b", positive=True)
    top_threshold = sp.Symbol("mu_t", positive=True)
    casimir_fundamental = sp.Symbol("C_F", positive=True)
    quantum_chromodynamics = {casimir_fundamental: sp.Rational(4, 3)}
    # alpha_s/(4 pi) as in Larin:1993tq, at mu_0 in the three-flavor theory.
    coupling_reference = sp.Symbol("coupling_at_reference_scale", positive=True)
    # Independent symbols for alpha_s/(4 pi) at mu_c, mu_b and mu.
    coupling_charm = sp.Symbol("coupling_at_charm_threshold", positive=True)
    coupling_bottom = sp.Symbol("coupling_at_bottom_threshold", positive=True)
    coupling = sp.Symbol("coupling_at_scale", positive=True)
    strong_coupling_names = {
        coupling_reference: r"\alpha_s(\mu_0)", coupling_charm: r"\alpha_s(\mu_c)",
        coupling_bottom: r"\alpha_s(\mu_b)", coupling: r"\alpha_s(\mu)"}
    symbol_names = {symbol: name.replace(r"\alpha_s", r"\frac{\alpha_s}{4\pi}")
                    for symbol, name in strong_coupling_names.items()}

    def beta_zero(number):
        """Coefficient in mu^2 d(alpha_s/(4 pi))/dmu^2 = -beta_0 (alpha_s/(4 pi))^2."""
        return 11 - sp.Rational(2, 3) * number

    def one_loop_coupling(initial_coupling, initial_scale, number):
        return initial_coupling / (1 + beta_zero(number) * initial_coupling
                                   * sp.log(scale**2 / initial_scale**2))

    def interval_logarithm(number, upper_coupling, lower_coupling):
        """Logarithm of the ratio accumulated in one interval with fixed flavor number."""
        return 6 * casimir_fundamental * number / beta_zero(number) * (upper_coupling - lower_coupling)

    # Couplings in the three intervals, continuous at each threshold.
    coupling_three = one_loop_coupling(coupling_reference, reference_scale, 3)
    coupling_three_at_charm = coupling_three.subs(scale, charm_threshold)
    coupling_four = one_loop_coupling(coupling_three_at_charm, charm_threshold, 4)
    coupling_four_at_bottom = coupling_four.subs(scale, bottom_threshold)
    coupling_five = one_loop_coupling(coupling_four_at_bottom, bottom_threshold, 5)
    couplings = {3: coupling_three, 4: coupling_four, 5: coupling_five}
    thresholds = {4: charm_threshold, 5: bottom_threshold}

    # Logarithm of the ratio of the matrix element at mu to that at mu_0, in
    # terms of the couplings at the interval ends (product of exponentials).
    log_ratio_by_couplings = {
        3: interval_logarithm(3, coupling, coupling_reference),
        4: (interval_logarithm(3, coupling_charm, coupling_reference)
            + interval_logarithm(4, coupling, coupling_charm)),
        5: (interval_logarithm(3, coupling_charm, coupling_reference)
            + interval_logarithm(4, coupling_bottom, coupling_charm)
            + interval_logarithm(5, coupling, coupling_bottom))}
    explicit_couplings = {coupling_charm: coupling_three_at_charm,
                          coupling_bottom: coupling_four_at_bottom}
    log_ratio = {number: expression.subs(explicit_couplings).subs(coupling, couplings[number])
                 for number, expression in log_ratio_by_couplings.items()}
    ratio = {number: sp.exp(expression) for number, expression in log_ratio.items()}
    to_strong_coupling = {symbol: sp.Symbol(symbol.name.replace("coupling", "strong_coupling"), positive=True)
                          / (4 * sp.pi) for symbol in strong_coupling_names}
    strong_coupling_symbol_names = {
        sp.Symbol(symbol.name.replace("coupling", "strong_coupling"), positive=True): name
        for symbol, name in strong_coupling_names.items()}
    log_ratio_strong_coupling = {
        number: sp.expand(expression.subs(to_strong_coupling).subs(quantum_chromodynamics))
        for number, expression in log_ratio_by_couplings.items()}

    # Checks: one-loop equation and evolution equation in each interval.
    coupling_residuals = {number: sp.simplify(
        scale_derivative(couplings[number], scale) + beta_zero(number) * couplings[number]**2)
        for number in couplings}
    evolution_residuals = {number: sp.simplify(
        scale_derivative(log_ratio[number], scale)
        + 6 * casimir_fundamental * number * couplings[number]**2) for number in couplings}
    ratio_at_reference_scale = sp.simplify(ratio[3].subs(scale, reference_scale))

    # Continuity and jumps at each threshold (above minus below).
    coupling_discontinuities = {}
    ratio_discontinuities = {}
    log_derivative_jumps = {}
    expected_log_derivative_jumps = {}
    coupling_derivative_jumps = {}
    for number, threshold in thresholds.items():
        coupling_at_threshold = couplings[number - 1].subs(scale, threshold)
        coupling_discontinuities[number] = sp.simplify(
            couplings[number].subs(scale, threshold) - coupling_at_threshold)
        ratio_discontinuities[number] = sp.simplify(
            ratio[number].subs(scale, threshold) / ratio[number - 1].subs(scale, threshold) - 1)
        log_derivative_jumps[number] = sp.simplify(
            (scale_derivative(log_ratio[number], scale)
             - scale_derivative(log_ratio[number - 1], scale)).subs(scale, threshold))
        expected_log_derivative_jumps[number] = -6 * casimir_fundamental * coupling_at_threshold**2
        coupling_derivative_jumps[number] = sp.simplify(
            (scale_derivative(couplings[number], scale)
             - scale_derivative(couplings[number - 1], scale)).subs(scale, threshold)
            / coupling_at_threshold**2)
    log_derivative_jump_residuals = {number: sp.simplify(
        log_derivative_jumps[number] - expected_log_derivative_jumps[number]) for number in thresholds}
    # Size of the jump of mu^2 d ln(ratio)/dmu^2 in terms of alpha_s at the threshold.
    threshold_strong_coupling = sp.Symbol("strong_coupling_at_threshold", positive=True)
    log_derivative_jump_strong_coupling = (-6 * casimir_fundamental
                                           * (threshold_strong_coupling / (4 * sp.pi))**2
                                           ).subs(quantum_chromodynamics)

    # Active flavors as a function of mu = sqrt(|q^2|); thresholds symbolic.
    active_flavor_table = (
        {"interval": scale < charm_threshold, "active_flavors": 3, "supported": True,
         "comment": "up, down, strange"},
        {"interval": sp.And(charm_threshold < scale, scale < bottom_threshold),
         "active_flavors": 4, "supported": True, "comment": "charm added at mu_c = m_c(m_c)"},
        {"interval": sp.And(bottom_threshold < scale, scale < top_threshold),
         "active_flavors": 5, "supported": True, "comment": "bottom added at mu_b = m_b(m_b)"},
        {"interval": scale > top_threshold, "active_flavors": None, "supported": False,
         "comment": "top-quark region: outside the supported domain, electroweak effects not treated"})
    return locals()


class MassiveCoefficientTests(unittest.TestCase):
    """Part A: Clifford decomposition of the massive free-quark contractions."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_massive_coefficient()

    def test_clifford_basis_is_orthogonal_and_oriented(self):
        result = self.result
        overlaps = result["basis_overlaps"]
        self.assertEqual(overlaps, sp.diag(*result["basis_norms"]))
        self.assertTrue(all(norm != 0 for norm in result["basis_norms"]))
        gamma = result["gamma"]
        self.assertEqual(sp.trace(gamma[0] * gamma[1] * gamma[2] * gamma[3] * result["gamma5"]), -4 * sp.I)

    def test_only_axial_and_scalar_structures_appear(self):
        result = self.result
        for first in range(4):
            for second in range(4):
                difference = (-result["contraction_sum"][first][second] / result["denominator"]
                              - result["reconstructed"][first][second]).applyfunc(sp.simplify)
                self.assertEqual(difference, sp.zeros(4))
                for position, label in enumerate(result["basis_labels"]):
                    coefficient = result["clifford_coefficients"][first][second][position]
                    self.assertEqual(sp.simplify(
                        coefficient - result["expected_coefficients"][first][second][position]), 0)
                    if not (label == "scalar" or label.startswith("axial")):
                        self.assertEqual(coefficient, 0)
                    if label.startswith("axial"):
                        # The axial coefficient contains the mass only through q^2 - m^2.
                        self.assertEqual(sp.simplify(coefficient * result["denominator"]).has(
                            result["quark_mass"]), False)

    def test_mass_linear_terms_are_symmetric_and_scalar(self):
        result = self.result
        for first in range(4):
            for second in range(4):
                self.assertEqual(result["mass_linear_antisymmetric"][first][second], sp.zeros(4))
                self.assertEqual(result["mass_linear"][first][second],
                                 2 * result["metric"][first, second] * sp.eye(4))

    def test_massless_limit_reproduces_massless_script(self):
        from derive_alp_hard_current_ope import derive as derive_massless
        massless = derive_massless()
        result = self.result
        replacement = dict(zip(result["hard_momentum"], massless["hard_momentum"]))
        axial_positions = [position for position, label in enumerate(result["basis_labels"])
                           if label.startswith("axial")]
        metric = result["metric"]
        for first in range(4):
            for second in range(4):
                for axial, position in enumerate(axial_positions):
                    # The massless script stores the coefficient with an upper axial index.
                    coefficient = (result["clifford_coefficients"][first][second][position]
                                   .subs(result["quark_mass"], 0).subs(replacement) * metric[axial, axial])
                    self.assertEqual(sp.simplify(
                        coefficient - massless["fourier_coefficients"][first][second][axial]), 0)

    def test_heavy_flavor_factor_and_limits(self):
        result = self.result
        mass, magnitude = result["quark_mass"], result["spacelike_magnitude"]
        self.assertEqual(sp.simplify(result["heavy_flavor_factor_spacelike"]
                                     - magnitude / (magnitude + mass**2)), 0)
        self.assertEqual(result["axial_coefficient"].subs(mass, 0), result["massless_axial_coefficient"])
        self.assertEqual(result["large_momentum_limit"], 1)
        self.assertEqual(result["small_momentum_limit"], 0)
        self.assertEqual(sp.simplify(result["large_momentum_series"]
                                     - (1 - mass**2 / magnitude + mass**4 / magnitude**2)), 0)
        self.assertEqual(sp.simplify(result["small_momentum_series"]
                                     - (magnitude / mass**2 - magnitude**2 / mass**4)), 0)
        self.assertEqual(result["heavy_flavor_factor_spacelike"].subs({magnitude: 3, mass: 1}),
                         sp.Rational(3, 4))


class CouplingMatchingTests(unittest.TestCase):
    """Part B: matching of alpha_s at a heavy-quark threshold."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_coupling_matching()

    def test_logarithmic_coefficients_follow_from_the_beta_functions(self):
        result = self.result
        self.assertEqual(result["derived_first_order_logarithm"], -sp.Rational(1, 6))
        self.assertEqual(result["derived_second_order_squared_logarithm"], sp.Rational(1, 36))
        self.assertEqual(result["derived_second_order_logarithm"], -sp.Rational(19, 24))
        self.assertEqual(result["coefficient_differences"], (0, 0, 0))
        # In general the single-logarithm coefficient is beta_1 difference plus
        # a term proportional to the first-order constant.
        self.assertEqual(sp.simplify(result["derived_second_order_logarithm_general"]
                                     - (-sp.Rational(19, 24) - result["first_order_constant"] / 3)), 0)

    def test_constant_is_not_fixed_by_consistency(self):
        result = self.result
        self.assertNotIn(result["second_order_constant"], result["consistency_solution"])
        for equation in result["consistency_equations"]:
            self.assertFalse(sp.sympify(equation).has(result["second_order_constant"]))
        self.assertEqual(result["literature_residual_low_orders"], (0, 0, 0, 0))
        # The first residual of the truncated series is of fourth order.
        self.assertNotEqual(sp.expand(result["literature_residual"].coeff(result["coupling"], 4)), 0)

    def test_running_mass_form_is_equivalent(self):
        self.assertEqual(self.result["reconciling_mass_anomalous_dimension"], [1])

    def test_jump_at_threshold(self):
        result = self.result
        coupling, constant = result["coupling"], result["second_order_constant"]
        self.assertEqual(result["first_order_term_at_threshold"], 0)
        self.assertEqual(result["coupling_jump"], constant * coupling**3)
        self.assertEqual(result["coefficient_jump"], -constant * coupling**3)
        self.assertEqual(result["coefficient_jump_numerical_constant"],
                         -sp.Rational(11, 72) * coupling**3)
        self.assertEqual(result["coefficient_jump_lowest_power"], 3)


class SingletEvolutionTests(unittest.TestCase):
    """Part C: singlet matrix-element ratio across the charm and bottom thresholds."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_singlet_evolution()

    def test_evolution_equations_in_each_interval(self):
        result = self.result
        for number in (3, 4, 5):
            self.assertEqual(result["coupling_residuals"][number], 0)
            self.assertEqual(result["evolution_residuals"][number], 0)
        self.assertEqual(result["ratio_at_reference_scale"], 1)

    def test_continuity_at_thresholds(self):
        result = self.result
        for number in (4, 5):
            self.assertEqual(result["coupling_discontinuities"][number], 0)
            self.assertEqual(result["ratio_discontinuities"][number], 0)

    def test_derivative_jump_at_thresholds(self):
        result = self.result
        for number in (4, 5):
            self.assertEqual(result["log_derivative_jump_residuals"][number], 0)
            self.assertNotEqual(result["log_derivative_jumps"][number], 0)
            # beta_0 decreases by 2/3 per added flavor.
            self.assertEqual(result["coupling_derivative_jumps"][number], sp.Rational(2, 3))
        strong_coupling = result["threshold_strong_coupling"]
        self.assertEqual(sp.simplify(result["log_derivative_jump_strong_coupling"]
                                     + strong_coupling**2 / (2 * sp.pi**2)), 0)

    def test_fixed_flavor_form_and_strong_coupling_form(self):
        result = self.result
        # Below the charm threshold the result is the fixed n_f = 3 ratio of
        # derive_alp_singlet_axial_rg.py.
        from derive_alp_singlet_axial_rg import derive as derive_fixed_flavor
        fixed = derive_fixed_flavor()
        fixed_log_ratio = fixed["leading_log_ratio"].subs(fixed["beta_values"]).subs(
            {fixed["flavor_number"]: 3, fixed["coupling"]: result["coupling"],
             fixed["coupling_reference"]: result["coupling_reference"],
             fixed["casimir_fundamental"]: result["casimir_fundamental"]})
        self.assertEqual(sp.simplify(result["log_ratio_by_couplings"][3] - fixed_log_ratio), 0)
        names = {symbol.name: symbol for symbol in result["strong_coupling_symbol_names"]}
        expected = 2 / sp.pi * (
            sp.Rational(1, 3) * (names["strong_coupling_at_charm_threshold"]
                                 - names["strong_coupling_at_reference_scale"])
            + sp.Rational(12, 25) * (names["strong_coupling_at_bottom_threshold"]
                                     - names["strong_coupling_at_charm_threshold"])
            + sp.Rational(15, 23) * (names["strong_coupling_at_scale"]
                                     - names["strong_coupling_at_bottom_threshold"]))
        self.assertEqual(sp.simplify(result["log_ratio_strong_coupling"][5] - expected), 0)

    def test_active_flavor_table(self):
        table = self.result["active_flavor_table"]
        self.assertEqual([row["active_flavors"] for row in table], [3, 4, 5, None])
        self.assertEqual([row["supported"] for row in table], [True, True, True, False])
        for row in table:
            for symbol in row["interval"].free_symbols:
                self.assertTrue(symbol.is_Symbol)
            self.assertFalse(any(atom.is_Float for atom in row["interval"].atoms(sp.Number)))


def main():
    massive = derive_massive_coefficient()
    names = massive["symbol_names"]
    print("PART A. One free quark of mass m, S_F(k)=i(kslash+m)/(k^2-m^2+i0), leading order in alpha_s,")
    print("local term without derivatives; metric (+---), epsilon0123=+1; +i0 implied in denominators.")
    show("sum of the two contractions before the factor -1/(q^2-m^2), mu=1 nu=2", massive["contraction_sum"][1][2])
    show("the same, mu=1 nu=1", massive["contraction_sum"][1][1])
    show("Clifford basis labels", sp.Matrix([sp.Symbol(label.replace(" ", "_")) for label in massive["basis_labels"]]).T)
    show("Clifford coefficients of the full contraction, mu=1 nu=2 (order of the basis above)",
         sp.Matrix(massive["clifford_coefficients"][1][2]).T)
    show("Clifford coefficients of the full contraction, mu=1 nu=1", sp.Matrix(massive["clifford_coefficients"][1][1]).T)
    show("antisymmetric part in mu, nu of the terms linear in m, mu=1 nu=2", massive["mass_linear_antisymmetric"][1][2])
    show("coefficient of epsilon^{mu nu rho sigma} q_rho (qbar gamma_sigma gamma5 q)", massive["axial_coefficient"])
    show("coefficient of g^{mu nu} (qbar q)", massive["scalar_density_coefficient"])
    show("massive axial coefficient divided by the massless one", massive["heavy_flavor_factor"])
    show("the same for spacelike q^2", massive["heavy_flavor_factor_spacelike"], names)
    show("expansion for |q^2| >> m^2", massive["large_momentum_series"], names)
    show("expansion for |q^2| << m^2", massive["small_momentum_series"], names)

    matching = derive_coupling_matching()
    names = matching["symbol_names"]
    print("PART B. alpha_s matching at a heavy-quark threshold, mu_h = m_h(m_h), series in alpha_s/pi.")
    show("coefficient of (alpha_s/pi) ln(mu^2/mu_h^2) required by the two beta functions",
         matching["derived_first_order_logarithm"])
    show("coefficient of (alpha_s/pi)^2 ln^2(mu^2/mu_h^2) required by the two beta functions",
         matching["derived_second_order_squared_logarithm"])
    show("coefficient of (alpha_s/pi)^2 ln(mu^2/mu_h^2) required by the two beta functions, for a general first-order constant",
         matching["derived_second_order_logarithm_general"],
         {matching["first_order_constant"]: r"\text{(first-order constant)}"})
    show("the same with vanishing first-order constant", matching["derived_second_order_logarithm"])
    show("derived minus Chetyrkin:1997sg for the three logarithmic coefficients",
         sp.Matrix(matching["coefficient_differences"]).T)
    show("alpha_s^{(n_f-1)}/alpha_s^{(n_f)} of Chetyrkin:1997sg", matching["literature_ratio"], names)
    show("one-loop quark-mass anomalous-dimension coefficient that makes the m_h(mu) form of Chetyrkin:1997un equal to it",
         sp.Matrix(matching["reconciling_mass_anomalous_dimension"]))
    show("alpha_s^{(n_f-1)}/pi - alpha_s^{(n_f)}/pi at mu = mu_h",
         matching["coupling_jump"].subs(matching["second_order_constant"], sp.Rational(11, 72)), names)
    show("jump of 1 - alpha_s/pi at mu = mu_h (below minus above)",
         matching["coefficient_jump_numerical_constant"], names)

    singlet = derive_singlet_evolution()
    names = singlet["symbol_names"]
    print("PART C. Leading-order singlet axial-current matrix element with 3, 4, 5 active flavors.")
    print("One-loop alpha_s continuous at mu_c and mu_b (the constant 11/72 is neglected).")
    print("ASSUMPTION: matrix element continuous across each threshold (no matching coefficient available).")
    for number, region in ((3, "mu < mu_c"), (4, "mu_c < mu < mu_b"), (5, "mu > mu_b")):
        show(f"logarithm of the ratio of the matrix element at mu to that at mu_0, {region}",
             singlet["log_ratio_by_couplings"][number], names)
    show("the same for mu > mu_b with C_F = 4/3 in terms of alpha_s",
         singlet["log_ratio_strong_coupling"][5], singlet["strong_coupling_symbol_names"])
    show("mu^2 d ln(ratio)/dmu^2 + 6 C_F n (alpha_s/(4 pi))^2 for n = 3, 4, 5",
         sp.Matrix([singlet["evolution_residuals"][number] for number in (3, 4, 5)]).T)
    show("discontinuity of alpha_s/(4 pi) at mu_c, mu_b",
         sp.Matrix([singlet["coupling_discontinuities"][number] for number in (4, 5)]).T)
    show("relative discontinuity of the ratio at mu_c, mu_b",
         sp.Matrix([singlet["ratio_discontinuities"][number] for number in (4, 5)]).T)
    show("jump of mu^2 d ln(ratio)/dmu^2 at a threshold (above minus below), C_F = 4/3",
         singlet["log_derivative_jump_strong_coupling"],
         {singlet["threshold_strong_coupling"]: r"\alpha_s(\mu_{\text{threshold}})"})
    show("jump of mu^2 d(alpha_s/(4 pi))/dmu^2 at a threshold divided by (alpha_s/(4 pi))^2",
         sp.Matrix([singlet["coupling_derivative_jumps"][number] for number in (4, 5)]).T)
    print("Active flavors versus mu = sqrt(|q^2|); mu_c = m_c(m_c), mu_b = m_b(m_b), mu_t = m_t(m_t):")
    for row in singlet["active_flavor_table"]:
        flavors = row["active_flavors"] if row["supported"] else "not supported"
        print(f"  {sp.latex(row['interval'])} | active flavors: {flavors} | {row['comment']}")


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
