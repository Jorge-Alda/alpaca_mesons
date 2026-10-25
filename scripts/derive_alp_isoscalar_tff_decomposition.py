"""Derive current-labelled no-mixing ALP tree TFF bookkeeping.

Run .venv/bin/python scripts/derive_alp_isoscalar_tff_decomposition.py
or append --test for its focused unit tests. Natural units, GeV,
metric (+---), epsilon^{0123}=+1, i M=i e^2 R F. Stable neutral
vectors; isospin-diagonal even sector; arbitrary diagonal U_a,Y_a.
The supplied total electromagnetic contact is kept outside this split.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp

import derive_alp_doubly_virtual_tff as tree


def show(label, expression):
    """Print derived expressions through SymPy's standard LaTeX printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


@lru_cache(maxsize=1)
def derive():
    """Expand the existing ordered vertices in two independent currents.

    The explicit diagonal matrices are supplied in the indexed tree/rho
    records; no general matrix is assumed commuting. Formal lambda variables
    only select current labels and are not new physical couplings.
    """
    charge = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    interchange = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    interchanged_charge = interchange * charge * interchange.T
    currents = ((charge + interchanged_charge) / 2,
                (charge - interchanged_charge) / 2)
    theta = sp.Symbol("theta_V", real=True)
    rho = sp.diag(1, -1, 0) / sp.sqrt(2)
    ideal_omega = sp.diag(1, 1, 0) / sp.sqrt(2)
    ideal_phi = sp.diag(0, 0, 1)
    generators = (rho, sp.cos(theta) * ideal_omega - sp.sin(theta) * ideal_phi,
                  sp.sin(theta) * ideal_omega + sp.cos(theta) * ideal_phi)
    decay_constant = sp.Symbol("f_a", positive=True)
    axial_coefficients = sp.symbols("c_u^A c_d^A c_s^A", real=True)
    pseudoscalar_coefficients = sp.symbols("c_u^P c_d^P c_s^P", real=True)
    chi_values = sp.symbols("chi_u chi_d chi_s", real=True)
    axial_source = 2 * sp.diag(*axial_coefficients) / decay_constant
    pseudoscalar_source = sp.diag(*pseudoscalar_coefficients)
    chi = sp.diag(*chi_values)
    mass_source = -2 * (pseudoscalar_source * chi + chi * pseudoscalar_source) / decay_constant
    source_constant, single_trace, double_trace, light_spurion, strange_spurion = sp.symbols(
        "F_V z_V y_V x_l x_s", real=True)
    vacuum_spurion = sp.diag(light_spurion, light_spurion, strange_spurion)

    def residues(current):
        photon_source = (source_constant * current
                         + single_trace * (vacuum_spurion * current + current * vacuum_spurion)
                         + double_trace * sp.trace(vacuum_spurion * current) * sp.eye(3))
        return tuple(sp.simplify(sp.trace(generator * photon_source)) for generator in generators)

    current_residues = tuple(residues(current) for current in currents)
    full_residues = residues(charge)
    alp_mass = sp.Symbol("m_a", positive=True)

    def vjp(generator, current):
        axial_trace = sp.trace(generator * current * axial_source)
        mass_trace = sp.trace(generator * current * mass_source)
        result, _ = tree.ordered_numerators(axial_trace, mass_trace)
        return sp.expand(result.subs(tree.tensor.m, alp_mass))

    def vvp(first_generator, second_generator):
        _, result = tree.ordered_numerators(
            sp.trace(first_generator * second_generator * axial_source),
            sp.trace(first_generator * second_generator * mass_source))
        return sp.expand(result.subs(tree.tensor.m, alp_mass))

    vjp_numerators = tuple(tuple(vjp(generator, current) for generator in generators)
                           for current in currents)
    full_vjp = tuple(vjp(generator, charge) for generator in generators)
    pair_numerators = tuple(tuple(vvp(first, second) for second in generators)
                            for first in generators)
    first_virtuality, second_virtuality = tree.tensor.s1, tree.tensor.s2
    swap = {first_virtuality: second_virtuality, second_virtuality: first_virtuality}
    masses = sp.symbols("M_rho M_omega M_phi", positive=True)
    names = ("rho", "omega", "phi")
    # Abstract W functions keep the assembly readable. Their stable-pole
    # definitions below use the residues derived from the stated even source.
    response_functions = tuple(tuple(sp.Function(f"W_{name}^{{({label})}}")
                                     for name in names) for label in range(2))
    first_responses = tuple(tuple(function(first_virtuality) for function in row)
                            for row in response_functions)
    second_responses = tuple(tuple(function(second_virtuality) for function in row)
                             for row in response_functions)
    formal_labels = sp.symbols("lambda_10 lambda_11 lambda_20 lambda_21")
    first_labels, second_labels = formal_labels[:2], formal_labels[2:]
    assembled = 2 * sp.sqrt(2) * sum(
        sum(first_labels[label] * first_responses[label][index] for label in range(2))
        * sum(second_labels[label] * vjp_numerators[label][index] for label in range(2))
        + sum(second_labels[label] * second_responses[label][index] for label in range(2))
        * sum(first_labels[label] * vjp_numerators[label][index].xreplace(swap) for label in range(2))
        for index in range(3))
    assembled -= 2 * sum(
        sum(first_labels[label] * first_responses[label][first] for label in range(2))
        * sum(second_labels[label] * second_responses[label][second] for label in range(2))
        * pair_numerators[first][second]
        for first in range(3) for second in range(3))
    current_terms = tuple(tuple(sp.diff(assembled, first_labels[first], second_labels[second]).subs(
        dict.fromkeys(formal_labels, 0)) for second in range(2)) for first in range(2))
    full_abstract_responses = tuple(sp.Function(f"W_{name}") for name in names)
    recombination = {response_functions[1][index](virtuality):
                     full_abstract_responses[index](virtuality) - response_functions[0][index](virtuality)
                     for index in range(3) for virtuality in (first_virtuality, second_virtuality)}
    recombined = sp.expand(sum(sum(row) for row in current_terms).xreplace(recombination))
    full_first = tuple(function(first_virtuality) for function in full_abstract_responses)
    full_second = tuple(function(second_virtuality) for function in full_abstract_responses)
    full_vjp_scalar = 2 * sp.sqrt(2) * sum(
        full_first[index] * full_vjp[index]
        + full_second[index] * full_vjp[index].xreplace(swap) for index in range(3))
    pair_terms = tuple(tuple(-2 * full_first[first] * full_second[second] * pair_numerators[first][second]
                            for second in range(3)) for first in range(3))
    full_pair_scalar = sum(sum(row) for row in pair_terms)
    full_scalar = full_vjp_scalar + full_pair_scalar
    mixed_pair_scalar = sum(pair_terms[0][index] + pair_terms[index][0] for index in (1, 2))
    rho_union = (2 * sp.sqrt(2) * (full_first[0] * full_vjp[0]
                 + full_second[0] * full_vjp[0].xreplace(swap))
                 + pair_terms[0][0] + mixed_pair_scalar)
    isoscalar_union = (full_vjp_scalar - 2 * sp.sqrt(2) * (
        full_first[0] * full_vjp[0] + full_second[0] * full_vjp[0].xreplace(swap))
        + sum(pair_terms[first][second] for first in range(3) for second in range(3)
              if first in (1, 2) or second in (1, 2)))
    stable_substitutions = {full_abstract_responses[index](virtuality):
                            full_residues[index] / (masses[index]**2 - virtuality)
                            for index in range(3) for virtuality in (first_virtuality, second_virtuality)}
    stored_vjp, stored_pair = tree.assemble_scalar(masses, full_residues, full_vjp, pair_numerators)
    stable_difference = sp.expand(full_scalar.xreplace(stable_substitutions) - stored_vjp - stored_pair)
    first_momentum = sp.Matrix(sp.symbols("q_10:14"))
    second_momentum = sp.Matrix(sp.symbols("q_20:24"))
    # epsilon_lower=-epsilon_upper for (+---). This is the actual coefficient
    # of epsilon1^mu epsilon2^nu; no division by a tensor is used.
    photon_tensor = sp.Matrix(4, 4, lambda first, second: -sum(
        sp.LeviCivita(first, second, third, fourth)
        * first_momentum[third] * second_momentum[fourth]
        for third in range(4) for fourth in range(4)))
    ward_first = (first_momentum.T * photon_tensor).applyfunc(sp.expand)
    ward_second = (photon_tensor * second_momentum).applyfunc(sp.expand)
    return locals()


class CurrentDecompositionTests(unittest.TestCase):
    """Check independent labels, source dependence and stored normalization."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_charge_and_photon_source_recombination(self):
        result = self.result
        self.assertEqual(sum(result["currents"], sp.zeros(3)), result["charge"])
        for index in range(3):
            self.assertEqual(sp.simplify(sum(row[index] for row in result["current_residues"])
                                        - result["full_residues"][index]), 0)
        self.assertEqual(result["current_residues"][0][0], 0)
        self.assertEqual(result["current_residues"][1][1:], (0, 0))

    def test_current_recombination_and_bose_exchange(self):
        result = self.result
        self.assertEqual(sp.expand(result["recombined"] - result["full_scalar"]), 0)
        for first in range(2):
            for second in range(2):
                difference = (result["current_terms"][first][second]
                              - result["current_terms"][second][first].xreplace(result["swap"]))
                self.assertEqual(sp.expand(difference), 0)

    def test_source_linearity_and_mixed_overlap(self):
        result = self.result
        coefficients = result["axial_coefficients"] + result["pseudoscalar_coefficients"]
        polynomial = sp.Poly(sp.expand(result["full_scalar"]), *coefficients)
        self.assertEqual(polynomial.total_degree(), 1)
        self.assertEqual(polynomial.TC(), 0)
        self.assertEqual(sp.expand(result["rho_union"] + result["isoscalar_union"]
                                  - result["mixed_pair_scalar"] - result["full_scalar"]), 0)
        self.assertNotEqual(sp.expand(result["mixed_pair_scalar"]), 0)
        exact_isospin = {
            result["axial_coefficients"][1]: result["axial_coefficients"][0],
            result["pseudoscalar_coefficients"][1]: result["pseudoscalar_coefficients"][0],
            result["chi_values"][1]: result["chi_values"][0],
        }
        self.assertEqual(sp.expand(result["mixed_pair_scalar"].subs(exact_isospin)), 0)

    def test_stable_pole_normalization_and_ward_identities(self):
        result = self.result
        self.assertEqual(result["stable_difference"], 0)
        self.assertEqual(result["ward_first"], sp.zeros(1, 4))
        self.assertEqual(result["ward_second"], sp.zeros(4, 1))


def main():
    result = derive()
    print("Stable tree poles; independent diagonal no-mixing sources; GeV; (+---).")
    print("Current labels refer to three-flavor isospin projections, not cut labels.")
    for label, current in enumerate(result["currents"]):
        show(f"Q^({label}) from up/down interchange", current)
        for index, name in enumerate(result["names"]):
            show(f"R_{name}^({label}) from the linear photon source", result["current_residues"][label][index])
    for index, name in enumerate(result["names"]):
        show(f"V_{name}, current (0)", sp.factor(result["vjp_numerators"][0][index]))
        show(f"V_{name}, current (1)", sp.factor(result["vjp_numerators"][1][index]))
    print("Each current coefficient is obtained by differentiating the labelled ordered assembly.")
    for first in range(2):
        for second in range(2):
            show(f"Resonance scalar for current labels ({first},{second})",
                 result["current_terms"][first][second])
    show("Current recombination minus full retained resonance scalar",
         sp.expand(result["recombined"] - result["full_scalar"]))
    show("Mixed rho-isoscalar VVP terms present in both vector unions", result["mixed_pair_scalar"])
    show("rho union + isoscalar union - mixed overlap - full resonance scalar",
         sp.expand(result["rho_union"] + result["isoscalar_union"]
                   - result["mixed_pair_scalar"] - result["full_scalar"]))
    show("Stable-pole scalar minus stored assembly", result["stable_difference"])
    show("Photon-1 Ward contraction", result["ward_first"])
    show("Photon-2 Ward contraction", result["ward_second"])
    print("Keep +widehat(c_gamma)/(4 pi^2 f_a) outside the current split.")
    print("No widths, hadronic cuts, numerical source choice or high-energy matching.")


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
