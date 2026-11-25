"""Calculate an implicit ALP-to-three-meson loop crossing from JSON inputs.

Input keys B_0, f_pi, m_0_squared, m_a_squared, s12, s23,
loop_scale_squared and R_max are numbers in GeV units. S and P each contain
seven 3x3 matrices, theta contains seven real coefficients, all ordered by
derivative order n=0..6 without factorials. The daughter_indices are three
zero-based indices in the sorted unmixed meson spectrum. bracket gives two
positive values of f_pi/f_a spanning a root on a continuous pole branch.
Use JSON strings such as "1/3" for exact rational input.

The same positive reference scale is used for non-tadpole scalar moments.
The actual internal mass is used for standalone tadpoles.
"""

import argparse
import json
import math
import unittest

import sympy as sp
from scipy.optimize import brentq

from alp_loop_crossing_poles import invariant_quadratic_blocks, exact_tree_amplitude
from derive_alp_finite_mixing_three_meson_tree import GENERATORS
from derive_alp_three_meson_onshell_driver import (
    ScalarMoment, crossing_condition, dressed_tree_vertex,
    formal_scalar_moment, on_shell_cubic_estimate,
)


def pole_expansions(meson_mass_matrix, mixing_source, alp_mass_squared, epsilon):
    """Return the canonical eigenbasis and pole expansions at fixed ALP mass.

    The meson mass matrix is real symmetric, the mixing source is the column
    appearing with positive sign in D_{na}/epsilon, and masses are in GeV².
    The calculation assumes distinct coupled meson eigenvalues and an ALP mass
    separated from them. All parameters may be exact SymPy expressions.
    """
    mass_matrix = sp.Matrix(meson_mass_matrix)
    source = sp.Matrix(mixing_source)
    meson_vectors = []
    meson_masses = []
    for mass, multiplicity, vectors in mass_matrix.eigenvects():
        if len(vectors) != multiplicity:
            raise ValueError("The meson mass matrix must be diagonalizable.")
        for vector in sp.GramSchmidt(vectors, orthonormal=True):
            meson_vectors.append(vector)
            meson_masses.append(mass)
    if all(not mass.free_symbols for mass in meson_masses):
        ordering = sorted(range(len(meson_masses)), key=lambda index: float(sp.N(meson_masses[index])))
        meson_vectors = [meson_vectors[index] for index in ordering]
        meson_masses = [meson_masses[index] for index in ordering]
    count = len(meson_masses)
    couplings = [(vector.T * source)[0] for vector in meson_vectors]
    for index, mass in enumerate(meson_masses):
        if sp.simplify(mass - alp_mass_squared) == 0 and couplings[index] != 0:
            raise ValueError("The ALP and a coupled meson are degenerate.")
        for other in range(index):
            if sp.simplify(mass - meson_masses[other]) == 0 and (
                    couplings[index] != 0 or couplings[other] != 0):
                raise ValueError("Diagonalize mixing within the degenerate meson subspace.")
    fixed_mass = sum(couplings[j]**2 / (alp_mass_squared - meson_masses[j])
                     for j in range(count) if couplings[j] != 0)
    overlap = sp.Matrix([couplings[j] / (meson_masses[j] - alp_mass_squared)
                         if couplings[j] != 0 else 0 for j in range(count)])
    alp_norm = sum(value**2 for value in overlap)
    alp_vector = sp.Matrix.vstack(
        epsilon * overlap - epsilon**3 * alp_norm * overlap / 2,
        sp.Matrix([1 - epsilon**2 * alp_norm / 2]),
    )
    vectors = []
    masses = []
    for index, mass in enumerate(meson_masses):
        coupling = couplings[index]
        if coupling == 0:
            vector = sp.zeros(count + 1, 1)
            vector[index] = 1
            vectors.append(vector)
            masses.append(mass)
            continue
        gap = mass - alp_mass_squared
        shift = coupling**2 / gap
        vector = sp.zeros(count + 1, 1)
        vector[index] = 1 - epsilon**2 * coupling**2 / (2 * gap**2)
        mixing_sum = 0
        for other in range(count):
            if other == index or couplings[other] == 0:
                continue
            value = couplings[other] * coupling / ((mass - meson_masses[other]) * gap)
            vector[other] = epsilon**2 * value
            mixing_sum += couplings[other] * value
        vector[count] = (-epsilon * coupling / gap
                         + epsilon**3 * (coupling * (shift + fixed_mass) / gap**2
                         - mixing_sum / gap + coupling**3 / (2 * gap**3)))
        vectors.append(vector)
        masses.append(mass + epsilon**2 * shift)
    return tuple(meson_vectors), tuple(masses) + (alp_mass_squared,), tuple(vectors) + (alp_vector,)


def solve_crossing(loop_coefficient, tree_amplitude, maximum_ratio, bracket):
    """Return epsilon and the signed residual for a bracketed positive crossing."""
    lower, upper = bracket
    if not (0 < lower < upper and loop_coefficient >= 0 and maximum_ratio > 0):
        raise ValueError("Require a positive ordered bracket and nonnegative loop coefficient.")

    def residual(epsilon):
        value = loop_coefficient * epsilon**3 - maximum_ratio * abs(tree_amplitude(epsilon))
        if not math.isfinite(value):
            raise ValueError("The crossing equation is not finite inside the bracket.")
        return value

    root = brentq(residual, lower, upper)
    return root, residual(root)


def calculate_loop_coefficient(invariants, alp_mass_squared, daughters, s12, s23,
                               loop_scale_squared):
    """Compute the leading mixing-dependent loop magnitude from dressed sources."""
    epsilon = sp.Symbol("epsilon", positive=True)
    meson_matrix, source = invariant_quadratic_blocks(invariants)
    meson_vectors, masses, eigenvectors = pole_expansions(
        meson_matrix, source, alp_mass_squared, epsilon)
    transform = sp.diag(sp.Matrix.hstack(*meson_vectors), sp.ones(1))
    eigenvectors = tuple(transform * vector for vector in eigenvectors)
    zero = sp.zeros(3)
    basis = tuple((generator, 0) for generator in GENERATORS) + ((zero, 1),)
    pole_directions = {}
    for index, vector in enumerate(eigenvectors):
        meson = sum((vector[j] * GENERATORS[j] for j in range(9)), sp.zeros(3))
        pole_directions[index] = (meson, vector[9])
    pole_masses = dict(enumerate(masses))
    vertex = dressed_tree_vertex(basis, invariants, epsilon)
    result = on_shell_cubic_estimate(
        masses, eigenvectors, (9, *daughters), pole_directions, pole_masses,
        basis, invariants, epsilon, 4, formal_scalar_moment(ScalarMoment),
        vertex, vertex, vertex, s12, s23)
    def moment_scale(term):
        if isinstance(term, ScalarMoment):
            return term.args[1]
        return loop_scale_squared

    size, _ = crossing_condition(result[0], epsilon, epsilon, sp.Integer(1), moment_scale)
    return sp.simplify(size / epsilon**3)


def main():
    """Read input JSON and solve the crossing within the supplied bracket."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="JSON file with invariant source arrays and crossing inputs")
    arguments = parser.parse_args()
    with open(arguments.input, encoding="utf-8") as file:
        data = json.load(file)
    invariants = {
        "B_0": sp.Rational(str(data["B_0"])),
        "f_pi": sp.Rational(str(data["f_pi"])),
        "anomaly_mass_squared": sp.Rational(str(data["m_0_squared"])),
        "S": [sp.Matrix(matrix).applyfunc(lambda value: sp.Rational(str(value)))
              for matrix in data["S"]],
        "P": [sp.Matrix(matrix).applyfunc(lambda value: sp.Rational(str(value)))
              for matrix in data["P"]],
        "theta": [sp.Rational(str(value)) for value in data["theta"]],
    }
    if min(len(invariants[name]) for name in ("S", "P", "theta")) < 7:
        raise ValueError("S, P and theta must specify derivatives through order six.")
    if any(matrix.shape != (3, 3) for name in ("S", "P") for matrix in invariants[name]):
        raise ValueError("Every S and P derivative must be a 3 by 3 matrix.")
    if any(matrix != matrix.T for name in ("S", "P") for matrix in invariants[name]):
        raise ValueError("This implementation requires real symmetric S and P matrices.")
    if invariants["S"][0] != sp.diag(*invariants["S"][0].diagonal()):
        raise ValueError("The exact tree amplitude requires a diagonal S[0].")
    alp_mass_squared = sp.Rational(str(data["m_a_squared"]))
    s12, s23 = (sp.Rational(str(value)) for value in (data["s12"], data["s23"]))
    loop_scale_squared = sp.Rational(str(data["loop_scale_squared"]))
    daughters = tuple(data["daughter_indices"])
    if len(daughters) != 3 or any(type(index) is not int or index not in range(9)
                                   for index in daughters):
        raise ValueError("daughter_indices must contain three integers between 0 and 8.")
    if min(invariants["f_pi"], invariants["B_0"], alp_mass_squared,
           loop_scale_squared, sp.Rational(str(data["R_max"]))) <= 0:
        raise ValueError("f_pi, B_0, m_a_squared, loop_scale_squared and R_max must be positive.")
    if len(data["bracket"]) != 2 or not 0 < data["bracket"][0] < data["bracket"][1]:
        raise ValueError("bracket must contain two increasing positive numbers.")
    if all(matrix == sp.zeros(3) for name in ("S", "P") for matrix in invariants[name][1:]) and all(
            value == 0 for value in invariants["theta"][1:]):
        raise ValueError("The ALP is decoupled and the tree amplitude vanishes.")
    coefficient = calculate_loop_coefficient(
        invariants, alp_mass_squared, daughters, s12, s23, loop_scale_squared)
    matrix, source = invariant_quadratic_blocks(invariants)

    def tree(epsilon):
        return exact_tree_amplitude(invariants, epsilon, float(alp_mass_squared),
                                    daughters, float(s12), float(s23), matrix, source)[0]

    root, residual = solve_crossing(float(coefficient), tree,
                                    float(data["R_max"]), tuple(data["bracket"]))
    print(json.dumps({"f_pi_over_f_a": root, "loop_coefficient": float(coefficient),
                      "residual": residual}, indent=2))


class CrossingTest(unittest.TestCase):
    def test_bracketed_crossing(self):
        root, residual = solve_crossing(2.0, lambda epsilon: epsilon, 1.0, (0.1, 1.0))
        self.assertAlmostEqual(root, 1 / math.sqrt(2))
        self.assertAlmostEqual(residual, 0)

    def test_symbolic_pole_equations(self):
        epsilon = sp.Symbol("epsilon", real=True)
        first_mass, second_mass, alp_mass = sp.symbols("m_1 m_2 m_a", real=True)
        first_source, second_source = sp.symbols("c_1 c_2", real=True)
        source = sp.Matrix([first_source, second_source])
        _, masses, vectors = pole_expansions(
            sp.diag(first_mass, second_mass), source, alp_mass, epsilon)
        squared_momentum = sp.Symbol("s")
        fixed_mass = first_source**2 / (alp_mass - first_mass) + second_source**2 / (alp_mass - second_mass)
        inverse = sp.diag(squared_momentum - first_mass, squared_momentum - second_mass,
                          squared_momentum - alp_mass + epsilon**2 * fixed_mass)
        inverse += epsilon * sp.Matrix([[0, 0, first_source], [0, 0, second_source],
                                        [first_source, second_source, 0]])
        for mass, vector in zip(masses, vectors):
            residual = inverse.subs(squared_momentum, mass) * vector
            self.assertEqual(residual.applyfunc(
                lambda entry: sp.simplify(sp.series(entry, epsilon, 0, 4).removeO())),
                sp.zeros(3, 1))


if __name__ == "__main__":
    main()


