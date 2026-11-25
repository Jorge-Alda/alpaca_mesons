"""Dressed-source mass matrix and exact tree pole projection.

Masses are in GeV, and squared masses and Mandelstam invariants in GeV².
The source dictionary follows derive_alp_finite_mixing_three_meson_tree.
"""
import numpy as np
import sympy as sp
import unittest
from scipy.optimize import linear_sum_assignment

from derive_alp_finite_mixing_three_meson_tree import (
    GENERATORS, matrix_C, pole_operator, three_meson_amplitude,
)


def invariant_quadratic_blocks(invariants):
    """Return the normalized nine-meson mass matrix and ALP mixing source."""
    meson_mass_matrix = sp.Matrix(9, 9, lambda i, j: sp.simplify(
        (GENERATORS[i] * pole_operator(GENERATORS[j], invariants, 0)).trace() / 2))
    source = matrix_C(invariants)
    mixing_source = sp.Matrix([sp.simplify((generator * source).trace())
                               for generator in GENERATORS])
    return meson_mass_matrix, mixing_source


def exact_tree_amplitude(invariants, epsilon, alp_mass_squared, daughter_indices,
                         s12, s23, meson_mass_matrix, mixing_source):
    """Evaluate exact tree amplitude using maximum-overlap pole assignment.

    Daughter indices refer to the sorted unmixed mass eigenvectors. Use an
    epsilon bracket on which this assignment is continuous. The fixed ALP pole
    determines the ALP diagonal mass entry through the Schur equation.
    """
    matrix = np.asarray(meson_mass_matrix, dtype=float)
    source = np.asarray(mixing_source, dtype=float).reshape(-1)
    unmixed_masses, unmixed_vectors = np.linalg.eigh(matrix)
    alp_diagonal = alp_mass_squared + epsilon**2 * source @ np.linalg.solve(
        matrix - alp_mass_squared * np.eye(9), source)
    mass_matrix = np.zeros((10, 10))
    mass_matrix[:9, :9] = matrix
    mass_matrix[:9, 9] = -epsilon * source
    mass_matrix[9, :9] = -epsilon * source
    mass_matrix[9, 9] = alp_diagonal
    masses, vectors = np.linalg.eigh(mass_matrix)
    references = np.eye(10)
    references[:9, :9] = unmixed_vectors
    overlap = references.T @ vectors
    rows, columns = linear_sum_assignment(-np.abs(overlap)**2)
    assignment = dict(zip(rows, columns))
    indices = [9, *daughter_indices]
    legs = []
    external_masses = []
    for index in indices:
        column = assignment[index]
        direction = vectors[:, column]
        if overlap[index, column] < 0:
            direction = -direction
        meson_matrix = sum((sp.Float(direction[i]) * GENERATORS[i] for i in range(9)), sp.zeros(3))
        legs.append((meson_matrix, sp.Float(direction[9])))
        external_masses.append(float(masses[column]))
    bare_mass_squared = (alp_diagonal + epsilon**2 * float(invariants['B_0'] * invariants['S'][2].trace())
                         - 2 * epsilon**2 * float(invariants['anomaly_mass_squared'] * invariants['theta'][1]**2 / 3))
    mandelstam = {(1, 2): s12, (2, 3): s23,
                  (1, 3): sum(external_masses) - s12 - s23}
    result = three_meson_amplitude(invariants, epsilon, bare_mass_squared,
                                   legs, external_masses, mandelstam)['total']
    return complex(sp.N(result)), external_masses, assignment


class InvariantQuadraticBlocksTest(unittest.TestCase):
    def test_symbolic_sources(self):
        quark_up, quark_down, quark_strange, b0, anomaly = sp.symbols(
            'm_u m_d m_s B_0 m_0_squared', real=True)
        pseudoscalar, theta = sp.symbols('p theta', real=True)
        invariants = {
            'B_0': b0,
            'anomaly_mass_squared': anomaly,
            'S': [sp.diag(quark_up, quark_down, quark_strange)],
            'P': [sp.zeros(3), sp.diag(pseudoscalar, 0, 0)],
            'theta': [0, theta],
        }
        meson_mass_matrix, mixing_source = invariant_quadratic_blocks(invariants)
        self.assertEqual(meson_mass_matrix.shape, (9, 9))
        self.assertEqual(meson_mass_matrix, meson_mass_matrix.T)
        self.assertEqual(sp.simplify(sum(value**2 for value in mixing_source)
                                     - 2 * (matrix_C(invariants)**2).trace()), 0)

    def test_fixed_alp_pole_schur_complement(self):
        meson_mass_squared, alp_mass_squared, mixing, epsilon = sp.symbols(
            'm_n_squared m_a_squared coupling epsilon', real=True, nonzero=True)
        alp_diagonal = alp_mass_squared + epsilon**2 * mixing**2 / (
            meson_mass_squared - alp_mass_squared)
        mass_matrix = sp.Matrix([[meson_mass_squared, -epsilon * mixing],
                                 [-epsilon * mixing, alp_diagonal]])
        self.assertEqual(sp.factor((alp_mass_squared * sp.eye(2) - mass_matrix).det()), 0)


if __name__ == '__main__':
    unittest.main()
