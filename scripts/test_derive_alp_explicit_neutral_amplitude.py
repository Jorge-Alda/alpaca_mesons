"""Independent check of the explicit cubic formula against the joint CCWZ expansion."""
import unittest
import sympy as sp
from derive_alp_explicit_neutral_amplitude import (
    aligned_three_pole_amplitude, aligned_cubic_terms, cubic_amplitude, pole_overlap,
    overlap_flavor_matrix, decoupled_kaon_component_amplitude,
    aligned_component_terms,
    aligned_sources,
    aligned_vector_scalar_sources,
)
from derive_alp_lo_joint_vertices import vertex, quadratic_blocks


class ExplicitCubicTest(unittest.TestCase):
    def test_same_source_u3_spectrum_and_aligned_components(self):
        import numpy as np
        from scipy.linalg import eigh
        from symbolic_utils import u3_generators
        basis = [u3_generators()[i] for i in (0, 3, 6, 7, 8)]
        zero = sp.zeros(3)
        diagonal = sp.diag(sp.Rational(1, 3), sp.Rational(1, 4), sp.Rational(1, 5))
        ds = sp.Matrix([[0, 0, 0], [0, 0, 1+sp.I], [0, 1-sp.I, 0]])
        chi = [sp.diag(sp.Rational(3, 10), sp.Rational(14, 25), sp.Rational(8, 5)),
               sp.I*(diagonal+ds/5), diagonal/7, sp.I*diagonal/11]
        vector = [ds/8, zero, zero, zero]
        axial = [diagonal/5+ds/20, zero, zero, zero]
        right = [v+a for v, a in zip(vector, axial)]
        left = [v-a for v, a in zip(vector, axial)]
        theta = [0, sp.Rational(1, 10), sp.Rational(1, 15), 0]
        f_pi, f_a, anomaly, bare = sp.Integer(1), sp.Integer(8), sp.Rational(7, 10), sp.Rational(23, 10)
        kinetic, mass = quadratic_blocks(basis, chi, right, left, theta, f_pi, f_a, anomaly, bare)
        values, vectors = eigh(np.array(mass, dtype=float), np.array(kinetic, dtype=float))
        parent = int(np.argmax(np.abs(vectors[-1, :])))
        others = [i for i in range(len(values)) if i != parent and abs(vectors[-1, i]) > 1e-7]
        self.assertGreaterEqual(len(others), 2)
        daughters = sorted(others, key=lambda i: abs(vectors[2, i])+abs(vectors[3, i]), reverse=True)[:2]
        indices = [parent]+daughters
        directions = []
        for index in indices:
            direction = vectors[:, index]
            if direction[-1] < 0:
                direction = -direction
            directions.append(sp.Matrix(direction))
        matrices = [sum((direction[i]*basis[i] for i in range(len(basis))), zero)
                    for direction in directions]
        s_a, s_1, s_2 = [sp.Float(values[i], 16) for i in indices]
        for squared_mass, direction, matrix in zip((s_a, s_1, s_2), directions, matrices):
            inverse = overlap_flavor_matrix(squared_mass, list(chi[0].diagonal()), chi[1],
                                            right[0], left[0], theta[1], f_pi, f_a, anomaly)
            difference = inverse-matrix/direction[-1]
            self.assertLess(max(abs(complex(sp.N(value))) for value in difference), 1e-8)
        gram = sp.Matrix([[s_a, (s_2-s_a-s_1)/2, (s_1-s_a-s_2)/2],
                          [(s_2-s_a-s_1)/2, s_1, (s_a-s_1-s_2)/2],
                          [(s_1-s_a-s_2)/2, (s_a-s_1-s_2)/2, s_2]])
        original = cubic_amplitude(matrices, [d[-1] for d in directions], gram,
                                    chi, right, left, theta, f_pi, f_a, anomaly)
        parent_overlap = matrices[0]/directions[0][-1]
        aligned_matrices = [matrices[i]-directions[i][-1]*parent_overlap for i in (1, 2)]
        chi_p, right_p, left_p, theta_p = aligned_sources(parent_overlap, chi, right, left, theta, f_pi, f_a)
        vector_p = [(r+l)/2 for r, l in zip(right_p, left_p)]
        axial_p = [(r-l)/2 for r, l in zip(right_p, left_p)]
        scalar_p = [(c+c.conjugate().T)/4 for c in chi_p]
        pseudoscalar_p = [(c-c.conjugate().T)/(4*sp.I) for c in chi_p]
        # Divide each daughter direction by its ALP component and restore all
        # three K-normalized external ALP components after evaluation.
        actual = sum(aligned_component_terms(aligned_matrices[0]/directions[1][-1],
                     aligned_matrices[1]/directions[2][-1], gram, vector_p, axial_p,
                     scalar_p, pseudoscalar_p, theta_p, 1, f_pi, f_a, anomaly))
        actual *= sp.prod(direction[-1] for direction in directions)
        self.assertLess(abs(complex(sp.N(actual-original))), 1e-10)

    def test_vector_axial_scalar_pseudoscalar_components(self):
        matrices = [sp.Matrix([[1, 2, 0], [2, 0, 1+sp.I], [0, 1-sp.I, -1]]),
                    sp.Matrix([[0, 1, sp.I], [1, 2, 0], [-sp.I, 0, 0]]),
                    sp.Matrix([[1, 0, 2], [0, -1, sp.I], [2, -sp.I, 1]])]
        scalar = [matrices[0], matrices[1], matrices[2], matrices[0]]
        pseudoscalar = [matrices[1], matrices[2], matrices[0], matrices[1]]
        vector = [matrices[0], matrices[1]]
        axial = [matrices[2], matrices[0]]
        B_0, f_pi, f_a, anomaly = map(sp.Integer, (7, 2, 3, 5))
        chi = [2*B_0*(s+sp.I*p) for s, p in zip(scalar, pseudoscalar)]
        right = [v+a for v, a in zip(vector, axial)]
        left = [v-a for v, a in zip(vector, axial)]
        theta = [0, sp.Rational(2, 3), sp.Rational(3, 5)]
        momenta = sp.Matrix([[3, 2, 1], [2, 5, 4], [1, 4, 6]])
        expected = aligned_cubic_terms(matrices[1], matrices[2], 1, 1,
                    momenta, chi, right, left, theta, f_pi, f_a, anomaly)
        actual = aligned_component_terms(matrices[1], matrices[2], momenta,
                    vector, axial, scalar, pseudoscalar, theta, B_0,
                    f_pi, f_a, anomaly)
        self.assertEqual(tuple(sp.simplify(a-b) for a, b in zip(actual, expected)), (0, 0, 0))
        chi_p, right_p, left_p, _ = aligned_sources(matrices[0], chi, right, left,
                                  theta, f_pi, f_a)
        vector_p, axial_p, scalar_p, pseudoscalar_p = aligned_vector_scalar_sources(
                                  matrices[0], vector, axial, scalar, pseudoscalar, f_pi, f_a)
        for n in range(2):
            self.assertEqual(sp.simplify(vector_p[n]-(right_p[n]+left_p[n])/2), sp.zeros(3))
            self.assertEqual(sp.simplify(axial_p[n]-(right_p[n]-left_p[n])/2), sp.zeros(3))
        for n in range(4):
            self.assertEqual(sp.simplify(scalar_p[n]-(chi_p[n]+chi_p[n].conjugate().T)/(4*B_0)), sp.zeros(3))
            self.assertEqual(sp.simplify(pseudoscalar_p[n]-(chi_p[n]-chi_p[n].conjugate().T)/(4*sp.I*B_0)), sp.zeros(3))
        # Stationary CP-even mass source: P_0=0. Check the compact derivatives
        # used for the physical pole-aligned three-leg amplitude.
        pseudoscalar[0] = sp.zeros(3)
        vector_p, axial_p, scalar_p, pseudoscalar_p = aligned_vector_scalar_sources(
                                  matrices[0], vector, axial, scalar, pseudoscalar, f_pi, f_a)
        anticommutator = lambda a, b: a*b+b*a
        commutator = lambda a, b: a*b-b*a
        coefficient = f_a/(2*f_pi)
        checks = [vector_p[0]-vector[0], axial_p[0]-axial[0]+coefficient*matrices[0],
                  axial_p[1]-axial[1]+sp.I*coefficient*commutator(matrices[0], vector[0]),
                  scalar_p[1]-scalar[1],
                  pseudoscalar_p[2]-pseudoscalar[2]+2*coefficient*anticommutator(matrices[0], scalar[1]),
                  scalar_p[3]-scalar[3]-3*coefficient*anticommutator(matrices[0], pseudoscalar[2])
                    +3*coefficient**2*anticommutator(matrices[0], anticommutator(matrices[0], scalar[1]))]
        self.assertTrue(all(sp.simplify(check) == sp.zeros(3) for check in checks))

    def test_restricted_complex_ds_vector_and_bose_symmetry(self):
        zero = sp.zeros(3)
        first = sp.diag(1, -1, 0)
        second = sp.Matrix([[0, 0, 0], [0, 0, 1+sp.I], [0, 1-sp.I, 0]])
        vector = sp.Matrix([[0, 0, 0], [0, 0, 2+3*sp.I], [0, 2-3*sp.I, 0]])
        B0, mu, md, ms = sp.symbols('B0 m_u m_d m_s', real=True)
        chi = [2*B0*sp.diag(mu, md, ms), vector, sp.I*vector, vector]
        right = [vector, zero, zero, zero]
        left = [vector, zero, zero, zero]
        theta = [0, 1, 2, 0]
        q = sp.symbols('q0:3', real=True)
        gram = sp.Matrix([[q[i]*q[j] for j in range(3)] for i in range(3)])
        first_alp, second_alp = sp.Rational(2, 3), sp.Rational(3, 5)
        parameters = tuple(map(sp.Integer, (2, 3, 5)))
        explicit = sum(aligned_cubic_terms(first, second, first_alp, second_alp,
                       gram, chi, right, left, theta, *parameters))
        direct = vertex([zero, first, second], [1, first_alp, second_alp], q,
                        chi, right, left, theta, *parameters, sp.Integer(7))
        self.assertEqual(sp.expand(explicit-direct), 0)
        swapped = sum(aligned_cubic_terms(second, first, second_alp, first_alp,
                      gram.extract([0, 2, 1], [0, 2, 1]), chi, right, left,
                      theta, *parameters))
        self.assertEqual(sp.expand(explicit-swapped), 0)
        mixing_6, mixing_7 = sp.Rational(3, 5), sp.Rational(4, 5)
        kaon = sp.Matrix([[0, 0, 0], [0, 0, -mixing_7-sp.I*mixing_6],
                          [0, -mixing_7+sp.I*mixing_6, 0]])
        actual = decoupled_kaon_component_amplitude(first, mixing_6, mixing_7,
                                                   gram, chi, right, left, *parameters[:2])
        expected = sum(aligned_cubic_terms(kaon, first, 0, 1, gram, chi, right,
                                          left, theta, *parameters))
        self.assertEqual(sp.expand(actual-expected), 0)

    def test_general_source_component_inverse(self):
        basis = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]),
                 sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
        zero = sp.zeros(2)
        chi = [sp.diag(2, 5), sp.Matrix([[1+sp.I, 2+3*sp.I], [1-sp.I, 2-sp.I]]),
               sp.diag(1, 2)]
        right = [sp.Matrix([[1, 2+sp.I], [2-sp.I, 3]]), zero, zero]
        left = [sp.Matrix([[2, sp.I], [-sp.I, 1]]), zero, zero]
        theta = [0, sp.Rational(2, 3), 0]
        f_pi, f_a, anomaly, bare = map(sp.Integer, (2, 3, 7, 11))
        kinetic, mass = quadratic_blocks(basis, chi, right, left, theta,
                                         f_pi, f_a, anomaly, bare)
        squared_mass = sp.Rational(1, 3)
        components, _ = pole_overlap(kinetic, mass, squared_mass)
        from_blocks = sum((components[i]*basis[i] for i in range(4)), zero)
        explicit = overlap_flavor_matrix(squared_mass, [sp.Integer(2), sp.Integer(5)],
                                         chi[1], right[0], left[0], theta[1],
                                         f_pi, f_a, anomaly)
        self.assertEqual(sp.simplify(from_blocks-explicit), zero)

    def test_aligned_poles_and_zero_alp_daughter(self):
        kinetic = sp.eye(3)
        mass_squared = sp.Matrix([[2, 0, 1], [0, 5, 0], [1, 0, 2]])
        basis = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]),
                 sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
        kinetic = sp.eye(5)
        mass_squared = sp.eye(5)
        zero = sp.zeros(2)
        chi = [sp.diag(2, 3), sp.Matrix([[0, sp.I], [sp.I, 0]]),
               sp.diag(1, -1), sp.diag(2, 1)]
        right = [sp.Matrix([[0, 1], [1, 0]]), sp.eye(2), zero, zero]
        left = [zero]*4
        theta = [0, 1, 2, 0]
        q = sp.symbols('q0:3', real=True)
        gram = sp.Matrix([[q[i]*q[j] for j in range(3)] for i in range(3)])
        parameters = tuple(map(sp.Integer, (2, 3, 5)))
        got = aligned_three_pole_amplitude(
            basis, kinetic, mass_squared, [0, 1, 5], gram, chi, right, left,
            theta, *parameters, zero_alp_directions={2: sp.Matrix([0, 1, 0, 0])})
        matrices = [overlap_flavor_matrix(s, [sp.Integer(2), sp.Integer(3)],
                    chi[1], right[0], left[0], theta[1], *parameters)
                    for s in (0, 1)] + [basis[1]]
        norms = [sp.sqrt(1+sum((sp.trace(gen*matrix)/2)**2 for gen in basis))
                 for matrix in matrices[:2]]
        expected = vertex(matrices, [1, 1, 0], q, chi, right, left, theta,
                          *parameters, sp.Integer(7))/sp.prod(norms)
        self.assertEqual(sp.expand(got-expected), 0)

    def test_general_complex_sources_and_noncommuting_directions(self):
        matrices = [sp.Matrix([[1, 2, 0], [2, 0, 1], [0, 1, -1]]),
                    sp.Matrix([[0, 1, sp.I], [1, 2, 0], [-sp.I, 0, 0]]),
                    sp.Matrix([[1, 0, 2], [0, -1, sp.I], [2, -sp.I, 1]])]
        zero = sp.zeros(3)
        chi = [sp.diag(2, 3, 5), matrices[0]+sp.I*matrices[1],
               matrices[2]+sp.I*matrices[0], matrices[1]+sp.I*matrices[2]]
        right = [matrices[0], matrices[1], zero, zero]
        left = [matrices[2], matrices[0], zero, zero]
        theta = [0, sp.Rational(2, 3), sp.Rational(3, 5), 0]
        alps = [sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(-1, 4)]
        q = sp.symbols('q0:3', real=True)
        f_pi, f_a, m_0_squared, bare_mass_squared = map(sp.Integer, (2, 3, 5, 7))
        expected = vertex(matrices, alps, q, chi, right, left, theta, f_pi, f_a, m_0_squared, bare_mass_squared)
        actual = cubic_amplitude(matrices, alps, q, chi, right, left, theta, f_pi, f_a, m_0_squared)
        self.assertEqual(sp.expand(actual-expected), 0)
        gram = sp.Matrix([[q[i]*q[j] for j in range(3)] for i in range(3)])
        terms = aligned_cubic_terms(matrices[1], matrices[2], alps[1], alps[2],
                                    gram, chi, right, left, theta, f_pi, f_a,
                                    m_0_squared)
        aligned = cubic_amplitude([zero, matrices[1], matrices[2]],
                                  [1, alps[1], alps[2]], gram, chi, right, left,
                                  theta, f_pi, f_a, m_0_squared)
        self.assertEqual(sp.expand(sum(terms)-aligned), 0)


if __name__ == '__main__':
    unittest.main()
