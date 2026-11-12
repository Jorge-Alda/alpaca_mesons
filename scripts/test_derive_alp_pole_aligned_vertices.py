import unittest
import random

import sympy as sp

from derive_alp_lo_joint_vertices import quadratic_blocks
from derive_alp_pole_aligned_vertices import alignment_matrix, aligned_amplitude, cubic_map, nonlinear_contact_parts, transformed_vertex_parts
from derive_alp_pole_projected_amplitudes import ccwz_vertex, projected_amplitude
from solve_alp_meson_spectrum import solve_spectrum_exact
from symbolic_utils import u3_generators


class PoleAlignedVerticesTest(unittest.TestCase):
    def setUp(self):
        self.basis = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]),
                      sp.diag(1, -1), sp.eye(2)]
        zero = sp.zeros(2)
        r = self.basis[2]/2
        # A separate linear ALP counterterm cancels the chi_1 tadpole.
        self.sources = ([sp.eye(2), sp.eye(2), self.basis[2], self.basis[0], zero],
                        [r, self.basis[0]/3, zero, zero, zero],
                        [-r, self.basis[0]/3, zero, zero, zero],
                        [sp.S.Zero]*5, sp.S.One, sp.Integer(2), sp.S.Zero,
                        sp.Rational(236, 15))
        self.kinetic, self.mass = quadratic_blocks(self.basis, *self.sources)
        self.spectrum = solve_spectrum_exact(self.kinetic, self.mass)
        self.alignment = alignment_matrix(self.spectrum.pole_vector(4))

    def test_cubic_map_and_exceptional_pole(self):
        directions = [sp.eye(5)[:, i] for i in (4, 0, 0)]
        result = cubic_map(self.basis, self.alignment, directions, 1)
        self.assertNotEqual(result, sp.zeros(5, 1))
        self.assertEqual(result[-1], 0)
        a, b, c = sp.symbols('a b c', real=True)
        meson = b*self.basis[0] + c*self.basis[0]
        shift = sum((self.alignment[i, -1]*self.basis[i] for i in range(4)), sp.zeros(2))
        commutator = lambda x, y: x*y-y*x
        coordinate = (a*a*commutator(shift, commutator(shift, meson))/24
                      + a*commutator(meson, commutator(shift, meson))/12)
        derivative = coordinate.diff(a, b, c)
        expected = sp.Matrix([sp.trace(generator*derivative)/2 for generator in self.basis] + [0])
        self.assertEqual(sp.simplify(result-expected), sp.zeros(5, 1))
        directions = [sp.eye(5)[:, i] for i in (4, 4, 0)]
        result = cubic_map(self.basis, self.alignment, directions, 1)
        meson = c*self.basis[0]
        alp = a+b
        coordinate = (alp**2*commutator(shift, commutator(shift, meson))/24
                      + alp*commutator(meson, commutator(shift, meson))/12)
        derivative = coordinate.diff(a, b, c)
        expected = sp.Matrix([sp.trace(generator*derivative)/2 for generator in self.basis] + [0])
        self.assertNotEqual(result, sp.zeros(5, 1))
        self.assertEqual(sp.simplify(result-expected), sp.zeros(5, 1))
        with self.assertRaises(ValueError):
            alignment_matrix(sp.Matrix([1, 0, 0, 0, 0]))

    def test_same_source_on_shell_and_off_shell_contact(self):
        # Parent pole 16 and a daughter at 59/60 both have ALP components.
        indices = [4, 1, 1, 0]
        poles = [sp.Matrix(self.spectrum.pole_vector(i)) for i in indices]
        self.assertNotEqual(poles[-1][-1], 0)
        directions = [self.alignment.inv()*v for v in poles]
        momenta = [sp.eye(4)[:, i] for i in range(4)]
        masses = [self.spectrum.pole_squared_masses[i] for i in indices]
        # Three independent scalar products define an on-shell four-leg Gram matrix.
        gram = sp.zeros(4)
        for i in range(3):
            gram[i, i] = masses[i]
        gram[0, 1] = gram[1, 0] = -sp.Rational(2)
        gram[0, 2] = gram[2, 0] = -sp.Rational(3)
        gram[1, 2] = gram[2, 1] = (masses[3] - sum(masses[:3]))/2 - gram[0, 1] - gram[0, 2]
        for i in range(3):
            gram[i, 3] = gram[3, i] = -sum(gram[i, j] for j in range(3))
        gram[3, 3] = masses[3]
        self.assertEqual(gram*sp.ones(4, 1), sp.zeros(4, 1))
        linear, nonlinear = transformed_vertex_parts(self.basis, self.sources, self.alignment, directions, momenta, gram)
        original = ccwz_vertex(self.basis, self.sources, poles, momenta, gram)
        self.assertEqual(sp.simplify(linear-original), 0)
        self.assertEqual(sp.simplify(nonlinear), 0)
        kinetic_part, mass_part = nonlinear_contact_parts(self.basis, self.sources, self.alignment, directions, momenta, gram)
        self.assertNotEqual(kinetic_part, 0)
        self.assertEqual(sp.simplify(kinetic_part+mass_part), 0)
        off_shell = gram.copy()
        off_shell[0, 0] += 1
        for i in range(1, 4):
            off_shell[0, i] -= sp.Rational(1, 3)
            off_shell[i, 0] -= sp.Rational(1, 3)
        for i in range(1, 4):
            for j in range(i+1, 4):
                off_shell[i, j] += sp.Rational(1, 6)
                off_shell[j, i] += sp.Rational(1, 6)
        _, nonlinear_off = transformed_vertex_parts(self.basis, self.sources, self.alignment, directions, momenta, off_shell)
        self.assertNotEqual(sp.simplify(nonlinear_off), 0)

    def test_two_and_three_daughters(self):
        def original(directions, momenta, gram):
            return ccwz_vertex(self.basis, self.sources, directions, momenta, gram)
        indices3 = [4, 1, 1]
        squared = [self.spectrum.pole_squared_masses[i] for i in indices3]
        gram3 = sp.zeros(3)
        for i in range(3):
            gram3[i, i] = squared[i]
        for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            gram3[i, j] = gram3[j, i] = (squared[k]-squared[i]-squared[j])/2
        amplitude3 = projected_amplitude(self.spectrum, indices3, gram3, original)
        self.assertNotEqual(amplitude3, 0)
        self.assertEqual(sp.simplify(aligned_amplitude(self.basis, self.sources, self.spectrum, indices3, gram3)
                                     - amplitude3), 0)
        # Two deterministic, nonexceptional kinematic points; all legs on shell.
        randomizer = random.Random(31415)
        indices4 = [4, 1, 1, 0]
        squared = [self.spectrum.pole_squared_masses[i] for i in indices4]
        for _ in range(2):
            gram = sp.zeros(4)
            for i in range(3):
                gram[i, i] = squared[i]
            gram[0, 1] = gram[1, 0] = -sp.Rational(randomizer.randrange(2, 5))
            gram[0, 2] = gram[2, 0] = -sp.Rational(randomizer.randrange(2, 5))
            gram[1, 2] = gram[2, 1] = (squared[3]-sum(squared[:3]))/2-gram[0, 1]-gram[0, 2]
            for i in range(3):
                gram[i, 3] = gram[3, i] = -sum(gram[i, j] for j in range(3))
            gram[3, 3] = squared[3]
            result = projected_amplitude(self.spectrum, indices4, gram, original)
            directions = [sp.Matrix(self.spectrum.pole_vector(i)) for i in indices4]
            momenta = [sp.eye(4)[:, i] for i in range(4)]
            self.assertNotEqual(sp.simplify(result-original(directions, momenta, gram)), 0)
            self.assertEqual(sp.simplify(aligned_amplitude(self.basis, self.sources, self.spectrum, indices4, gram)
                                         - result), 0)

    def test_u3_strangeness_changing_vector_source(self):
        basis = u3_generators()
        zero = sp.zeros(3)
        real_ds, imaginary_ds = sp.symbols('real_ds imaginary_ds', real=True)
        vector_source = real_ds*basis[6]+imaginary_ds*basis[7]
        sources = ([sp.diag(1, 2, 3), zero, zero, zero],
                   [vector_source, zero, zero, zero],
                   [vector_source, zero, zero, zero], [sp.S.Zero]*4,
                   sp.S.One, sp.Integer(2), sp.S.Zero, sp.Integer(4))
        alignment = sp.eye(10)
        alignment[3, 9] = sp.Rational(1, 3)
        directions = [sp.eye(10)[:, i] for i in (9, 1, 4)]
        momenta = [sp.eye(3)[:, i] for i in range(3)]
        gram = sp.Matrix([[9, -3, -6], [-3, 1, 2], [-6, 2, 4]])
        linear, nonlinear = transformed_vertex_parts(basis, sources, alignment, directions, momenta, gram)
        direct = ccwz_vertex(basis, sources, [alignment*d for d in directions], momenta, gram)
        self.assertEqual(sp.simplify(linear-direct), 0)
        self.assertEqual(nonlinear, 0)
        self.assertNotEqual(sp.diff(linear, imaginary_ds), 0)


if __name__ == '__main__':
    unittest.main()
