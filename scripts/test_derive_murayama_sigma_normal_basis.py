import unittest
import numpy as np
import sympy as sp
from scripts.derive_murayama_sigma_normal_basis import derive


class NormalBasisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_shift_removes_mixed_inverse_propagator_at_selected_pole(self):
        r = self.result
        self.assertTrue(all(sp.simplify(x) == 0 for x in
            (r['meson'] - r['s'] * sp.eye(3))*r['shift'] + r['mixing']))

    def test_residues_match_eigenvectors(self):
        r = self.result
        numeric = {r['B0']: 2, r['mu']: sp.Rational(1, 5), r['md']: sp.Rational(1, 3),
                   r['ms']: 2, r['m1']: 1, r['fpi']: sp.Rational(1, 4), r['fphi']: sp.Rational(2, 3)}
        matrix = np.array(r['mass'].subs(numeric), dtype=float)
        masses, vectors = np.linalg.eigh(matrix)
        polynomial = sp.lambdify(r['s'], r['meson_polynomial'].subs(numeric), 'numpy')
        for index, mass in enumerate(masses):
            residue = polynomial(mass) / np.prod(np.delete(mass-masses, index))
            self.assertAlmostEqual(residue, vectors[3, index]**2, places=10)


if __name__ == '__main__':
    unittest.main()
