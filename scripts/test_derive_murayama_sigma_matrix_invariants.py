"""Checks for the radial neutral vertex and matrix field covariance."""
import unittest
import sympy as sp
from scripts.derive_murayama_sigma_matrix_invariants import derive


class MatrixVertexTests(unittest.TestCase):
    def test_radial_source_and_polar_kinetic_term(self):
        mass_matrix, radial_mass, radial_kinetic, vertex = derive()
        fa = next(symbol for symbol in mass_matrix.free_symbols if str(symbol) == 'f_a')
        self.assertEqual(sp.simplify(radial_mass - mass_matrix/fa), sp.zeros(4))
        self.assertEqual(sp.simplify(radial_kinetic - sp.diag(0, 0, 0, 2/fa)), sp.zeros(4))
        self.assertEqual(vertex.shape, (4, 4))

    def test_pole_residue_contraction_under_nonorthogonal_change(self):
        # Two distinct poles with rank-one residues, and a nonorthogonal shift.
        kinetic = sp.eye(2)
        mass = sp.diag(2, 5)
        vertex = sp.Matrix([[3, 2], [2, 7]])
        transform = sp.Matrix([[1, 2], [0, 1]])
        residue_i, residue_j = sp.diag(1, 0), sp.diag(0, 1)
        inv = transform.inv()
        transformed_i = inv * residue_i * inv.T
        transformed_j = inv * residue_j * inv.T
        transformed_vertex = transform.T * vertex * transform
        self.assertEqual(sp.trace(residue_i*vertex*residue_j*vertex),
                         sp.trace(transformed_i*transformed_vertex*transformed_j*transformed_vertex))
        for pole, residue in ((2, transformed_i), (5, transformed_j)):
            inverse = (sp.Symbol('s')*(transform.T*kinetic*transform) - transform.T*mass*transform).inv()
            self.assertEqual(inverse.applyfunc(lambda entry: sp.limit((sp.Symbol('s')-pole)*entry, sp.Symbol('s'), pole)), residue)


if __name__ == '__main__':
    unittest.main()
