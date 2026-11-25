"""First-order on-shell three-meson amplitude from 1PI vertex callbacks.

Run: PYTHONPATH=scripts .venv/bin/python scripts/derive_alp_three_meson_onshell_assembly.py
All inputs are symbolic.  The inverse propagator uses D(s)=s I-M² and the
one-loop correction Σ(s), with incoming momentum conventions for vertices.
"""
import unittest
import sympy as sp

from derive_alp_three_meson_onshell_kinematics import daughter_gram, variation


def assemble(masses_squared, pole_vectors, external_poles, tree_three,
             tree_four, tree_five, self_energy, one_point, loop_three, loop_four,
             s12, s23, squared_momentum):
    """Compute the complete linear correction from projected 1PI functions.

    pole_vectors lists a full orthonormal basis, external_poles selects the ALP
    followed by three daughters, and masses_squared lists their full spectrum.
    All 1PI callbacks use Minkowski signs. Vertex callbacks receive directions
    followed by incoming momentum tuples in the basis (p1,p2,p3). self_energy
    receives symbolic s and returns the full matrix Σ(s). tree_three/four/five
    return multilinear projected vertices. The loop and tree functions carry
    formal scalar moments as needed. The derivative in s acts on those moments.
    """
    dimension = len(pole_vectors)
    basis = [sp.eye(dimension)[:, index] for index in range(dimension)]
    external_vectors = [pole_vectors[index] for index in external_poles]
    external_masses = [masses_squared[index] for index in external_poles]
    momenta = [sp.Matrix([-1, -1, -1]), sp.Matrix([1, 0, 0]),
               sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
    zero = sp.zeros(3, 1)
    gram = daughter_gram(external_masses, s12, s23)
    s = sp.Symbol('s', real=True)
    identity = sp.eye(dimension)
    mass_matrix = sum((mass * vector * vector.T for mass, vector in zip(masses_squared, pole_vectors)), sp.zeros(dimension))
    inverse = lambda value: value * identity - mass_matrix
    vacuum = -inverse(0).inv() * one_point
    # The caller supplies the canonically normalized ALP field direction.
    alp_direction = basis[-1]
    def sigma(value):
        loop = self_energy(s).subs(s, value)
        vacuum_term = sp.Matrix(dimension, dimension, lambda i, j:
            tree_three((basis[i], basis[j], vacuum), (zero, zero, zero), gram))
        return loop + vacuum_term
    alp_vector = external_vectors[0]
    alp_mass = external_masses[0]
    overlap = (alp_vector.T * alp_direction)[0]
    bare_mass_shift = (alp_vector.T * sigma(alp_mass) * alp_vector)[0] / overlap**2
    def fixed_sigma(value):
        return sigma(value) - bare_mass_shift * alp_direction * alp_direction.T
    mass_shifts = []
    vector_shifts = []
    for pole in external_poles:
        mass, vector = masses_squared[pole], pole_vectors[pole]
        projected = fixed_sigma(mass)
        mass_shifts.append(-(vector.T * projected * vector)[0])
        norm_shift = -(vector.T * self_energy(s).diff(s).subs(s, mass) * vector)[0] / 2
        shifted = norm_shift * vector
        for index in range(dimension):
            if index != pole:
                projection = sp.simplify((pole_vectors[index].T * projected * vector)[0])
                difference = sp.simplify(mass - masses_squared[index])
                if difference == 0:
                    if projection != 0:
                        raise ValueError('Nonzero loop mixing in a degenerate pole subspace requires diagonalization.')
                    continue
                shifted -= pole_vectors[index] * projection / difference
        vector_shifts.append(shifted)
    directions = tuple(external_vectors)
    correction = loop_four(directions, tuple(momenta)) + tree_five(directions + (vacuum,), tuple(momenta) + (zero,), gram)
    def connected_tree(vectors, gram_matrix):
        amplitude = tree_four(tuple(vectors), tuple(momenta), gram_matrix)
        for daughter in range(1, 4):
            other = [index for index in range(1, 4) if index != daughter]
            q = momenta[other[0]] + momenta[other[1]]
            left = sp.Matrix([[tree_three((vectors[0], vectors[daughter], direction), (momenta[0], momenta[daughter], q), gram_matrix) for direction in basis]])
            right = sp.Matrix([tree_three((vectors[other[0]], vectors[other[1]], direction), (momenta[other[0]], momenta[other[1]], -q), gram_matrix) for direction in basis])
            if left != sp.zeros(1, dimension) and right != sp.zeros(dimension, 1):
                amplitude -= (left * inverse(squared_momentum(q, gram_matrix)).inv() * right)[0]
        return amplitude
    for daughter in range(1, 4):
        other = [index for index in range(1, 4) if index != daughter]
        q = momenta[other[0]] + momenta[other[1]]
        left_momenta = (momenta[0], momenta[daughter], q)
        right_momenta = (momenta[other[0]], momenta[other[1]], -q)
        left_directions = (external_vectors[0], external_vectors[daughter])
        right_directions = (external_vectors[other[0]], external_vectors[other[1]])
        left = sp.Matrix([[tree_three(left_directions + (direction,), left_momenta, gram) for direction in basis]])
        right = sp.Matrix([tree_three(right_directions + (direction,), right_momenta, gram) for direction in basis])
        delta_left = sp.Matrix([[loop_three(left_directions + (direction,), left_momenta) + tree_four(left_directions + (direction, vacuum), left_momenta + (zero,), gram) for direction in basis]])
        delta_right = sp.Matrix([loop_three(right_directions + (direction,), right_momenta) + tree_four(right_directions + (direction, vacuum), right_momenta + (zero,), gram) for direction in basis])
        if left == sp.zeros(1, dimension) and right == sp.zeros(dimension, 1):
            continue
        q_squared = squared_momentum(q, gram)
        propagator = inverse(q_squared).inv()
        correction += (-delta_left*propagator*right-left*propagator*delta_right+left*propagator*fixed_sigma(q_squared)*propagator*right)[0]
    parameter = sp.Symbol('t', real=True)
    varied_vectors = [vector + parameter * shift for vector, shift in zip(external_vectors, vector_shifts)]
    varied_gram = gram + parameter * variation(external_masses, s12, s23, mass_shifts)
    correction += sp.diff(connected_tree(varied_vectors, varied_gram), parameter).subs(parameter, 0)
    return correction, vacuum, bare_mass_shift, tuple(mass_shifts), tuple(vector_shifts)


class OnShellAssemblyTest(unittest.TestCase):
    def test_pole_basis_larger_than_external_and_zero_loop(self):
        masses = sp.symbols('m_0:5', nonzero=True, real=True)
        vectors = [sp.eye(5)[:, index] for index in range(5)]
        external = (4, 0, 1, 2)
        result = assemble(masses, vectors, external,
            lambda *args: 0, lambda *args: 0, lambda *args: 0,
            lambda s: sp.zeros(5), sp.zeros(5, 1), lambda *args: 0,
            lambda *args: 0, *sp.symbols('s12 s23'), lambda q, gram: (q.T*gram*q)[0])
        self.assertEqual(result[0], 0)
        self.assertEqual(result[1], sp.zeros(5, 1))
        self.assertEqual(result[3], (0, 0, 0, 0))

    def test_residues_fixed_mass_and_kinematic_shift(self):
        masses = sp.symbols('m_0:5', nonzero=True, real=True)
        slopes = sp.symbols('b_0:5', real=True)
        coupling, vertex_loop = sp.symbols('g h', real=True)
        vectors = [sp.eye(5)[:, index] for index in range(5)]
        external = (4, 0, 1, 2)
        def tree_four(directions, momenta, gram):
            return (coupling + gram[0, 0]) * directions[0][4] * directions[1][0] * directions[2][1] * directions[3][2]
        def loop_four(directions, momenta):
            return vertex_loop * directions[0][4] * directions[1][0] * directions[2][1] * directions[3][2]
        result = assemble(masses, vectors, external,
            lambda *args: 0, tree_four, lambda *args: 0,
            lambda s: sp.diag(*(s * slope for slope in slopes)), sp.zeros(5, 1),
            lambda *args: 0, loop_four, *sp.symbols('s12 s23'),
            lambda q, gram: (q.T*gram*q)[0])
        expected = vertex_loop - (coupling + masses[0]) * (
            slopes[4] + slopes[0] + slopes[1] + slopes[2]) / 2 - slopes[0] * masses[0]
        self.assertEqual(sp.simplify(result[0] - expected), 0)
        self.assertEqual(result[2], slopes[4] * masses[4])
        self.assertEqual(result[3], (0, -slopes[0]*masses[0], -slopes[1]*masses[1], -slopes[2]*masses[2]))

    def test_cubic_exchange_vertex_corrections(self):
        masses = sp.symbols('m_0:4', nonzero=True, real=True)
        s12, s23, coupling, loop_coupling = sp.symbols('s12 s23 g h', real=True)
        vectors = [sp.eye(4)[:, index] for index in range(4)]
        def product(directions):
            return sp.prod(sum(direction) for direction in directions)
        def tree_three(directions, momenta, gram):
            return coupling * product(directions)
        def loop_three(directions, momenta):
            return loop_coupling * product(directions)
        result = assemble(masses, vectors, (3, 0, 1, 2), tree_three,
            lambda *args: 0, lambda *args: 0, lambda s: sp.zeros(4),
            sp.zeros(4, 1), loop_three, lambda *args: 0, s12, s23,
            lambda q, gram: (q.T*gram*q)[0])
        s13 = sum(masses) - s12 - s23
        expected = -2 * coupling * loop_coupling * sum(
            1/(q_squared-mass) for q_squared in (s12, s23, s13) for mass in masses)
        self.assertEqual(sp.simplify(result[0] - expected), 0)

    def test_degenerate_external_state(self):
        meson_mass, other_mass, alp_mass = sp.symbols('m m_other m_a', nonzero=True, real=True)
        masses = (meson_mass, meson_mass, other_mass, alp_mass)
        vectors = [sp.eye(4)[:, index] for index in range(4)]
        s12, s23 = sp.symbols('s12 s23', real=True)
        mixing = sp.Symbol('mixing', nonzero=True, real=True)
        def run(matrix):
            return assemble(masses, vectors, (3, 0, 1, 2), lambda *args: 0,
                lambda *args: 0, lambda *args: 0, lambda s: matrix,
                sp.zeros(4, 1), lambda *args: 0, lambda *args: 0,
                s12, s23, lambda q, gram: (q.T*gram*q)[0])
        self.assertEqual(run(sp.zeros(4))[0], 0)
        with self.assertRaisesRegex(ValueError, 'degenerate pole subspace'):
            run(sp.Matrix([[0, mixing, 0, 0], [mixing, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]))

if __name__ == '__main__':
    unittest.main()
