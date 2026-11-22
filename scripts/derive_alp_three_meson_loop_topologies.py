"""Symbolic one-loop graph counting and variation of the four-point tree.

Run: .venv/bin/python scripts/derive_alp_three_meson_loop_topologies.py
Test: PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_loop_topologies

All vertices and propagators are symbolic.  No masses or loop functions are
assigned numerical values.  The parameter t marks one-loop corrections.
"""
from itertools import combinations, combinations_with_replacement, permutations
import unittest
import sympy as sp


def one_loop_valences():
    """Enumerate 1PI polygons with four labeled external fields."""
    return [tuple(item + 2 for item in partition)
            for number in range(1, 5)
            for partition in combinations_with_replacement(range(1, 5), number)
            if sum(partition) == 4]


def variation():
    """Differentiate a cubic-exchange contribution and a four-field vertex."""
    t = sp.Symbol('t', real=True)
    vertex_four, vertex_four_one = sp.symbols('V_4 delta_V_4')
    vertex_left, vertex_left_one = sp.symbols('V_3L delta_V_3L', commutative=False)
    vertex_right, vertex_right_one = sp.symbols('V_3R delta_V_3R', commutative=False)
    inverse_propagator, self_energy = sp.symbols('D Pi', commutative=False)
    # First-order expansion of (D+t Pi)^(-1), obtained by differentiating
    # the standard noncommutative inverse rather than assuming commutativity.
    propagator = (inverse_propagator + t * self_energy)**(-1)
    exchange = (vertex_left + t * vertex_left_one) * propagator * (vertex_right + t * vertex_right_one)
    correction = sp.diff(vertex_four + t * vertex_four_one - exchange, t).subs(t, 0).expand()
    return correction


def bubble_wick_factors():
    """Count the pairings of two internal legs at each end of a bubble."""
    field_a, field_p = sp.symbols('a P')
    # Differentiation extracts the two propagating fields from each vertex.
    identical = sp.diff(field_a**2 / sp.factorial(2), field_a, 2)**2 / sp.factorial(2)
    distinct = sp.diff(field_a * field_p, field_a, field_p)**2
    return identical, distinct


def bubble_partitions():
    """Count four-point partitions without exchanging the ends twice."""
    fields = tuple(range(4))
    three_five = tuple(combinations(fields, 1))
    four_four = tuple(pair for pair in combinations(fields, 2) if 0 in pair)
    return three_five, four_four


def pole_and_vacuum_variation():
    """Solve the linearized pole, residue and one-point equations."""
    mass_squared, other_mass_squared = sp.symbols('m_r^2 m_n^2', real=True)
    self_energy_diagonal, self_energy_off_diagonal, self_energy_derivative = sp.symbols(
        r'\Pi_{rr} \Pi_{nr} \Pi_{rr}^{\prime}', real=True
    )
    pole_shift = sp.Symbol(r'\delta m_r^2', real=True)
    overlap_shift = sp.Symbol(r'\delta v_n', real=True)
    normalization_shift = sp.Symbol(r'\delta v_r', real=True)
    vacuum_inverse_propagator = sp.Symbol('D(0)', nonzero=True)
    one_point = sp.Symbol(r'\delta V_1')
    vacuum_shift = sp.Symbol(r'\delta\Phi_0')
    equations = (
        pole_shift + self_energy_diagonal,
        (mass_squared - other_mass_squared) * overlap_shift + self_energy_off_diagonal,
        2 * normalization_shift + self_energy_derivative,
        vacuum_inverse_propagator * vacuum_shift + one_point,
    )
    solution = sp.solve(equations, (pole_shift, overlap_shift, normalization_shift, vacuum_shift),
                        dict=True)[0]
    return solution, equations


def projected_quartic_potential():
    """Derive quartic pole contractions from the source Taylor expansion.

    Leg matrices are arbitrary and the inverse pole operator is kept as a
    matrix.  The common external ALP factors are stripped.  Kinetic vertices
    retain independent momenta and must be supplied separately.
    """
    dimension = 3
    quark_mass = sp.MatrixSymbol('m_q', dimension, dimension)
    scalar_second = sp.MatrixSymbol(r'\widehat{S}_2', dimension, dimension)
    pseudoscalar_first = sp.MatrixSymbol(r'\widehat{P}_1', dimension, dimension)
    inverse_pole_operator = sp.MatrixSymbol(
        r'\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})', dimension, dimension
    )
    meson_first, meson_second, meson_third = [
        sp.MatrixSymbol(name, dimension, dimension) for name in ('p_1', 'p_2', 'p_3')
    ]
    f_pi, b0 = sp.symbols(r'f_\pi B_0', nonzero=True)

    def symmetrized_trace(source, matrices):
        return sum(
            sp.Trace(source * sp.MatMul(*ordering))
            for ordering in permutations(matrices)
        ) / sp.factorial(len(matrices))

    # This is the potential part of the already-derived V_{n_a;n_P}.
    zero_alp = lambda matrices: b0 * symmetrized_trace(quark_mass, matrices) / f_pi**2
    one_alp = lambda matrices: -b0 * symmetrized_trace(pseudoscalar_first, matrices) / f_pi**2
    two_alps = lambda matrices: -b0 * symmetrized_trace(scalar_second, matrices) / f_pi**2
    pole_direction = 2 * inverse_pole_operator
    single = (zero_alp((pole_direction, meson_first, meson_second, meson_third))
              + one_alp((meson_first, meson_second, meson_third)))
    double = (zero_alp((pole_direction, pole_direction, meson_first, meson_second))
              + 2 * one_alp((pole_direction, meson_first, meson_second))
              + two_alps((meson_first, meson_second)))
    return sp.expand(single), sp.expand(double)


def main():
    print('One-loop valences:')
    for graph in one_loop_valences():
        print(sp.latex(sp.Tuple(*graph)))
    print('First-order variation of one cubic-exchange channel:')
    print(sp.latex(variation()))
    print('Bubble factors for identical and distinct internal fields:')
    print(sp.latex(sp.Tuple(*bubble_wick_factors())))
    print('Numbers of 3+5 and 4+4 external-leg partitions:')
    print(sp.latex(sp.Tuple(*(len(channels) for channels in bubble_partitions()))))
    print('Pole and vacuum shifts from the linear equations:')
    shifts, _ = pole_and_vacuum_variation()
    for quantity, expression in shifts.items():
        print(sp.latex(sp.Eq(quantity, expression)))
    print('Single and double ALP projected potential quartic vertices:')
    for expression in projected_quartic_potential():
        print(sp.latex(expression))


class OneLoopTopologyTest(unittest.TestCase):
    def test_valences(self):
        self.assertEqual(one_loop_valences(), [(6,), (3, 5), (4, 4), (3, 3, 4), (3, 3, 3, 3)])

    def test_exchange_and_bubble(self):
        self.assertEqual(bubble_wick_factors(), (sp.Rational(1, 2), 1))
        self.assertEqual(len(sp.Add.make_args(variation())), 4)
        self.assertEqual(tuple(map(len, bubble_partitions())), (4, 3))

    def test_projected_source_content(self):
        single, double = projected_quartic_potential()
        self.assertEqual(len(sp.Add.make_args(single)), 30)
        self.assertEqual(len(sp.Add.make_args(double)), 20)

    def test_pole_and_vacuum_shifts(self):
        solution, equations = pole_and_vacuum_variation()
        self.assertEqual([sp.simplify(equation.subs(solution)) for equation in equations], [0] * 4)


if __name__ == '__main__':
    main()
