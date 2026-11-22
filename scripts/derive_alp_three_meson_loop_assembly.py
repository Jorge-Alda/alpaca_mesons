"""Combinatorics of one-loop projected vertices from the Hessian expansion.

Run: .venv/bin/python scripts/derive_alp_three_meson_loop_assembly.py
Tests: PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_loop_assembly

The external labels are formal leg labels. Internal labels identify exact pole
propagators and carry no numerical masses. A cycle begins with a propagator,
followed by a vertex differentiated with respect to the associated block.
The weights use the Euclidean one-loop Hessian convention of the memory record.
"""
from collections import defaultdict
from itertools import permutations, product
import unittest
from unittest.mock import patch

import sympy as sp
from derive_alp_dressed_vertex_power_counting import contract_potential_vertex


def kinetic_vertex(legs, dot, decay_constant):
    """Return the all-incoming four- or six-meson kinetic vertex.

    Each leg is a pair of a flavour matrix and an arbitrary momentum label.
    dot accepts two momentum labels and returns their Lorentz scalar product.
    The minus sign follows from the two incoming Fourier derivatives.
    """
    number_of_legs = len(legs)
    if number_of_legs not in (4, 6):
        return sp.S.Zero
    total = sp.S.Zero
    for ordering in permutations(legs):
        if number_of_legs == 4:
            first, first_derivative, second, second_derivative = ordering
            first_commutator = first[0] * first_derivative[0] - first_derivative[0] * first[0]
            second_commutator = second[0] * second_derivative[0] - second_derivative[0] * second[0]
            coefficient = sp.Rational(1, 48) / decay_constant**2
        else:
            first, second, first_derivative, third, fourth, second_derivative = ordering
            first_commutator = second[0] * first_derivative[0] - first_derivative[0] * second[0]
            first_commutator = first[0] * first_commutator - first_commutator * first[0]
            second_commutator = fourth[0] * second_derivative[0] - second_derivative[0] * fourth[0]
            second_commutator = third[0] * second_commutator - second_commutator * third[0]
            coefficient = sp.Rational(1, 1440) / decay_constant**4
        total -= coefficient * dot(first_derivative[1], second_derivative[1]) * (
            first_commutator * second_commutator).trace()
    return sp.expand(total)


def projected_vertex(legs, momenta, invariants, epsilon, dot, euclidean=False):
    """Contract dressed potential and kinetic vertices with pole directions.

    Each leg is (meson matrix, ALP component), as in contract_potential_vertex.
    Momenta are incoming and labels are passed unchanged to dot. In Euclidean
    mode, dot uses the Euclidean metric and the potential changes sign.
    """
    potential = contract_potential_vertex(legs, invariants, epsilon)
    kinetic = kinetic_vertex(tuple((leg[0], momentum) for leg, momentum in zip(legs, momenta)),
                             dot, invariants['f_pi'])
    return sp.expand(kinetic - potential if euclidean else kinetic + potential)


def routed_cycle_numerator(cycle, external_directions, external_momenta, pole_directions,
                           loop_momentum, dot, invariants, epsilon, euclidean=False):
    """Construct projected vertices and propagator shifts for a labeled cycle.

    Momenta are tuples of coefficients in a common basis. The sum of all
    external momenta must vanish. Internal momentum entering the first vertex
    is loop_momentum, and all momenta passed to each vertex are incoming.
    Returns the product of vertices and the shifts of the ordered propagators.
    """
    zero = (sp.S.Zero,) * len(loop_momentum)

    def add(first, second):
        return tuple(sp.expand(left + right) for left, right in zip(first, second))

    shifts = []
    vertex_product = sp.S.One
    shift = zero
    for position, (block, incoming_pole) in enumerate(cycle):
        shifts.append(shift)
        next_shift = shift
        for label in block:
            next_shift = add(next_shift, external_momenta[label])
        outgoing_pole = cycle[(position + 1) % len(cycle)][1]
        incoming_momentum = add(loop_momentum, shift)
        outgoing_momentum = tuple(-entry for entry in add(loop_momentum, next_shift))
        legs = tuple(external_directions[label] for label in block) + (
            pole_directions[incoming_pole], pole_directions[outgoing_pole])
        momenta = tuple(external_momenta[label] for label in block) + (
            incoming_momentum, outgoing_momentum)
        vertex_product *= projected_vertex(legs, momenta, invariants, epsilon, dot, euclidean=euclidean)
        shift = next_shift
    if any(sp.simplify(component) != 0 for component in shift):
        raise ValueError('External momenta must sum to zero.')
    return sp.expand(vertex_product), tuple(shifts)


def on_shell_one_loop_amplitude(four_point_correction, exchanges,
                               external_pole_corrections, kinematic_correction,
                               vacuum_correction):
    """Assemble the first-order four-leg amplitude in a fixed field basis.

    Each exchange is (propagator, left vertex, right vertex, left correction,
    right correction, inverse-propagator correction). Row/column matrices or
    commuting scalar inputs can be supplied, consistently within each tuple.
    The final three inputs are separately projected on-shell corrections.
    """
    amplitude = four_point_correction
    for propagator, left, right, left_correction, right_correction, inverse_correction in exchanges:
        terms = (-left_correction * propagator * right,
                 -left * propagator * right_correction,
                 left * propagator * inverse_correction * propagator * right)
        for term in terms:
            amplitude += term[0, 0] if isinstance(term, sp.MatrixBase) else term
    amplitude += sum(external_pole_corrections, sp.S.Zero)
    amplitude += kinematic_correction + vacuum_correction
    if isinstance(amplitude, sp.MatrixBase):
        return amplitude[0, 0]
    return sp.expand(amplitude)


def fixed_alp_pole_correction(inverse_propagator_correction, pole_vector, alp_direction):
    """Retune the bare ALP mass so the selected pole stays fixed at one loop.

    The canonically normalized ALP direction is a column. The pole must have a
    nonzero overlap with it. The correction uses D(s) = s I - M².
    """
    pole_overlap = (pole_vector.T * alp_direction)[0]
    pole_shift = (pole_vector.T * inverse_propagator_correction * pole_vector)[0]
    return inverse_propagator_correction - (
        pole_shift / pole_overlap**2) * alp_direction * alp_direction.T


def vacuum_displacement(stationary_hessian, one_point_correction):
    """Solve the linearized stationarity equation using its stated Hessian."""
    return -stationary_hessian.inv() * one_point_correction


def ordered_partitions(external_legs, number_of_vertices):
    """Return the ordered partitions into nonempty external-leg blocks."""
    for assignment in product(range(number_of_vertices), repeat=len(external_legs)):
        if set(assignment) != set(range(number_of_vertices)):
            continue
        yield tuple(tuple(leg for leg, block in zip(external_legs, assignment)
                          if block == vertex)
                    for vertex in range(number_of_vertices))


def canonical_cycle(blocks, propagators):
    """Identify cyclic rotations while retaining momentum and pole ordering."""
    cycle = tuple(zip(blocks, propagators))
    return min(cycle[index:] + cycle[:index] for index in range(len(cycle)))


def projected_vertices(external_legs, poles, alp_pole=None):
    """Generate weighted cycles for derivatives of one to four external legs.

    Set alp_pole to select cycles containing exactly one ALP-like propagator.
    Each entry is (weight, ((external block, incoming pole), ...)). The vertex
    in a block connects its incoming pole to the next block's incoming pole.
    """
    coefficients = defaultdict(lambda: sp.S.Zero)
    for number_of_vertices in range(1, len(external_legs) + 1):
        coefficient = sp.Rational((-1)**(number_of_vertices + 1), 2 * number_of_vertices)
        for blocks in ordered_partitions(external_legs, number_of_vertices):
            for propagators in product(poles, repeat=number_of_vertices):
                if alp_pole is not None and propagators.count(alp_pole) != 1:
                    continue
                coefficients[canonical_cycle(blocks, propagators)] += coefficient
    return tuple((weight, cycle) for cycle, weight in sorted(coefficients.items()) if weight)


def evaluate_projected_vertices(external_legs, poles, vertex, propagator, trace, alp_pole=None):
    """Evaluate the sum with user-supplied vertex, propagator and trace callables.

    vertex(block, pole_in, pole_out) returns the differentiated Hessian vertex.
    propagator(pole, block) supplies its routed propagator or formal placeholder.
    trace(factors) performs the momentum integration and internal contraction,
    or constructs a symbolic representation of them. No loop is evaluated here.
    """
    result = 0
    for coefficient, cycle in projected_vertices(external_legs, poles, alp_pole):
        factors = []
        for position, (block, incoming_pole) in enumerate(cycle):
            outgoing_pole = cycle[(position + 1) % len(cycle)][1]
            factors.extend((propagator(incoming_pole, block),
                            vertex(block, incoming_pole, outgoing_pole)))
        result += coefficient * trace(tuple(factors))
    return result


class HessianExpansionTest(unittest.TestCase):
    def test_euclidean_potential_sign(self):
        potential, kinetic = sp.symbols('potential kinetic')
        with patch.dict(globals(), {'contract_potential_vertex': lambda *args: potential,
                                    'kinetic_vertex': lambda *args: kinetic}):
            legs = ((sp.eye(3), 0),)
            self.assertEqual(projected_vertex(legs, (0,), {'f_pi': 1}, 1, None),
                             kinetic + potential)
            self.assertEqual(projected_vertex(legs, (0,), {'f_pi': 1}, 1, None, euclidean=True),
                             kinetic - potential)

    def test_on_shell_exchange_and_pole_fixing(self):
        propagator, left, right, left_shift, right_shift, self_energy, direct = sp.symbols(
            'g v_l v_r dv_l dv_r sigma v4')
        calculated = on_shell_one_loop_amplitude(
            direct, [(propagator, left, right, left_shift, right_shift, self_energy)],
            [left_shift], right_shift, 0)
        parameter = sp.Symbol('t')
        exact = (direct * parameter + left_shift * parameter + right_shift * parameter
                 - (left + parameter * left_shift) *
                 (1 / propagator + parameter * self_energy)**-1 *
                 (right + parameter * right_shift))
        self.assertEqual(sp.simplify(calculated - sp.diff(exact, parameter).subs(parameter, 0)), 0)
        matrix_self_energy = sp.Matrix([[sp.Symbol('s11'), sp.Symbol('s12')],
                                        [sp.Symbol('s12'), sp.Symbol('s22')]])
        pole = sp.Matrix([sp.Symbol('v1'), sp.Symbol('v2', nonzero=True)])
        alp_direction = sp.Matrix([0, 1])
        fixed = fixed_alp_pole_correction(matrix_self_energy, pole, alp_direction)
        self.assertEqual(sp.simplify((pole.T * fixed * pole)[0]), 0)
        stationary_hessian = sp.diag(sp.Symbol('h1', nonzero=True),
                                     sp.Symbol('h2', nonzero=True))
        source = sp.Matrix(sp.symbols('j1 j2'))
        self.assertEqual(stationary_hessian * vacuum_displacement(stationary_hessian, source)
                         + source, sp.zeros(2, 1))

    def test_cycle_routing_and_cyclic_rotation(self):
        loop = sp.Symbol('k')
        momentum = sp.Symbol('p')
        directions = {1: 'first', 2: 'second'}
        momenta = {1: (momentum,), 2: (-momentum,)}
        pole_directions = {'a': 'alp', 'n': 'meson'}

        def symbolic_vertex(legs, incoming_momenta, invariants, epsilon, dot, euclidean=False):
            self.assertEqual(sp.simplify(sum(item[0] for item in incoming_momenta)), 0)
            return sp.Function(legs[0])(*(item[0] for item in incoming_momenta))

        cycle = (((1,), 'a'), ((2,), 'n'))
        with patch.dict(globals(), {'projected_vertex': symbolic_vertex}):
            numerator, shifts = routed_cycle_numerator(
                cycle, directions, momenta, pole_directions, (loop,), None, None, None)
            rotated_numerator, rotated_shifts = routed_cycle_numerator(
                cycle[1:] + cycle[:1], directions, momenta, pole_directions,
                (loop + momentum,), None, None, None)
            self.assertEqual(shifts, ((0,), (momentum,)))
            self.assertEqual(rotated_shifts, ((0,), (-momentum,)))
            self.assertEqual(sp.simplify(numerator - rotated_numerator), 0)
            with self.assertRaises(ValueError):
                routed_cycle_numerator(cycle, directions, {1: (momentum,), 2: (momentum,)},
                                       pole_directions, (loop,), None, None, None)

    def test_kinetic_vertex_against_direct_exponential(self):
        first_matrix = sp.Matrix([[0, 1], [1, 0]])
        second_matrix = sp.Matrix([[1, 0], [0, -1]])
        first_field, second_field = sp.symbols('p_1 p_2')
        field = first_field * first_matrix + second_field * second_matrix
        derivative = first_matrix + second_matrix
        zero = sp.zeros(2)
        for number_of_legs in (4, 6):
            left = []
            right = []
            for degree in range(1, number_of_legs):
                derivative_power = sum(
                    (field**position * derivative * field**(degree - position - 1)
                     for position in range(degree)), zero)
                left.append(sp.I**degree * derivative_power / sp.factorial(degree))
                right.append((-sp.I)**degree * derivative_power / sp.factorial(degree))
            lagrangian = sum(
                (left[degree - 1] * right[number_of_legs - degree - 1]
                 for degree in range(1, number_of_legs)), zero).trace() / 4
            half = number_of_legs // 2
            expected = -sp.diff(lagrangian, first_field, half, second_field, half).subs(
                {first_field: 0, second_field: 0})
            legs = tuple((first_matrix, position) for position in range(half)) + tuple(
                (second_matrix, position + half) for position in range(half))
            self.assertEqual(sp.simplify(kinetic_vertex(legs, lambda left, right: 1, 1)
                                         - expected), 0)

    def test_two_point_weights(self):
        cycles = projected_vertices((1, 2), ('a',))
        self.assertEqual(len(cycles), 2)
        self.assertEqual({tuple(len(block) for block, _ in cycle): weight
                          for weight, cycle in cycles}, {(2,): sp.Rational(1, 2),
                                                         (1, 1): -sp.Rational(1, 2)})

    def test_four_point_topologies_and_weights(self):
        cycles = projected_vertices((1, 2, 3, 4), ('a',))
        totals = defaultdict(lambda: sp.S.Zero)
        for weight, cycle in cycles:
            totals[tuple(sorted(len(block) + 2 for block, _ in cycle))] += weight
        self.assertEqual(dict(totals), {(6,): sp.Rational(1, 2), (3, 5): -2,
                                        (4, 4): -sp.Rational(3, 2),
                                        (3, 3, 4): 6, (3, 3, 3, 3): -3})

    def test_distinct_propagators_and_selection(self):
        cycles = projected_vertices((1, 2), ('a', 'n'), 'a')
        self.assertEqual(sorted(weight for weight, _ in cycles),
                         [-sp.Rational(1, 2), -sp.Rational(1, 2), sp.Rational(1, 2)])
        self.assertTrue(all(sum(pole == 'a' for _, pole in cycle) == 1 for _, cycle in cycles))

    def test_callback_receives_adjacent_poles(self):
        result = evaluate_projected_vertices(
            (1, 2), ('a',), lambda block, left, right: sp.Symbol('h' + ''.join(map(str, block))),
            lambda pole, block: sp.Symbol('g'), lambda factors: sp.Mul(*factors))
        g, h1, h2, h12 = sp.symbols('g h1 h2 h12')
        self.assertEqual(result, (g*h12 - g**2*h1*h2)/2)


if __name__ == '__main__':
    for number_of_external_legs in range(1, 5):
        cycles = projected_vertices(tuple(range(1, number_of_external_legs + 1)), ('a',))
        print(f'External legs {number_of_external_legs}:')
        print(sp.latex(sp.Tuple(*(sp.Tuple(weight, sp.Tuple(*(len(block) + 2 for block, _ in cycle)))
                                  for weight, cycle in cycles))))
