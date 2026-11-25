"""Connect the projected cycle expansion to the on-shell assembly.

Run: PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_onshell_driver
Inputs remain symbolic. Scalar moments are supplied as a callable whose
combined-mass dependence can be differentiated. No loop integral is evaluated.
"""
import unittest
import sympy as sp

from derive_alp_loop_1pi_callback import minkowski_one_loop
from derive_alp_three_meson_loop_assembly import projected_vertex
from derive_alp_three_meson_onshell_assembly import assemble
from derive_alp_three_meson_onshell_kinematics import daughter_gram


class ScalarMoment(sp.Function):
    """Unevaluated Euclidean scalar moment with its mass derivative."""
    nargs = 2

    def fdiff(self, argindex=1):
        if argindex == 2:
            power, combined_mass = self.args
            return -power * ScalarMoment(power + 1, combined_mass)
        raise sp.ArgumentIndexError(self, argindex)


def formal_scalar_moment(moment_function):
    """Build a normalized, unevaluated simplex integral of scalar moments.

    moment_function(power, combined_mass) supplies a differentiable formal
    scalar moment. Keeping the integral unevaluated preserves its dependence
    on external invariants for pole-residue differentiation.
    """
    def integrate(power, combined_mass, parameters, exponents, coefficient):
        normalization = sp.factorial(sum(exponents) - 1)
        normalization /= sp.prod(sp.factorial(exponent - 1) for exponent in exponents)
        weights = parameters + (1 - sum(parameters),)
        integrand = coefficient * moment_function(power, combined_mass)
        integrand *= normalization * sp.prod(weight**(exponent - 1)
                                              for weight, exponent in zip(weights, exponents))
        if not parameters:
            return integrand
        limits = [(parameter, 0, 1 - sum(parameters[index + 1:]))
                  for index, parameter in enumerate(parameters)]
        return sp.Integral(integrand, *limits)
    return integrate


def dressed_tree_vertex(basis_directions, invariants, epsilon):
    """Build the Minkowski tree vertex for arbitrary field directions."""
    def vertex(directions, momenta, gram):
        legs = []
        for direction in directions:
            meson = sum((direction[index] * basis_directions[index][0]
                         for index in range(len(basis_directions))), sp.zeros(3))
            alp = sum(direction[index] * basis_directions[index][1]
                      for index in range(len(basis_directions)))
            legs.append((meson, alp))

        def dot(first, second):
            return (sp.Matrix(first).T * gram * sp.Matrix(second))[0]

        return projected_vertex(tuple(legs), tuple(momenta), invariants, epsilon, dot)
    return vertex


def crossing_condition(one_loop_amplitude, epsilon, exact_tree_amplitude,
                       maximum_ratio, moment_scale):
    """Construct an estimate after full amplitude and pole differentiation.

    moment_scale receives a formal scalar-moment term and returns its squared
    scale. Distinct integrals remain distinct unless algebraically identical.
    Absolute values are applied after identical terms have been combined.
    """
    cubic = sp.expand(one_loop_amplitude).coeff(epsilon, 3)
    integral_atoms = cubic.atoms(sp.Integral)
    standalone = {item for item in cubic.atoms(ScalarMoment)
                  if cubic.coeff(item) != 0}
    loop_terms = integral_atoms | standalone
    if sp.simplify(cubic.subs({term: 0 for term in loop_terms})) != 0:
        raise ValueError('The loop amplitude contains unrecognized terms.')
    grouped_coefficients = {}
    representatives = {}
    for term in loop_terms:
        coefficient = sp.expand(cubic).coeff(term)
        if isinstance(term, sp.Integral):
            moments = term.function.atoms(ScalarMoment)
            if len(moments) != 1:
                raise ValueError('Each integral must contain one scalar moment.')
            moment = next(iter(moments))
            numerator = sp.cancel(term.function / moment)
            parameter_integral = sp.integrate(numerator, *term.limits)
        else:
            moment = term
            parameter_integral = sp.S.One
        power = moment.args[0]
        if power <= 0:
            continue
        limits = term.limits if isinstance(term, sp.Integral) else ()
        key = (power, sp.expand(moment.args[1]), limits)
        grouped_coefficients[key] = grouped_coefficients.get(key, sp.S.Zero) + coefficient * parameter_integral
        representatives[key] = term
    size = sp.S.Zero
    for key, coefficient in grouped_coefficients.items():
        power = key[0]
        scale_squared = moment_scale(representatives[key])
        size += sp.Abs(sp.simplify(coefficient)) * scale_squared**(2-power) / (16*sp.pi**2)
    size *= sp.Abs(epsilon)**3
    return size, sp.Eq(size, maximum_ratio * sp.Abs(exact_tree_amplitude))


def on_shell_cubic_estimate(masses_squared, pole_vectors, external_poles,
                           pole_directions, pole_masses_squared, basis_directions,
                           invariants, epsilon, dimension, scalar_moment,
                           tree_three, tree_four, tree_five, s12, s23,
                           alp_pole=None, max_order=3):
    """Assemble the truncated Minkowski correction from projected cycles.

    The scalar moment callable receives (power, mass, parameters, exponents,
    coefficient) and must retain the full parameter dependence and normalized
    simplex measure. It must allow symbolic differentiation. The returned
    expression retains terms through max_order in epsilon. Pole vectors,
    masses and tree callbacks can themselves contain epsilon.
    """
    field_dimension = len(basis_directions)
    basis = tuple(sp.eye(field_dimension)[:, index] for index in range(field_dimension))
    external_masses = [masses_squared[index] for index in external_poles]
    on_shell_gram = -daughter_gram(external_masses, s12, s23)

    def loop_series(directions, momenta, gram):
        return sum((epsilon**order * minkowski_one_loop(
            directions, momenta, basis_directions, pole_directions,
            pole_masses_squared, invariants, epsilon, gram, dimension,
            scalar_moment, order, alp_pole) for order in range(max_order + 1)), sp.S.Zero)

    zero = sp.zeros(3, 1)
    one_point = sp.Matrix([loop_series((direction,), (zero,), on_shell_gram)
                           for direction in basis])

    def self_energy(squared_momentum):
        # A two-point function depends on its single external squared momentum.
        gram = sp.diag(-squared_momentum, 0, 0)
        momentum = sp.Matrix([1, 0, 0])
        return sp.Matrix(field_dimension, field_dimension, lambda i, j:
                         loop_series((basis[i], basis[j]), (momentum, -momentum), gram))

    def loop_vertex(directions, momenta):
        return loop_series(directions, momenta, on_shell_gram)

    result = assemble(
        masses_squared, pole_vectors, external_poles, tree_three, tree_four,
        tree_five, self_energy, one_point, loop_vertex, loop_vertex, s12, s23,
        lambda momentum, gram: (momentum.T * gram * momentum)[0])
    return tuple(value.applyfunc(lambda entry: sp.series(entry, epsilon, 0, max_order + 1).removeO())
                 if isinstance(value, sp.MatrixBase) else
                 tuple(sp.series(entry, epsilon, 0, max_order + 1).removeO() for entry in value)
                 if isinstance(value, tuple) and (not value or not isinstance(value[0], sp.MatrixBase)) else
                 value if isinstance(value, tuple) else
                 sp.series(value, epsilon, 0, max_order + 1).removeO() for value in result)


class DriverTest(unittest.TestCase):
    def test_equal_moments_combined_before_magnitudes(self):
        epsilon, x, mass, q, tree, maximum = sp.symbols('epsilon x m2 q T R', positive=True)
        moment = ScalarMoment(2, mass + x*q)
        first = sp.Integral(x*moment, (x, 0, 1))
        second = sp.Integral((1-x)*moment, (x, 0, 1))
        size, _ = crossing_condition(epsilon**3*(first-second), epsilon,
                                     epsilon*tree, maximum, lambda term: mass)
        self.assertEqual(size, 0)

    def test_crossing_after_moment_derivative(self):
        epsilon, coefficient, tree, maximum_ratio, mass, invariant = sp.symbols(
            'epsilon c T R m2 s', positive=True)
        moment = ScalarMoment(1, mass + invariant)
        loop_amplitude = epsilon**3 * coefficient * sp.diff(moment, invariant)
        size, condition = crossing_condition(loop_amplitude, epsilon, epsilon*tree,
                                             maximum_ratio, lambda term: mass)
        self.assertEqual(size, coefficient*epsilon**3/(16*sp.pi**2))
        self.assertEqual(condition, sp.Eq(size, maximum_ratio*epsilon*tree))

    def test_nonzero_projected_tadpole(self):
        epsilon, b0, decay_constant, source = sp.symbols('epsilon B_0 f_pi c')
        masses = sp.symbols('m_1 m_2 m_3 m_a', positive=True)
        basis = tuple(sp.eye(4)[:, index] for index in range(4))
        zero = sp.zeros(3)
        meson = sp.diag(1, 0, 0)
        pseudoscalar = [zero] * 7
        pseudoscalar[3] = source * meson
        invariants = {'S': [zero]*7, 'P': pseudoscalar, 'theta': [0]*7,
                      'f_pi': decay_constant, 'B_0': b0,
                      'anomaly_mass_squared': sp.Symbol('m_0^2')}
        result = on_shell_cubic_estimate(
            masses, basis, (3, 0, 1, 2), {'a': (zero, 1)}, {'a': masses[3]},
            ((meson, 0), (meson, 0), (meson, 0), (zero, 1)),
            invariants, epsilon, sp.Symbol('d'),
            formal_scalar_moment(ScalarMoment),
            lambda *args: 0, lambda *args: 0, lambda *args: 0,
            *sp.symbols('s12 s23'), alp_pole='a')
        self.assertEqual(result[0], -b0*source*epsilon**3*ScalarMoment(1, masses[3]) /
                         (2*decay_constant**4))
        tree, maximum_ratio = sp.symbols('T R', positive=True)
        size, condition = crossing_condition(result[0], epsilon, epsilon*tree,
                                             maximum_ratio, lambda term: masses[3])
        self.assertEqual(size, sp.Abs(b0*source/decay_constant**4) * sp.Abs(epsilon)**3 *
                         masses[3] / (32*sp.pi**2))
        self.assertEqual(condition, sp.Eq(size, maximum_ratio*sp.Abs(epsilon)*tree))

    def test_normalized_simplex_and_derivative(self):
        x, invariant, mass = sp.symbols('x s m')
        moment = sp.Function('J')
        integrate = formal_scalar_moment(lambda power, squared_mass: moment(squared_mass))
        result = integrate(2, mass + x * invariant, (x,), (2, 1), 1 + x)
        expected = sp.Integral(2*x*(1+x)*moment(mass+x*invariant), (x, 0, 1))
        self.assertEqual(result, expected)
        self.assertEqual(sp.diff(result, invariant), sp.diff(expected, invariant))

    def test_zero_loop_assembly(self):
        masses = sp.symbols('m_0:4', positive=True)
        basis = tuple(sp.eye(4)[:, index] for index in range(4))
        zero = sp.zeros(3)
        field_directions = ((zero, 0),) * 3 + ((zero, 1),)
        invariants = {'S': [zero]*7, 'P': [zero]*7, 'theta': [0]*7,
                      'f_pi': sp.Symbol('f_pi', nonzero=True), 'B_0': sp.Symbol('B_0'),
                      'anomaly_mass_squared': sp.Symbol('m_0^2')}
        result = on_shell_cubic_estimate(
            masses, basis, (3, 0, 1, 2), {'a': (zero, 1)}, {'a': masses[3]},
            field_directions, invariants, sp.Symbol('epsilon'), sp.Symbol('d'),
            lambda *args: sp.S.Zero, lambda *args: 0, lambda *args: 0,
            lambda *args: 0, *sp.symbols('s12 s23'), max_order=0)
        self.assertEqual(result[0], 0)
        self.assertEqual(result[1], sp.zeros(4, 1))


if __name__ == '__main__':
    unittest.main()
