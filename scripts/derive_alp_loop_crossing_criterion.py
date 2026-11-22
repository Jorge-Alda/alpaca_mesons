"""Symbolic finite-mixing crossing algebra.

Run: .venv/bin/python scripts/derive_alp_loop_crossing_criterion.py
Tests: PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_crossing_criterion

All physical inputs are symbolic. The mass dimensions are those of the linked
tree and vertex records. The coefficient used below is a local placeholder for
the dimensionless, process-dependent ratio of reduced contractions, including
the loop factor. Its dependence on epsilon is kept until explicitly frozen.
"""
import unittest

import sympy as sp


def derive_crossing():
    epsilon = sp.Symbol(r'\epsilon', positive=True)
    maximum_ratio = sp.Symbol(r'R_{\max}', positive=True)
    overlap_trace = sp.Symbol(
        r'\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})]^2\rangle',
        nonnegative=True,
    )
    coefficient = sp.Function('dimensionless_loop_coefficient')(epsilon**2)
    normalization = 1 + 2 * epsilon**2 * overlap_trace
    external_factor = epsilon / sp.sqrt(normalization)
    loop_factor = epsilon**2 / normalization
    ratio = sp.cancel(external_factor * loop_factor / external_factor) * coefficient
    criterion = sp.Eq(ratio, maximum_ratio)

    frozen_coefficient = sp.Symbol('dimensionless_loop_coefficient', positive=True)
    crossing_squared = sp.solve(
        criterion.subs(coefficient, frozen_coefficient).subs(epsilon**2, sp.Symbol('z')),
        sp.Symbol('z'),
    )[0]
    crossing = sp.sqrt(crossing_squared)
    leading_ratio = sp.series(ratio, epsilon, 0, 4)

    squared_momentum = sp.Symbol('q^2', real=True)
    alp_pole_mass_squared = sp.Symbol('m_a^2', real=True)
    bare_alp_mass_squared = sp.Symbol('m_{a,0}^2', real=True)
    mass_correction = sp.Symbol('dressed_mass_correction', real=True)
    mixing_at_momentum = sp.Symbol(
        r'\langle\mathcal{C}\mathbb{M}_{q^2}^{-1}(\mathcal{C})\rangle', real=True
    )
    mixing_at_pole = sp.Symbol(
        r'\langle\mathcal{C}\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})\rangle', real=True
    )
    inverse_propagator = (
        squared_momentum - bare_alp_mass_squared - epsilon**2 * mass_correction
        + 2 * epsilon**2 * mixing_at_momentum
    )
    pole_condition = inverse_propagator.subs({
        squared_momentum: alp_pole_mass_squared,
        mixing_at_momentum: mixing_at_pole,
    })
    tuned_bare_mass = sp.solve(pole_condition, bare_alp_mass_squared)[0]
    on_shell_inverse_propagator = sp.factor(inverse_propagator.subs(
        bare_alp_mass_squared, tuned_bare_mass
    ))
    return (criterion, crossing_squared, crossing, leading_ratio,
            on_shell_inverse_propagator, ratio, epsilon, maximum_ratio,
            frozen_coefficient, overlap_trace, squared_momentum,
            alp_pole_mass_squared, mixing_at_momentum, mixing_at_pole)


def main():
    (criterion, crossing_squared, crossing, leading_ratio, propagator, ratio,
     epsilon, maximum_ratio, coefficient, overlap, *_) = derive_crossing()
    variable_coefficient = sp.Function('dimensionless_loop_coefficient')(epsilon**2)
    universal_factor = sp.cancel(ratio / variable_coefficient)
    leading_factor = sp.series(universal_factor, epsilon, 0, 4)
    root_residual = sp.simplify(
        criterion.lhs.subs(variable_coefficient, coefficient)
        .subs(epsilon, crossing) - maximum_ratio
    )
    expressions = [
        ('Common ALP-line factor in the ratio', universal_factor),
        ('Small epsilon expansion of that factor', leading_factor),
        ('Residual after substituting the positive frozen-coefficient root', root_residual),
        ('Exact inverse propagator at fixed ALP pole', propagator),
    ]
    for label, expression in expressions:
        print(f'{label}:\n{sp.latex(expression)}\n')
    print('The frozen-coefficient solution requires its denominator to be positive.')
    print('The implicit criterion retains the epsilon dependence of the reduced amplitudes.')


class CrossingCriterionTest(unittest.TestCase):
    def test_frozen_solution_and_power(self):
        (criterion, crossing_squared, _, _, _, ratio, epsilon, maximum_ratio,
         coefficient, overlap, *_) = derive_crossing()
        self.assertEqual(sp.simplify(crossing_squared - maximum_ratio / (
            coefficient - 2 * maximum_ratio * overlap)), 0)
        self.assertEqual(sp.simplify(criterion.lhs.subs(
            sp.Function('dimensionless_loop_coefficient')(epsilon**2), coefficient
        ).subs(epsilon, sp.sqrt(crossing_squared)) - maximum_ratio), 0)
        self.assertEqual(sp.limit((ratio / epsilon**2).subs(
            sp.Function('dimensionless_loop_coefficient')(epsilon**2), coefficient
        ), epsilon, 0), coefficient)

    def test_pole_subtraction(self):
        (_, _, _, _, propagator, _, epsilon, _, _, _, momentum, pole,
         mixing_momentum, mixing_pole) = derive_crossing()
        self.assertEqual(sp.simplify(propagator - (
            momentum - pole + 2 * epsilon**2 * (mixing_momentum - mixing_pole)
        )), 0)
        self.assertEqual(sp.simplify(propagator.subs({momentum: pole,
                                                    mixing_momentum: mixing_pole})), 0)


if __name__ == '__main__':
    main()
