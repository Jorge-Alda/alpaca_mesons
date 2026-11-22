"""Symbolic limits of an isolated ALP-meson pole pair.

Run: .venv/bin/python scripts/derive_alp_loop_crossing_limits.py
Tests: PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_crossing_limits
Masses and mixing entries are symbolic and have units of GeV squared.
"""
import unittest

import sympy as sp
from solve_alp_meson_spectrum import inverse_propagator


def derive_limits():
    epsilon = sp.Symbol(r'\epsilon', positive=True)
    meson_mass_squared = sp.Symbol(r'm_{n,0}^2', real=True)
    mass_gap = sp.Symbol(r'\left(M^2_{\hat a\hat a,\mathrm{can}}-m_{n,0}^2\right)', positive=True)
    mixing_entry = sp.Symbol(r'\langle\lambda_n\mathcal{C}\rangle', positive=True)
    momentum_squared = sp.Symbol('s', real=True)
    mass_matrix = sp.Matrix([
        [meson_mass_squared, -epsilon * mixing_entry],
        [-epsilon * mixing_entry, meson_mass_squared + mass_gap],
    ])
    inverse = inverse_propagator(sp.eye(2), mass_matrix, momentum_squared)
    characteristic = sp.factor(inverse.det())
    poles = sorted(sp.solve(characteristic, momentum_squared), key=sp.default_sort_key)
    lower = next(pole for pole in poles if sp.simplify(sp.limit(pole, epsilon, 0) - meson_mass_squared) == 0)
    upper = next(pole for pole in poles if pole != lower)
    residue_lower = (inverse.adjugate().subs(momentum_squared, lower) / (lower-upper)).applyfunc(sp.simplify)
    residue_upper = (inverse.adjugate().subs(momentum_squared, upper) / (upper-lower)).applyfunc(sp.simplify)
    return (epsilon, meson_mass_squared, mass_gap, mixing_entry, momentum_squared,
            inverse, characteristic, lower, upper, residue_lower, residue_upper)


def spectator_equation():
    """Derive the exact Schur equation with a second unmixed meson."""
    epsilon = sp.Symbol(r'\epsilon', positive=True)
    meson_mass_squared = sp.Symbol(r'm_{n,0}^2', real=True)
    other_meson_mass_squared = sp.Symbol(r'm_{n^{\prime},0}^2', real=True)
    alp_mass_entry = sp.Symbol(r'M^2_{\hat a\hat a,\mathrm{can}}', real=True)
    mixing = sp.Symbol(r'\langle\lambda_n\mathcal{C}\rangle', real=True)
    other_mixing = sp.Symbol(r'\langle\lambda_{n^{\prime}}\mathcal{C}\rangle', real=True)
    momentum_squared = sp.Symbol('s', real=True)
    mass_matrix = sp.Matrix([
        [meson_mass_squared, 0, -epsilon * mixing],
        [0, other_meson_mass_squared, -epsilon * other_mixing],
        [-epsilon * mixing, -epsilon * other_mixing, alp_mass_entry],
    ])
    inverse = inverse_propagator(sp.eye(3), mass_matrix, momentum_squared)
    determinant = sp.factor(inverse.det())
    # Integrating out the second meson retains its full momentum dependence.
    schur = sp.factor((momentum_squared - meson_mass_squared) * (
        momentum_squared - alp_mass_entry
        - epsilon**2 * other_mixing**2 / (momentum_squared - other_meson_mass_squared)
    ) - epsilon**2 * mixing**2)
    return determinant, schur, momentum_squared, other_meson_mass_squared


def main():
    (epsilon, _, mass_gap, mixing_entry, _, inverse, characteristic, lower,
     upper, residue_lower, residue_upper) = derive_limits()
    results = [
        ('Inverse propagator', inverse), ('Pole equation', characteristic),
        ('Lower pole', lower), ('Upper pole', upper),
        ('Lower pole at fixed nonzero gap', sp.series(lower, epsilon, 0, 3)),
        ('Upper pole at fixed nonzero gap', sp.series(upper, epsilon, 0, 3)),
        ('Pole splitting at degeneracy', sp.simplify((upper-lower).subs(mass_gap, 0))),
        ('Lower residue at degeneracy', residue_lower.subs(mass_gap, 0)),
        ('Upper residue at degeneracy', residue_upper.subs(mass_gap, 0)),
        ('ALP weight at the lower pole in decoupling', sp.limit(residue_lower[1, 1], epsilon, 0)),
        ('ALP weight at the upper pole in decoupling', sp.limit(residue_upper[1, 1], epsilon, 0)),
        ('Mixing-to-gap ratio', epsilon * mixing_entry / mass_gap),
    ]
    for label, expression in results:
        print(f'{label}:\n{sp.latex(expression)}\n')
    determinant, schur, momentum_squared, other_meson_mass_squared = spectator_equation()
    print('Exact Schur equation with the other meson:')
    print(sp.latex(schur))
    print('Characteristic determinant after restoring the other meson denominator:')
    print(sp.latex(sp.factor((momentum_squared-other_meson_mass_squared)*schur-determinant)))


class PoleLimitsTest(unittest.TestCase):
    def test_poles_and_residues(self):
        (epsilon, meson_mass_squared, mass_gap, mixing_entry, momentum_squared, inverse, _,
         lower, upper, residue_lower, residue_upper) = derive_limits()
        for pole, residue in ((lower, residue_lower), (upper, residue_upper)):
            self.assertEqual((inverse.subs(momentum_squared, pole) * residue).applyfunc(sp.simplify), sp.zeros(2))
            self.assertEqual((residue * residue - residue).applyfunc(sp.simplify), sp.zeros(2))
        self.assertEqual((residue_lower + residue_upper).applyfunc(sp.simplify), sp.eye(2))
        self.assertEqual(sp.simplify((upper-lower).subs(mass_gap, 0)), 2*epsilon*mixing_entry)
        self.assertEqual(residue_lower[1, 1].subs(mass_gap, 0), sp.Rational(1, 2))
        self.assertEqual(sp.limit(lower, epsilon, 0), meson_mass_squared)
        self.assertEqual(sp.limit(upper, epsilon, 0), meson_mass_squared + mass_gap)
        self.assertEqual(sp.limit(residue_lower[1, 1], epsilon, 0), 0)

    def test_fixed_gap_shift(self):
        epsilon, _, mass_gap, mixing_entry, _, _, _, lower, upper, _, _ = derive_limits()
        self.assertEqual(sp.simplify(sp.diff(lower, epsilon, 2).subs(epsilon, 0)), -2*mixing_entry**2/mass_gap)
        self.assertEqual(sp.simplify(sp.diff(upper, epsilon, 2).subs(epsilon, 0)), 2*mixing_entry**2/mass_gap)

    def test_spectator_and_joint_limit(self):
        determinant, schur, momentum_squared, other_meson_mass_squared = spectator_equation()
        self.assertEqual(sp.factor(determinant-(momentum_squared-other_meson_mass_squared)*schur), 0)
        epsilon, _, mass_gap, mixing_entry, _, _, _, lower, upper, residue_lower, _ = derive_limits()
        # This path approaches degeneracy while the mixing-to-gap ratio is fixed.
        self.assertEqual(sp.simplify((upper-lower).subs(mass_gap, epsilon*mixing_entry)),
                         sp.sqrt(5)*epsilon*mixing_entry)
        self.assertEqual(sp.simplify(residue_lower[1, 1].subs(mass_gap, epsilon*mixing_entry)),
                         (1-1/sp.sqrt(5))/2)
        self.assertEqual(sp.limit(epsilon**3/epsilon, epsilon, 0), 0)
        self.assertEqual(sp.limit(epsilon**2/(1+2*epsilon**2*mixing_entry**2/mass_gap**2),
                                  epsilon, 0), 0)


if __name__ == '__main__':
    main()
