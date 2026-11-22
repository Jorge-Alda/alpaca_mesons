"""Symbolic fixed-pole expansion of the mixed propagator and meson poles.

Run: .venv/bin/python scripts/derive_alp_three_meson_loop_power.py
Masses are squared masses in GeV^2. All physical inputs remain symbolic.
"""
import unittest

import sympy as sp


def derive_expansion():
    epsilon = sp.Symbol(r'\epsilon', real=True)
    s, meson_mass_one, meson_mass_two, alp_mass = sp.symbols(
        's m_1^2 m_2^2 m_a^2', real=True
    )
    mixing_one, mixing_two = sp.symbols('C_1 C_2', real=True)
    unmixed = sp.diag(s - meson_mass_one, s - meson_mass_two, s - alp_mass)
    mixing = sp.Matrix([[0, 0, mixing_one], [0, 0, mixing_two],
                        [mixing_one, mixing_two, 0]])
    mass_fixing = sp.diag(0, 0, mixing_one**2 / (alp_mass - meson_mass_one)
                         + mixing_two**2 / (alp_mass - meson_mass_two))
    inverse_unmixed = unmixed.inv()
    propagator = (inverse_unmixed - epsilon * inverse_unmixed * mixing * inverse_unmixed
                  + epsilon**2 * (inverse_unmixed * mixing * inverse_unmixed * mixing * inverse_unmixed
                                  - inverse_unmixed * mass_fixing * inverse_unmixed))
    mass_shift = mixing_one**2 / (meson_mass_one - alp_mass)
    alp_component = -epsilon * mixing_one / (meson_mass_one - alp_mass)
    other_component = epsilon**2 * mixing_one * mixing_two / (
        (meson_mass_one - meson_mass_two) * (meson_mass_one - alp_mass))
    normalization = 1 - epsilon**2 * mixing_one**2 / (2 * (meson_mass_one - alp_mass)**2)
    vector = sp.Matrix([normalization, other_component, alp_component])
    pole_s = meson_mass_one + epsilon**2 * mass_shift
    inverse = unmixed + epsilon * mixing + epsilon**2 * mass_fixing
    residual = (inverse.subs(s, pole_s) * vector).applyfunc(
        lambda entry: sp.series(entry, epsilon, 0, 3).removeO().factor())
    norm = sp.series((vector.T * vector)[0], epsilon, 0, 3).removeO().factor()
    inverse_residual = (inverse * propagator - sp.eye(3)).applyfunc(
        lambda entry: sp.series(entry, epsilon, 0, 3).removeO().factor())
    return (propagator.applyfunc(sp.factor), mass_shift, vector, residual,
            norm, inverse_residual, epsilon, s, alp_mass, meson_mass_one,
            meson_mass_two, mixing_one, mixing_two, mass_fixing)


def main():
    propagator, mass_shift, vector, residual, norm, inverse_residual, *_ = derive_expansion()
    for label, value in [('Propagator through second order', propagator),
                         ('Meson pole shift coefficient', mass_shift),
                         ('Normalized meson pole vector', vector),
                         ('Pole equation residual', residual),
                         ('Vector norm', norm), ('Inverse residual', inverse_residual)]:
        print(f'{label}:\n{sp.latex(value)}\n')


class PropagatorExpansionTest(unittest.TestCase):
    def test_inverse_and_pole(self):
        _, _, _, residual, norm, inverse_residual, *_ = derive_expansion()
        self.assertEqual(residual, sp.zeros(3, 1))
        self.assertEqual(norm, 1)
        self.assertEqual(inverse_residual, sp.zeros(3))

    def test_fixed_alp_mass(self):
        (_, _, _, _, _, _, epsilon, s, alp_mass, meson_mass_one,
         meson_mass_two, mixing_one, mixing_two, mass_fixing) = derive_expansion()
        self.assertEqual(sp.simplify(mass_fixing[2, 2] - mixing_one**2 / (alp_mass - meson_mass_one)
                                     - mixing_two**2 / (alp_mass - meson_mass_two)), 0)


if __name__ == '__main__':
    main()
