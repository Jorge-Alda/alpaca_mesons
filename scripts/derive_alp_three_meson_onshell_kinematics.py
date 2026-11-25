"""On-shell kinematics for the one-loop three-meson amplitude.

Run: PYTHONPATH=scripts .venv/bin/python scripts/derive_alp_three_meson_onshell_kinematics.py
All squared masses and invariants have units GeV² and remain symbolic.
"""
import unittest
import sympy as sp


def daughter_gram(masses_squared, s12, s23):
    """Return the Gram matrix of incoming daughter momenta at fixed s12,s23.

    The order of masses is the ALP pole followed by the three daughter poles.
    Momentum conservation fixes the third daughter-pair invariant.
    """
    alp_mass, first_mass, second_mass, third_mass = masses_squared
    s13 = alp_mass + first_mass + second_mass + third_mass - s12 - s23
    return sp.Matrix([
        [first_mass, (s12-first_mass-second_mass)/2, (s13-first_mass-third_mass)/2],
        [(s12-first_mass-second_mass)/2, second_mass, (s23-second_mass-third_mass)/2],
        [(s13-first_mass-third_mass)/2, (s23-second_mass-third_mass)/2, third_mass],
    ])


def variation(masses_squared, s12, s23, mass_shifts):
    """Return the first variation of the Gram matrix at fixed s12 and s23."""
    gram = daughter_gram(masses_squared, s12, s23)
    return sum((sp.diff(gram, mass) * shift for mass, shift in zip(masses_squared, mass_shifts)), sp.zeros(3))


class KinematicTest(unittest.TestCase):
    def test_invariants_and_momentum_conservation(self):
        masses = sp.symbols('m_a2 m_12 m_22 m_32', real=True)
        s12, s23 = sp.symbols('s_12 s_23', real=True)
        gram = daughter_gram(masses, s12, s23)
        self.assertEqual(sp.simplify(gram[0, 0] + gram[1, 1] + 2*gram[0, 1] - s12), 0)
        self.assertEqual(sp.simplify(gram[1, 1] + gram[2, 2] + 2*gram[1, 2] - s23), 0)
        self.assertEqual(sp.simplify(sum(gram) - masses[0]), 0)

    def test_variation_holds_two_invariants_fixed(self):
        masses = sp.symbols('m_a2 m_12 m_22 m_32', real=True)
        shifts = sp.symbols('delta_m_a2 delta_m_12 delta_m_22 delta_m_32', real=True)
        s12, s23, t = sp.symbols('s_12 s_23 t', real=True)
        varied = daughter_gram(tuple(mass+t*shift for mass, shift in zip(masses, shifts)), s12, s23)
        self.assertEqual(variation(masses, s12, s23, shifts), varied.diff(t).subs(t, 0))
        self.assertEqual(sp.simplify(sum(variation(masses, s12, s23, shifts)) - shifts[0]), 0)


if __name__ == '__main__':
    masses = sp.symbols('m_a2 m_12 m_22 m_32', real=True)
    shifts = sp.symbols('delta_m_a2 delta_m_12 delta_m_22 delta_m_32', real=True)
    s12, s23 = sp.symbols('s_12 s_23', real=True)
    print('Daughter momentum Gram matrix:')
    print(sp.latex(daughter_gram(masses, s12, s23)))
    print('First variation at fixed independent invariants:')
    print(sp.latex(variation(masses, s12, s23, shifts)))
