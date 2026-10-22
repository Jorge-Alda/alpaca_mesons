"""Derive the odd contribution to a -> gamma P1 P2 decay widths.

Invocation: .venv/bin/python scripts/derive_alp_photon_two_pseudoscalar_width.py
Unit tests: append --test to the invocation.
Inputs: memory/alp/three-meson/alp-photon-pp-amplitude.md.
Natural units (GeV), metric (+---), epsilon^0123=+1. The coefficient
A_{gamma PP} includes coherent exchange terms in the recorded amplitude.
Local anomalous terms and even radiative amplitudes are not added here.
"""
import argparse
import unittest

import sympy as sp


def derive_width():
    """Obtain the polarization contraction and exact three-body phase space."""
    alp_mass, first_mass, second_mass = sp.symbols('m_a m_1 m_2', positive=True)
    pair_mass_squared = sp.Symbol('s', positive=True)
    second_invariant, third_invariant, cosine = sp.symbols('t u z', real=True)
    symmetry_factor = sp.Symbol('S_PP', positive=True)
    amplitude = sp.Function('A_gamma_PP')(pair_mass_squared, second_invariant)
    third_replacement = sp.solve(sp.Eq(
        pair_mass_squared + second_invariant + third_invariant,
        alp_mass**2 + first_mass**2 + second_mass**2), third_invariant)[0]
    photon_first_product = (third_replacement-first_mass**2)/2
    photon_second_product = (second_invariant-second_mass**2)/2
    first_second_product = (pair_mass_squared-first_mass**2-second_mass**2)/2
    gram_matrix = sp.Matrix([
        [0, photon_first_product, photon_second_product],
        [photon_first_product, first_mass**2, first_second_product],
        [photon_second_product, first_second_product, second_mass**2]])
    # Levi-Civita contraction gives -J_mu J^mu=det Gram(k,p1,p2).
    # J.k=0 ensures the real-photon projector may be replaced by -g.
    polarization_contraction = sp.factor(gram_matrix.det())
    pair_kallen = sp.expand(pair_mass_squared**2 + first_mass**4 + second_mass**4
        - 2*pair_mass_squared*first_mass**2 - 2*pair_mass_squared*second_mass**2
        - 2*first_mass**2*second_mass**2)
    photon_energy = (alp_mass**2-pair_mass_squared)/(2*sp.sqrt(pair_mass_squared))
    second_energy = (pair_mass_squared+second_mass**2-first_mass**2)/(2*sp.sqrt(pair_mass_squared))
    pair_momentum = sp.sqrt(pair_kallen)/(2*sp.sqrt(pair_mass_squared))
    second_invariant_of_cosine = sp.simplify(second_mass**2 +
        2*photon_energy*second_energy - 2*photon_energy*pair_momentum*cosine)
    lower_boundary = second_invariant_of_cosine.subs(cosine, 1)
    upper_boundary = second_invariant_of_cosine.subs(cosine, -1)
    jacobian = sp.factor(-sp.diff(second_invariant_of_cosine, cosine))
    angular_contraction = sp.factor(polarization_contraction.subs(
        second_invariant, second_invariant_of_cosine))
    # dPhi3 = ds/(2pi) dPhi2(a;gamma,Q) dPhi2(Q;P1,P2),
    # after integrating overall orientation and the relative azimuth.
    phase_space_cosine = ((alp_mass**2-pair_mass_squared)/(8*sp.pi*alp_mass**2)
        * sp.sqrt(pair_kallen)/(16*sp.pi*pair_mass_squared)/(2*sp.pi))
    phase_space_mandelstam = sp.simplify(phase_space_cosine/jacobian)
    width_multiplier = sp.simplify(phase_space_mandelstam/(2*alp_mass*symmetry_factor))
    squared_amplitude = polarization_contraction*sp.Abs(amplitude)**2
    mandelstam_width_density = width_multiplier*squared_amplitude
    integrated_width = sp.Integral(mandelstam_width_density,
        (second_invariant, lower_boundary, upper_boundary),
        (pair_mass_squared, (first_mass+second_mass)**2, alp_mass**2))
    angular_width_multiplier = sp.factor(width_multiplier*jacobian*angular_contraction)
    # This angular integral is applicable only when the full coefficient is
    # independent of t, such as a retained s-channel contribution alone.
    single_invariant_multiplier = sp.factor(sp.integrate(
        angular_width_multiplier, (cosine, -1, 1)))
    return locals()


class PhotonWidthTests(unittest.TestCase):
    """Independent kinematic and phase-space normalization checks."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_width()

    def test_photon_ward_identity_and_physical_polarizations(self):
        """Compare the invariant result with explicit transverse photons."""
        energy, momentum, cosine = sp.symbols('E p z', real=True)
        first_energy, second_energy = sp.symbols('E_1 E_2', real=True)
        photon = sp.Matrix([energy, 0, 0, energy])
        first = sp.Matrix([first_energy, -momentum*sp.sqrt(1-cosine**2), 0, -momentum*cosine])
        second = sp.Matrix([second_energy, momentum*sp.sqrt(1-cosine**2), 0, momentum*cosine])
        current = sp.Matrix([sum(sp.LeviCivita(mu, nu, rho, sigma)*first[nu]*second[rho]*photon[sigma]
            for nu in range(4) for rho in range(4) for sigma in range(4)) for mu in range(4)])
        self.assertEqual(sp.simplify((current.T*photon)[0]), 0)
        transverse_sum = sp.simplify(current[1]**2+current[2]**2)
        self.assertEqual(sp.simplify(transverse_sum-
            energy**2*momentum**2*(first_energy+second_energy)**2*(1-cosine**2)), 0)
        metric = sp.diag(1, -1, -1, -1)
        vectors = sp.Matrix.hstack(photon, first, second)
        self.assertEqual(sp.simplify((vectors.T*metric*vectors).det()-transverse_sum), 0)

    def test_massless_phase_space_volume(self):
        """Integrate the triangular massless Dalitz domain exactly."""
        result = self.result
        mass, invariant = result['alp_mass'], result['pair_mass_squared']
        volume = sp.integrate(result['phase_space_mandelstam']*(mass**2-invariant),
            (invariant, 0, mass**2))
        self.assertEqual(sp.simplify(volume-mass**2/(256*sp.pi**3)), 0)

    def test_angular_endpoints(self):
        result = self.result
        for cosine in [-1, 1]:
            self.assertEqual(sp.simplify(result['angular_contraction'].subs(result['cosine'], cosine)), 0)


def main():
    """Print the derivation or run its minimal unit tests."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', nargs='?', const=True, default=False, help='Run symbolic unit tests.')
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=['photon-width'], exit=True)
        return
    result = derive_width()
    print('Natural units; metric (+---); real photon; no initial spin average. '
        'S_PP=1 for every listed nonidentical pair. Odd amplitude only; no line-shape approximation.')
    for label, name in [
        ('u', 'third_replacement'), ('Gram matrix', 'gram_matrix'),
        ('Delta_G(s,t), with u eliminated', 'polarization_contraction'),
        ('Polarization-summed squared amplitude', 'squared_amplitude'),
        ('t_-(s)', 'lower_boundary'), ('t_+(s)', 'upper_boundary'),
        ('dPhi3/(ds dt)', 'phase_space_mandelstam'),
        ('dGamma_odd/(ds dt)', 'mandelstam_width_density'),
        ('Gamma_odd as a double Mandelstam integral', 'integrated_width'),
        ('t(s,z), auxiliary angular substitution', 'second_invariant_of_cosine'),
        ('Absolute dt/dz', 'jacobian'),
        ('Polarization contraction in PP rest frame', 'angular_contraction'),
        ('Multiplier of squared amplitude in dGamma/(ds dt)', 'width_multiplier'),
        ('Multiplier of |A(s,t(z))|^2 in dGamma/(ds dz)', 'angular_width_multiplier'),
        ('Multiplier of |A(s)|^2 in dGamma/ds, only for t-independent coefficient', 'single_invariant_multiplier')]:
        print(f'{label}:\n{sp.latex(result[name])}\n')
    print('Integrate (m_1+m_2)^2 <= s <= m_a^2 and t minimum <= t <= t maximum. '
        'Below threshold the width is zero. All exchanges are summed before squaring; '
        'a pole inside phase space requires a specified finite line shape.')


if __name__ == '__main__':
    main()
