"""Derive preliminary-fit null projections for a gamma* -> pi+ pi-.

Run: .venv/bin/python scripts/derive_alp_virtual_photon_pion_nulls.py
Test: .venv/bin/python scripts/derive_alp_virtual_photon_pion_nulls.py --test

Natural units (GeV), metric (+---), canonical epsilon^0123=+1.
No-mixing sources; original seven VJP/four VVP and leading G_V only.
Use crossed decay invariants s=(p+ + p-)^2, t=(k+p-)^2,
u=(k+p+)^2, k=-q_gamma, k^2=y. The common overall crossing
phase does not affect whether a null projection vanishes. Strong isospin
sets charged and neutral rho masses equal. Local anomalous terms, loops,
rescattering, and further odd operators are outside this derivation.
"""

import argparse
import unittest

import sympy as sp

from derive_alp_vector_two_pseudoscalar import derive as derive_vertices
from derive_radiative_tff_amplitude import m, s1, s2


def show(label, expression):
    """Render a derived expression with the project's LaTeX convention."""
    print(label + ':\n' + sp.latex(expression) + '\n')


def derive():
    """Project operator contractions before doing the charged-pion traces."""
    vertices = derive_vertices()
    coefficients = sp.symbols('c_1:8 d_1:5', real=True)
    electric_charge, even_coupling, decay_constant, alp_scale = sp.symbols(
        'e G_V f_pi f_a', nonzero=True)
    photon_squared, pair_squared, first_squared, second_squared = sp.symbols(
        'y s t u', real=True)
    alp_mass, pion_mass, rho_mass = sp.symbols('m_a m_pi M_rho', positive=True)
    axial_up, axial_down, axial_strange = sp.symbols('c_u^A c_d^A c_s^A', real=True)
    axial = sp.diag(axial_up, axial_down, axial_strange)
    charge = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    first_outgoing, second_outgoing = sp.zeros(3), sp.zeros(3)
    first_outgoing[1, 0] = sp.sqrt(2)
    second_outgoing[0, 1] = sp.sqrt(2)

    ordered_row = sp.Matrix([[sp.diff(vertices['ordered'], coefficient)
                              for coefficient in coefficients]])
    strong_row = sp.Matrix([[sp.diff(vertices['strong'], coefficient)
                             for coefficient in coefficients]])
    direct_projection = (ordered_row * vertices['nulls']).applyfunc(sp.factor)
    strong_projection = (strong_row * vertices['nulls']).applyfunc(sp.factor)

    def commutator(first, second):
        return first * second - second * first

    def derivative_trace(vector, source):
        return sp.trace((vector * charge + charge * vector) * source) / 2

    # A complete orthonormal Hermitian U(3) vector basis. The nonzero
    # flavour traces select only the rho triplet for these pion currents.
    vector_generators = [sp.diag(1, -1, 0) / sp.sqrt(2),
                         sp.diag(1, 1, -2) / sp.sqrt(6),
                         sp.eye(3) / sp.sqrt(3)]
    for first in range(3):
        for second in range(first + 1, 3):
            symmetric, antisymmetric = sp.zeros(3), sp.zeros(3)
            symmetric[first, second] = symmetric[second, first] = 1 / sp.sqrt(2)
            antisymmetric[first, second] = -sp.I / sp.sqrt(2)
            antisymmetric[second, first] = sp.I / sp.sqrt(2)
            vector_generators.extend((symmetric, antisymmetric))

    sources = (2 * axial / alp_scale, -second_outgoing / decay_constant,
               -first_outgoing / decay_constant)
    currents = (commutator(first_outgoing, second_outgoing),
                commutator(first_outgoing, axial),
                commutator(second_outgoing, axial))
    denominators = (pair_squared - rho_mass**2, first_squared - rho_mass**2,
                    second_squared - rho_mass**2)
    prefactors = (-even_coupling / (sp.sqrt(2) * decay_constant**2),
                 sp.sqrt(2) * even_coupling / (decay_constant * alp_scale),
                 -sp.sqrt(2) * even_coupling / (decay_constant * alp_scale))
    exchange_projections = []
    for source, current, denominator, prefactor in zip(
            sources, currents, denominators, prefactors):
        flavour_sum = sum(derivative_trace(vector, source) * sp.trace(vector * current)
                          for vector in vector_generators)
        projection = direct_projection.subs({vertices['D']: flavour_sum,
                                              s2: photon_squared})
        exchange_projections.append((projection * (-2 * electric_charge)
                                     * prefactor / denominator).applyfunc(sp.factor))
    total_projection = sum(exchange_projections, sp.zeros(1, 4)).applyfunc(sp.factor)

    # A scalar coefficient multiplying epsilon(eps,p+,p-,k) expands in
    # P'_ell(z); f1=(3/4) integral (1-z^2) A dz. A constant projects to itself.
    cosine = sp.Symbol('z', real=True)
    projector_normalization = sp.integrate(sp.Rational(3, 4) * (1 - cosine**2),
                                           (cosine, -1, 1))
    # The equal-u,d axial source is a physical allowed profile for which
    # both crossed contributions vanish, leaving an exact nonzero P wave.
    isoscalar_projection = total_projection.subs(axial_down, axial_up)
    p_wave_isoscalar_projection = (projector_normalization * isoscalar_projection
                                   ).applyfunc(sp.factor)
    return locals()


class VirtualPhotonNullTests(unittest.TestCase):
    """Independent limiting and tensor-identity checks."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_three_structural_nulls(self):
        self.assertEqual(self.result['total_projection'][:, 1:], sp.zeros(1, 3))
        self.assertEqual(self.result['strong_projection'], sp.zeros(1, 4))

    def test_real_photon_limit(self):
        result = self.result
        self.assertEqual(result['total_projection'].subs(result['photon_squared'], 0),
                         sp.zeros(1, 4))

    def test_surviving_p_wave(self):
        result = self.result
        self.assertEqual(result['projector_normalization'], 1)
        self.assertNotEqual(result['p_wave_isoscalar_projection'][0], 0)
        self.assertEqual(result['exchange_projections'][1].subs(
            result['axial_down'], result['axial_up']), sp.zeros(1, 4))
        self.assertEqual(result['exchange_projections'][2].subs(
            result['axial_down'], result['axial_up']), sp.zeros(1, 4))

    def test_basis_normalization(self):
        generators = self.result['vector_generators']
        self.assertEqual(sp.Matrix([[sp.trace(first * second)
            for second in generators] for first in generators]), sp.eye(9))


def main():
    """Print the derivation or execute its focused unit tests."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    if parser.parse_args().test:
        unittest.main(argv=[__file__])
        return
    result = derive()
    print(__doc__)
    show('Exact fitted null columns in c1,...,c7,d1,...,d4 order',
         result['vertices']['nulls'])
    show('Ordered VJP tensor row times null columns', result['direct_projection'])
    show('VVP tensor row times null columns', result['strong_projection'])
    for label, projection in zip(('s', 't', 'u'), result['exchange_projections']):
        show(label + ' channel pion-amplitude null projection', projection)
    show('Complete charged-pion scalar null projection', result['total_projection'])
    show('P-wave null projection for c_u^A=c_d^A',
         result['p_wave_isoscalar_projection'])


if __name__ == '__main__':
    main()
