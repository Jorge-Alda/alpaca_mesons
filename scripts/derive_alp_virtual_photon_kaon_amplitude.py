"""Derive retained odd a gamma* -> K Kbar currents from original operators.

Run .venv/bin/python scripts/derive_alp_virtual_photon_kaon_amplitude.py
Tests: append --test. Natural units GeV, (+---), epsilon^0123=+1.
Sources use Alda:2026nrd no-mixing convention. This is the original
seven-VJP/four-VVP plus leading G_V contribution, not a full kaon amplitude.
Independent anomalous contacts, loops and rescattering remain separate.
"""

import argparse
from functools import lru_cache
import unittest

import sympy as sp

from derive_alp_vector_two_pseudoscalar import derive as derive_vertices
from derive_radiative_tff_amplitude import contractions, q1, q2, e1, m, s1, s2


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def commutator(first, second):
    return first * second - second * first


def anticommutator(first, second):
    return first * second + second * first


@lru_cache(maxsize=1)
def derive():
    """Re-evaluate original contractions at nonzero photon virtuality."""
    vertices = derive_vertices()
    electric_charge, even_coupling, pion_decay_constant, alp_scale = sp.symbols(
        'e G_V f_pi f_a', nonzero=True)
    photon_squared, pair_squared, first_squared, second_squared = sp.symbols('y s t u', real=True)
    alp_mass, first_mass, second_mass, internal_mass = sp.symbols('m_a m_1 m_2 M_R', positive=True)
    normalization_mass = vertices['M_V']
    c4 = sp.Symbol('c_4', real=True)
    internal_vector, response, charge, axial, pseudoscalar, chi = (
        sp.MatrixSymbol(name, 3, 3) for name in ('T_R', 'W_gamma', 'Q', 'c_q^A', 'c_q^P', 'chi'))
    first_outgoing, second_outgoing = (
        sp.MatrixSymbol(name, 3, 3) for name in ('lambda_1^dagger', 'lambda_2^dagger'))
    sources = ((2 * axial / alp_scale, -2 * anticommutator(pseudoscalar, chi) / alp_scale,
                alp_mass, pair_squared, None),
               (-second_outgoing / pion_decay_constant,
                anticommutator(second_outgoing, chi) / pion_decay_constant,
                second_mass, first_squared, second_outgoing),
               (-first_outgoing / pion_decay_constant,
                anticommutator(first_outgoing, chi) / pion_decay_constant,
                first_mass, second_squared, first_outgoing))
    currents = (commutator(first_outgoing, second_outgoing),
                commutator(first_outgoing, axial), commutator(second_outgoing, axial))
    prefactors = (-even_coupling / (sp.sqrt(2) * pion_decay_constant**2),
                  sp.sqrt(2) * even_coupling / (pion_decay_constant * alp_scale),
                  -sp.sqrt(2) * even_coupling / (pion_decay_constant * alp_scale))
    pieces, exchanges = [], []
    for source, current, prefactor in zip(sources, currents, prefactors):
        derivative_source, mass_source, source_mass, vector_squared, outgoing = source
        converted = sp.sqrt(2) * electric_charge * vertices['strong'].subs({
            vertices['D']: sp.Trace(anticommutator(internal_vector, response) * derivative_source) / 2,
            vertices['P']: sp.Trace(anticommutator(internal_vector, response) * mass_source) / 2,
            m: source_mass, s1: photon_squared, s2: vector_squared})
        direct = -2 * electric_charge * vertices['ordered'].subs({
            vertices['D']: sp.Trace(anticommutator(internal_vector, charge) * derivative_source) / 2,
            vertices['P']: sp.Trace(anticommutator(internal_vector, charge) * mass_source) / 2,
            m: source_mass, s1: vector_squared, s2: photon_squared})
        charged_term = sp.S.Zero
        if outgoing is not None:
            lagrangian_coefficient = sp.I * (-sp.I * electric_charge / pion_decay_constant) * 2 * c4 / normalization_mass
            charged_term = sp.simplify(-vertices['first'][2] * lagrangian_coefficient) * sp.Trace(
                internal_vector * commutator(commutator(outgoing, charge), chi))
        pieces.append((converted, direct, charged_term))
        exchanges.append(prefactor * sp.Trace(internal_vector * current)
                         * sum(pieces[-1]) / (vector_squared - internal_mass**2))

    # Explicit matrices are the supplied physical flavor generators, not a
    # component replacement of an unspecified abstract matrix calculation.
    charge_entries = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    axial_entries = sp.diag(*sp.symbols('c_u^A c_d^A c_s^A', real=True))
    pseudoscalar_entries = sp.diag(*sp.symbols('c_u^P c_d^P c_s^P', real=True))
    chi_entries = sp.diag(*sp.symbols('chi_u chi_d chi_s', real=True))
    theta = sp.Symbol('theta_V', real=True)
    rho = sp.diag(1, -1, 0) / sp.sqrt(2)
    ideal_omega, ideal_phi = sp.diag(1, 1, 0) / sp.sqrt(2), sp.diag(0, 0, 1)
    generators = [rho, sp.cos(theta) * ideal_omega - sp.sin(theta) * ideal_phi,
                  sp.sin(theta) * ideal_omega + sp.cos(theta) * ideal_phi]
    neutral_masses = sp.symbols('M_rho M_omega M_phi', positive=True)
    charged_rho_mass, charged_kstar_mass, neutral_kstar_mass = sp.symbols(
        'M_rho_ch M_Kstar_ch M_Kstar_0', positive=True)
    vector_masses = list(neutral_masses)
    for row, column, mass_value in ((0, 1, charged_rho_mass), (0, 2, charged_kstar_mass), (1, 2, neutral_kstar_mass)):
        symmetric, antisymmetric = sp.zeros(3), sp.zeros(3)
        symmetric[row, column] = symmetric[column, row] = 1 / sp.sqrt(2)
        antisymmetric[row, column], antisymmetric[column, row] = -sp.I / sp.sqrt(2), sp.I / sp.sqrt(2)
        generators.extend((symmetric, antisymmetric))
        vector_masses.extend((mass_value, mass_value))
    residue_values = sp.symbols('R_rho R_omega R_phi', real=True)
    response_entries = sum((residue * generator / (mass_value**2 - photon_squared)
                            for residue, generator, mass_value in zip(residue_values, generators[:3], neutral_masses)), sp.zeros(3))
    generic_response_entries = sp.diag(*sp.symbols('w_u w_d w_s'))
    substitutions = {charge: charge_entries, axial: axial_entries, pseudoscalar: pseudoscalar_entries,
                     chi: chi_entries, response: generic_response_entries}
    kaon_masses = sp.symbols('M_Kplus M_K0', positive=True)
    channel_terms, allowed_vectors, c4_flavor, channel_substitutions = [], [], [], []
    for light_index, kaon_mass in zip((0, 1), kaon_masses):
        outgoing_first, outgoing_second = sp.zeros(3), sp.zeros(3)
        outgoing_first[2, light_index] = sp.sqrt(2)
        outgoing_second[light_index, 2] = sp.sqrt(2)
        channel_substitution = dict(substitutions)
        channel_substitution.update({first_outgoing: outgoing_first, second_outgoing: outgoing_second,
                                     first_mass: kaon_mass, second_mass: kaon_mass})
        channel_substitutions.append(channel_substitution)
        terms, nonzero = [], []
        for exchange in exchanges:
            vector_terms = [sp.factor(exchange.subs(channel_substitution).subs(
                {internal_vector: generator, internal_mass: mass_value}).doit())
                for generator, mass_value in zip(generators, vector_masses)]
            terms.append(sp.factor(sum(vector_terms)))
            nonzero.append(tuple(index for index, term in enumerate(vector_terms) if term != 0))
        channel_terms.append(terms)
        allowed_vectors.append(nonzero)
        c4_flavor.append([sp.simplify(piece[2].subs(channel_substitution).subs(internal_vector, generator).doit())
                          for piece in pieces for generator in generators])

    # The tensor is epsilon(eps,p1,p2,k); replacing eps by k gives zero
    # at arbitrary k^2, independently of flavor and scalar denominators.
    first_momentum, second_momentum, photon_momentum = (
        sp.Matrix(sp.symbols(name + '_0:4')) for name in ('p_1', 'p_2', 'k'))
    ward = sp.Matrix.hstack(photon_momentum, first_momentum, second_momentum, photon_momentum).det()
    phase_space_squared, kaon_mass = sp.symbols('s_K M_K', positive=True)
    kaon_velocity = sp.sqrt(1 - 4 * kaon_mass**2 / phase_space_squared)
    polar_cosine = sp.Symbol('z', real=True)
    angular_weight = sp.integrate(1 - polar_cosine**2, (polar_cosine, -1, 1))
    # Same charged distinct-particle two-body cut derivation as pion record:
    # phi integral of cos^2(phi) and 2 Im T give respectively pi and 1/2.
    azimuth = sp.Symbol('phi', real=True)
    azimuth_weight = sp.integrate(sp.cos(azimuth)**2, (azimuth, 0, 2 * sp.pi))
    cut_coefficient = sp.simplify(sp.Rational(1, 2) * phase_space_squared * kaon_velocity**3
                                 * azimuth_weight * angular_weight / (64 * sp.pi**2))
    return locals()


class VirtualKaonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_real_photon_operator_limit(self):
        from derive_alp_photon_two_pseudoscalar import derive as derive_real
        result, real = self.result, derive_real()
        replacements = {real['internal_vector']: result['internal_vector'],
                        real['photon_response']: result['response'], real['charge']: result['charge'],
                        real['axial']: result['axial'], real['pseudoscalar']: result['pseudoscalar'],
                        real['mass']: result['chi'], real['first_generator']: sp.Adjoint(result['first_outgoing']),
                        real['second_generator']: sp.Adjoint(result['second_outgoing'])}
        for virtual, real_piece in zip(result['pieces'], real['radiative_pieces']):
            for derived_piece, reference_piece in zip(virtual, real_piece[1:]):
                self.assertEqual(sp.expand(derived_piece.subs(result['photon_squared'], 0).doit()
                                           - reference_piece.subs(replacements).doit()), 0)

    def test_ward_and_original_operator_contractions(self):
        self.assertEqual(self.result['ward'], 0)
        self.assertEqual(contractions(q1, q2, e1, q2, False), [0] * 7)

    def test_flavor_channels_and_c4(self):
        result = self.result
        self.assertEqual(result['allowed_vectors'][0], [(0, 1, 2), (5, 6), (5, 6)])
        self.assertEqual(result['allowed_vectors'][1], [(0, 1, 2), (7, 8), (7, 8)])
        self.assertTrue(any(term != 0 for term in result['c4_flavor'][0]))
        self.assertTrue(all(term == 0 for term in result['c4_flavor'][1]))

    def test_sources_and_cut_normalization(self):
        result = self.result
        coefficients = list(result['axial_entries'].diagonal()) + list(result['pseudoscalar_entries'].diagonal())
        for terms in result['channel_terms']:
            scalar = sum(terms)
            self.assertEqual(sp.expand(scalar.subs(dict.fromkeys(coefficients, 0))), 0)
            self.assertEqual(sp.hessian(scalar, coefficients).applyfunc(sp.simplify), sp.zeros(6))
        self.assertEqual(sp.simplify(result['cut_coefficient']
            - result['phase_space_squared'] * result['kaon_velocity']**3 / (96 * sp.pi)), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    if parser.parse_args().test:
        unittest.main(argv=['derive_alp_virtual_photon_kaon_amplitude.py'])
        return
    result = derive()
    print('GeV; (+---); epsilon^0123=+1. k^2=y, s+t+u=m_a^2+2 M_K^2+y; '
          'full neutral W_gamma(y), stable tree poles with analytic i0. '
          'No independent local anomalous term, rescattering, or high-energy completion.')
    show('Original ordered tensor contractions', result['vertices']['first'])
    show('Original VJP numerator', result['vertices']['ordered'])
    show('Original VVP numerator', result['vertices']['strong'])
    show('Full signed neutral photon response', result['response_entries'])
    for routing, pieces, exchange in zip(('s', 't', 'u'), result['pieces'], result['exchanges']):
        for name, piece in zip(('converted VVP', 'direct VJP', 'charged c4'), pieces):
            show(routing + ' ' + name, piece)
        show(routing + ' exchange per Hermitian nonet generator', exchange)
    for label, terms, allowed in zip(('K+ K-', 'K0 Kbar0'), result['channel_terms'], result['allowed_vectors']):
        print(label + ' nonzero vector indices (rho,omega,phi,rho1,rho2,K*1,K*2,K*01,K*02): ' + str(allowed))
        for routing, term in zip(('s', 't', 'u'), terms):
            show(label + ' ' + routing + ' scalar with W_gamma diagonal entries unevaluated', term)
    show('Photon Ward contraction', result['ward'])
    show('Reduced P-wave angular weight integral', result['angular_weight'])
    show('Distinct-particle single-cut coefficient per kaon channel', result['cut_coefficient'])


if __name__ == '__main__':
    main()
