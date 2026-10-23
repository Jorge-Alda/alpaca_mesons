"""Derive the retained odd-resonance a gamma* -> pi+ pi- amplitude.

Run: .venv/bin/python scripts/derive_alp_virtual_photon_pion_amplitude.py
Tests: .venv/bin/python scripts/derive_alp_virtual_photon_pion_amplitude.py --test

Natural units (GeV), metric (+---), epsilon^0123=+1. Sources already use
the Alda no-mixing convention. The scalar multiplies the existing decay
tensor epsilon(epsilon_gamma,p_plus,p_minus,k), with k=-q_gamma and k^2=y.
The retained terms are seven VJP, four VVP, and leading G_V, to order e.
Diagonal axial and pseudoscalar sources and strong isospin are explicit
assumptions. Independent anomalous contacts, loops and rescattering are
excluded. Stable pole denominators mean their analytic i0 continuation.
The P wave is (3/4) integral dz (1-z^2) times this scalar; conversion to a
current convention with an explicit i or a factor e must be done separately.
"""

import argparse
import unittest

import sympy as sp

from derive_alp_vector_two_pseudoscalar import derive as derive_vertices
from derive_radiative_tff_amplitude import m, s1, s2


def show(label, expression):
    """Print a deterministic labelled equation with SymPy's renderer."""
    print(label + ':\n' + sp.latex(expression) + '\n')


def derive():
    """Contract original operators, sum flavours, and integrate the P wave."""
    vertices = derive_vertices()
    electric_charge, even_coupling, pion_decay_constant, alp_scale = sp.symbols(
        'e G_V f_pi f_a', nonzero=True)
    photon_squared, pair_squared, first_squared, second_squared = sp.symbols(
        'y s t u', real=True)
    alp_mass, pion_mass, rho_mass = sp.symbols('m_a m_pi M_rho', positive=True)
    axial_up, axial_down, axial_strange = sp.symbols('c_u^A c_d^A c_s^A', real=True)
    pseudoscalar_up, pseudoscalar_down, pseudoscalar_strange = sp.symbols(
        'c_u^P c_d^P c_s^P', real=True)
    photon_up, photon_down, photon_strange = sp.symbols('w_u w_d w_s', real=True)
    strange_mass_source = sp.Symbol('chi_s', positive=True)
    axial = sp.diag(axial_up, axial_down, axial_strange)
    pseudoscalar = sp.diag(pseudoscalar_up, pseudoscalar_down, pseudoscalar_strange)
    mass_source = sp.diag(pion_mass**2, pion_mass**2, strange_mass_source)
    charge = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    photon_response = sp.diag(photon_up, photon_down, photon_strange)
    first_outgoing, second_outgoing = sp.zeros(3), sp.zeros(3)
    first_outgoing[1, 0] = sp.sqrt(2)
    second_outgoing[0, 1] = sp.sqrt(2)

    def commutator(first, second):
        return first * second - second * first

    def anticommutator(first, second):
        return first * second + second * first

    vector_generators = [sp.diag(1, -1, 0) / sp.sqrt(2),
                         sp.diag(1, 1, -2) / sp.sqrt(6), sp.eye(3) / sp.sqrt(3)]
    for first in range(3):
        for second in range(first + 1, 3):
            symmetric, antisymmetric = sp.zeros(3), sp.zeros(3)
            symmetric[first, second] = symmetric[second, first] = 1 / sp.sqrt(2)
            antisymmetric[first, second] = -sp.I / sp.sqrt(2)
            antisymmetric[second, first] = sp.I / sp.sqrt(2)
            vector_generators.extend((symmetric, antisymmetric))

    derivative_sources = (2 * axial / alp_scale, -second_outgoing / pion_decay_constant,
                          -first_outgoing / pion_decay_constant)
    mass_sources = (-2 * anticommutator(pseudoscalar, mass_source) / alp_scale,
                    anticommutator(second_outgoing, mass_source) / pion_decay_constant,
                    anticommutator(first_outgoing, mass_source) / pion_decay_constant)
    currents = (commutator(first_outgoing, second_outgoing),
                commutator(first_outgoing, axial), commutator(second_outgoing, axial))
    source_masses = (alp_mass, pion_mass, pion_mass)
    vector_squared_values = (pair_squared, first_squared, second_squared)
    prefactors = (-even_coupling / (sp.sqrt(2) * pion_decay_constant**2),
                  sp.sqrt(2) * even_coupling / (pion_decay_constant * alp_scale),
                  -sp.sqrt(2) * even_coupling / (pion_decay_constant * alp_scale))
    flavour_traces, channel_numerators, exchanges = [], [], []
    for derivative_source, pseudoscalar_source, current, source_mass, vector_squared, prefactor in zip(
            derivative_sources, mass_sources, currents, source_masses,
            vector_squared_values, prefactors):
        trace_sums = []
        for photon_matrix, source in ((photon_response, derivative_source),
                                      (photon_response, pseudoscalar_source),
                                      (charge, derivative_source), (charge, pseudoscalar_source)):
            trace_sums.append(sp.simplify(sum(
                sp.trace(anticommutator(vector, photon_matrix) * source)
                * sp.trace(vector * current) / 2 for vector in vector_generators)))
        flavour_traces.append(sp.Matrix(trace_sums))
        converted = sp.sqrt(2) * electric_charge * vertices['strong'].subs({
            vertices['D']: trace_sums[0], vertices['P']: trace_sums[1],
            m: source_mass, s1: photon_squared, s2: vector_squared})
        direct = -2 * electric_charge * vertices['ordered'].subs({
            vertices['D']: trace_sums[2], vertices['P']: trace_sums[3],
            m: source_mass, s1: vector_squared, s2: photon_squared})
        numerator = sp.factor(prefactor * (converted + direct))
        channel_numerators.append(numerator)
        exchanges.append(numerator / (vector_squared - rho_mass**2))
    # c4 contains a light-quark mass commutator, hence vanishes before
    # projection in strong isospin, for either charged-pion source.
    c4_commutators = [commutator(commutator(generator, charge), mass_source)
                     for generator in (first_outgoing, second_outgoing)]
    amplitude = sum(exchanges)

    # U(3) completeness, Tr(T_R T_S)=delta_RS, contracts each flavour
    # sum to Tr({C,W} U)/2. Keep matrices abstract in the displayed result.
    axial_matrix, pseudoscalar_matrix, mass_matrix, charge_matrix, response_matrix = (
        sp.MatrixSymbol(name, 3, 3) for name in
        ('c_q^A', 'c_q^P', 'chi', 'Q', 'W_gamma'))
    first_pion_matrix, second_pion_matrix = (
        sp.MatrixSymbol(name, 3, 3) for name in
        ('lambda_1^dagger', 'lambda_2^dagger'))
    abstract_derivative_sources = (
        2 * axial_matrix / alp_scale, -second_pion_matrix / pion_decay_constant,
        -first_pion_matrix / pion_decay_constant)
    abstract_mass_sources = (
        -2 * anticommutator(pseudoscalar_matrix, mass_matrix) / alp_scale,
        anticommutator(second_pion_matrix, mass_matrix) / pion_decay_constant,
        anticommutator(first_pion_matrix, mass_matrix) / pion_decay_constant)
    abstract_currents = (commutator(first_pion_matrix, second_pion_matrix),
                         commutator(first_pion_matrix, axial_matrix),
                         commutator(second_pion_matrix, axial_matrix))
    matrix_substitutions = dict(zip(
        (axial_matrix, pseudoscalar_matrix, mass_matrix, charge_matrix,
         response_matrix, first_pion_matrix, second_pion_matrix),
        (axial, pseudoscalar, mass_source, charge, photon_response,
         first_outgoing, second_outgoing)))
    matrix_flavour_traces, matrix_channel_numerators = [], []
    for derivative_source, mass_insertion, current, source_mass, vector_squared, prefactor in zip(
            abstract_derivative_sources, abstract_mass_sources, abstract_currents,
            source_masses, vector_squared_values, prefactors):
        traces = [sp.Trace(anticommutator(current, photon_matrix) * source) / 2
                  for photon_matrix, source in
                  ((response_matrix, derivative_source), (response_matrix, mass_insertion),
                   (charge_matrix, derivative_source), (charge_matrix, mass_insertion))]
        matrix_flavour_traces.append(traces)
        converted = sp.sqrt(2) * electric_charge * vertices['strong'].subs({
            vertices['D']: traces[0], vertices['P']: traces[1],
            m: source_mass, s1: photon_squared, s2: vector_squared})
        direct = -2 * electric_charge * vertices['ordered'].subs({
            vertices['D']: traces[2], vertices['P']: traces[3],
            m: source_mass, s1: vector_squared, s2: photon_squared})
        matrix_channel_numerators.append(prefactor * (converted + direct))

    # Integrate each affine numerator by polynomial division. The three
    # flavour sums make the t and u numerators the same function of r^2.
    vector_squared = sp.Symbol('r_squared', real=True)
    crossed_numerator = channel_numerators[1].subs(first_squared, vector_squared)
    quotient, remainder = map(sp.factor, sp.div(
        crossed_numerator, vector_squared - rho_mass**2, vector_squared))
    cosine = sp.Symbol('z', real=True)
    angular_center, angular_slope = sp.symbols('a b', real=True, nonzero=True)
    primitive = sp.integrate((1 - cosine**2) / (angular_center + angular_slope * cosine), cosine)
    angular_integral = sp.simplify(sp.Rational(3, 4)
        * (primitive.subs(cosine, 1) - primitive.subs(cosine, -1)))
    triangle = (pair_squared - alp_mass**2 - photon_squared)**2 - 4 * alp_mass**2 * photon_squared
    center = (alp_mass**2 + photon_squared + 2 * pion_mass**2 - pair_squared) / 2
    slope = sp.sqrt(1 - 4 * pion_mass**2 / pair_squared) * sp.sqrt(triangle) / 2
    # a=t0-Mrho^2, b=kappa; retaining log(a+b)-log(a-b) specifies
    # analytic continuation from the same pole prescription in each log.
    projected_crossed = 2 * quotient + 2 * remainder * angular_integral
    p_wave = exchanges[0] + projected_crossed.subs({angular_center: center - rho_mass**2,
                                                   angular_slope: slope})
    matrix_amplitude = sum(numerator / (squared - rho_mass**2)
                           for numerator, squared in
                           zip(matrix_channel_numerators, vector_squared_values))
    matrix_crossed_numerator = matrix_channel_numerators[1].subs(
        first_squared, vector_squared)
    matrix_quotient = sp.diff(matrix_crossed_numerator, vector_squared)
    matrix_remainder = matrix_crossed_numerator.subs(vector_squared, rho_mass**2)
    matrix_p_wave = (matrix_channel_numerators[0] / (pair_squared - rho_mass**2)
                     + 2 * matrix_quotient + 2 * matrix_remainder * angular_integral.subs(
                         {angular_center: center - rho_mass**2, angular_slope: slope}))
    return locals()


class VirtualPhotonAmplitudeTests(unittest.TestCase):
    """Check flavour, fitted nulls and the analytic angular integration."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_charged_pion_flavour_and_isospin(self):
        result = self.result
        self.assertTrue(all(matrix == sp.zeros(3) for matrix in result['c4_commutators']))
        self.assertEqual(sp.simplify(result['channel_numerators'][1].subs(
            result['first_squared'], result['second_squared']) - result['channel_numerators'][2]), 0)
        self.assertEqual(result['channel_numerators'][1].subs(
            result['axial_down'], result['axial_up']), 0)
        for expression in result['channel_numerators']:
            self.assertFalse(expression.has(result['axial_strange'], result['pseudoscalar_strange'],
                                            result['photon_strange']))

    def test_matrix_trace_numerators_against_component_derivation(self):
        result = self.result
        for traces, component_traces in zip(result['matrix_flavour_traces'],
                                             result['flavour_traces']):
            for trace, component in zip(traces, component_traces):
                self.assertEqual(sp.simplify(
                    trace.subs(result['matrix_substitutions']).doit() - component), 0)
        for numerator, component in zip(result['matrix_channel_numerators'],
                                         result['channel_numerators']):
            self.assertEqual(sp.simplify(
                numerator.subs(result['matrix_substitutions']).doit() - component), 0)
        for matrix_name, component_name in (('matrix_amplitude', 'amplitude'),
                                            ('matrix_quotient', 'quotient'),
                                            ('matrix_remainder', 'remainder'),
                                            ('matrix_p_wave', 'p_wave')):
            self.assertEqual(sp.simplify(result[matrix_name].subs(
                result['matrix_substitutions']).doit() - result[component_name]), 0)

    def test_existing_null_derivation(self):
        from derive_alp_virtual_photon_pion_nulls import derive as derive_nulls
        result = self.result
        null_result = derive_nulls()
        coefficients = sp.symbols('c_1:8 d_1:5', real=True)
        projection = sp.Matrix([[sp.diff(result['amplitude'], coefficient)
                                 for coefficient in coefficients]]) * result['vertices']['nulls']
        self.assertEqual((projection - null_result['total_projection']).applyfunc(sp.simplify),
                         sp.zeros(1, 4))

    def test_real_photon_against_existing_amplitude(self):
        from derive_alp_photon_two_pseudoscalar import derive as derive_real_photon
        result, real_result = self.result, derive_real_photon()
        replacements = {
            real_result['photon_response']: result['photon_response'],
            real_result['charge']: result['charge'],
            real_result['first_generator']: result['first_outgoing'].H,
            real_result['second_generator']: result['second_outgoing'].H,
            real_result['axial']: result['axial'],
            real_result['pseudoscalar']: result['pseudoscalar'],
            real_result['mass']: result['mass_source'],
            real_result['first_mass']: result['pion_mass'],
            real_result['second_mass']: result['pion_mass'],
            real_result['internal_mass']: result['rho_mass'],
        }
        for index, exchange in enumerate(real_result['exchanges']):
            flavour_sum = sum(exchange.subs(replacements).subs(
                real_result['internal_vector'], vector).doit()
                for vector in result['vector_generators'])
            self.assertEqual(sp.simplify(flavour_sum - result['exchanges'][index].subs(
                result['photon_squared'], 0)), 0)

    def test_projection_against_independent_quadrature(self):
        from scipy.integrate import quad
        result = self.result
        for center, slope in ((5., 1.), (3., 0.2)):
            symbolic = float(result['angular_integral'].subs({result['angular_center']: center,
                                                              result['angular_slope']: slope}))
            numerical = quad(lambda z: .75 * (1-z*z) / (center+slope*z), -1, 1)[0]
            self.assertAlmostEqual(symbolic, numerical, places=10)

    def test_zero_angular_slope_limit(self):
        result = self.result
        center, slope = result['angular_center'], result['angular_slope']
        # Positive center fixes the real branch for this elementary limit.
        positive_center = sp.Symbol('positive_center', positive=True)
        limit = sp.limit(result['angular_integral'].subs(center, positive_center), slope, 0)
        self.assertEqual(sp.simplify(limit - 1 / positive_center), 0)


def main():
    """Print derivation steps, or run the focused unit tests."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    if parser.parse_args().test:
        unittest.main(argv=[__file__])
        return
    result = derive()
    print(__doc__)
    show('Original VVP contraction', result['vertices']['strong'])
    show('Original ordered VJP contraction', result['vertices']['ordered'])
    for name, traces, numerator in zip(('s', 't', 'u'), result['matrix_flavour_traces'],
                                       result['matrix_channel_numerators']):
        show(name + ' U(3) traces (D_W,P_W,D_Q,P_Q)', sp.Matrix(traces))
        show(name + ' U(3) trace numerator', numerator)
    show('Full retained U(3) scalar amplitude', result['matrix_amplitude'])
    show('U(3) crossed polynomial quotient', result['matrix_quotient'])
    show('U(3) crossed pole residue', result['matrix_remainder'])
    show('Normalized angular pole integral', result['angular_integral'])
    show('Angular center t0', result['center'])
    show('Angular slope kappa', result['slope'])
    show('U(3) P wave of retained amplitude', result['matrix_p_wave'])


if __name__ == '__main__':
    main()
