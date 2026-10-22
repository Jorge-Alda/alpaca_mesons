"""Derive no-mixing a -> omega gamma and phi gamma at leading e.

Run: .venv/bin/python scripts/derive_physical_alp_vgamma.py
Tests: .venv/bin/python scripts/derive_physical_alp_vgamma.py --test

Natural units (GeV), metric +---, canonical epsilon^0123=+1. The
original eleven odd operators and the declared even electromagnetic terms
are retained, with narrow external vectors and no extra ALP meson poles.
The imported indexed contractions use the opposite epsilon orientation;
the -2 factors below include that adapter and the diagonal anticommutator.
"""
import argparse
import json
from pathlib import Path
import unittest

import sympy as sp

from derive_radiative_tff_amplitude import (
    contractions, e1, e2, m, q1, q2, reduce, s1, s2)


def show(label, expression):
    """Print one derived expression with SymPy's LaTeX renderer."""
    print(label + ':\n' + sp.latex(expression) + '\n')


def derive():
    """Build the Lorentz, flavor, mass, response and null calculations."""
    c = sp.symbols('c_1:8', real=True)
    d = sp.symbols('d_1:5', real=True)
    derivative_trace, mass_trace = sp.symbols('D P', real=True)
    operator_mass, external_mass, electric_charge = sp.symbols(
        'M_V M_ext e', positive=True)
    alp_mass = sp.Symbol('m_a', positive=True)
    first = contractions(q1, q2, e1, e2)
    second = contractions(q2, q1, e2, e1)
    ordered_vjp = reduce(-2/operator_mass * sum(
        c[index]*(mass_trace if index == 2 else derivative_trace)*first[index]
        for index in range(7))).subs(m, alp_mass)
    ordered_vvp = reduce(-2*(
        derivative_trace*(d[0]*(first[0]+second[0])
                          + d[2]*(first[4]+second[4])
                          + d[3]*(first[6]+second[6]))
        + mass_trace*d[1]*(first[2]+second[2]))).subs(m, alp_mass)

    charge = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    light_spurion, strange_spurion = sp.symbols('x_l x_s', real=True)
    vacuum_spurion = sp.diag(light_spurion, light_spurion, strange_spurion)
    vector_source, single_source, double_source = sp.symbols('F_V z_V y_V', real=True)
    source = (vector_source*charge
              + single_source*(vacuum_spurion*charge+charge*vacuum_spurion)
              + double_source*sp.trace(vacuum_spurion*charge)*sp.eye(3))
    ideal_generators = (
        sp.diag(1, -1, 0)/sp.sqrt(2),
        sp.diag(1, 1, 0)/sp.sqrt(2), sp.diag(0, 0, 1))
    angle, decay_constant = sp.symbols('theta_V f_a', real=True)
    physical_generators = (
        ideal_generators[0],
        sp.cos(angle)*ideal_generators[1]-sp.sin(angle)*ideal_generators[2],
        sp.sin(angle)*ideal_generators[1]+sp.cos(angle)*ideal_generators[2])
    # Physical pole masses and mixing angle are independent inputs here.
    # No vector-mass fit relation is imposed by this derivation.
    pole_masses = sp.symbols('M_rho M_omega M_phi', positive=True)
    mass_squared = sp.diag(*(pole_mass**2 for pole_mass in pole_masses))
    residues = sp.Matrix([sp.trace(generator*source)
                          for generator in physical_generators])
    response_fields = sp.Matrix(sp.symbols('r_rho r_omega r_phi'))
    # Factor F_gamma^2/4. Differentiating gives M^2 r=sqrt(2)e R.
    response_lagrangian = ((response_fields.T*mass_squared*response_fields)[0]
                       - 2*sp.sqrt(2)*electric_charge*(residues.T*response_fields)[0])
    response_equations = sp.Matrix([sp.diff(response_lagrangian, field)/2
                                   for field in response_fields])
    response_solution = sp.solve(list(response_equations), list(response_fields))
    response_components = sp.Matrix([
        sp.simplify(response_solution[field]/(sp.sqrt(2)*electric_charge))
        for field in response_fields])
    photon_response = sum((entry*generator for entry, generator
                           in zip(response_components, physical_generators)), sp.zeros(3))
    axial = sp.diag(*sp.symbols('c_u_A c_d_A c_s_A', real=True))
    pseudoscalar = sp.diag(*sp.symbols('c_u_P c_d_P c_s_P', real=True))
    chi = sp.diag(*sp.symbols('chi_u chi_d chi_s', real=True))
    derivative_source = 2*axial/decay_constant
    mass_source = -2*(pseudoscalar*chi+chi*pseudoscalar)/decay_constant
    generators = {
        'omega': physical_generators[1], 'phi': physical_generators[2]}
    # Keep the response entries independent for readable final expressions;
    # the exact mass/source expressions above define each entry.
    response_diagonal = sp.diag(*sp.symbols('W_u W_d W_s', real=True))
    amplitudes = {}
    contributions = {}
    for name, generator in generators.items():
        direct = (-2*electric_charge/external_mass*ordered_vjp).subs({
            derivative_trace: sp.trace(generator*charge*derivative_source),
            mass_trace: sp.trace(generator*charge*mass_source),
            s1: external_mass**2, s2: 0})
        converted = (sp.sqrt(2)*electric_charge/external_mass*ordered_vvp).subs({
            derivative_trace: sp.trace(generator*response_diagonal*derivative_source),
            mass_trace: sp.trace(generator*response_diagonal*mass_source),
            s1: external_mass**2, s2: 0})
        contributions[name] = (sp.factor(direct), sp.factor(converted))
        amplitudes[name] = sp.factor(direct+converted)
    parameters = c+d
    nulls = sp.zeros(11, 4)
    for column, row in enumerate((6, 10)):
        nulls[row, column] = 1
    for row in (1, 4, 5):
        nulls[row, 2] = 1
    for row in (0, 4):
        nulls[row, 3] = 1
    gradient = sp.Matrix([[sp.diff(amplitudes[name], coefficient)
                           for coefficient in parameters] for name in generators])
    projection = (gradient*nulls).applyfunc(sp.simplify)
    from summarize_odd_preliminary_fit import main as reconstruct_central_fit
    _, fitted_columns = reconstruct_central_fit(write_output=False)
    exact_projector = nulls*(nulls.T*nulls).inv()*nulls.T
    residual = (sp.eye(11)-exact_projector)*sp.Matrix(fitted_columns)
    stored_null_residuals = {'central_single_trace_negative': float(residual.norm())}
    assert projection == sp.zeros(2, 4)
    assert max(stored_null_residuals.values()) < 1e-12
    return locals()


class AmplitudeTests(unittest.TestCase):
    """Independent identities relevant to the physical amplitude."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_four_fit_nulls(self):
        self.assertEqual(self.result['projection'], sp.zeros(2, 4))

    def test_photon_ward_identity_before_scalar_projection(self):
        self.assertEqual(contractions(q1, q2, e1, q2, False), [0]*7)
        self.assertEqual(contractions(q2, q1, q2, e1, False), [0]*7)

    def test_response_stationarity(self):
        result = self.result
        equations = result['response_equations'].subs(dict(zip(
            result['response_fields'], sp.sqrt(2)*result['electric_charge']
            * result['response_components'])), simultaneous=True)
        self.assertEqual(equations.applyfunc(sp.simplify), sp.zeros(3, 1))

    def test_universal_common_mass_response(self):
        result = self.result
        common_mass = sp.Symbol('M_common', positive=True)
        response = result['photon_response'].subs({
            **{pole_mass: common_mass for pole_mass in result['pole_masses']},
            result['single_source']: 0, result['double_source']: 0})
        difference = response-result['vector_source']*result['charge']/common_mass**2
        self.assertEqual(difference.applyfunc(sp.trigsimp), sp.zeros(3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=['derive_physical_alp_vgamma.py'])
        return
    result = derive()
    print('Tree level, leading electromagnetic order; masses in GeV. '
          'M = G epsilon(eps_V,eps_gamma,q_V,q_gamma).')
    for label, key in (
            ('Indexed reference VJP contractions', 'first'),
            ('Canonical ordered VJP numerator', 'ordered_vjp'),
            ('Canonical VVP numerator, both assignments', 'ordered_vvp'),
            ('Even electromagnetic source S_gamma', 'source'),
            ('Input physical diagonal pole mass-squared matrix', 'mass_squared'),
            ('Signed physical rho, omega, phi source residues', 'residues'),
            ('Zero-virtuality response equations', 'response_equations'),
            ('Exact zero-virtuality W matrix', 'photon_response'),
            ('No-mixing derivative source', 'derivative_source'),
            ('No-mixing mass source', 'mass_source'),
            ('Fit coefficient null columns c7,d4,common shift,c1=c5 shift', 'nulls'),
            ('Two physical amplitude gradients times null columns', 'projection')):
        expression = result[key]
        show(label, sp.Matrix(expression) if isinstance(expression, list) else expression)
    for label, residual in result['stored_null_residuals'].items():
        show(label+' stored-null distance to the exact null subspace', sp.Float(residual))
    for name, pieces in result['contributions'].items():
        physical_mass = result['pole_masses'][1 if name == 'omega' else 2]
        show(name+' gamma external mass prescription',
             sp.Eq(result['external_mass'], physical_mass))
        substitution = {result['external_mass']: physical_mass}
        show(name+' direct G in no-mixing quark couplings', pieces[0].subs(substitution))
        show(name+' VVP photon-conversion G in no-mixing quark couplings',
             pieces[1].subs(substitution))
        show(name+' complete G in no-mixing quark couplings',
             result['amplitudes'][name].subs(substitution))


if __name__ == '__main__':
    main()
