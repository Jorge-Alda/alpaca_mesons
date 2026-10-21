#!/usr/bin/env python3
"""Derive flavour support of the recorded strong a -> (omega,phi) PP trees.

Run with ``.venv/bin/python scripts/derive_alp_vpp_selection_rules.py``.
All traces are dimensionless. Masses, pole denominators and the independent
d1,d2,d3 structures are kept distinct: a channel is allowed when at least one
coefficient is nonzero. Degenerate real components of charged intermediate
fields are summed before testing a coefficient. No numerical tolerance is used.
The external generator diag(v_n,v_n,v_s) covers both fitted isoscalar states;
ideal mixing is a separate specialization. Overall generator normalizations
cannot affect the resulting support.
"""

from itertools import combinations_with_replacement

import sympy as sp


def matrix_unit(row, column):
    """Return the elementary three-flavour matrix."""
    result = sp.zeros(3)
    result[row, column] = 1
    return result


def commutator(first, second):
    """Return the matrix commutator."""
    return first * second - second * first


def anticommutator(first, second):
    """Return the matrix anticommutator."""
    return first * second + second * first


def trace(expression):
    """Expand a finite flavour trace exactly."""
    return sp.expand(sp.trace(expression))


cosine, sine = sp.symbols('cos_theta_P sin_theta_P', real=True)
vector_nonstrange, vector_strange = sp.symbols('v_n v_s', real=True)
external_vector = sp.diag(vector_nonstrange, vector_nonstrange, vector_strange)
quark_mass = sp.diag(*sp.symbols('chi_u chi_d chi_s', real=True))
pion = sp.diag(1, -1, 0) / sp.sqrt(2)
octet = sp.diag(1, 1, -2) / sp.sqrt(6)
singlet = sp.eye(3) / sp.sqrt(3)
eta = cosine * octet - sine * singlet
eta_prime = sine * octet + cosine * singlet
generators = {'pi0': pion, 'eta': eta, 'eta_prime': eta_prime}
charges = {'pi0': 0, 'eta': 0, 'eta_prime': 0}
for positive, negative, row, column, charge in (
        ('pi+', 'pi-', 0, 1, 1), ('K+', 'K-', 0, 2, 1),
        ('K0', 'K0bar', 1, 2, 0)):
    generators[positive] = matrix_unit(row, column)
    generators[negative] = matrix_unit(column, row)
    charges[positive], charges[negative] = charge, -charge


def real_pair(row, column):
    """Return the two normalized Hermitian components of one charged field."""
    first, second = matrix_unit(row, column), matrix_unit(column, row)
    return ((first + second) / sp.sqrt(2),
            -sp.I * (first - second) / sp.sqrt(2))


vector_groups = ((pion,), (sp.diag(*sp.symbols('omega_n omega_n omega_s')),),
                 (sp.diag(*sp.symbols('phi_n phi_n phi_s')),),
                 real_pair(0, 1), real_pair(0, 2), real_pair(1, 2))
pseudoscalar_groups = ((pion,), (eta,), (eta_prime,),
                       real_pair(0, 1), real_pair(0, 2), real_pair(1, 2))
diagonal_source = sp.diag(*sp.symbols('c_u c_d c_s'))
neutral_strange_source = (sp.Symbol('c_ds') * matrix_unit(1, 2)
                          + sp.Symbol('c_sd') * matrix_unit(2, 1))
# Neutral ALP and neutral external vector require an electrically neutral source.
# The us/su entries allowed by |Delta S| alone carry electric charge.
strange_source = neutral_strange_source


def odd_coefficients(first, second, axial_source, pseudoscalar_source):
    """Derive independent s,t,u trace coefficients from the amplitude record."""
    result = []
    for group in vector_groups:
        for source in (axial_source,
                       anticommutator(pseudoscalar_source, quark_mass)):
            result.append(sum(trace(anticommutator(external_vector, internal) * source)
                              * trace(internal * commutator(first, second))
                              for internal in group))
        for emitted, spectator in ((second, first), (first, second)):
            for source in (emitted, anticommutator(emitted, quark_mass)):
                result.append(sum(trace(anticommutator(external_vector, internal) * source)
                                  * trace(internal * commutator(spectator, axial_source))
                                  for internal in group))
    return tuple(sp.expand(value) for value in result)


def even_coefficients(first, second):
    """Derive both scalar-exchange trace coefficients for each pole group."""
    result = []
    for group in pseudoscalar_groups:
        for spectator, emitted in ((second, first), (first, second)):
            result.append(sp.expand(sum(
                trace(anticommutator(spectator, internal)
                      * commutator(strange_source, quark_mass))
                * trace(external_vector * commutator(internal, emitted))
                for internal in group)))
    return tuple(result)


def derive_channels():
    """Enumerate every unordered pair, including charge-forbidden controls."""
    result = {'even': {}, 'odd_diagonal': {}, 'odd_neutral_strange': {},
              'odd_neutral_strange_pseudoscalar_only': {}}
    for first_name, second_name in combinations_with_replacement(generators, 2):
        first, second = generators[first_name].H, generators[second_name].H
        pair = (first_name, second_name)
        result['even'][pair] = even_coefficients(first, second)
        result['odd_diagonal'][pair] = odd_coefficients(
            first, second, diagonal_source, diagonal_source)
        result['odd_neutral_strange'][pair] = odd_coefficients(
            first, second, neutral_strange_source, sp.zeros(3))
        result['odd_neutral_strange_pseudoscalar_only'][pair] = odd_coefficients(
            first, second, sp.zeros(3), neutral_strange_source)
    return result


def allowed_channels(coefficients, substitutions=None):
    """Return channels with at least one nonzero independent coefficient."""
    substitutions = substitutions or {}
    return {pair for pair, values in coefficients.items()
            if any(sp.expand(value.subs(substitutions)) != 0 for value in values)}


if __name__ == '__main__':
    DERIVED = derive_channels()
    print('Exact flavour coefficients; generic eta mixing and diagonal chi; natural units.')
    print('External neutral vector generator:', sp.latex(external_vector))
    for sector, coefficients in DERIVED.items():
        print(sector, sorted(allowed_channels(coefficients)))
        print('ideal omega', sorted(allowed_channels(coefficients,
              {vector_nonstrange: 1 / sp.sqrt(2), vector_strange: 0})))
        print('ideal phi', sorted(allowed_channels(coefficients,
              {vector_nonstrange: 0, vector_strange: 1})))
        print('exact isospin', sorted(allowed_channels(coefficients,
              {quark_mass[0, 0]: quark_mass[1, 1]})))
        if sector != 'odd_diagonal':
            print('c_ds only', sorted(allowed_channels(coefficients,
                  {sp.Symbol('c_sd'): 0})))
            print('c_sd only', sorted(allowed_channels(coefficients,
                  {sp.Symbol('c_ds'): 0})))
        for pair in sorted(allowed_channels(coefficients)):
            nonzero = [sp.factor(value) for value in coefficients[pair] if value != 0]
            print(pair, 'representative trace =', sp.latex(nonzero[0]))
