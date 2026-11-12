#!/usr/bin/env python3
"""Which source Taylor coefficients enter the exact LO quadratic blocks K and M^2.

Run: .venv/bin/python scripts/check_alp_quadratic_taylor_dependence.py

Conventions (Alda:2026nrd, as in scripts/derive_alp_lo_joint_vertices.py):
natural units; U = u^2 = exp(i Phi/f_pi); x = a/f_a; chi(x) = 2 B_0 M(x),
M = S + i P with S, P Hermitian; r_mu = (d_mu x) C_R(x), l_mu = (d_mu x) C_L(x);
L = f_pi^2/4 tr(u_mu u^mu + chi_+) - f_pi^2 m_0^2/(4 N) (-tr Phi/f_pi - 2 Theta(x))^2
    + (d a d a - m_a0^2 a^2)/2,  N = 3.
Taylor coefficients are derivatives at x = 0: F(x) = sum_n F_n x^n / n!, n = 0..3,
with Theta_0 = 0 (aligned vacuum, chi(0) Hermitian is NOT imposed on the symbols
here; the ALP tadpole is not discussed).

All coefficient matrices are fully generic symbolic 3x3 Hermitian matrices
(S_n, P_n, C_{R,n}, C_{L,n}); the two meson generators are generic symbolic
Hermitian matrices. Two independent constructions are compared:
  A. quadratic_blocks from derive_alp_lo_joint_vertices.py (momentum space);
  B. an independent position-space expansion in a field-scaling parameter
     epsilon (Phi, d Phi, a, d a -> epsilon times themselves) keeping every
     source Taylor term n = 0..3, then extracting the epsilon^2 part and taking
     second derivatives in the field values and field derivatives. B also
     reports the mixed (field)(field derivative) terms, which would be total
     derivatives such as d_mu(a^2)/2.
"""
from __future__ import annotations

import os
import sys
from math import factorial

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from derive_alp_lo_joint_vertices import quadratic_blocks  # noqa: E402

SIZE = 3
MAX_TAYLOR = 3
f_pi, f_a, B_0, m_0_sq, m_a0_sq = sp.symbols('f_pi f_a B_0 m0sq ma0sq', positive=True)


def hermitian_matrix(label):
    """Generic Hermitian matrix with real symbols label_d{i}, label_r{ij}, label_i{ij}."""
    matrix = sp.zeros(SIZE)
    for i in range(SIZE):
        matrix[i, i] = sp.Symbol(f'{label}_d{i}', real=True)
        for j in range(i + 1, SIZE):
            real_part = sp.Symbol(f'{label}_r{i}{j}', real=True)
            imaginary_part = sp.Symbol(f'{label}_i{i}{j}', real=True)
            matrix[i, j] = real_part + sp.I * imaginary_part
            matrix[j, i] = real_part - sp.I * imaginary_part
    return matrix


scalar = [hermitian_matrix(f'S{n}') for n in range(MAX_TAYLOR + 1)]
pseudoscalar = [hermitian_matrix(f'P{n}') for n in range(MAX_TAYLOR + 1)]
chi_coefficients = [2 * B_0 * (scalar[n] + sp.I * pseudoscalar[n]) for n in range(MAX_TAYLOR + 1)]
right_coefficients = [hermitian_matrix(f'CR{n}') for n in range(MAX_TAYLOR + 1)]
left_coefficients = [hermitian_matrix(f'CL{n}') for n in range(MAX_TAYLOR + 1)]
theta_coefficients = [sp.S.Zero] + [sp.Symbol(f'Theta{n}', real=True) for n in range(1, MAX_TAYLOR + 1)]
generator_one = hermitian_matrix('GA')
generator_two = hermitian_matrix('GB')

FIELD_LABELS = ['P_A', 'P_B', 'a']


def coefficient_label(symbol):
    """Map an entry symbol to the coefficient it belongs to (M_n, C_R,n, ...)."""
    name = symbol.name
    head = name.split('_')[0]
    for prefix, label in (('Theta', 'Theta_'), ('CR', 'C_R,'), ('CL', 'C_L,'),
                          ('S', 'S_'), ('P', 'P_')):
        if head.startswith(prefix) and head[len(prefix):].isdigit():
            return f'{label}{head[len(prefix):]}'
    if head in ('GA', 'GB'):
        return None
    return name


def dependence(expression):
    labels = {coefficient_label(s) for s in sp.sympify(expression).free_symbols}
    return sorted(label for label in labels if label is not None)


# ---------------- construction B: independent epsilon expansion ----------------
def truncated_product(left, right, degree):
    return [sum((left[k] * right[n - k] for k in range(n + 1)), sp.zeros(SIZE))
            for n in range(degree + 1)]


def scalar_series_times_matrix(scalar_series, matrix, degree):
    return [sp.expand(scalar_series[n]) * matrix for n in range(degree + 1)]


def construction_b():
    """Return K, M^2 and the mixed field/derivative block from the epsilon expansion."""
    degree = 2
    epsilon = sp.Symbol('epsilon')
    field_values = sp.symbols('phiA phiB alp', real=True)
    field_derivatives = sp.symbols('dphiA dphiB dalp', real=True)
    phi = field_values[0] * generator_one + field_values[1] * generator_two
    dphi = field_derivatives[0] * generator_one + field_derivatives[1] * generator_two
    identity, zero = sp.eye(SIZE), sp.zeros(SIZE)

    def exponential_parts(sign):
        """Degree-graded u (sign=+1) or u^dagger (sign=-1) and its derivative."""
        step = sign * sp.I * phi / (2 * f_pi)
        dstep = sign * sp.I * dphi / (2 * f_pi)
        value = [identity] + [step**n / factorial(n) for n in range(1, degree + 1)]
        derivative = [zero] + [sum((step**j * dstep * step**(n - j - 1) for j in range(n)), zero)
                               / factorial(n) for n in range(1, degree + 1)]
        return value, derivative

    u, du = exponential_parts(+1)
    ud, dud = exponential_parts(-1)
    # x and d x as degree-graded scalars; x^n/n! has degree n.
    x_power = [[sp.S.Zero] * (degree + 1) for _ in range(MAX_TAYLOR + 1)]
    for n in range(MAX_TAYLOR + 1):
        if n <= degree:
            x_power[n][n] = (field_values[2] / f_a)**n / factorial(n)
    dx = [sp.S.Zero, field_derivatives[2] / f_a] + [sp.S.Zero] * (degree - 1)

    def source_function(coefficients):
        series = [zero] * (degree + 1)
        for n in range(MAX_TAYLOR + 1):  # n = 3 kept; it only feeds degree 3
            term = scalar_series_times_matrix(x_power[n], coefficients[n], degree)
            series = [series[k] + term[k] for k in range(degree + 1)]
        return series

    def times_dx(series):
        return [sum((dx[k] * series[n - k] for k in range(n + 1)), zero) for n in range(degree + 1)]

    right = times_dx(source_function(right_coefficients))
    left = times_dx(source_function(left_coefficients))
    chi = source_function(chi_coefficients)
    chi_dagger = [matrix.H for matrix in chi]
    first = truncated_product(ud, du, degree)
    second = truncated_product(u, dud, degree)
    third = truncated_product(truncated_product(ud, right, degree), u, degree)
    fourth = truncated_product(truncated_product(u, left, degree), ud, degree)
    u_mu = [sp.I * (first[n] - second[n]) + third[n] - fourth[n] for n in range(degree + 1)]
    chi_plus = [p + q for p, q in zip(truncated_product(truncated_product(ud, chi, degree), ud, degree),
                                      truncated_product(truncated_product(u, chi_dagger, degree), u, degree))]
    kinetic_trace = sum((sp.trace(u_mu[k] * u_mu[degree - k]) for k in range(degree + 1)), sp.S.Zero)
    theta = sum(theta_coefficients[n] * (epsilon * field_values[2] / f_a)**n / factorial(n)
                for n in range(1, MAX_TAYLOR + 1))
    topological = sp.expand((-epsilon * sp.trace(phi) / f_pi - 2 * theta)**2).coeff(epsilon, degree)
    bare = (field_derivatives[2]**2 - m_a0_sq * field_values[2]**2) / 2
    lagrangian = sp.expand(f_pi**2 / 4 * (kinetic_trace + sp.trace(chi_plus[degree]))
                           - f_pi**2 * m_0_sq / (4 * SIZE) * topological + bare)
    kinetic = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff(lagrangian, field_derivatives[i], field_derivatives[j])))
    mass_squared = sp.Matrix(3, 3, lambda i, j: sp.expand(-sp.diff(lagrangian, field_values[i], field_values[j])))
    mixed = sp.Matrix(3, 3, lambda i, j: sp.expand(sp.diff(lagrangian, field_values[i], field_derivatives[j])))
    return kinetic, mass_squared, mixed


def main():
    kinetic_b, mass_b, mixed_b = construction_b()
    kinetic_a, mass_a = quadratic_blocks([generator_one, generator_two], chi_coefficients,
                                         right_coefficients, left_coefficients, theta_coefficients,
                                         f_pi, f_a, m_0_sq, m_a0_sq)
    agree = all(sp.expand(kinetic_a[i, j] - kinetic_b[i, j]) == 0 and sp.expand(mass_a[i, j] - mass_b[i, j]) == 0
                for i in range(3) for j in range(3))
    print('Fully symbolic generic Hermitian 3x3 coefficients S_n, P_n, C_R,n, C_L,n (n=0..3), Theta_n (n=1..3),')
    print('generic Hermitian meson generators G_A, G_B. M_n = S_n + i P_n; S_n, P_n listed separately.')
    print(f'Construction A (derive_alp_lo_joint_vertices.quadratic_blocks) == construction B (epsilon expansion): {agree}')
    print(f'Mixed (field)(field-derivative) terms in B, i.e. total derivatives such as d(a^2)/2: '
          f'{"none (all zero)" if mixed_b == sp.zeros(3) else mixed_b}')
    print('\nDependence table (coefficient labels appearing in each entry):')
    entries = [('P-P', 0, 1), ('P-P diag', 0, 0), ('P-a', 0, 2), ('a-a', 2, 2)]
    for name, i, j in entries:
        print(f'  {name:9s} kinetic K : {dependence(kinetic_b[i, j])}')
        print(f'  {name:9s} mass    M2: {dependence(mass_b[i, j])}')

    # (iii) exact M_2 term in M^2_aa, and P_2 absent.
    second_order = {s: 0 for s in (scalar[2].free_symbols | pseudoscalar[2].free_symbols)}
    m2_part = sp.expand(mass_b[2, 2] - mass_b[2, 2].subs(second_order))
    claim = -f_pi**2 / (2 * f_a**2) * sp.re(sp.trace(chi_coefficients[2]))
    claim_s = -f_pi**2 / (2 * f_a**2) * 2 * B_0 * sp.trace(scalar[2])
    print(f'\nM2_aa(M_2 part) - (-(f_pi^2/(2 f_a^2)) Re tr chi_2) = {sp.simplify(m2_part - claim)}')
    print(f'M2_aa(M_2 part) - (-(f_pi^2/f_a^2) B_0 tr S_2)        = {sp.expand(m2_part - claim_s)}')
    print(f'P_2 symbols in M2_aa: {sorted(str(s) for s in mass_b[2, 2].free_symbols & pseudoscalar[2].free_symbols)}')
    # Theta_1 parts.
    theta_one = theta_coefficients[1]
    print(f'Theta_1 part of M2_Pa: {sp.factor(sp.expand(mass_b[0, 2] - mass_b[0, 2].subs(theta_one, 0)))}')
    print(f'Theta_1 part of M2_aa: {sp.factor(sp.expand(mass_b[2, 2] - mass_b[2, 2].subs(theta_one, 0)))}')
    # C_R,0 and C_L,0 only through C_A,0 = (C_R,0 - C_L,0)/2?
    vector = hermitian_matrix('CV0')
    substitution = {}
    for symbol in right_coefficients[0].free_symbols:
        tail = symbol.name[len('CR0'):]
        substitution[symbol] = sp.Symbol('CV0' + tail, real=True) + sp.Symbol('CA0' + tail, real=True)
    for symbol in left_coefficients[0].free_symbols:
        tail = symbol.name[len('CL0'):]
        substitution[symbol] = sp.Symbol('CV0' + tail, real=True) - sp.Symbol('CA0' + tail, real=True)
    vector_symbols = vector.free_symbols
    vector_free = all(not (sp.expand(entry.subs(substitution)).free_symbols & vector_symbols)
                      for entry in list(kinetic_b) + list(mass_b))
    print(f'C_R,0 and C_L,0 enter only via C_A,0 = (C_R,0 - C_L,0)/2 (C_V,0 drops out): {vector_free}')


if __name__ == '__main__':
    main()
