#!/usr/bin/env python3
"""Derive finite CCWZ source covariance; run with .venv/bin/python."""

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    """Use the product rule for arbitrary invertible, noncommuting matrices."""
    u, g, h, source = sp.symbols("u g h r", commutative=False)
    du, dg, dh = sp.symbols("du dg dh", commutative=False)
    mass = sp.Symbol("chi", commutative=False)
    transformed_u = g * u / h
    derivative = dg * u / h + g * du / h - g * u / h * dh / h
    transformed_source = g * source / g - sp.I * dg / g
    transformed_one_form = sp.expand(
        (h / u / g) * (derivative - sp.I * transformed_source * transformed_u)
    )
    original_one_form = u ** -1 * (du - sp.I * source * u)
    expected_one_form = sp.expand(h * original_one_form / h - dh / h)
    assert sp.expand(transformed_one_form - expected_one_form) == 0
    right_form, left_form = sp.symbols("R L", commutative=False)
    transformed_right = h * right_form / h - dh / h
    transformed_left = h * left_form / h - dh / h
    transformed_vielbein_abstract = sp.expand(sp.I * (transformed_right - transformed_left))
    transformed_connection = sp.expand((transformed_right + transformed_left) / 2)
    assert sp.expand(transformed_vielbein_abstract - h * sp.I * (right_form - left_form) / h) == 0
    assert sp.expand(transformed_connection - h * (right_form + left_form) / (2 * h) + dh / h) == 0

    # u' = g_R u h^-1 = h u g_L^-1 implies these two g matrices.
    uprime = sp.Symbol("uprime", commutative=False)
    gr = uprime * h / u
    gl = uprime ** -1 * h * u
    transformed_U = sp.expand(gr * u ** 2 * gl ** -1)
    assert transformed_U == uprime ** 2
    transformed_mass = sp.expand(gr * mass * gl ** -1)
    transformed_plus_first = sp.expand(uprime ** -1 * transformed_mass * uprime ** -1)
    expected_mass = sp.expand(h / u * mass / u / h)
    assert sp.expand(transformed_plus_first - expected_mass) == 0
    mass_dagger = sp.Symbol("chi_dagger", commutative=False)
    transformed_mass_dagger = sp.expand(gl * mass_dagger * gr ** -1)
    transformed_plus_second = sp.expand(uprime * transformed_mass_dagger * uprime)
    expected_mass_dagger = sp.expand(h * u * mass_dagger * u / h)
    assert sp.expand(transformed_plus_second - expected_mass_dagger) == 0

    # Symbolic U(2) check of cyclic trace, using an arbitrary invertible matrix.
    entries_h = sp.symbols("h_11 h_12 h_21 h_22")
    entries_vielbein = sp.symbols("u_11 u_12 u_21 u_22")
    entries_mass = sp.symbols("chi_11 chi_12 chi_21 chi_22")
    h_matrix = sp.Matrix(2, 2, entries_h)
    vielbein = sp.Matrix(2, 2, entries_vielbein)
    mass_block = sp.Matrix(2, 2, entries_mass)
    transformed_vielbein = h_matrix * vielbein * h_matrix.inv()
    transformed_mass_block = h_matrix * mass_block * h_matrix.inv()
    trace_vielbein = sp.factor(sp.trace(transformed_vielbein ** 2) - sp.trace(vielbein ** 2))
    trace_mass = sp.factor(sp.trace(transformed_mass_block) - sp.trace(mass_block))
    assert trace_vielbein == trace_mass == 0
    return original_one_form, transformed_one_form, transformed_U, transformed_vielbein_abstract, transformed_connection, transformed_plus_first, transformed_plus_second, trace_vielbein, trace_mass


def axial_map_expansion():
    """Compare the symmetric axial map with its linear coordinate shift."""
    a, meson_order = sp.symbols("a meson_order", commutative=True)
    phi, generator = sp.symbols("Phi W", commutative=False)
    group = 1 - sp.I * a * generator / 2
    meson_matrix = sum((sp.I * meson_order * phi) ** n / sp.factorial(n) for n in range(4))
    shifted_exponent = sp.I * (meson_order * phi - a * generator)
    shifted_matrix = sum(shifted_exponent ** n / sp.factorial(n) for n in range(4))
    difference = sp.expand(group * meson_matrix * group - shifted_matrix)
    return (sp.expand(difference.coeff(a, 1).coeff(meson_order, n)) for n in (1, 2))


def main():
    print("Natural units; arbitrary local chiral matrices; no meson expansion.")
    original, transformed, transformed_U, vielbein, connection, mass, mass_dagger, trace_vielbein, trace_mass = derive()
    show("Original right one-form", original)
    show("Transformed right one-form", transformed)
    show("Transformed meson matrix", transformed_U)
    show("Transformed vielbein", vielbein)
    show("Transformed connection", connection)
    show("Transformed first mass building block", mass)
    show("Transformed second mass building block", mass_dagger)
    show("Trace difference for two vielbeins", trace_vielbein)
    show("Trace difference for the mass building block", trace_mass)
    axial_quadratic, axial_cubic = axial_map_expansion()
    assert axial_quadratic == 0
    assert axial_cubic != 0
    show("Axial map difference at one ALP and one meson", axial_quadratic)
    show("Axial map difference at one ALP and two mesons", axial_cubic)
    print("The left one-form follows by u -> u^-1 and r -> ell.")
    print("Their difference and sum give u_mu and Gamma_mu; the adjoint mass term gives chi_+/-.")


if __name__ == "__main__":
    main()
