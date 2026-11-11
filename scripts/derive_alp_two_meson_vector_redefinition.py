#!/usr/bin/env python3
"""Derive a U(2) vector redefinition of an ALP two-meson vertex."""

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    """Compare the vertex before and after a nonlinear vector rotation."""
    alp, first_meson, second_meson = sp.symbols("x P_1 P_2", real=True)
    alp_derivative, first_derivative, second_derivative = sp.symbols(
        "dx dP_1 dP_2", real=True
    )
    first_mass_squared, second_mass_squared = sp.symbols("m_1_sq m_2_sq", real=True)
    coupling = sp.Symbol("c_q_VS", real=True, nonzero=True)

    sigma_one = sp.Matrix([[0, 1], [1, 0]])
    sigma_two = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma_three = sp.diag(1, -1)
    meson_matrix = first_meson * sigma_three + second_meson * sigma_one
    generator = sigma_two / 2
    variation = sp.expand(sp.I * (generator * meson_matrix - meson_matrix * generator))
    first_variation = sp.simplify(sp.trace(sigma_three * variation) / 2)
    second_variation = sp.simplify(sp.trace(sigma_one * variation) / 2)
    assert first_variation == second_meson
    assert second_variation == -first_meson

    # Invert P'_1=P_1+xP_2, P'_2=P_2-xP_1 to first order in x.
    old_first = first_meson - alp * second_meson
    old_second = second_meson + alp * first_meson
    old_first_derivative = first_derivative - alp_derivative * second_meson - alp * second_derivative
    old_second_derivative = second_derivative + alp_derivative * first_meson + alp * first_derivative
    transformed_free = sp.expand(
        (old_first_derivative**2 + old_second_derivative**2
         - first_mass_squared * old_first**2 - second_mass_squared * old_second**2) / 2
    )
    contact = alp * sp.expand(transformed_free).coeff(alp, 1).subs(alp_derivative, 0)
    contact += alp_derivative * sp.expand(transformed_free).coeff(alp_derivative, 1).subs(alp, 0)
    # A vector-source coupling is already linear in a and is unchanged by
    # the coordinate map at this order.  Lorentz products are kept symbolic.
    alp_dot_first, alp_dot_second = sp.symbols("p_a_dot_p_1 p_a_dot_p_2", real=True)
    first_momentum_squared, second_momentum_squared = sp.symbols("p_1_sq p_2_sq", real=True)
    contact_vertex = sp.expand(
        -sp.diff(contact, alp_derivative, first_meson, second_derivative) * alp_dot_second
        -sp.diff(contact, alp_derivative, second_meson, first_derivative) * alp_dot_first
        +sp.diff(contact, alp, first_meson, second_meson)
    )
    momentum_conservation = {alp_dot_first: -first_momentum_squared - sp.Symbol("p_1_dot_p_2"),
                             alp_dot_second: -second_momentum_squared - sp.Symbol("p_1_dot_p_2")}
    # Both scalar products with p_a=-p_1-p_2 have the same cross term.
    contact_with_conservation = sp.expand(contact_vertex.subs(momentum_conservation))
    on_shell = {first_momentum_squared: first_mass_squared, second_momentum_squared: second_mass_squared}
    physical_interaction = coupling * alp_derivative * (
        first_meson * second_derivative - second_meson * first_derivative
    )
    physical_vertex = sp.expand(
        -sp.diff(physical_interaction, alp_derivative, first_meson, second_derivative) * alp_dot_second
        -sp.diff(physical_interaction, alp_derivative, second_meson, first_derivative) * alp_dot_first
    )
    physical_on_shell = sp.expand(physical_vertex.subs(momentum_conservation).subs(on_shell))
    assert sp.simplify(contact_with_conservation.subs(on_shell)) == 0
    assert sp.simplify(physical_on_shell - coupling * (second_mass_squared - first_mass_squared)) == 0
    return variation, contact, contact_with_conservation, physical_on_shell


def main():
    print("Natural units; x=a/f_a; tree level; first order in x; Hermitian U(2) generators.")
    variation, contact, contact_vertex, physical_on_shell = derive()
    show("Vector change of the meson matrix", variation)
    show("Contact terms from the quadratic Lagrangian", contact)
    show("Contact vertex with momentum conservation", contact_vertex)
    show("On-shell vector-source vertex", physical_on_shell)


if __name__ == "__main__":
    main()
