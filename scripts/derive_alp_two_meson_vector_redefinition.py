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


def derive_mixed_pole_relation():
    """Check the linear pole shift and the nonlinear off-shell contact."""
    invariant_mass_squared, first_squared, second_squared = sp.symbols(
        "s p_1_sq p_2_sq", real=True
    )
    first_mass_squared, second_mass_squared = sp.symbols("m_1_sq m_2_sq", real=True)
    kinetic_alp, kinetic_meson, kinetic_mixing = sp.symbols("K_aa K_nn K_an", real=True)
    mass_alp, mass_meson, mass_mixing = sp.symbols("M_aa M_nn M_an", real=True)
    pole_meson_component, rotation = sp.symbols("w_n v", real=True)
    alp_vertex, meson_vertex = sp.symbols("V_a12 V_n12", real=True)

    inverse_propagator = sp.Matrix([
        [invariant_mass_squared * kinetic_alp - mass_alp,
         invariant_mass_squared * kinetic_mixing - mass_mixing],
        [invariant_mass_squared * kinetic_mixing - mass_mixing,
         invariant_mass_squared * kinetic_meson - mass_meson],
    ])
    change_of_fields = sp.Matrix([[1, 0], [pole_meson_component, 1]])
    transformed_inverse_propagator = change_of_fields.T * inverse_propagator * change_of_fields
    off_diagonal = sp.expand(transformed_inverse_propagator[1, 0])
    assert sp.simplify(off_diagonal - (inverse_propagator[1, 0]
                       + pole_meson_component * inverse_propagator[1, 1])) == 0

    first_inverse = first_squared - first_mass_squared
    second_inverse = second_squared - second_mass_squared
    pole_vector_piece = alp_vertex + pole_meson_component * meson_vertex
    contact = rotation * (second_inverse - first_inverse)
    original_contact = derive()[2]
    assert sp.simplify(contact - rotation * original_contact) == 0
    transformed_vertex = pole_vector_piece + contact
    assert sp.simplify((transformed_vertex - pole_vector_piece) - contact) == 0
    assert sp.simplify(contact.subs({first_squared: first_mass_squared,
                                     second_squared: second_mass_squared})) == 0
    return off_diagonal, pole_vector_piece, contact, transformed_vertex


def main():
    print("Natural units; x=a/f_a; tree level; first order in x; Hermitian U(2) generators.")
    variation, contact, contact_vertex, physical_on_shell = derive()
    show("Vector change of the meson matrix", variation)
    show("Contact terms from the quadratic Lagrangian", contact)
    show("Contact vertex with momentum conservation", contact_vertex)
    show("On-shell vector-source vertex", physical_on_shell)
    off_diagonal, pole_piece, contact, transformed_vertex = derive_mixed_pole_relation()
    show("Transformed inverse-propagator mixing", off_diagonal)
    show("Coherent pole-vector contribution", pole_piece)
    show("Formal off-shell vector-rotation contact", contact)
    show("Transformed three-point vertex", transformed_vertex)


if __name__ == "__main__":
    main()
