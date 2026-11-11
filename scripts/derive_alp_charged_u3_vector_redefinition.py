#!/usr/bin/env python3
"""Check a charge-conserving U(3) a K pi vertex after a vector rotation.

Run: .venv/bin/python scripts/derive_alp_charged_u3_vector_redefinition.py
Natural units; tree level and LO chiral order. The mixing overlap is exact.
"""
import sympy as sp
from derive_alp_flavour_violating_two_pseudoscalars import (
    B_0, E, cS_ds, cV_ds, extracted_channel, f_a, m_d, m_s, m_u, masses,
)


def derive():
    """Extract the charged rotation, contact, and invariant on-shell vertex."""
    pion, pion_conjugate, kaon, kaon_conjugate = sp.symbols('pi_plus pi_minus K_plus K_minus')
    alp, alp_derivative = sp.symbols('a da')
    pion_derivative, pion_conjugate_derivative, kaon_derivative, kaon_conjugate_derivative = sp.symbols(
        'dpi_plus dpi_minus dK_plus dK_minus')
    rotation = sp.Symbol('v', real=True)
    pion_mass_squared, kaon_mass_squared = sp.symbols('m_pi_sq m_K_sq', real=True)
    pion_momentum_squared, kaon_momentum_squared, pair_dot = sp.symbols('p_pi_sq p_K_sq p_pi_dot_p_K', real=True)
    alp_overlap = sp.Symbol('w_a', real=True, nonzero=True)
    meson_matrix = sp.sqrt(2)*(pion*E(0,1)+pion_conjugate*E(1,0)+kaon*E(0,2)+kaon_conjugate*E(2,0))
    generator = (E(1,2)+E(2,1))/2
    charge = sp.diag(sp.Rational(2,3), -sp.Rational(1,3), -sp.Rational(1,3))
    variation = sp.expand(sp.I*(generator*meson_matrix-meson_matrix*generator))
    assert charge*generator-generator*charge == sp.zeros(3)
    assert sp.simplify(variation[0,1]/sp.sqrt(2)+sp.I*kaon/2) == 0
    assert sp.simplify(variation[2,0]/sp.sqrt(2)-sp.I*pion_conjugate/2) == 0

    # Invert Phi'=Phi + v a i[T,Phi] and differentiate by the product rule.
    old_pion = pion + sp.I*rotation*alp*kaon/2
    old_pion_conjugate = pion_conjugate-sp.I*rotation*alp*kaon_conjugate/2
    old_kaon = kaon+sp.I*rotation*alp*pion/2
    old_kaon_conjugate = kaon_conjugate-sp.I*rotation*alp*pion_conjugate/2
    variables = (alp,pion,pion_conjugate,kaon,kaon_conjugate)
    derivatives = (alp_derivative,pion_derivative,pion_conjugate_derivative,kaon_derivative,kaon_conjugate_derivative)
    def derivative(expression):
        return sum(sp.diff(expression, variable)*variation for variable,variation in zip(variables,derivatives))
    lagrangian = (derivative(old_pion)*derivative(old_pion_conjugate)
                  +derivative(old_kaon)*derivative(old_kaon_conjugate)
                  -pion_mass_squared*old_pion*old_pion_conjugate
                  -kaon_mass_squared*old_kaon*old_kaon_conjugate)
    contact = sp.expand(lagrangian).coeff(rotation,1)
    # all incoming: k_a=-k_pi-k_K, and each derivative contributes -i k.
    alp_dot_pion = -pion_momentum_squared-pair_dot
    alp_dot_kaon = -kaon_momentum_squared-pair_dot
    vertex = sp.expand(
        sp.diff(contact,alp,pion,kaon_conjugate)
        -alp_dot_pion*sp.diff(contact,alp_derivative,pion_derivative,kaon_conjugate)
        -alp_dot_kaon*sp.diff(contact,alp_derivative,pion,kaon_conjugate_derivative)
        -pair_dot*sp.diff(contact,alp,pion_derivative,kaon_conjugate_derivative))
    vertex_on_shell = sp.simplify(vertex.subs({pion_momentum_squared:pion_mass_squared,
                                               kaon_momentum_squared:kaon_mass_squared}))
    assert vertex != 0 and vertex_on_shell == 0
    # Exact ALP overlap of an arbitrarily mixed incoming neutral pole: pure
    # meson cubic interactions vanish in the even-parity LO chiral theory.
    physical_amplitude = sp.simplify(alp_overlap*extracted_channel('K+','pi-').subs({
        masses['K+']**2:B_0*(m_u+m_s), masses['pi-']**2:B_0*(m_u+m_d)}))
    invariant_amplitude = -alp_overlap*B_0*(m_s-m_d)*(cV_ds-cS_ds)/f_a
    assert sp.simplify(physical_amplitude-invariant_amplitude) == 0
    transformed_amplitude = physical_amplitude + sp.I*rotation*alp_overlap*vertex_on_shell
    assert sp.simplify(transformed_amplitude-invariant_amplitude) == 0
    return vertex, physical_amplitude, vertex_on_shell


def main():
    contact, amplitude, on_shell_contact = derive()
    print('Charge-conserving ds rotation; tree level; exact incoming ALP pole overlap.')
    print('Off-shell vertex contact:',sp.latex(contact))
    print('On-shell contact:',sp.latex(on_shell_contact))
    print('Physical K+ pi- amplitude:',sp.latex(amplitude))


if __name__ == '__main__':
    main()
