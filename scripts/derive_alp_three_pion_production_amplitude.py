"""Assemble the retained a + gamma* -> pi+ pi- pi0 production amplitude.

Run .venv/bin/python scripts/derive_alp_three_pion_production_amplitude.py.
Inputs are the indexed LO Born, vector pole factorization and pole phase
records under memory/alp/lepton-pair/. Natural units (GeV), metric (+---),
epsilon^{0123}=+1, S=1+iT, p+q=k=p_++p_-+p_0; q is incoming.
No finite width, fitted decay normalization, rescattering approximation or
regular off-shell completion is chosen. The photon coupling e is retained.
The full-charge Born term is NOT a derived three-pion I=0 projection.
"""

import sympy as sp

import derive_alp_virtual_photon_three_pion_born as born
import derive_alp_three_pion_vector_pole as pole
import derive_radiative_tff_amplitude as tensors


def show(label, expression):
    """Print an analytical expression using the standard SymPy renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    """Derive source coefficients, transverse structures and retained sum.

    Explicit Lorentz components implement standard epsilon contractions;
    SymPy LeviCivita and matrix multiplication perform the tensor algebra.
    Algebraic Ward identities are evaluated away from pole singularities;
    the causal prescriptions are restored in the physical expressions.
    """
    source = born.derive()
    metric = sp.diag(1, -1, -1, -1)
    photon = sp.Matrix(sp.symbols('q_0:4', real=True))
    positive = sp.Matrix(sp.symbols('p_plus_0:4', real=True))
    negative = sp.Matrix(sp.symbols('p_minus_0:4', real=True))
    neutral = sp.Matrix(sp.symbols('p_zero_0:4', real=True))
    momentum = positive + negative + neutral
    photon_squared = (photon.T * metric * photon)[0]
    photon_positive = (photon.T * metric * positive)[0]
    photon_negative = (photon.T * metric * negative)[0]
    positive_denominator = photon_squared - 2 * photon_positive
    negative_denominator = photon_squared - 2 * photon_negative
    born_transverse_numerator = (
        negative_denominator * (2 * positive - photon)
        - positive_denominator * (2 * negative - photon))
    born_ward = sp.expand((photon.T * metric * born_transverse_numerator)[0])
    electric_charge = sp.Symbol('e', real=True)
    born_current = (electric_charge * source['source_amplitude']
                    * born_transverse_numerator
                    / (positive_denominator * negative_denominator))

    lowered = [metric * vector for vector in (positive, negative, neutral)]
    decay_current = sp.Matrix([sum(
        sp.LeviCivita(mu, nu, alpha, beta)
        * lowered[0][nu] * lowered[1][alpha] * lowered[2][beta]
        for nu in range(4) for alpha in range(4) for beta in range(4))
        for mu in range(4)])
    lowered_momentum, lowered_photon = metric * momentum, metric * photon
    production_tensor = sp.Matrix(4, 4, lambda alpha, mu: sum(
        sp.LeviCivita(alpha, mu, rho, sigma)
        * lowered_momentum[rho] * lowered_photon[sigma]
        for rho in range(4) for sigma in range(4)))
    photon_tensor_ward = (production_tensor * metric * photon).applyfunc(sp.expand)
    vector_tensor_ward = (production_tensor.T * metric * momentum).applyfunc(sp.expand)
    decay_ward = sp.expand((momentum.T * metric * decay_current)[0])
    # The longitudinal polarization sum vanishes by the derived tensor Ward
    # identity, leaving E^T g Pi g J = -E^T g J.
    pole_transverse_vector = -production_tensor.T * metric * decay_current

    direct, pair, _ = pole.derive_production_numerator()
    virtuality, photon_virtuality = sp.symbols('s y', real=True)
    ordered_substitutions = {tensors.s1: virtuality, tensors.s2: photon_virtuality}
    direct = direct.subs(ordered_substitutions)
    pair = pair.subs(ordered_substitutions)
    other_generator = sp.MatrixSymbol('T_j', 3, 3)
    direct_generator = sp.MatrixSymbol('T_i', 3, 3)
    retained_pole_coefficient = 0
    production_numerators = {}
    pole_coefficients = {}
    for name in ('omega', 'phi'):
        generator = sp.MatrixSymbol('T_' + name, 3, 3)
        mass = sp.Symbol('M_' + name, positive=True)
        direct_source = direct.xreplace({direct_generator: generator})
        converted_source = sum(
            sp.Function('W_' + intermediate)(photon_virtuality)
            * pair.xreplace({direct_generator: generator,
                            other_generator: sp.MatrixSymbol('T_' + intermediate, 3, 3)})
            for intermediate in ('rho', 'omega', 'phi'))
        numerator = -2 * electric_charge * direct_source + sp.sqrt(2) * electric_charge * converted_source
        production_numerators[name] = numerator
        residue = numerator.subs(virtuality, mass**2)
        decay_function = sp.Function('F_' + name + '_on')(sp.Symbol('t'), sp.Symbol('u'))
        # zeta_i=-1 from the existing crossing/LSZ calculation. Any physical
        # decay-normalization phase remains in the complex F_i_on.
        coefficient = -sp.I * decay_function * residue / (mass * (mass**2 - virtuality))
        pole_coefficients[name] = coefficient
        retained_pole_coefficient += coefficient
    # Derive the pole Ward contraction in its factored tensor order rather
    # than expanding thousands of cancelling component monomials.
    pole_ward = -(photon_tensor_ward.T * metric * decay_current)[0]
    total_ward = sp.simplify(
        electric_charge * source['source_amplitude'] * born_ward
        / (positive_denominator * negative_denominator)
        + retained_pole_coefficient * pole_ward)
    return locals()


def main():
    """Display reproducible intermediate expressions and scope limits."""
    result = derive()
    print('GeV; (+---); epsilon^0123=+1; S=1+iT; incoming photon; arbitrary q^2.')
    print('U_a=2 c_q^A/f_a; Y_a=-2 {c_q^P,chi}/f_a; independent no-mixing sources.')
    print('W_j(y)=R_j/(M_j^2-y-i0); all three neutral photon responses retained.')
    for name in ('positive_derivative_weight', 'negative_derivative_weight',
                 'neutral_derivative_weight', 'pseudoscalar_mass_weight',
                 'source_amplitude', 'contact_current_coefficient'):
        show(name, result['source'][name])
    show('Born transverse numerator', result['born_transverse_numerator'])
    show('Photon contraction of Born transverse numerator', result['born_ward'])
    show('Full-charge Born current', result['born_current'])
    show('Three-pion decay-current transversality', result['decay_ward'])
    show('Three-pion decay current J', result['decay_current'])
    show('Production tensor E with vector and photon indices', result['production_tensor'])
    show('Photon contraction of production tensor', result['photon_tensor_ward'])
    show('Vector contraction of production tensor', result['vector_tensor_ward'])
    for name in ('omega', 'phi'):
        show('Source-resolved N_' + name + '(s,y)', result['production_numerators'][name])
        show('Strict pole coefficient multiplying -E^T g J: ' + name,
             result['pole_coefficients'][name])
    show('Total retained photon Ward contraction', result['total_ward'])
    print('Retained amplitude = displayed Born current + sum_i coefficient_i (-E^T g J).')
    print('Restore +i0 in charged-pion denominators and -i0 in M_i^2-s.')
    print('Strict pole residues only: normalized complex F_i_on and its physical phase remain inputs.')
    print('Regular source-dependent terms and left-hand cuts beyond LO Born remain undetermined.')
    print('The three-pion I=0, J^PC=1^{--} projection of the full-charge Born term remains open.')
    print('Its value cannot be inferred from vanishing minimal coupling to the light-quark isoscalar current.')


if __name__ == '__main__':
    main()
