"""Derive ALP + virtual photon -> pi+ pi- pi0 through a spin-one V_i pole.

Run: .venv/bin/python scripts/derive_alp_three_pion_vector_pole.py
This executable derivation has no test mode.

GeV, metric (+---), canonical epsilon^{0123}=+1. The free tensor Hessian
is the one used in derive_alp_vector_two_pseudoscalar.py. The decay source
is assumed to be a pure k wedge J continuation of the measured on-shell
vector decay. This assumption fixes only a pole representative, not the
regular off-shell amplitude or independent VPPP contact terms. No width,
Khuri--Treiman solution, fit, or double-spectral density is supplied.
Inputs: memory/alp/lepton-pair/alp-doubly-virtual-tff.md and
memory/rcht/odd-parity/dispersive-three-pion-treatment-plan.md.
"""
from itertools import combinations

import sympy as sp

from derive_alp_doubly_virtual_tff import ordered_numerators
import derive_radiative_tff_amplitude as tensors


def derive_tensor_response():
    """Invert the six-component Hessian before imposing the decay source.

    The independent components are contravariant V^{mu nu}. A source enters
    as V_{mu nu} B^{mu nu}; both antisymmetric orders are included. Source
    coefficients carry their own Fourier phases. The stationary response
    therefore has no additional Feynman propagator phase attached here.
    """
    metric = sp.diag(1, -1, -1, -1)
    energy, longitudinal, mass = sp.symbols('k_0 k_z M_i', real=True)
    momentum = sp.Matrix([energy, 0, 0, longitudinal])
    virtuality = energy**2 - longitudinal**2
    pairs = tuple(combinations(range(4), 2))
    components = sp.Matrix(sp.symbols('V_0:6'))
    source_components = sp.Matrix(sp.symbols('B_0:6'))
    field, source = sp.zeros(4), sp.zeros(4)
    for position, (mu, nu) in enumerate(pairs):
        field[mu, nu], field[nu, mu] = components[position], -components[position]
        source[mu, nu], source[nu, mu] = source_components[position], -source_components[position]
    divergence = field.T * metric * momentum
    lagrangian = -(divergence.T * metric * divergence)[0] / 2
    lagrangian += mass**2 * sum(metric[mu, mu] * metric[nu, nu]
                               * field[mu, nu]**2 for mu, nu in pairs) / 2
    linear = 2 * sum(metric[mu, mu] * metric[nu, nu]
                     * field[mu, nu] * source[mu, nu] for mu, nu in pairs)
    hessian = sp.hessian(lagrangian, components)
    inverse_response = (-hessian.inv() * sp.Matrix([
        sp.diff(linear, component) for component in components])).applyfunc(sp.factor)
    contraction = source.T * metric * momentum
    spin_one_source = (momentum * contraction.T - contraction * momentum.T) / virtuality
    response_formula = (2 * spin_one_source / (virtuality - mass**2)
                        - 2 * (source - spin_one_source) / mass**2)
    expected = sp.Matrix([response_formula[mu, nu] for mu, nu in pairs])
    residual = (inverse_response - expected).applyfunc(sp.factor)
    current = sp.Matrix(sp.symbols('J_0:4'))
    decay_source = momentum * current.T - current * momentum.T
    replacement = dict(zip(source_components, [decay_source[mu, nu] for mu, nu in pairs]))
    projected_decay = spin_one_source.subs(replacement).applyfunc(sp.factor)
    decay_response = inverse_response.subs(replacement).applyfunc(sp.factor)
    local_decay = (source - spin_one_source).subs(replacement).applyfunc(sp.factor)
    return dict(hessian=hessian, response=inverse_response, residual=residual,
                decay_response=decay_response, decay_source=decay_source,
                projected_decay=projected_decay, local_decay=local_decay,
                momentum=momentum, mass=mass, virtuality=virtuality, pairs=pairs)


def derive_production_numerator():
    """Contract the retained VJP/VVP tensors with independent trace sources."""
    generator = sp.MatrixSymbol('T_i', 3, 3)
    other_generator = sp.MatrixSymbol('T_j', 3, 3)
    charge = sp.MatrixSymbol('Q', 3, 3)
    axial = sp.MatrixSymbol('U_a', 3, 3)
    mass_source = sp.MatrixSymbol('Y_a', 3, 3)
    direct, _ = ordered_numerators(sp.Trace(generator * charge * axial),
                                  sp.Trace(generator * charge * mass_source))
    _, pair = ordered_numerators(sp.Trace(generator * other_generator * axial),
                                sp.Trace(generator * other_generator * mass_source))
    alp_mass = sp.Symbol('m_a', positive=True)
    direct = direct.subs(tensors.m, alp_mass)
    pair = pair.subs(tensors.m, alp_mass)
    electric_charge = sp.Symbol('e', real=True)
    response = sp.Function('W_j')(tensors.s2)
    return direct, pair, sp.expand(-2 * electric_charge * direct
                                  + sp.sqrt(2) * electric_charge * response * pair)


def derive_three_pion_transversality():
    """Build the epsilon current and derive k.J=0 with k=p+ + p- + p0."""
    metric = sp.diag(1, -1, -1, -1)
    positive, negative, neutral = [sp.Matrix(sp.symbols(name + '_0:4'))
                                   for name in ('p_plus', 'p_minus', 'p_zero')]
    # epsilon^{mu nu alpha beta} contracts three lower-index momenta.
    lowered = [metric * vector for vector in (positive, negative, neutral)]
    current = sp.Matrix([sum(sp.LeviCivita(mu, nu, alpha, beta)
                            * lowered[0][nu] * lowered[1][alpha] * lowered[2][beta]
                            for nu in range(4) for alpha in range(4) for beta in range(4))
                         for mu in range(4)])
    ward = sp.expand(((positive + negative + neutral).T * metric * current)[0])
    return current, ward


def derive_pole_factorization():
    """Contract the production tensor with the physical spin-one projector.

    Define M_prod=C_i epsilon(n*,epsilon,k,q) and M_dec=i n.J F_i.
    The propagator is i(-g+k k/M_i^2)/(s-M_i^2+i0). Then multiplying
    i M_prod, propagator, and i M_dec determines M_pole below. The map
    C_i=zeta_i N_i(M_i^2,y)/M_i retains zeta_i as a convention input:
    a common relative LSZ/Fourier phase is not fixed by the separately
    recorded on-shell decay amplitude and wedge-normalized numerator.
    """
    metric = sp.diag(1, -1, -1, -1)
    energy, longitudinal, mass = sp.symbols('k_0 k_z M_i', real=True)
    momentum = sp.Matrix([energy, 0, 0, longitudinal])
    photon = sp.Matrix(sp.symbols('q_0:4', real=True))
    current = sp.Matrix(sp.symbols('J_0:4'))
    lowered_momentum, lowered_photon = metric * momentum, metric * photon
    production_tensor = sp.Matrix(4, 4, lambda alpha, mu: sum(
        sp.LeviCivita(alpha, mu, rho, sigma)
        * lowered_momentum[rho] * lowered_photon[sigma]
        for rho in range(4) for sigma in range(4)))
    projector = -metric + momentum * momentum.T / mass**2
    contraction = (production_tensor.T * metric * projector * metric * current).applyfunc(sp.factor)
    transverse_result = -production_tensor.T * metric * current
    virtuality = sp.Symbol('s', real=True)
    coupling = sp.Function('C_i')(tensors.s2)
    decay_function = sp.Symbol('F_i')
    # i M_prod * [i Pi/(s-M^2)] * i M_dec, with M_dec=i J F_i.
    i_amplitude = (sp.I * coupling) * (sp.I / (virtuality - mass**2)) * (sp.I * sp.I * decay_function) * contraction
    amplitude = (i_amplitude / sp.I).applyfunc(sp.factor)
    return dict(tensor=production_tensor, projector=projector, contraction=contraction,
                reduction=(contraction-transverse_result).applyfunc(sp.factor),
                photon_ward=(production_tensor * metric * photon).applyfunc(sp.factor),
                vector_ward=(production_tensor.T * metric * momentum).applyfunc(sp.factor),
                amplitude=amplitude)


def show(label, expression):
    """Print every analytical expression with SymPy's LaTeX printer."""
    print(f'{label}:\n{sp.latex(expression)}\n')


def main():
    """Print the derivation, with its pole-only scope stated explicitly."""
    print('Natural units GeV; (+---); epsilon^{0123}=+1; stable poles.')
    print('Pure spin-one decay-source continuation assumed; regular contacts undetermined.')
    result = derive_tensor_response()
    show('Six-component free tensor Hessian', result['hessian'])
    show('Stationary response to an arbitrary antisymmetric source B', result['response'])
    show('Inverse-Hessian minus spin-one/local-projector decomposition', result['residual'])
    show('Pure spin-one three-pion source k wedge J', result['decay_source'])
    show('Local tensor-projector contribution on k wedge J', result['local_decay'])
    show('Full tensor response on k wedge J', result['decay_response'])
    current, ward = derive_three_pion_transversality()
    show('Three-pion epsilon current without its scalar decay function', current)
    show('Three-pion current transversality', ward)
    direct, pair, production = derive_production_numerator()
    show('Derived ordered VJP numerator V_i gamma', direct)
    show('Derived ordered VVP numerator P_ij', pair)
    show('Production numerator: displayed VVP summand is summed over j', production)
    virtuality, mass = sp.symbols('s M_i', real=True)
    numerator = sp.Function('N_i')(virtuality, tensors.s2)
    residue = numerator.subs(virtuality, mass**2)
    show('Strict production pole coefficient, multiplying the on-shell decay residue',
         residue / (mass**2 - virtuality))
    show('Regular difference of the polynomial continuation from that pole',
         (numerator - residue) / (mass**2 - virtuality))
    print('The decay residue is defined by M(V_i -> 3pi)=i epsilon(n,p+,p-,p0) F_i.')
    print('Its normalization and phase must be supplied in the same S-matrix convention.')
    print('The on-shell decay amplitude does not determine its off-shell source or regular terms.')
    pole = derive_pole_factorization()
    show('Production epsilon tensor with open vector and photon indices', pole['tensor'])
    show('Spin-one polarization sum', pole['projector'])
    show('Production tensor times photon momentum', pole['photon_ward'])
    show('Production tensor times vector momentum', pole['vector_ward'])
    show('Projector contraction minus the transverse reduction', pole['reduction'])
    show('ALP + virtual photon -> 3pi pole amplitude with M_prod=C_i epsilon(n*,epsilon,k,q), M_dec=i n.J F_i', pole['amplitude'])
    phase = sp.Symbol('zeta_i')
    coupling = sp.Function('C_i')(tensors.s2)
    show('Production coupling with convention input zeta_i',
         sp.Eq(coupling, phase * residue / mass))
    print('J^mu=epsilon^{mu nu alpha beta} p+_nu p-_alpha p0_beta; F_i is the on-shell decay function.')
    print('C_i includes production crossing and relative LSZ phase; zeta_i is not fixed here.')
    print('Only strict s=M_i^2 residue is claimed. No regular off-shell terms are inferred.')


if __name__ == '__main__':
    main()
