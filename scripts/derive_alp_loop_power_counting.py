#!/usr/bin/env python3
"""One-loop power counting for ALP lines in the LO ALP-meson theory at finite f_pi/f_a.

Conventions (Alda:2026nrd, Kaiser:1998ds):
  * Units: GeV.  x = a/f_a,  U = u^2 = exp(i Phi/f_pi),  Phi = sum_i lambda_i P_i,
    tr(lambda_i lambda_j) = 2 delta_ij.
  * Dressed frame (input hypothesis, see memory/alp/loop-corrections/finite-mixing/power-counting-and-validity/loop-counting/alp-loop-power-counting.md):
    the ALP currents are removed by an x-dependent chiral rotation, so the LO
    Lagrangian is
        L = f^2/4 <u_mu u^mu + chi_+> - f^2 m_0^2/(4 N_f) X_TC^2 + (da)^2/2 - m_a0^2 a^2/2,
    with chi = 2 B_0 hat M(x), u_mu free of sources and X_TC = -tr Phi/f_pi - 2 hat Theta(x).
  * Linear-spurion convention with constant diagonal axial current:
        M(x) = m_q + i x {c_q^P, m_q},  C_A = c_q^A (diagonal),  C_V = 0,
    so W_R = exp(i x c_A), W_L = exp(-i x c_A) and
        hat M(x) = exp(-2 i x c_A) (1 + 2 i x c_P) m_q,
        hat M_n = (-2i)^n m_q c_A^(n-1) (c_A - n c_P)   (n >= 1),
        hat Theta(x) = c_GA x,  c_GA = c_G + tr c_A.
  * Loop function (dimensional regularization, quadratic divergences vanish since
    int d^dk 1 = 0):
        A_0(m^2) = -mu^(4-d) int d^dk/(2 pi)^d  i/(k^2 - m^2 + i0)
                 = m^2/(16 pi^2) [1 - log(m^2/mu^2)]   after MSbar subtraction of
        m^2/(16 pi^2) [2/(4-d) - gamma_E + log 4 pi].
    Then <phi^2> = -A_0, and with V^L_{IJKL} the fourth derivative of the Lagrangian
    (amputated vertex without the factor i)
        delta M^2_IJ = (1/2) sum_n V^L_{IJKL} v_K^(n) v_L^(n) A_0(m_n^2).
    For momentum-dependent vertices: legs (p, -p, k, -k), k.k -> m_n^2, p.k -> 0,
    and Sigma_IJ(s) = -(1/2) sum V^L(s, m_n^2) v v A_0, Gamma(s) = s K - M^2 + Sigma(s).
  * Unmixed approximation: ALP and meson propagators taken diagonal in the
    (P_u, P_d, P_s, a) flavour basis, with masses 2 B_0 m_q (m_0 neglected in
    internal meson lines) and M^2_aa for the ALP.  It is labelled wherever used.

Minimal invocation:
    .venv/bin/python scripts/derive_alp_loop_power_counting.py
Tests:
    .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_loop_power_counting.py
"""
from __future__ import annotations

import itertools

import numpy as np
import sympy as sp

from derive_alp_lo_joint_vertices import vertex
from solve_alp_meson_spectrum import solve_spectrum_numerical

FLAVOURS = ('u', 'd', 's')
NEUTRAL_FIELD_LABELS = ('P_u', 'P_d', 'P_s', 'a')
CHARGED_PAIRS = ((0, 1), (0, 2), (1, 2))
NUMBER_OF_FLAVOURS = 3

# Illustrative inputs (GeV), isospin limit.  m_0^2 from the LO U(3) trace relation
# m_pi^2 + m_eta^2 + m_eta'^2 = sum_q 2 B_0 m_q + m_0^2 (assumption, illustrative).
ILLUSTRATIVE = {
    'pion_mass': 0.13957,
    'kaon_mass': 0.49368,
    'eta_mass': 0.547862,
    'eta_prime_mass': 0.95778,
    'pion_decay_constant': 0.09207,
    'renormalization_scale': 1.0,
}


# ----------------------------------------------------------------------------- L1

def dressed_mass_derivative_closed_form(order, quark_masses, axial_coefficients,
                                        pseudoscalar_coefficients):
    """hat M_n = d^n hat M/dx^n at x=0 for diagonal inputs (sympy diagonal matrix)."""
    entries = []
    for mass, axial, pseudoscalar in zip(quark_masses, axial_coefficients,
                                         pseudoscalar_coefficients):
        if order == 0:
            entries.append(mass)
        else:
            entries.append((-2*sp.I)**order*mass*axial**(order-1)*(axial-order*pseudoscalar))
    return sp.diag(*entries)


def dressed_mass_function(x, quark_masses, axial_coefficients, pseudoscalar_coefficients):
    """hat M(x) = W_R^dag M(x) W_L with W_R = exp(i x c_A), W_L = exp(-i x c_A)."""
    entries = [sp.exp(-sp.I*x*axial)*(mass + 2*sp.I*x*pseudoscalar*mass)*sp.exp(-sp.I*x*axial)
               for mass, axial, pseudoscalar in zip(quark_masses, axial_coefficients,
                                                    pseudoscalar_coefficients)]
    return sp.diag(*entries)


def mass_term_vertex_closed_form(meson_generators, number_of_alps, dressed_mass_derivative,
                                 pion_decay_constant, alp_decay_constant, B_0):
    """Amputated vertex (no factor i) from f^2/4 <chi_+> with n_a ALPs and n_P mesons.

    V = B_0 f^(2-n_P) f_a^(-n_a) 2^(-1-n_P) sum_sigma sum_k
        [(-i)^n_P <lam..lam hatM_{n_a} lam..lam> + i^n_P <lam..lam hatM_{n_a}^dag lam..lam>]/(k!(n_P-k)!)
    with lambda_sigma(1..k) to the left of hat M.  Momentum independent.
    """
    number_of_mesons = len(meson_generators)
    size = dressed_mass_derivative.rows
    hermitian_conjugate = dressed_mass_derivative.conjugate().T
    total = sp.S.Zero
    for ordering in itertools.permutations(range(number_of_mesons)):
        for split in range(number_of_mesons+1):
            left = sp.eye(size)
            for index in ordering[:split]:
                left = left*meson_generators[index]
            right = sp.eye(size)
            for index in ordering[split:]:
                right = right*meson_generators[index]
            weight = sp.factorial(split)*sp.factorial(number_of_mesons-split)
            total += ((-sp.I)**number_of_mesons*(left*dressed_mass_derivative*right).trace()
                      + sp.I**number_of_mesons*(left*hermitian_conjugate*right).trace())/weight
    return sp.expand(B_0*pion_decay_constant**(2-number_of_mesons)
                     * alp_decay_constant**(-number_of_alps)*total/2**(1+number_of_mesons))


def dressed_frame_vertex(meson_generators, alp_components, momentum_symbols,
                         dressed_mass_derivatives, B_0, singlet_coefficient,
                         pion_decay_constant, alp_decay_constant,
                         anomaly_mass_squared, bare_alp_mass_squared):
    """Vertex of the dressed-frame LO Lagrangian using the general joint-vertex code.

    Currents are zero, chi_n = 2 B_0 hat M_n, hat Theta = c_GA x (linear).
    """
    order = len(meson_generators)
    size = dressed_mass_derivatives[0].rows
    zero = sp.zeros(size)
    chi = [2*B_0*matrix for matrix in dressed_mass_derivatives]
    chi = chi + [zero]*(order+1-len(chi))
    currents = [zero]*(order+1)
    theta = [0, singlet_coefficient] + [0]*order
    return vertex(meson_generators, alp_components, momentum_symbols, chi, currents, currents,
                  theta, pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                  bare_alp_mass_squared)


def flavour_basis_generators():
    """Neutral sqrt(2) E_qq and charged pairs (sqrt(2) E_qq', sqrt(2) E_q'q)."""
    def unit(row, column):
        matrix = sp.zeros(3)
        matrix[row, column] = sp.sqrt(2)
        return matrix
    neutral = [unit(q, q) for q in range(3)]
    charged = [(unit(q, r), unit(r, q)) for q, r in CHARGED_PAIRS]
    return neutral, charged


def neutral_dressed_lagrangian(neutral_fields, alp, quark_masses, axial_coefficients,
                               pseudoscalar_coefficients, B_0, singlet_coefficient,
                               pion_decay_constant, alp_decay_constant,
                               anomaly_mass_squared, bare_alp_mass_squared):
    """Exact dressed LO potential part of L for diagonal fields (P_u, P_d, P_s, a).

    L = f^2 B_0 sum_q m_q [cos th_q + 2 c_Pq x sin th_q] - f^2 m_0^2/(4N)(sqrt2 sum P/f + 2 c_GA x)^2
        - m_a0^2 a^2/2,   th_q = sqrt2 P_q/f + 2 c_Aq x.
    The kinetic term of commuting fields is exactly quadratic and is omitted.
    """
    x = alp/alp_decay_constant
    angle = [sp.sqrt(2)*field/pion_decay_constant + 2*axial*x
             for field, axial in zip(neutral_fields, axial_coefficients)]
    mass_part = pion_decay_constant**2*B_0*sum(
        mass*(sp.cos(theta) + 2*pseudoscalar*x*sp.sin(theta))
        for mass, theta, pseudoscalar in zip(quark_masses, angle, pseudoscalar_coefficients))
    topological = (sp.sqrt(2)*sum(neutral_fields)/pion_decay_constant + 2*singlet_coefficient*x)
    return (mass_part - pion_decay_constant**2*anomaly_mass_squared/(4*NUMBER_OF_FLAVOURS)
            * topological**2 - bare_alp_mass_squared*alp**2/2)


def neutral_derivative_tensor(lagrangian, fields, order):
    """All order-th derivatives of the Lagrangian at zero fields (sympy Array)."""
    zero = {field: 0 for field in fields}
    shape = (len(fields),)*order
    values = {}
    for indices in itertools.product(range(len(fields)), repeat=order):
        key = tuple(sorted(indices))
        if key not in values:
            values[key] = sp.simplify(sp.diff(lagrangian, *[fields[i] for i in key]).subs(zero))
    return sp.MutableDenseNDimArray([values[tuple(sorted(i))] for i in
                                     itertools.product(range(len(fields)), repeat=order)], shape)


# ----------------------------------------------------------------------------- L2

def loop_function_A0(squared_mass, squared_scale):
    """Renormalized A_0(m^2) = m^2/(16 pi^2) [1 - log(m^2/mu^2)] (sympy or float)."""
    if isinstance(squared_mass, (int, float)):
        if squared_mass < 0:
            raise ValueError('negative squared mass in a loop: tachyonic tree spectrum')
        if squared_mass == 0:
            return 0.0
        return squared_mass/(16*np.pi**2)*(1-np.log(squared_mass/squared_scale))
    return squared_mass/(16*sp.pi**2)*(1-sp.log(squared_mass/squared_scale))


def tadpole_mass_correction(quartic_tensor, pole_squared_masses, pole_vectors, squared_scale,
                            pole_indices=None):
    """delta M^2_IJ = (1/2) sum_n V_IJKL v_K^n v_L^n A_0(m_n^2) for a numpy quartic tensor."""
    indices = range(len(pole_squared_masses)) if pole_indices is None else pole_indices
    correction = np.zeros(quartic_tensor.shape[:2])
    for n in indices:
        vector = pole_vectors[:, n]
        correction += 0.5*np.einsum('ijkl,k,l->ij', quartic_tensor, vector, vector) \
            * loop_function_A0(float(pole_squared_masses[n]), squared_scale)
    return correction


def coleman_weinberg_potential(field_values, mass_squared, quartic_potential_tensor, squared_scale):
    """V_1 = sum_n m_n^4(phi) [log(m_n^2/mu^2) - 3/2]/(64 pi^2) with
    M^2(phi)_KL = M^2_KL + (1/2) V^pot_KLIJ phi_I phi_J (no cubic terms)."""
    field_mass = mass_squared + 0.5*np.einsum('klij,i,j->kl', quartic_potential_tensor,
                                              field_values, field_values)
    eigenvalues = np.linalg.eigvalsh(field_mass)
    return float(np.sum(eigenvalues**2*(np.log(eigenvalues/squared_scale)-1.5))/(64*np.pi**2))


def loop_kinematics(expression, momentum_symbols, squared_momentum, internal_squared_mass):
    """Legs (p,-p,k,-k): q1 q2 -> -s, q3 q4 -> -k^2 = -m^2, terms linear in k -> 0."""
    q = momentum_symbols
    polynomial = sp.Poly(sp.expand(expression), *q)
    result = sp.S.Zero
    for exponents, coefficient in polynomial.terms():
        if sum(exponents) == 0:
            result += coefficient
        elif exponents == (1, 1, 0, 0):
            result += -coefficient*squared_momentum
        elif exponents == (0, 0, 1, 1):
            result += -coefficient*internal_squared_mass
        elif sum(exponents) == 2 and max(exponents) == 1:
            continue  # p.k terms integrate to zero
        else:
            raise ValueError(f'unexpected momentum structure {exponents}')
    return sp.expand(result)


def charged_loop_self_energy(external_generators, external_alp_components, quark_mass_values,
                             axial_values, pseudoscalar_values, B_0, singlet_coefficient,
                             pion_decay_constant, alp_decay_constant, anomaly_mass_squared,
                             squared_momentum, squared_scale, symbolic_loop=False):
    """Sigma_IJ(s) from charged-meson tadpoles (momentum-dependent kinetic vertices included).

    Charged pair (q,q') propagates with m^2 = B_0 (m_q + m_q'); a complex field contributes
    -V_{IJ phi phi*}(s, m^2) A_0(m^2).
    """
    _, charged = flavour_basis_generators()
    momenta = sp.symbols('q1:5', real=True)
    derivatives = [dressed_mass_derivative_closed_form(n, quark_mass_values, axial_values,
                                                       pseudoscalar_values) for n in range(5)]
    size = len(external_generators)
    sigma = sp.zeros(size)
    for (pair, (generator, conjugate)) in zip(CHARGED_PAIRS, charged):
        squared_mass = B_0*(quark_mass_values[pair[0]]+quark_mass_values[pair[1]])
        loop = (sp.Symbol(f'A0_{FLAVOURS[pair[0]]}{FLAVOURS[pair[1]]}') if symbolic_loop
                else loop_function_A0(squared_mass, squared_scale))
        for i in range(size):
            for j in range(i, size):
                expression = dressed_frame_vertex(
                    [external_generators[i], external_generators[j], generator, conjugate],
                    [external_alp_components[i], external_alp_components[j], 0, 0], momenta,
                    derivatives, B_0, singlet_coefficient, pion_decay_constant,
                    alp_decay_constant, anomaly_mass_squared, 0)
                value = loop_kinematics(expression, momenta, squared_momentum, squared_mass)
                sigma[i, j] = sigma[j, i] = sp.expand(sigma[i, j] - value*loop)
    return sigma


def unmixed_closed_forms():
    """Labelled unmixed-approximation closed forms (sympy) for classes (a), (b), (c)."""
    f, f_a, B_0 = sp.symbols('f_pi f_a B_0', positive=True)
    masses = sp.symbols('m_u m_d m_s', positive=True)
    axial = sp.symbols('c_Au c_Ad c_As', real=True)
    pseudoscalar = sp.symbols('c_Pu c_Pd c_Ps', real=True)
    A0 = sp.Function('A_0')
    alp_mass_squared = sp.Symbol('m_a^2', positive=True)
    second_order = [axial[q]*(axial[q]-2*pseudoscalar[q]) for q in range(3)]
    fourth_order = [axial[q]**3*(axial[q]-4*pseudoscalar[q]) for q in range(3)]
    meson_squared = {(q, r): B_0*(masses[q]+masses[r]) for q in range(3) for r in range(3)}
    tree_alp_quark = 4*B_0*f**2/f_a**2*sum(masses[q]*second_order[q] for q in range(3))
    results = {
        # (a) ALP line in meson entry (q, q'):  ratio to tree B_0 (m_q + m_q')
        'alp_line_in_meson_entry_ratio': {
            (q, r): 2*(masses[q]*second_order[q]+masses[r]*second_order[r])/(masses[q]+masses[r])
            * A0(alp_mass_squared)/f_a**2 for q in range(3) for r in range(q, 3)},
        # (b) meson lines in M^2_aa (flavour basis, m_0 neglected internally)
        'meson_lines_in_alp_entry': 4*B_0/f_a**2*sum(
            masses[q]*second_order[q]*A0(meson_squared[(q, r)]) for q in range(3) for r in range(3)),
        # (c) ALP line in M^2_aa
        'alp_line_in_alp_entry': 8*f**2*B_0/f_a**4*sum(
            masses[q]*fourth_order[q] for q in range(3))*A0(alp_mass_squared),
        'tree_alp_entry_quark_part': tree_alp_quark,
        'symbols': dict(f=f, f_a=f_a, B_0=B_0, masses=masses, axial=axial,
                        pseudoscalar=pseudoscalar, A0=A0, alp_mass_squared=alp_mass_squared),
    }
    return results


def illustrative_neutral_theory(ratio_f_pi_over_f_a, axial=1.0, pseudoscalar=0.0,
                                gluon_coefficient=0.0, bare_alp_mass_squared=0.0, inputs=None):
    """Numerical neutral (K=1, M^2) and quartic tensor in the dressed frame (floats)."""
    inputs = ILLUSTRATIVE if inputs is None else inputs
    f = inputs['pion_decay_constant']
    two_B0_m = [inputs['pion_mass']**2, inputs['pion_mass']**2,
                2*inputs['kaon_mass']**2-inputs['pion_mass']**2]
    anomaly = inputs['eta_mass']**2+inputs['eta_prime_mass']**2-2*inputs['kaon_mass']**2
    f_a = f/ratio_f_pi_over_f_a
    fields = sp.symbols('P_u P_d P_s a', real=True)
    axial_values = [axial]*3
    singlet = gluon_coefficient+sum(axial_values)
    lagrangian = neutral_dressed_lagrangian(
        fields[:3], fields[3], [value/2 for value in two_B0_m], axial_values, [pseudoscalar]*3,
        1.0, singlet, f, f_a, anomaly, bare_alp_mass_squared)
    hessian = np.array(neutral_derivative_tensor(lagrangian, fields, 2).tolist(), dtype=float)
    quartic = np.array(neutral_derivative_tensor(lagrangian, fields, 4).tolist(), dtype=float)
    return np.eye(4), -hessian, quartic, dict(f_a=f_a, singlet=singlet, two_B0_m=two_B0_m,
                                              anomaly=anomaly)


def pole_sum_decomposition(ratio_f_pi_over_f_a, **options):
    """Neutral-line tadpole delta M^2 exact in mixing, split into ALP-like and meson-like poles."""
    scale = ILLUSTRATIVE['renormalization_scale']**2
    kinetic, mass_squared, quartic, info = illustrative_neutral_theory(ratio_f_pi_over_f_a, **options)
    spectrum = solve_spectrum_numerical(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
    alp_like = int(np.argmax(np.abs(spectrum.pole_vectors[3, :])))
    meson_like = [n for n in range(4) if n != alp_like]
    return dict(
        spectrum=spectrum, alp_like_index=alp_like, mass_squared=mass_squared, info=info,
        alp_like=tadpole_mass_correction(quartic, spectrum.pole_squared_masses,
                                         spectrum.pole_vectors, scale, [alp_like]),
        meson_like=tadpole_mass_correction(quartic, spectrum.pole_squared_masses,
                                           spectrum.pole_vectors, scale, meson_like))


# ----------------------------------------------------------------------------- L3

def imaginary_part_bubble_function(squared_momentum, first_squared_mass, second_squared_mass):
    """Im B_0 = pi * (length of x in [0,1] with x m1^2 + (1-x) m2^2 - x(1-x) s < 0)."""
    x = sp.Symbol('x', real=True)
    roots = sp.solve(sp.Eq(x*first_squared_mass+(1-x)*second_squared_mass
                           - x*(1-x)*squared_momentum, 0), x)
    return sp.simplify(sp.pi*sp.Abs(roots[1]-roots[0]))


def bubble_unitarity_check():
    """Im a_0(i->i) from the s-channel bubble = sigma_f S_f |a_0(i->f)|^2 (a_0 = M/(16 pi))."""
    above_threshold, alp_mass_sq, coupling = sp.symbols('w m_a^2 V_PPaa', positive=True)
    s = 4*alp_mass_sq+above_threshold  # s > 4 m_a^2
    identical_factor = sp.Rational(1, 2)
    im_B0 = sp.simplify(imaginary_part_bubble_function(s, alp_mass_sq, alp_mass_sq))
    im_loop_amplitude = identical_factor*coupling**2*im_B0/(16*sp.pi**2)
    phase_space = sp.sqrt(1-4*alp_mass_sq/s)
    lhs = im_loop_amplitude/(16*sp.pi)
    rhs = phase_space*identical_factor*(coupling/(16*sp.pi))**2
    return sp.simplify(lhs-rhs), im_B0


def bubble_to_tree_ratios():
    """|loop/tree| for P_q P_q -> P_q P_q (flavour basis, unmixed) by intermediate state."""
    f, f_a, B_0, m_q, c_AP, kappa, B0_loop = sp.symbols(
        'f_pi f_a B_0 m_q c_AP kappa_q B0_loop', positive=True)
    meson_four = 4*B_0*m_q/f**2
    meson_three_alp = 4*sp.sqrt(2)*B_0*m_q*c_AP/(f*f_a)
    meson_two_alp_two = 8*B_0*m_q*kappa/f_a**2
    loop = B0_loop/(16*sp.pi**2)
    return dict(
        PP=sp.simplify(sp.Rational(1, 2)*meson_four**2*loop/meson_four),
        aP=sp.simplify(meson_three_alp**2*loop/meson_four),
        aa=sp.simplify(sp.Rational(1, 2)*meson_two_alp_two**2*loop/meson_four),
        vertices=dict(PPPP=meson_four, PPPa=meson_three_alp, PPaa=meson_two_alp_two),
        symbols=dict(f=f, f_a=f_a, B_0=B_0, m_q=m_q, c_AP=c_AP, kappa=kappa, B0_loop=B0_loop))


# ----------------------------------------------------------------------------- L4

def regime_map():
    """Closed-form crossing values r_* = f_pi/f_a where loop/tree = 1 (|A_0| -> m^2/16pi^2)."""
    r, f, kappa, sum_kappa_m2, c_GA, m0sq, m_a, m_qq, c_AP, B0_loop = sp.symbols(
        'r f_pi kappa_q Sigma_kappa_m2 c_GA m_0^2 m_a m_qq c_AP B0_loop', positive=True)
    loop = 16*sp.pi**2*f**2
    ratios = {
        'alp_line_meson_entry_quark_mass': 4*kappa*r**4*sum_kappa_m2/loop,
        'alp_line_meson_entry_topological': sp.Rational(4, 3)*kappa*c_GA**2*m0sq*r**4/loop,
        'alp_line_meson_entry_fixed_mass': 2*kappa*r**2*m_a**2/loop,
        'alp_line_alp_entry_any_mass': 4*r**4*sum_kappa_m2/loop,
        'bubble_aP': 4*c_AP**2*m_qq**2*r**2*B0_loop/loop,
        'bubble_aa': 4*kappa**2*m_qq**2*r**4*B0_loop/loop,
    }
    crossings = {name: sp.solve(sp.Eq(value, 1), r)[0] for name, value in ratios.items()}
    return ratios, crossings, dict(r=r, f=f, kappa=kappa, sum_kappa_m2=sum_kappa_m2, c_GA=c_GA,
                                   m0sq=m0sq, m_a=m_a, m_qq=m_qq, c_AP=c_AP, B0_loop=B0_loop)


def illustrative_regime_numbers(axial=1, singlet=3, fixed_alp_mass=1.0):
    """Crossing values for c_A = c (c_P = 0), c_GA = singlet, |B_0| -> 1 (illustrative)."""
    _, crossings, sym = regime_map()
    m_pi, m_K, f = (ILLUSTRATIVE['pion_mass'], ILLUSTRATIVE['kaon_mass'],
                    ILLUSTRATIVE['pion_decay_constant'])
    m_ss_sq = 2*m_K**2-m_pi**2
    anomaly = ILLUSTRATIVE['eta_mass']**2+ILLUSTRATIVE['eta_prime_mass']**2-2*m_K**2
    common = {sym['f']: f, sym['kappa']: axial**2, sym['sum_kappa_m2']: axial**2*(2*m_pi**2+m_ss_sq),
              sym['c_GA']: singlet, sym['m0sq']: anomaly, sym['m_a']: fixed_alp_mass,
              sym['c_AP']: axial, sym['B0_loop']: 1}
    numbers = {}
    for name, value in crossings.items():
        if name.startswith('bubble'):
            for label, mass_sq in (('uu', m_pi**2), ('ss', m_ss_sq)):
                numbers[f'{name}_{label}'] = float(value.subs(common).subs(sym['m_qq'], np.sqrt(mass_sq)))
        else:
            numbers[name] = float(value.subs(common))
    return numbers


# ----------------------------------------------------------------------------- L5

def contact_amplitude(coefficient, factor_species, derivative_pairs, leg_species, leg_dot):
    """Amplitude (no factor i) of coefficient * prod(fields) with derivative pairs contracted.

    factor_species: species of each field factor; derivative_pairs: index pairs of factors
    carrying contracted derivatives; leg_dot(i, j) = p_i.p_j for incoming momenta.
    Each derivative gives -i p, so a contracted pair gives -p_i.p_j.
    """
    total = sp.S.Zero
    for assignment in itertools.permutations(range(len(leg_species))):
        if any(factor_species[f] != leg_species[assignment[f]] for f in range(len(factor_species))):
            continue
        term = sp.S.One
        for first, second in derivative_pairs:
            term *= -leg_dot(assignment[first], assignment[second])
        total += term
    return sp.expand(coefficient*total)


def mandelstam_dot(squared_masses, s, t):
    """p_i.p_j (all incoming) for 2->2 with u = sum m^2 - s - t."""
    u = sum(squared_masses)-s-t
    m = squared_masses
    table = {(0, 1): (s-m[0]-m[1])/2, (2, 3): (s-m[2]-m[3])/2,
             (0, 2): (t-m[0]-m[2])/2, (1, 3): (t-m[1]-m[3])/2,
             (0, 3): (u-m[0]-m[3])/2, (1, 2): (u-m[1]-m[2])/2}
    def dot(i, j):
        if i == j:
            return m[i]
        return table[(min(i, j), max(i, j))]
    return dot


def singlet_derivative_terms():
    """Dressed-frame expansion of the Kaiser:1998ds psi-tilde dependent terms (assumption).

    psi_tilde = tr Phi/f + 2 c_GA a/f_a, D_mu theta -> 2 c_GA d_mu a/f_a (currents removed).
    uu_second_derivative_ratio: second psi_tilde derivative of the <u_mu u^mu> coefficient at 0,
    divided by its value f^2/4.  V5_at_zero, V5_second_derivative: V_5 and V_5'' at 0.
    """
    f, f_a, c_GA, s, t, m_P_sq, m_a_sq = sp.symbols('f_pi f_a c_GA s t m_P^2 m_a^2', positive=True)
    uu_ratio, V5_zero, V5_second = sp.symbols(
        'uu_second_derivative_ratio V5_at_zero V5_second_derivative', real=True)
    # (i) (f^2/4)(uu_ratio/2) psi^2 <u u>, psi^2 -> 4 c_GA^2 a^2/f_a^2, <u u> -> (2/f^2) sum (dP)^2
    coefficient_aa_dPdP = f**2/4*uu_ratio/2*4*c_GA**2/f_a**2*2/f**2
    dot = mandelstam_dot([m_P_sq, m_P_sq, m_a_sq, m_a_sq], s, t)
    pp_to_aa_kinetic = contact_amplitude(coefficient_aa_dPdP, ['a', 'a', 'P', 'P'], [(2, 3)],
                                         ['P', 'P', 'a', 'a'], dot)
    # cross term psi^2 -> 2 (sqrt6 eta_0/f)(2 c_GA a/f_a): P_i P_i -> a eta_0
    coefficient_a_eta0_dPdP = f**2/4*uu_ratio/2*2*sp.sqrt(6)/f*2*c_GA/f_a*2/f**2
    eta0_sq = sp.Symbol('m_eta0^2', positive=True)
    dot_mixed = mandelstam_dot([m_P_sq, m_P_sq, m_a_sq, eta0_sq], s, t)
    pp_to_a_eta0 = contact_amplitude(coefficient_a_eta0_dPdP, ['a', 'e', 'P', 'P'], [(2, 3)],
                                     ['P', 'P', 'a', 'e'], dot_mixed)
    # V_4 D psi D theta (absent in the V_4 = 0 convention): eta_0 - a kinetic mixing K_{eta0 a}
    V4_zero = sp.Symbol('V4_at_zero', real=True)
    eta0_alp_kinetic_mixing = V4_zero*sp.sqrt(6)/f*2*c_GA/f_a
    # (ii) V_5(psi)(D theta)^2: V_5(0) term -> ALP kinetic renormalization (L = K_aa (da)^2/2)
    alp_kinetic_shift = 2*V5_zero*4*c_GA**2/f_a**2
    # V_5''/2 psi^2 (D theta)^2 with psi^2 -> (tr Phi)^2/f^2 = 6 eta_0^2/f^2
    coefficient_eta0_eta0_dada = V5_second/2*6/f**2*4*c_GA**2/f_a**2
    eta0_to_aa = contact_amplitude(coefficient_eta0_eta0_dada, ['e', 'e', 'a', 'a'], [(2, 3)],
                                   ['e', 'e', 'a', 'a'], dot)
    # a^2 (da)^2 term: aa -> aa is not energy growing on shell
    coefficient_aa_dada = V5_second/2*4*c_GA**2/f_a**2*4*c_GA**2/f_a**2
    dot_aa = mandelstam_dot([m_a_sq]*4, s, t)
    aa_to_aa = sp.simplify(contact_amplitude(coefficient_aa_dada, ['a', 'a', 'a', 'a'], [(2, 3)],
                                             ['a', 'a', 'a', 'a'], dot_aa))
    # S-wave, identical pairs in both states: a_0 = M/(32 pi) at high energy (constant in angle)
    def unitarity_scale(amplitude, multiplicity=1):
        leading = sp.limit(amplitude/s, s, sp.oo)
        return sp.solve(sp.Eq(sp.sqrt(multiplicity)*sp.Abs(leading)*s/(32*sp.pi),
                              sp.Rational(1, 2)), s)[0]
    return dict(
        pp_to_aa_kinetic=pp_to_aa_kinetic, alp_kinetic_shift=alp_kinetic_shift,
        pp_to_a_eta0=pp_to_a_eta0, eta0_alp_kinetic_mixing=eta0_alp_kinetic_mixing,
        eta0_to_aa=eta0_to_aa, aa_to_aa=aa_to_aa,
        scale_pp_to_aa_single=unitarity_scale(pp_to_aa_kinetic),
        scale_pp_to_aa_nonet=unitarity_scale(pp_to_aa_kinetic, 9),
        scale_eta0_to_aa=unitarity_scale(eta0_to_aa),
        symbols=dict(f=f, f_a=f_a, c_GA=c_GA, s=s, t=t, m_P_sq=m_P_sq, m_a_sq=m_a_sq,
                     uu_ratio=uu_ratio, V5_zero=V5_zero, V5_second=V5_second, V4_zero=V4_zero,
                     eta0_sq=eta0_sq))


# ----------------------------------------------------------------------------- L6

def toy_frame_comparison(values):
    """One meson + ALP: derivative frame vs dressed frame one-loop tadpoles.

    Derivative frame: L = (dP - v da)^2/2 + (da)^2/2 + f^2 B_0 m cos(sqrt2 P/f) - m_a0^2 a^2/2,
    v = sqrt2 c f/f_a.  Dressed frame: P = P_hat + v a.  Returns pole shifts
    delta m_n^2 = v_n^T delta M^2 v_n in both frames and the off-shell aa entries.
    """
    f, m_P_sq, c, f_a, m_a0_sq, scale = (values[k] for k in
                                         ('f', 'm_P_sq', 'c', 'f_a', 'm_a0_sq', 'scale'))
    mixing = np.sqrt(2)*c*f/f_a
    quartic_meson = 2*m_P_sq/f**2  # d^4/dP^4 of f^2 B_0 m cos(sqrt2 P/f)
    derivative_quartic = np.zeros((2, 2, 2, 2))
    derivative_quartic[0, 0, 0, 0] = quartic_meson
    transformation = np.array([[1.0, mixing], [0.0, 1.0]])  # (P, a) = T (P_hat, a)
    dressed_quartic = np.einsum('abcd,ai,bj,ck,dl->ijkl', derivative_quartic,
                                transformation, transformation, transformation, transformation)
    frames = {
        'derivative': (np.array([[1.0, -mixing], [-mixing, 1+mixing**2]]),
                       np.diag([m_P_sq, m_a0_sq]), derivative_quartic),
        'dressed': (np.eye(2), transformation.T @ np.diag([m_P_sq, m_a0_sq]) @ transformation,
                    dressed_quartic),
    }
    output = {}
    for name, (kinetic, mass_sq, quartic) in frames.items():
        spectrum = solve_spectrum_numerical(kinetic, mass_sq, field_labels=('P', 'a'))
        correction = tadpole_mass_correction(quartic, spectrum.pole_squared_masses,
                                             spectrum.pole_vectors, scale)
        shifts = [float(spectrum.pole_vectors[:, n] @ correction @ spectrum.pole_vectors[:, n])
                  for n in range(2)]
        output[name] = dict(poles=spectrum.pole_squared_masses, shifts=shifts,
                            off_shell=correction)
    return output


def toy_frame_series():
    """Symbolic O((f/f_a)^2) check of the on-shell meson-pole shift in both frames."""
    eps, m_P_sq, gap, quartic, A_P, A_a = sp.symbols('epsilon m_P^2 Delta V_PPPP A_P A_a',
                                                     positive=True)
    m_a0_sq = m_P_sq+gap  # bare ALP entry above the meson entry (nondegenerate)
    def first_order_pole_data(kinetic, mass_sq):
        # meson-like pole via perturbation theory in epsilon to second order (exact sympy series)
        s = sp.Symbol('s')
        polynomial = sp.expand((s*kinetic-mass_sq).det())
        roots = sp.solve(polynomial, s)
        root = [rt for rt in roots if sp.simplify(sp.limit(rt, eps, 0)-m_P_sq) == 0][0]
        null = (root*kinetic-mass_sq).nullspace()[0]
        norm = sp.sqrt((null.T*kinetic*null)[0])
        other = [rt for rt in roots if rt != root][0]
        null_other = (other*kinetic-mass_sq).nullspace()[0]
        norm_other = sp.sqrt((null_other.T*kinetic*null_other)[0])
        return root, null/norm, other, null_other/norm_other
    T = sp.Matrix([[1, eps], [0, 1]])
    results = {}
    for name, kinetic, mass_sq in (
            ('derivative', sp.Matrix([[1, -eps], [-eps, 1+eps**2]]), sp.diag(m_P_sq, m_a0_sq)),
            ('dressed', sp.eye(2), T.T*sp.diag(m_P_sq, m_a0_sq)*T)):
        root, vector, other, vector_other = first_order_pole_data(kinetic, mass_sq)
        tensor = sp.MutableDenseNDimArray.zeros(2, 2, 2, 2)
        base = sp.MutableDenseNDimArray.zeros(2, 2, 2, 2)
        base[0, 0, 0, 0] = quartic
        Tm = sp.eye(2) if name == 'derivative' else T
        for i, j, k, l in itertools.product(range(2), repeat=4):
            tensor[i, j, k, l] = Tm[0, i]*Tm[0, j]*Tm[0, k]*Tm[0, l]*quartic
        def contract(v_ext, v_int):
            return sum(tensor[i, j, k, l]*v_ext[i]*v_ext[j]*v_int[k]*v_int[l]
                       for i, j, k, l in itertools.product(range(2), repeat=4))/2
        shift = contract(vector, vector)*A_P+contract(vector, vector_other)*A_a
        results[name] = sp.factor(sp.simplify(sp.series(sp.simplify(shift), eps, 0, 3).removeO()))
    return results


# ----------------------------------------------------------------------------- main

def main():
    print('L1  hat M_n (linear spurion, constant diagonal c_A):  (-2i)^n m_q c_A^(n-1)(c_A - n c_P)')
    print('L1  V_{n_a,n_P} = B_0 f^(2-n_P) f_a^(-n_a) x [trace of hat M_{n_a} with generators]')
    forms = unmixed_closed_forms()
    print('L2 (a) ALP line / meson entry (u,u):', forms['alp_line_in_meson_entry_ratio'][(0, 0)])
    print('L2 (b) meson lines in M^2_aa:', forms['meson_lines_in_alp_entry'])
    print('L2 (c) ALP line in M^2_aa:', forms['alp_line_in_alp_entry'])
    print('L2 tree M^2_aa (quark part):', forms['tree_alp_entry_quark_part'])
    for ratio in (0.3, 1.0, 3.0):
        data = pole_sum_decomposition(ratio, gluon_coefficient=-3.0)
        print(f'L2 pole sum r={ratio}: poles {np.round(data["spectrum"].pole_squared_masses, 5)}'
              f' ALP-like index {data["alp_like_index"]}')
        print('   ALP-like-pole delta M^2 diag', np.round(np.diag(data['alp_like']), 6))
        print('   meson-like-pole delta M^2 diag', np.round(np.diag(data['meson_like']), 6))
    residual, im_B0 = bubble_unitarity_check()
    print('L3 Im B_0 =', im_B0, ' unitarity residual =', residual)
    print('L3 loop/tree:', {k: v for k, v in bubble_to_tree_ratios().items()
                            if k in ('PP', 'aP', 'aa')})
    _, crossings, _ = regime_map()
    print('L4 crossings:', crossings)
    print('L4 illustrative (c=1, c_P=0, c_GA=3, m_a=1 GeV, |B_0|=1):',
          {k: round(v, 3) for k, v in illustrative_regime_numbers().items()})
    terms = singlet_derivative_terms()
    for key in ('pp_to_aa_kinetic', 'pp_to_a_eta0', 'eta0_alp_kinetic_mixing', 'alp_kinetic_shift', 'eta0_to_aa', 'aa_to_aa',
                'scale_pp_to_aa_single', 'scale_pp_to_aa_nonet', 'scale_eta0_to_aa'):
        print('L5', key, '=', terms[key])
    print('L6 toy series:', toy_frame_series())


if __name__ == '__main__':
    main()
