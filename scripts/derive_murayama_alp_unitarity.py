#!/usr/bin/env python3
"""Tree-level partial-wave unitarity of ALP scattering in the only-up model.

Model: Murayama:2026ioh / DiLuzio:2026poh with charges (n_u, n_d, n_s) = (1, 0, 0),
f_a = f_varphi, theta = 0, CP-conserving vacuum, LO U(3) chiral Lagrangian, tree level.

Conventions (natural units, every dimensionful quantity in GeV):
  * Polar fields:      Phi = (f_a + sigma) exp(i a / f_a) / sqrt(2).
  * Cartesian fields:  Phi = (f_a + s + i p) / sqrt(2).
  * Potential at U = 1 (the quark-mass term of the chiral Lagrangian for n_u = 1):
        V = m_varphi^2 |Phi|^2 + lambda |Phi|^4 - B_0 m_u f_pi^2 (1 + sigma/f_a) cos(a/f_a)
          (+ an optional even potential V_ext(a) of the phase only, Part D).
    Here lambda |Phi|^4 = (lambda/4) [(f_a+s)^2 + p^2]^2.
  * Feynman rule of an n-point vertex with all momenta incoming: i * vertex, with
    vertex = d^n L / d fields (derivatives replaced by -i p).  Amplitudes are defined
    by i M = sum of diagrams, so an exchange of a scalar of mass m in the channel x gives
    vertex_1 * vertex_2 / (m^2 - x).
  * Partial waves: (a_J)_fi = N_i N_f / (32 pi) int_{-1}^{1} dcos P_J M_fi,
    N = 1/sqrt(2) for two identical particles, else 1.  Primary condition
    |Re a_0| <= 1/2, secondary |a_0| <= 1.

Approximation used everywhere unless stated: ALP-meson mass mixing is neglected and the
external ALP is the field a (polar) or p (Cartesian) with its unmixed mass
d^2V/da^2 = B_0 m_u f_pi^2/f_a^2.  The expression B_0 m_u f_pi^2/f_a^2 is written out:
the symbol m_0 already denotes the topological (U(1)_A) mass in this repository.

Minimal invocation example:
    .venv/bin/python scripts/derive_murayama_alp_unitarity.py
Tests:
    .venv/bin/python -m unittest scripts.test_derive_murayama_alp_unitarity
"""
import sympy as sp

B_0, m_u, f_pi, f_a, lam = sp.symbols('B_0 m_u f_pi f_a lambda', positive=True)
mass_phi_squared = sp.symbols('m_varphi_squared', real=True)
mandelstam_s, mandelstam_t, mandelstam_u = sp.symbols('s t u', real=True)
cos_theta = sp.symbols('cos_theta', real=True)
external_mass_squared = sp.symbols('m_ext_squared', positive=True)
radial_field, phase_field = sp.symbols('sigma a', real=True)
cartesian_radial, cartesian_angular = sp.symbols('s_field p_field', real=True)
external_potential_second, external_potential_fourth = sp.symbols('V_ext_2 V_ext_4', real=True)

UP_SOURCE = B_0*m_u*f_pi**2          # coefficient of the quark-mass term at U = 1, GeV^4
UNMIXED_ALP_MASS_SQUARED = UP_SOURCE/f_a**2


# ---------------------------------------------------------------- potentials and vacuum
def polar_potential(include_external=False):
    """Potential at U = 1 in polar fields; optional even phase-only potential V_ext(a).

    V_ext is represented by its Taylor coefficients V_ext''(0) a^2/2 + V_ext''''(0) a^4/24.
    """
    rho_squared = (f_a + radial_field)**2/2                  # |Phi|^2, independent of a
    potential = (mass_phi_squared*rho_squared + lam*rho_squared**2
                 - UP_SOURCE*(1 + radial_field/f_a)*sp.cos(phase_field/f_a))
    if include_external:
        potential += (external_potential_second*phase_field**2/2
                      + external_potential_fourth*phase_field**4/24)
    return potential


def cartesian_potential():
    """Potential at U = 1 in Cartesian fields (quark source linear in the fields)."""
    rho_squared = ((f_a + cartesian_radial)**2 + cartesian_angular**2)/2
    return (mass_phi_squared*rho_squared + lam*rho_squared**2
            - UP_SOURCE*(1 + cartesian_radial/f_a))


def vacuum_solution():
    """Stationarity in sigma fixes m_varphi^2; return it and the radial mass."""
    potential = polar_potential()
    at_vacuum = {radial_field: 0, phase_field: 0}
    vacuum_mass = sp.solve(sp.diff(potential, radial_field).subs(at_vacuum), mass_phi_squared)[0]
    radial_mass_squared = sp.simplify(sp.diff(potential, radial_field, 2).subs(at_vacuum)
                                      .subs(mass_phi_squared, vacuum_mass))
    return vacuum_mass, radial_mass_squared


VACUUM_PHI_MASS_SQUARED, RADIAL_MASS_SQUARED = vacuum_solution()
m_sigma_squared = sp.symbols('m_sigma_squared', positive=True)


def lambda_from_radial_mass():
    """lambda = (m_sigma^2 - B_0 m_u f_pi^2/f_a^2)/(2 f_a^2)."""
    return sp.solve(sp.Eq(RADIAL_MASS_SQUARED, m_sigma_squared), lam)[0]


# ---------------------------------------------------------------- Part A: polar vertices
def polar_vertices(include_external=False):
    """Vertices of the polar Lagrangian needed for a a -> a a.

    Returns (unmixed ALP mass squared, contact a^4 vertex, momentum-independent part of
    the sigma a a vertex, coefficient c of sigma (d a)^2 in the kinetic term).
    The kinetic term (1/2)(1 + sigma/f_a)^2 (d a)^2 gives the sigma a a vertex
    -2 c p_1.p_2 for incoming ALP momenta p_1, p_2.
    """
    potential = polar_potential(include_external).subs(mass_phi_squared, VACUUM_PHI_MASS_SQUARED)
    at_vacuum = {radial_field: 0, phase_field: 0}
    alp_mass_squared = sp.simplify(sp.diff(potential, phase_field, 2).subs(at_vacuum))
    contact = sp.simplify(-sp.diff(potential, phase_field, 4).subs(at_vacuum))
    radial_alp_alp_potential_vertex = sp.simplify(
        -sp.diff(potential, radial_field, phase_field, phase_field).subs(at_vacuum))
    kinetic_coefficient = sp.Rational(1, 2)*(1 + radial_field/f_a)**2
    sigma_da_da_coefficient = sp.diff(kinetic_coefficient, radial_field).subs(radial_field, 0)
    return alp_mass_squared, contact, radial_alp_alp_potential_vertex, sigma_da_da_coefficient


def polar_sigma_alp_alp_vertex(exchanged_momentum_squared, alp_mass_squared_on_shell,
                               include_external=False):
    """Off-shell sigma a a vertex with on-shell ALPs of mass squared m_ext^2.

    2 p_1.p_2 = q^2 - 2 m_ext^2 for incoming p_1, p_2 with (p_1 + p_2)^2 = q^2.
    """
    _, _, potential_vertex, kinetic_coefficient = polar_vertices(include_external)
    two_p1_dot_p2 = exchanged_momentum_squared - 2*alp_mass_squared_on_shell
    # c sigma (d a)^2 with two identical ALPs: vertex 2 c (-i p_1).(-i p_2) = -c (2 p_1.p_2).
    return -kinetic_coefficient*two_p1_dot_p2 + potential_vertex


def polar_alp_alp_amplitude(alp_mass_squared_on_shell, include_external=False):
    """a a -> a a in polar fields, sigma exchange in s, t, u plus contact.

    The Mandelstam variables obey s + t + u = 4 m_ext^2 (eliminated through u).
    """
    _, contact, _, _ = polar_vertices(include_external)
    u_value = 4*alp_mass_squared_on_shell - mandelstam_s - mandelstam_t
    amplitude = contact
    for channel in (mandelstam_s, mandelstam_t, u_value):
        vertex = polar_sigma_alp_alp_vertex(channel, alp_mass_squared_on_shell, include_external)
        amplitude += vertex**2/(RADIAL_MASS_SQUARED - channel)
    return amplitude


# ------------------------------------------------------------ Cartesian vertices
def cartesian_vertices():
    """Cubic and quartic Cartesian vertices (constants) and the two masses."""
    potential = cartesian_potential().subs(mass_phi_squared, VACUUM_PHI_MASS_SQUARED)
    at_vacuum = {cartesian_radial: 0, cartesian_angular: 0}
    derivative = lambda *fields: sp.simplify(-sp.diff(potential, *fields).subs(at_vacuum))
    S, P = cartesian_radial, cartesian_angular
    return {
        'mass_s': -derivative(S, S), 'mass_p': -derivative(P, P),
        'sss': derivative(S, S, S), 'spp': derivative(S, P, P),
        'ssss': derivative(S, S, S, S), 'pppp': derivative(P, P, P, P),
        'sspp': derivative(S, S, P, P), 'ppp': derivative(P, P, P), 'spp_s': derivative(S, S, P),
    }


def cartesian_alp_alp_amplitude():
    """p p -> p p: -6 lambda contact plus sigma exchange with vertex -2 lambda f_a."""
    vertices = cartesian_vertices()
    u_value = 4*vertices['mass_p'] - mandelstam_s - mandelstam_t
    amplitude = vertices['pppp']
    for channel in (mandelstam_s, mandelstam_t, u_value):
        amplitude += vertices['spp']**2/(vertices['mass_s'] - channel)
    return amplitude


def closed_form_alp_alp_amplitude():
    """Expected closed form -6 lambda + 4 lambda^2 f_a^2 sum_x 1/(m_sigma^2 - x)."""
    u_value = 4*UNMIXED_ALP_MASS_SQUARED - mandelstam_s - mandelstam_t
    return -6*lam + 4*lam**2*f_a**2*sum(1/(RADIAL_MASS_SQUARED - x)
                                        for x in (mandelstam_s, mandelstam_t, u_value))


def hypothesis_polar_form():
    """Hypothesis A in the polar form, external mass B_0 m_u f_pi^2/f_a^2."""
    u_value = 4*UNMIXED_ALP_MASS_SQUARED - mandelstam_s - mandelstam_t
    return UP_SOURCE/f_a**4 + sum((x - UNMIXED_ALP_MASS_SQUARED)**2/(f_a**2*(RADIAL_MASS_SQUARED - x))
                                  for x in (mandelstam_s, mandelstam_t, u_value))


def polar_minus_cartesian_with_external_mass():
    """Polar minus Cartesian amplitude if the ALP mass squared in the kinematics and in
    2 p_1.p_2 is B_0 m_u f_pi^2/f_a^2 + delta instead of the unmixed value.

    Only a diagnostic: with mixing neglected the pole of a is at B_0 m_u f_pi^2/f_a^2.
    Returns (delta, difference, its fixed-angle s -> infinity limit, residual of
    difference - [8 delta/f_a^2 + sum_x ((2 lambda f_a^2 - 2 delta)^2 - 4 lambda^2 f_a^4)
    /(f_a^2 (m_sigma^2 - x))]).
    """
    mass_shift = sp.symbols('delta_m_squared', real=True)
    shifted_mass = UNMIXED_ALP_MASS_SQUARED + mass_shift
    polar = polar_alp_alp_amplitude(shifted_mass)
    # The Cartesian amplitude has no derivative couplings: it depends on s, t, u only.
    u_value = 4*shifted_mass - mandelstam_s - mandelstam_t
    cartesian = -6*lam + 4*lam**2*f_a**2*sum(1/(RADIAL_MASS_SQUARED - x)
                                             for x in (mandelstam_s, mandelstam_t, u_value))
    difference = sp.simplify(polar - cartesian)
    scale = sp.symbols('scale', positive=True)
    asymptotic = sp.limit(difference.subs({mandelstam_s: scale, mandelstam_t: -scale/3}), scale, sp.oo)
    expected = sp.expand(8*mass_shift/f_a**2) + sum(
        ((2*lam*f_a**2 - 2*mass_shift)**2 - 4*lam**2*f_a**4)/(f_a**2*(RADIAL_MASS_SQUARED - x))
        for x in (mandelstam_s, mandelstam_t, u_value))
    return mass_shift, difference, sp.simplify(asymptotic), sp.simplify(difference - expected)


# ------------------------------------------------------------ Part B: S wave
def pole_term_s_wave(pole_mass_squared, t_at_cos_minus_one, t_at_cos_plus_one):
    """int_{-1}^{1} dcos 1/(M^2 - t(cos)) for t linear in cos (exact)."""
    t_offset = (t_at_cos_plus_one + t_at_cos_minus_one)/2
    t_slope = (t_at_cos_plus_one - t_at_cos_minus_one)/2
    return sp.log((pole_mass_squared - t_offset + t_slope)
                  / (pole_mass_squared - t_offset - t_slope))/t_slope


def alp_alp_s_wave():
    """Exact a_0(a a -> a a) with N_i N_f = 1/2, using the closed form amplitude.

    t = -(s - 4 m_a^2)(1 - cos)/2 and u = -(s - 4 m_a^2)(1 + cos)/2; the t- and
    u-channel terms give equal logarithms.
    """
    momentum_squared_times_four = mandelstam_s - 4*UNMIXED_ALP_MASS_SQUARED
    t_of_cos = -momentum_squared_times_four*(1 - cos_theta)/2
    u_of_cos = -momentum_squared_times_four*(1 + cos_theta)/2
    amplitude = -6*lam + 4*lam**2*f_a**2*(1/(RADIAL_MASS_SQUARED - mandelstam_s)
                                          + 1/(RADIAL_MASS_SQUARED - t_of_cos)
                                          + 1/(RADIAL_MASS_SQUARED - u_of_cos))
    integral = sp.integrate(amplitude, (cos_theta, -1, 1), conds='none')
    return sp.simplify(integral/(2*32*sp.pi)), amplitude


def alp_alp_s_wave_closed_form():
    """a_0 = [-12 lambda + 8 lambda^2 f_a^2/(m_sigma^2 - s)
              + 16 lambda^2 f_a^2 log(1 + (s - 4 m_a^2)/m_sigma^2)/(s - 4 m_a^2)]/(64 pi)."""
    k = mandelstam_s - 4*UNMIXED_ALP_MASS_SQUARED
    return (-12*lam + 8*lam**2*f_a**2/(RADIAL_MASS_SQUARED - mandelstam_s)
            + 16*lam**2*f_a**2*sp.log(1 + k/RADIAL_MASS_SQUARED)/k)/(64*sp.pi)


def intermediate_energy_s_wave():
    """Leading term of a_0 for 4 m_a^2 << s << m_sigma^2, from the polar form.

    Massless external kinematics; keeps the term of order s^2/m_sigma^2 at fixed f_a.
    """
    t_of_cos = -mandelstam_s*(1 - cos_theta)/2
    u_of_cos = -mandelstam_s*(1 + cos_theta)/2
    leading = sum(x**2 for x in (mandelstam_s, t_of_cos, u_of_cos))/(f_a**2*m_sigma_squared)
    return sp.simplify(sp.integrate(leading, (cos_theta, -1, 1))/(2*32*sp.pi))


def intermediate_energy_scale():
    """s^2 at which the intermediate-energy a_0 reaches 1/2."""
    s_squared = sp.symbols('s_squared', positive=True)
    a0 = intermediate_energy_s_wave().subs(mandelstam_s**2, s_squared)
    return sp.solve(sp.Eq(a0, sp.Rational(1, 2)), s_squared)[0]


# ------------------------------------------------------------ Part C: coupled channels
def two_body_t_range(s_value, masses_in, masses_out):
    """t at cos = -1 and cos = +1 for 1 2 -> 3 4 with t = (p_1 - p_3)^2."""
    m1, m2 = masses_in
    m3, m4 = masses_out
    root_s = sp.sqrt(s_value)
    energy_1 = (s_value + m1**2 - m2**2)/(2*root_s)
    energy_3 = (s_value + m3**2 - m4**2)/(2*root_s)
    momentum_in = sp.sqrt(energy_1**2 - m1**2)
    momentum_out = sp.sqrt(energy_3**2 - m3**2)
    t_base = m1**2 + m3**2 - 2*energy_1*energy_3
    return t_base - 2*momentum_in*momentum_out, t_base + 2*momentum_in*momentum_out


def cartesian_channel_s_waves(s_value):
    """Exact tree S waves in the Cartesian theory for the channels
    {p p/sqrt2, s s/sqrt2} (2x2 matrix) and s p -> s p.

    Each amplitude is contact + exchanges; constant terms are integrated directly and
    pole terms with pole_term_s_wave.  Returns (matrix, s p wave, dict of exchange parts).
    """
    v = cartesian_vertices()
    m_s = sp.sqrt(v['mass_s'])
    m_p = sp.sqrt(v['mass_p'])

    def wave(contact, s_channel, t_terms, u_terms, masses_in, masses_out, normalization):
        t_minus, t_plus = two_body_t_range(s_value, masses_in, masses_out)
        total_mass_squared = sum(m**2 for m in masses_in + masses_out)
        u_minus, u_plus = total_mass_squared - s_value - t_minus, total_mass_squared - s_value - t_plus
        exchange = 2*s_channel
        exchange += sum(c*pole_term_s_wave(M2, t_minus, t_plus) for c, M2 in t_terms)
        exchange += sum(c*pole_term_s_wave(M2, u_minus, u_plus) for c, M2 in u_terms)
        return normalization*(2*contact + exchange)/(32*sp.pi), normalization*exchange/(32*sp.pi)

    half = sp.Rational(1, 2)
    spp, sss = v['spp'], v['sss']
    pp_pp = wave(v['pppp'], spp**2/(v['mass_s'] - s_value),
                 [(spp**2, v['mass_s'])], [(spp**2, v['mass_s'])], (m_p, m_p), (m_p, m_p), half)
    ss_ss = wave(v['ssss'], sss**2/(v['mass_s'] - s_value),
                 [(sss**2, v['mass_s'])], [(sss**2, v['mass_s'])], (m_s, m_s), (m_s, m_s), half)
    pp_ss = wave(v['sspp'], spp*sss/(v['mass_s'] - s_value),
                 [(spp**2, v['mass_p'])], [(spp**2, v['mass_p'])], (m_p, m_p), (m_s, m_s), half)
    sp_sp = wave(v['sspp'], spp**2/(v['mass_p'] - s_value),
                 [(sss*spp, v['mass_s'])], [(spp**2, v['mass_p'])], (m_s, m_p), (m_s, m_p), 1)
    matrix = sp.Matrix([[pp_pp[0], pp_ss[0]], [pp_ss[0], ss_ss[0]]])
    exchanges = {'pp_pp': pp_pp[1], 'ss_ss': ss_ss[1], 'pp_ss': pp_ss[1], 'sp_sp': sp_sp[1]}
    return matrix, sp_sp[0], exchanges


def asymptotic_contact_matrix():
    """Contact-only S-wave matrices: (lambda/16pi)[[-3,-1],[-1,-3]] and -lambda/(8 pi)."""
    v = cartesian_vertices()
    matrix = sp.Matrix([[v['pppp'], v['sspp']], [v['sspp'], v['ssss']]])/(32*sp.pi)
    return matrix, 2*v['sspp']/(32*sp.pi)


def coupling_bounds():
    """Upper bounds on lambda from |Re a_0| <= 1/2 and |a_0| <= 1 at s >> m_sigma^2."""
    matrix, sp_wave = asymptotic_contact_matrix()
    largest = max(abs(e) for e in matrix.eigenvals())
    out = {}
    for name, value in (('alp_alp', abs(matrix[0, 0])), ('coupled', largest), ('sigma_alp', abs(sp_wave))):
        out[name] = (sp.solve(sp.Eq(value, sp.Rational(1, 2)), lam)[0],
                     sp.solve(sp.Eq(value, 1), lam)[0])
    return out


def bound_on_radial_mass(lambda_bound):
    """Translate lambda <= bound into m_sigma^2 - B_0 m_u f_pi^2/f_a^2 <= 2 bound f_a^2."""
    return sp.expand(2*lambda_bound*f_a**2)


def intermediate_window_condition():
    """lambda above which s^2 = (48 pi/5) f_a^2 m_sigma^2 lies below s = m_sigma^2.

    The intermediate-energy scale is reached inside 4 m_a^2 << s << m_sigma^2 only if
    (48 pi/5) f_a^2 < m_sigma^2 = B_0 m_u f_pi^2/f_a^2 + 2 lambda f_a^2.
    """
    return sp.solve(sp.Eq(sp.Rational(48, 5)*sp.pi*f_a**2, RADIAL_MASS_SQUARED), lam)[0]


# ------------------------------------------------------------ Part D: phase-only potential
def asymptotic_constant_with_external_potential():
    """s -> infinity limit (at fixed angle) of the polar a a -> a a amplitude with V_ext(a).

    The ALP mass entering the kinematics is the unmixed one including V_ext''(0).
    """
    alp_mass_squared, _, _, _ = polar_vertices(include_external=True)
    amplitude = polar_alp_alp_amplitude(alp_mass_squared, include_external=True)
    scale = sp.symbols('scale', positive=True)
    fixed_angle = amplitude.subs({mandelstam_s: scale, mandelstam_t: -scale/3})
    return sp.simplify(sp.limit(fixed_angle, scale, sp.oo))


# ------------------------------------------------------------ Part E: mesons
def gell_mann_matrices():
    """lambda_1, lambda_2, lambda_3 embedded in U(3) (only the ones used here)."""
    l1 = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    l2 = sp.Matrix([[0, -sp.I, 0], [sp.I, 0, 0], [0, 0, 0]])
    l3 = sp.diag(1, -1, 0)
    return {1: l1, 2: l2, 3: l3}


def alp_three_meson_vertex(indices):
    """Contact vertex a P_i P_j P_k from f_pi^2 B_0 Re<M_hat U^dagger>, U = exp(i Phi/f_pi).

    Polar connection-free representative: M_hat = diag(m_u (1+sigma/f_a) e^{-i a/f_a}, m_d, m_s).
    Computed by direct expansion with symbolic fields.
    """
    lam_matrices = gell_mann_matrices()
    field_symbols = sp.symbols('P1 P2 P3', real=True)
    mesons = sum((field_symbols[k - 1]*lam_matrices[k] for k in (1, 2, 3)), sp.zeros(3))
    epsilon = sp.symbols('epsilon', real=True)
    scaled = epsilon*mesons/f_pi
    u_dagger = sp.eye(3)
    term = sp.eye(3)
    for order in range(1, 4):
        term = term*(-sp.I*scaled)/order
        u_dagger += term
    m_d, m_s = sp.symbols('m_d m_s', positive=True)
    mass_source = sp.diag(m_u*sp.exp(-sp.I*phase_field/f_a), m_d, m_s)
    trace = sp.expand((mass_source*u_dagger).trace()).coeff(epsilon, 3)
    lagrangian = f_pi**2*B_0*(trace + sp.conjugate(trace))/2
    derivative = sp.diff(lagrangian, phase_field, *[field_symbols[i - 1] for i in indices])
    return sp.simplify(derivative.subs({phase_field: 0, **{x: 0 for x in field_symbols}}))


def alp_meson_contact_bounds():
    """a P_3 -> P_3 P_3 and a P_3 -> pi^+ pi^- S waves from the constant contact vertex.

    a P_3 -> P_3 P_3: N_i N_f = 1/sqrt(2).  pi^+ pi^-: vertex V_{311} = V_{322}, N_i N_f = 1.
    """
    v333 = alp_three_meson_vertex((3, 3, 3))
    v311 = alp_three_meson_vertex((3, 1, 1))
    v322 = alp_three_meson_vertex((3, 2, 2))
    a0_333 = 2*v333/(sp.sqrt(2)*32*sp.pi)
    a0_charged = 2*v311/(32*sp.pi)
    ratio = sp.symbols('f_pi_over_f_a', positive=True)
    bounds = {}
    for name, a0 in (('P3P3P3', a0_333), ('P3pipi', a0_charged)):
        expr = sp.simplify(a0.subs(f_a, f_pi/ratio))
        bounds[name] = (a0, sp.solve(sp.Eq(expr, sp.Rational(1, 2)), ratio)[0])
    return v333, v311, v322, bounds


def alp_meson_scattering_polar_cartesian():
    """a P_3 -> a P_3 (mixing neglected): polar contact + t-channel sigma vs Cartesian.

    Polar: contact -B_0/(2 f_a^2) Re tr(M_hat_2 {lambda_3, lambda_3}) = B_0 m_u/f_a^2 and
    t-channel sigma exchange with vertices -(t - B_0 m_u f_pi^2/f_a^2)/f_a and -B_0 m_u/f_a.
    Cartesian: no contact, vertices -2 lambda f_a and -B_0 m_u/f_a.
    """
    contact = B_0*m_u/f_a**2
    radial_meson_vertex = -B_0*m_u/f_a
    polar = contact + polar_sigma_alp_alp_vertex(mandelstam_t, UNMIXED_ALP_MASS_SQUARED) \
        * radial_meson_vertex/(RADIAL_MASS_SQUARED - mandelstam_t)
    cartesian = cartesian_vertices()['spp']*radial_meson_vertex/(RADIAL_MASS_SQUARED - mandelstam_t)
    return sp.simplify(polar), sp.simplify(cartesian)


# ------------------------------------------------------------ numbers
NUMERICAL_INPUTS = {m_u: sp.Rational(216, 100000), f_pi: sp.Rational(9207, 100000),
                    B_0: sp.Rational(1349768, 10000000)**2/(sp.Rational(216 + 470, 100000))}


def main():
    vacuum_mass, radial_mass = VACUUM_PHI_MASS_SQUARED, RADIAL_MASS_SQUARED
    print('m_varphi^2 =', vacuum_mass, '; m_sigma^2 =', radial_mass)
    alp_mass, contact, pot_vertex, kin = polar_vertices()
    print('A: unmixed ALP mass^2 =', alp_mass, '; a^4 contact =', contact,
          '; sigma a a potential vertex =', pot_vertex, '; sigma (da)^2 coefficient =', kin)
    print('A: off-shell sigma a a vertex =',
          sp.factor(polar_sigma_alp_alp_vertex(sp.Symbol('q2'), alp_mass)))
    polar = polar_alp_alp_amplitude(alp_mass)
    print('A: polar - closed form =', sp.simplify(polar - closed_form_alp_alp_amplitude()))
    print('A: hypothesis polar - closed form =', sp.simplify(hypothesis_polar_form() - closed_form_alp_alp_amplitude()))
    print('A: Cartesian - closed form =', sp.simplify(cartesian_alp_alp_amplitude() - closed_form_alp_alp_amplitude()))
    _, _, asymptotic_difference, residual = polar_minus_cartesian_with_external_mass()
    print('A: shifted external mass: polar - Cartesian at s -> oo =', asymptotic_difference,
          '; residual of closed form =', residual)
    a0, _ = alp_alp_s_wave()
    sample = {B_0: sp.Rational(2656, 1000), m_u: sp.Rational(216, 100000),
              f_pi: sp.Rational(9207, 100000), f_a: sp.Rational(1, 2), lam: sp.Rational(13, 10)}
    for s_value in (sp.Rational(1, 5), sp.Rational(9, 10), 3, 50):
        exact = sp.N(a0.subs(sample).subs(mandelstam_s, s_value), 15)
        closed = sp.N(alp_alp_s_wave_closed_form().subs(sample).subs(mandelstam_s, s_value), 15)
        print(f'B: s = {s_value} GeV^2: integrated a_0 = {exact}, closed form = {closed}')
    print('B: a_0(s -> oo) =', sp.limit(alp_alp_s_wave_closed_form(), mandelstam_s, sp.oo))
    print('B: intermediate a_0 =', intermediate_energy_s_wave(), '; s^2 scale =', intermediate_energy_scale())
    print('B: intermediate scale inside s < m_sigma^2 only for lambda >', intermediate_window_condition())
    print('C: contact matrices =', asymptotic_contact_matrix())
    for s_value in (10**4, 10**8):
        matrix, sigma_alp_wave, exchanges = cartesian_channel_s_waves(sp.Integer(s_value))
        print(f'C: s = {s_value} GeV^2 (sample point): exchange parts =',
              {k: sp.N(v.subs(sample), 4) for k, v in exchanges.items()})
    print('C: bounds (primary, secondary) =', coupling_bounds())
    for name, (primary, _) in coupling_bounds().items():
        print(f'C: {name}: m_sigma^2 - B_0 m_u f_pi^2/f_a^2 <=', bound_on_radial_mass(primary))
    print('D: asymptotic constant with V_ext =', asymptotic_constant_with_external_potential())
    print('E: a P3 -> a P3 polar, Cartesian =', alp_meson_scattering_polar_cartesian())
    v333, v311, v322, bounds = alp_meson_contact_bounds()
    print('E: V_aP3P3P3 =', v333, '; V_aP3P1P1 =', v311, '; V_aP3P2P2 =', v322)
    for name, (a0_value, ratio_bound) in bounds.items():
        print(f'E: {name}: a_0 =', a0_value, '; f_pi/f_a <=', ratio_bound,
              '=', sp.N(ratio_bound.subs(NUMERICAL_INPUTS), 4))


if __name__ == '__main__':
    main()
