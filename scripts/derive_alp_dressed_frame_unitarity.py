#!/usr/bin/env python3
"""Dressed frame of the LO ALP-meson theory and tree-level partial-wave unitarity.

Run:   .venv/bin/python scripts/derive_alp_dressed_frame_unitarity.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_frame_unitarity.py

Conventions (Alda:2026nrd; memory/alp/chiral-transformations/source-invariants/finite-transformations/alp-general-chiral-redefinition-invariants.md):
  * natural units, masses and momenta in GeV, squared masses and s, t, u in GeV^2;
  * U = exp(i Phi / f_pi), Phi = P_i lambda_i, tr(lambda_i lambda_j) = 2 delta_ij,
    U(3) generators ordered (lambda_0, ..., lambda_8) as in symbolic_utils;
  * x = a / f_a; mass source chi(x) = 2 B_0 M(x); currents r_mu = (d_mu x) C_R(x),
    l_mu = (d_mu x) C_L(x), C_R = C_V + C_A, C_L = C_V - C_A;
  * transport matrices d_x W_eta = i C_eta W_eta, W_eta(0) = 1; dressed sources
    hat M = W_R^dag M W_L and hat Theta(x) = Theta(x) + int_0^x tr C_A;
  * vertices are derivatives of the LO Lagrangian (overall Feynman factor i
    omitted, all momenta incoming), as in derive_alp_lo_joint_vertices.py;
  * partial waves: (a_J)_{fi} = N_i N_f / (32 pi) int_{-1}^{1} d cos(theta)
    P_J(cos theta) T_{fi}, with N = 1/sqrt(2) for two identical particles and 1
    otherwise; elastic unitarity for s >> masses gives |Re a_J| <= 1/2;
  * H5: total ALP potential V_tot(a) = V_QCD(a) + V_ext(a), V_QCD = -L_LO at
    P = 0 in the dressed frame, V_ext an external function of a only (it contains
    the bare mass term m_{a,0}^2 a^2/2); a = 0 is assumed to be the minimum,
    V_tot'(0) = 0. Feynman rule i V with V = -d^n V_tot/da^n for pure-ALP vertices.

H2 numerics: the source Taylor arrays and all basis-field vertices are exact
(rational, sqrt(2/3), sqrt(3)); s and t are exact rationals. The pole masses
of the coupled (P_i, a) system are roots of a quartic or quintic, so the pole
masses, pole vectors and the on-shell Gram matrix are evaluated with
DIGITS = 50 significant digits (mpmath), and the two frames are compared to
40 digits. This replaces exact arithmetic only in the spectrum.

All numerical inputs below are illustrative, not phenomenological fits, except
the stated f_pi, quark-mass and pion-mass inputs used for the printed numbers.
"""
from __future__ import annotations

import itertools
import math

import mpmath as mp
import numpy as np
import sympy as sp

from derive_alp_lo_joint_vertices import quadratic_blocks, vertex
from derive_alp_pole_projected_amplitudes import ccwz_vertex, projected_amplitude
from solve_alp_meson_spectrum import PoleSpectrum
from symbolic_utils import u3_generators

DIGITS = 50

# --------------------------------------------------------------------------
# H1: removal of the derivative couplings by transport matrices
# --------------------------------------------------------------------------


def transport_derivatives(current_derivatives, order):
    """Taylor derivatives W_n of the transport matrix, d_x W = i C W, W(0) = 1.

    current_derivatives[r] is d^r C / dx^r at x = 0 (missing entries are zero).
    Uses W_{n+1} = i sum_r binom(n, r) C_r W_{n-r}.
    """
    size = current_derivatives[0].rows
    current = list(current_derivatives) + [sp.zeros(size)] * (order + 1)
    derivatives = [sp.eye(size)]
    for n in range(order):
        derivatives.append(sp.expand(sp.I * sum(
            (sp.binomial(n, r) * current[r] * derivatives[n - r] for r in range(n + 1)),
            sp.zeros(size))))
    return derivatives


def _series(derivatives, variable):
    return sum((value * variable**n / sp.factorial(n) for n, value in enumerate(derivatives)),
               sp.zeros(derivatives[0].rows))


def _taylor_derivatives(matrix_series, variable, order):
    expanded = matrix_series.applyfunc(sp.expand)
    return [expanded.applyfunc(lambda entry: entry.coeff(variable, n)) * sp.factorial(n)
            for n in range(order + 1)]


def dressed_mass_derivatives(mass_derivatives, right_derivatives, left_derivatives, order):
    """Taylor derivatives hat M_n of hat M = W_R^dag M W_L through `order`."""
    x = sp.Symbol('x', real=True)
    right = _series(transport_derivatives(right_derivatives, order), x)
    left = _series(transport_derivatives(left_derivatives, order), x)
    size = mass_derivatives[0].rows
    mass = _series(list(mass_derivatives) + [sp.zeros(size)] * (order + 1 - len(mass_derivatives)), x)
    product = (right.H * mass * left).applyfunc(lambda e: sp.expand(e).subs(sp.conjugate(x), x))
    return _taylor_derivatives(product, x, order)


def dressed_theta_derivatives(theta_derivatives, right_derivatives, left_derivatives, order):
    """Taylor derivatives of hat Theta(x) = Theta(x) + int_0^x tr C_A, C_A = (C_R - C_L)/2."""
    theta = list(theta_derivatives) + [0] * (order + 1)
    size = right_derivatives[0].rows
    right = list(right_derivatives) + [sp.zeros(size)] * (order + 1)
    left = list(left_derivatives) + [sp.zeros(size)] * (order + 1)
    return [theta[0]] + [sp.expand(theta[n] + sp.trace(right[n - 1] - left[n - 1]) / 2)
                         for n in range(1, order + 1)]


def transformed_current_derivatives(current_derivatives, order):
    """Taylor derivatives of C' = g C g^dag + i g d_x g^dag for g = W^dag.

    This is the connection law of the general-redefinition record evaluated on
    the redefinition that removes the current; H1 states that it vanishes.
    """
    x = sp.Symbol('x', real=True)
    size = current_derivatives[0].rows
    transport = _series(transport_derivatives(current_derivatives, order + 1), x)
    current = _series(list(current_derivatives) + [sp.zeros(size)] * (order + 1 - len(current_derivatives)), x)
    redefinition = transport.H.applyfunc(lambda e: e.subs(sp.conjugate(x), x))
    redefinition_dagger = transport
    transformed = (redefinition * current * redefinition_dagger
                   + sp.I * redefinition * redefinition_dagger.diff(x))
    return _taylor_derivatives(transformed, x, order)


def pullback_field_strength(current_derivatives):
    """F_{01} = d_0 r_1 - d_1 r_0 - i [r_0, r_1] for r_mu = (d_mu x) C(x(y)).

    The ALP field x depends on two coordinates (y0, y1); C(x) is the Taylor
    polynomial with the supplied derivatives. Returns the expanded matrix.
    """
    y0, y1 = sp.symbols('y0 y1', real=True)
    field = sp.Function('x', real=True)(y0, y1)
    current = _series(current_derivatives, field)
    r0 = field.diff(y0) * current
    r1 = field.diff(y1) * current
    strength = r1.diff(y0) - r0.diff(y1) - sp.I * (r0 * r1 - r1 * r0)
    return strength.applyfunc(sp.expand)


# --------------------------------------------------------------------------
# Sources of the two frames and spectrum at high precision
# --------------------------------------------------------------------------


def linear_spurion_sources(quark_mass, axial, pseudoscalar, gluon, B_0, order=4):
    """Original-frame Taylor arrays: chi = 2 B_0 (m_q + i x {c_P, m_q}), C_R = -C_L = c_A, Theta = c_G x."""
    zero = sp.zeros(3)
    chi = [2 * B_0 * quark_mass, 2 * B_0 * sp.I * (pseudoscalar * quark_mass + quark_mass * pseudoscalar)]
    chi += [zero] * (order + 1 - len(chi))
    right = [axial] + [zero] * order
    left = [-axial] + [zero] * order
    theta = [0, gluon] + [0] * (order - 1)
    return chi, right, left, theta


def dressed_sources(chi, right, left, theta, order=4):
    """Dressed-frame arrays: chi -> W_R^dag chi W_L, currents removed, Theta -> hat Theta."""
    zero = sp.zeros(chi[0].rows)
    return (dressed_mass_derivatives(chi, right, left, order), [zero] * (order + 1),
            [zero] * (order + 1), dressed_theta_derivatives(theta, right, left, order))


def _to_mp(value):
    return mp.mpmathify(str(sp.N(value, DIGITS + 10)))


def _to_sympy(value):
    return sp.Float(mp.nstr(value, DIGITS + 5), DIGITS)


def high_precision_spectrum(kinetic, mass_squared, labels=None):
    """Poles and K-orthonormal pole vectors with DIGITS-digit arithmetic.

    Cholesky reduction K = L L^T and symmetric eigenproblem (mpmath); the
    ALP component (last field) of every pole vector is made positive.
    """
    mp.mp.dps = DIGITS + 10
    size = kinetic.rows
    kin = mp.matrix([[_to_mp(kinetic[i, j]) for j in range(size)] for i in range(size)])
    mass = mp.matrix([[_to_mp(mass_squared[i, j]) for j in range(size)] for i in range(size)])
    lower_inverse = mp.cholesky(kin) ** -1
    reduced = lower_inverse * mass * lower_inverse.T
    eigenvalues, eigenvectors = mp.eigsy((reduced + reduced.T) / 2)
    vectors = lower_inverse.T * eigenvectors
    order = sorted(range(size), key=lambda n: eigenvalues[n])
    columns, poles = [], []
    for n in order:
        column = [vectors[i, n] for i in range(size)]
        if column[-1] < 0:
            column = [-entry for entry in column]
        columns.append(sp.Matrix([_to_sympy(entry) for entry in column]))
        poles.append(_to_sympy(eigenvalues[n]))
    pole_vectors = sp.Matrix.hstack(*columns)
    float_matrix = lambda m: m.applyfunc(lambda e: sp.Float(sp.N(e, DIGITS + 5), DIGITS))
    return PoleSpectrum(field_labels=tuple(labels or [f'F{i}' for i in range(size - 1)] + ['a']),
                        kinetic=float_matrix(kinetic), mass_squared=float_matrix(mass_squared),
                        pole_squared_masses=poles, pole_vectors=pole_vectors,
                        residues=[c * c.T for c in columns], exact=True)


def on_shell_gram(squared_masses, s, t):
    """Gram matrix p_i.p_j of all-incoming momenta with p_i^2 = m_i^2, s=(p1+p2)^2, t=(p1+p3)^2."""
    m1, m2, m3, m4 = squared_masses
    u = m1 + m2 + m3 + m4 - s - t
    return sp.Matrix([[m1, (s - m1 - m2) / 2, (t - m1 - m3) / 2, (u - m1 - m4) / 2],
                      [(s - m1 - m2) / 2, m2, (u - m2 - m3) / 2, (t - m2 - m4) / 2],
                      [(t - m1 - m3) / 2, (u - m2 - m3) / 2, m3, (s - m3 - m4) / 2],
                      [(u - m1 - m4) / 2, (t - m2 - m4) / 2, (s - m3 - m4) / 2, m4]])


_BASIS_VERTEX_CACHE = {}


def multilinear_vertex_function(exact_basis, sources, f_pi, f_a, anomaly_mass_squared):
    """Vertex function for projected_amplitude built from exact basis-field vertices.

    The CCWZ vertex is multilinear in the leg directions and symmetric under
    simultaneous permutations of legs and momenta. Exact vertices of the
    basis fields (mesons of exact_basis, then the ALP) are computed once with
    `vertex` for sorted field labels at f_a = 1 and contracted with the
    (numerical) pole directions; a vertex with n ALP legs scales as f_a^-n,
    since every ALP enters through x = a/f_a. Three- and four-leg vertices do
    not depend on the bare ALP mass. Equivalent to ccwz_vertex(exact_basis, ...).
    """
    size = len(exact_basis) + 1
    zero = sp.zeros(exact_basis[0].rows)
    frame_key = str((exact_basis, sources, f_pi, anomaly_mass_squared))
    cache = _BASIS_VERTEX_CACHE.setdefault(frame_key, {})

    def basis_terms(labels):
        if labels not in cache:
            symbols = sp.symbols(f'_r0:{len(labels)}', real=True)
            expression = vertex([exact_basis[i] if i < size - 1 else zero for i in labels],
                                [1 if i == size - 1 else 0 for i in labels], symbols, *sources,
                                f_pi, sp.Integer(1), anomaly_mass_squared, 0)
            terms = []
            for powers, coefficient in sp.Poly(expression, *symbols).terms():
                legs = tuple(i for i, power in enumerate(powers) for _ in range(power))
                terms.append((coefficient, legs))
            cache[labels] = terms
        if labels not in scaled:
            alp_legs = sum(1 for i in labels if i == size - 1)
            scaled[labels] = [(sp.N(coefficient / f_a**alp_legs, DIGITS + 5), legs)
                              for coefficient, legs in cache[labels]]
        return scaled[labels]

    scaled = {}

    def function(directions, momenta, gram):
        count = len(directions)
        products = [[(momenta[i].T * gram * momenta[j])[0] for j in range(count)] for i in range(count)]
        total = sp.S.Zero
        for labels in itertools.product(range(size), repeat=count):
            weight = sp.prod([directions[k][labels[k]] for k in range(count)])
            if weight == 0:
                continue
            order = sorted(range(count), key=lambda k: labels[k])
            for coefficient, legs in basis_terms(tuple(labels[k] for k in order)):
                factor = products[order[legs[0]]][order[legs[1]]] if legs else 1
                total += weight * coefficient * factor
        return sp.expand(total)

    return function


def frame_theory(generator_indices, sources, f_pi, f_a, anomaly_mass_squared, bare_alp_mass_squared):
    """Quadratic blocks, high-precision spectrum and a CCWZ vertex function for one frame."""
    generators = u3_generators()
    exact_basis = [sp.Matrix(generators[i]) for i in generator_indices]
    kinetic, mass_squared = quadratic_blocks(exact_basis, sources[0][:3], sources[1][:3],
                                             sources[2][:3], sources[3][:3], f_pi, f_a,
                                             anomaly_mass_squared, bare_alp_mass_squared)
    spectrum = high_precision_spectrum(kinetic, mass_squared)
    return kinetic, mass_squared, spectrum, multilinear_vertex_function(
        exact_basis, sources, f_pi, f_a, anomaly_mass_squared)


def pole_amplitude(spectrum, vertex_function, pole_indices, s, t):
    """On-shell 2->2 tree amplitude (contact + s, t, u exchanges, full D^-1) among poles."""
    gram = on_shell_gram([spectrum.pole_squared_masses[i] for i in pole_indices], s, t)
    return sp.N(projected_amplitude(spectrum, pole_indices, gram, vertex_function), DIGITS)


def alp_like_and_pion_like(spectrum, pion_field_index):
    """Indices of the pole with the largest ALP component and, among the others, the largest lambda_3 one."""
    vectors = spectrum.pole_vectors
    size = spectrum.size
    alp_like = max(range(size), key=lambda n: abs(vectors[size - 1, n]))
    pion_like = max((n for n in range(size) if n != alp_like),
                    key=lambda n: abs(vectors[pion_field_index, n]))
    return alp_like, pion_like


def illustrative_inputs(off_diagonal_axial=False):
    """Illustrative neutral-sector inputs in GeV (H2 frame comparison only).

    With real symmetric c_A, c_P and real diagonal m_q, the Lagrangian of both
    frames is invariant under U -> U^T with a fixed (c_A symmetric gives
    D_mu(U^T) = (D_mu U)^T; hat M(x) is symmetric). P_7 (lambda_7^T = -lambda_7)
    is odd and P_0, P_3, P_8, P_6, a are even, so P_7 neither mixes with nor is
    exchanged between states built from the even fields; it is omitted from
    the d-s basis (checked in the unit tests).
    """
    quark_mass = sp.diag(sp.Rational(1, 500), sp.Rational(1, 200), sp.Rational(1, 10))
    axial = sp.diag(1, sp.Rational(1, 2), -sp.Rational(1, 3))
    if off_diagonal_axial:
        axial[1, 2] = axial[2, 1] = sp.Rational(1, 4)
    return dict(quark_mass=quark_mass, axial=axial,
                pseudoscalar=sp.diag(sp.Rational(1, 5), -sp.Rational(1, 4), sp.Rational(1, 2)),
                gluon=sp.Rational(1, 3), B_0=sp.Rational(27, 10), f_pi=sp.Rational(23, 250),
                anomaly_mass_squared=sp.Rational(4, 5), bare_alp_mass_squared=sp.Integer(4),
                generator_indices=(3, 8, 0, 6) if off_diagonal_axial else (3, 8, 0))


def compare_frames(ratio, kinematic_points, off_diagonal_axial=False):
    """Return spectra and amplitudes of both frames for the three H2 channels.

    ratio = f_pi / f_a. Channels (all-incoming pole labels):
    (A, pi; A, pi), (pi, pi; A, A), (A, A; A, A), with A the ALP-like pole.
    """
    data = illustrative_inputs(off_diagonal_axial)
    f_a = data['f_pi'] / ratio
    original = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'],
                                      data['gluon'], data['B_0'])
    dressed = dressed_sources(*original)
    common = (data['f_pi'], f_a, data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
    results = {}
    for name, sources in (('original', original), ('dressed', dressed)):
        kinetic, mass_squared, spectrum, function = frame_theory(data['generator_indices'], sources, *common)
        alp, pion = alp_like_and_pion_like(spectrum, 0)
        channels = {'A pi -> A pi': [alp, pion, alp, pion], 'pi pi -> A A': [pion, pion, alp, alp],
                    'A A -> A A': [alp, alp, alp, alp]}
        amplitudes = {(channel, point): pole_amplitude(spectrum, function, indices, *point)
                      for channel, indices in channels.items() for point in kinematic_points}
        results[name] = dict(kinetic=kinetic, mass_squared=mass_squared, spectrum=spectrum,
                             amplitudes=amplitudes, alp_like=alp, pion_like=pion)
    return results


# --------------------------------------------------------------------------
# H3: O(s) coupled-channel S-wave matrix
# --------------------------------------------------------------------------


def kinetic_four_meson_vertex(generators, f_pi, legs, momenta_products):
    """Closed form of the four-meson vertex of (f_pi^2/4) <u_mu u^mu> without sources.

    V = -1/(24 f_pi^2) sum_sigma [tr(l1 l2 l3 l4) - tr(l1 l3 l2 l4)]_sigma p_sigma2 . p_sigma4,
    from L_4 = tr([Phi, d Phi][Phi, d Phi]) / (48 f_pi^2). legs are generator indices;
    momenta_products(i, j) returns p_i.p_j.
    """
    total = 0
    for perm in itertools.permutations(range(4)):
        g = [generators[legs[k]] for k in perm]
        total += ((g[0] * g[1] * g[2] * g[3]).trace() - (g[0] * g[2] * g[1] * g[3]).trace()) \
            * momenta_products(perm[1], perm[3])
    return sp.expand(-total / (24 * f_pi**2))


def s_wave_coefficient_tensor(generators_numeric):
    """Tensor P_IJKL (floats, units 1/f_pi^2) with a_0 O(s) part = N N s P / (16 pi).

    Massless kinematics: p1.p2 = p3.p4 = s/2, p1.p3 = p2.p4 = t/2, p1.p4 = p2.p3 = u/2;
    S-wave averages <t> = <u> = -s/2. Built from the closed form above with
    numpy traces; the coefficient of s is returned.
    """
    lam = np.array(generators_numeric, dtype=complex)
    four = np.einsum('aij,bjk,ckl,dli->abcd', lam, lam, lam, lam)
    # weight of p_i.p_j after S-wave averaging, in units of s: s/2 -> 1/2, t/2, u/2 -> -1/4
    pair_weight = np.array([[0, .5, -.25, -.25], [.5, 0, -.25, -.25],
                            [-.25, -.25, 0, .5], [-.25, -.25, .5, 0]])
    tensor = np.zeros((len(lam),) * 4, dtype=complex)
    letters = 'abcd'
    for perm in itertools.permutations(range(4)):
        # leg perm[k] sits at trace position k
        positions = ''.join(letters[perm.index(k)] for k in range(4))
        first = np.einsum('abcd->' + positions, four)
        swapped = np.einsum('acbd->' + positions, four)
        tensor += (first - swapped) * pair_weight[perm[1], perm[3]]
    return (-tensor / 24).real


def pair_matrix(tensor):
    """S-wave matrix over symmetric two-particle states, in units of s / (32 pi f_pi^2).

    (a_0)_{(IJ),(KL)} = N_IJ N_KL 2 s P_IJKL / (32 pi), N = 1/sqrt 2 for I = J.
    """
    size = tensor.shape[0]
    pairs = [(i, j) for i in range(size) for j in range(i, size)]
    norm = lambda p: 1 / math.sqrt(2) if p[0] == p[1] else 1.0
    return np.array([[2 * norm(p) * norm(q) * tensor[p[0], p[1], q[0], q[1]] for q in pairs] for p in pairs]), pairs


def pad_with_alp(tensor):
    """Embed the meson tensor in (P_i, a), with zero entries whenever an ALP index occurs."""
    size = tensor.shape[0]
    padded = np.zeros((size + 1,) * 4)
    padded[:size, :size, :size, :size] = tensor
    return padded


def rotate_tensor(tensor, orthogonal):
    """T'_{nmkl} = T_IJKL O_In O_Jm O_Kk O_Ll."""
    return np.einsum('ijkl,ia,jb,kc,ld->abcd', tensor, orthogonal, orthogonal, orthogonal, orthogonal)


def cayley_orthogonal(size, seed):
    """Exact rational orthogonal matrix (I - A)(I + A)^-1 from a random antisymmetric rational A."""
    import random
    generator = random.Random(seed)
    antisymmetric = sp.zeros(size)
    for i in range(size):
        for j in range(i + 1, size):
            value = sp.Rational(generator.randint(-5, 5), generator.randint(1, 4))
            antisymmetric[i, j], antisymmetric[j, i] = value, -value
    return (sp.eye(size) - antisymmetric) * (sp.eye(size) + antisymmetric).inv()


def singlet_s_wave(generator_indices, f_pi, mass_squared):
    """Exact S-wave of the flavour-singlet state sum_a |a a> / sqrt(n) for a common meson mass.

    Uses the CCWZ vertex generator with chi = m^2 * 1, no topological term and
    no sources; returns sum_{a,c} (1/(64 pi)) int dcos T_{aa,cc} / n with
    s + t + u = 4 m^2 and t = -(s - 4 m^2)(1 - cos)/2.
    """
    s, cos = sp.symbols('s c', real=True)
    t = -(s - 4 * mass_squared) * (1 - cos) / 2
    gram = on_shell_gram([mass_squared] * 4, s, t)
    generators = u3_generators()
    zero = sp.zeros(3)
    sources = ([mass_squared * sp.eye(3)] + [zero] * 4, [zero] * 5, [zero] * 5, [0] * 5,
               f_pi, sp.Integer(1), 0, 0)
    momenta = [sp.eye(4)[:, i] for i in range(4)]
    basis = [sp.Matrix(generators[i]) for i in generator_indices]
    total = 0
    for a in range(len(basis)):
        for c in range(a, len(basis)):
            directions = [sp.Matrix([1 if k == a else 0 for k in range(len(basis))] + [0])] * 2 \
                + [sp.Matrix([1 if k == c else 0 for k in range(len(basis))] + [0])] * 2
            amplitude = ccwz_vertex(basis, sources, directions, momenta, gram)
            total += (1 if a == c else 2) * amplitude
    projected = sp.integrate(total, (cos, -1, 1)) / (64 * sp.pi * len(basis))
    return sp.simplify(projected), s


def unitarity_scale_squared(flavours, f_pi, mass_squared, bound):
    """s at which N s/(32 pi f^2) - m^2/(16 pi N f^2) equals `bound` (1/2 or 1)."""
    s = sp.Symbol('s', positive=True)
    wave = flavours * s / (32 * sp.pi * f_pi**2) - mass_squared / (16 * sp.pi * flavours * f_pi**2)
    return sp.solve(sp.Eq(wave, bound), s)[0]


# --------------------------------------------------------------------------
# H4: constant contact vertices in the dressed frame
# --------------------------------------------------------------------------


def constant_vertex_closed_form(dressed_mass, meson_generators, alp_count, B_0, f_pi, f_a):
    """Vertex with alp_count ALPs and mesons lambda_{i1..ik} from (f_pi^2/4) <chi_+>.

    V = f_pi^2 B_0 / (f_pi^k f_a^n) Re[ (-i)^k / k! sum_sigma tr(hat M_n lambda_sigma1 ... lambda_sigmak) ],
    dressed_mass[n] = hat M_n.
    """
    k = len(meson_generators)
    summed = sum((sp.prod([meson_generators[i] for i in perm], sp.eye(3)) for perm in itertools.permutations(range(k))),
                 sp.zeros(3)) if k else sp.eye(3)
    value = (-sp.I)**k / sp.factorial(k) * (dressed_mass[alp_count] * summed).trace()
    return sp.expand(f_pi**2 * B_0 / (f_pi**k * f_a**alp_count) * sp.re(sp.expand(value)))


def topological_vertex_closed_form(theta_derivatives, meson_generators, alp_count, anomaly_mass_squared, f_pi, f_a):
    """Vertex from -(f_pi^2 m_0^2 / (4 N)) (tr Phi / f_pi + 2 hat Theta(x))^2, N = 3."""
    k = len(meson_generators)
    if k > 2:
        return sp.S.Zero
    # coefficient structure: (tr Phi)^2 / f^2 + 4 tr Phi Theta / f + 4 Theta^2
    traces = [g.trace() for g in meson_generators]
    if k == 2:
        value = 2 * traces[0] * traces[1] / f_pi**2 if alp_count == 0 else 0
    elif k == 1:
        value = 4 * traces[0] / f_pi * theta_derivatives[alp_count] if alp_count >= 1 else 0
    else:
        value = 4 * sum(sp.binomial(alp_count, r) * theta_derivatives[r] * theta_derivatives[alp_count - r]
                        for r in range(alp_count + 1))
    return sp.expand(-f_pi**2 * anomaly_mass_squared / 12 * value / f_a**alp_count)


def diagonal_linear_spurion_dressed_mass(quark_masses, axial, pseudoscalar, n):
    """hat M_n for M = m + i x {c_P, m}, constant c_A, all diagonal: m_q (-2i)^n a_q^(n-1) (a_q - n p_q)."""
    if n == 0:
        return sp.diag(*quark_masses)
    return sp.diag(*[m * (-2 * sp.I)**n * a**(n - 1) * (a - n * p)
                     for m, a, p in zip(quark_masses, axial, pseudoscalar)])


def contact_unitarity_bounds(dressed_mass, meson_generator, B_0, f_pi):
    """Upper bounds on f_pi / f_a from |Re a_0| <= 1/2, neglecting mass mixing.

    Returns dict with the constant amplitudes (in terms of r = f_pi/f_a) and the
    bounds for aa->aa, a P -> a P and P P -> a a (P = given generator).
    """
    r = sp.Symbol('r', positive=True)
    f_a = f_pi / r
    alp4 = constant_vertex_closed_form(dressed_mass, [], 4, B_0, f_pi, f_a)
    alp2_meson2 = constant_vertex_closed_form(dressed_mass, [meson_generator] * 2, 2, B_0, f_pi, f_a)
    waves = {'a a -> a a': alp4 / (32 * sp.pi), 'a P -> a P': alp2_meson2 / (16 * sp.pi),
             'P P -> a a': alp2_meson2 / (32 * sp.pi)}
    bounds = {}
    for name, wave in waves.items():
        power = 4 if name == 'a a -> a a' else 2
        coefficient = sp.simplify(wave / r**power)
        bounds[name] = (sp.Rational(1, 2) / sp.Abs(coefficient))**sp.Rational(1, power)
    return waves, bounds, r


# --------------------------------------------------------------------------
# H5: external ALP potential V_ext(a), not induced by QCD
# --------------------------------------------------------------------------


def qcd_alp_potential(dressed_mass, theta_derivatives, B_0, f_pi, f_a, anomaly_mass_squared, alp,
                      flavours=3):
    """QCD-induced ALP potential V_QCD(a) = -L_LO at P = 0, Taylor polynomial through a^4.

    V_QCD = -f_pi^2 B_0 Re tr hat M(a/f_a) + (f_pi^2 m_0^2 / N_f) hat Theta(a/f_a)^2,
    from (f_pi^2/4) <chi_+> = f_pi^2 B_0 Re <hat M U^dag> and
    -(f_pi^2 m_0^2/(4 N_f)) X_TC^2 with X_TC = -2 hat Theta at P = 0.
    dressed_mass[n] = hat M_n (not chi), theta_derivatives[n] = hat Theta_n.
    """
    order = len(dressed_mass) - 1
    x = alp / f_a
    mass = sum((sp.re(sp.expand(dressed_mass[n].trace())) * x**n / sp.factorial(n) for n in range(order + 1)))
    theta = sum((theta_derivatives[n] * x**n / sp.factorial(n) for n in range(min(len(theta_derivatives), order + 1))))
    potential = -f_pi**2 * B_0 * mass + f_pi**2 * anomaly_mass_squared * theta**2 / flavours
    return sp.expand(sp.series(sp.expand(potential), alp, 0, order + 1).removeO())


def potential_derivative_at_zero(potential, alp, n):
    """d^n V / da^n at a = 0."""
    return sp.expand(sp.diff(potential, alp, n).subs(alp, 0))


def alp_elastic_amplitude(quartic, cubic, alp_mass_squared, s, t):
    """Tree aa -> aa amplitude from V_tot: contact plus s, t, u ALP exchange.

    quartic = V_tot''''(0), cubic = V_tot'''(0). Feynman rule i V for a vertex
    V = d^n L / da^n = -d^n V_tot / da^n, propagator i/(q^2 - m_a^2), iM = sum of
    diagrams: M = -V'''' - V'''^2 [1/(s-m_a^2) + 1/(t-m_a^2) + 1/(u-m_a^2)],
    u = 4 m_a^2 - s - t. Mass mixing with mesons neglected.
    """
    u = 4 * alp_mass_squared - s - t
    return -quartic - cubic**2 * sum(1 / (channel - alp_mass_squared) for channel in (s, t, u))


def alp_elastic_s_wave(quartic, cubic, alp_mass_squared, s):
    """Closed-form S wave of aa -> aa, N_i = N_f = 1/sqrt 2, for s > 4 m_a^2.

    a_0 = [-V'''' - V'''^2/(s - m_a^2) + 2 V'''^2 ln((s - 3 m_a^2)/m_a^2)/(s - 4 m_a^2)] / (32 pi),
    using int_{-1}^{1} dcos / (t - m_a^2) = -2 ln((s - 3 m_a^2)/m_a^2)/(s - 4 m_a^2) for
    t = -(s - 4 m_a^2)(1 - cos)/2 (the u channel gives the same integral).
    """
    m2 = alp_mass_squared
    return (-quartic - cubic**2 / (s - m2) + 2 * cubic**2 * sp.log((s - 3 * m2) / m2) / (s - 4 * m2)) / (32 * sp.pi)


def quartic_unitarity_interval(qcd_quartic_coefficient, external_quartic, bound=sp.Rational(1, 2)):
    """Allowed interval of r = f_pi/f_a > 0 from |Re a_0(aa->aa)| <= bound at s >> m_a^2.

    With mass mixing neglected and V_tot'''(0) = 0 (or s large enough that the
    exchange terms are negligible), a_0 = (k r^4 - V_ext'''')/(32 pi), where the
    QCD contact vertex is V^QCD_aaaa = k r^4, k = B_0 Re tr hat M_4 / f_pi^2 (linear
    hat Theta). Returns (lower, upper) with upper possibly sp.oo, or None if no r > 0 works.
    """
    width = 32 * sp.pi * bound
    k, v = qcd_quartic_coefficient, external_quartic
    if k == 0:
        return (sp.Integer(0), sp.oo) if abs(v) <= width else None
    ends = sorted([(v - width) / k, (v + width) / k], key=lambda e: float(e))
    low, high = max(ends[0], 0, key=lambda e: float(e)), ends[1]
    if float(high) < 0:
        return None
    return (low**sp.Rational(1, 4), high**sp.Rational(1, 4))


def phase_model_quartic_check():
    """d^4/da^4 of -lambda |phi|^4 at fixed sigma, phi = (v + sigma) exp(i a / v)/sqrt(2).

    Returns (the potential term, its fourth a-derivative); the latter is zero.
    """
    alp, sigma = sp.symbols('a sigma', real=True)
    v, coupling = sp.symbols('v lambda', positive=True)
    field = (v + sigma) * sp.exp(sp.I * alp / v) / sp.sqrt(2)
    term = -coupling * sp.simplify(sp.expand(field * sp.conjugate(field)))**2
    return sp.factor(term), sp.simplify(sp.diff(term, alp, 4))


# --------------------------------------------------------------------------


def main():
    mp.mp.dps = 30
    print('H1: dressed-frame connection C\' through x^3 for generic 3x3 rational C(x):')
    import random
    gen = random.Random(7)

    def hermitian():
        m = sp.Matrix(3, 3, lambda i, j: sp.Rational(gen.randint(-4, 4), gen.randint(1, 3))
                      + sp.I * sp.Rational(gen.randint(-4, 4), gen.randint(1, 3)))
        return (m + m.H) / 2
    current = [hermitian() for _ in range(4)]
    print('  all zero:', all(m == sp.zeros(3) for m in transformed_current_derivatives(current, 3)))
    print('  pull-back field strength zero:', pullback_field_strength(current) == sp.zeros(3))

    print('\nH4: closed-form constant vertices (dressed frame), hat M_n = d^n hat M/dx^n at 0:')
    B_0, f_pi, f_a = sp.symbols('B_0 f_pi f_a', positive=True)
    mu, md, ms, cu, cd, cs, pu, pd, ps = sp.symbols('m_u m_d m_s c_u c_d c_s p_u p_d p_s', real=True)
    masses, axial, pseudo = (mu, md, ms), (cu, cd, cs), (pu, pd, ps)
    hat = [diagonal_linear_spurion_dressed_mass(masses, axial, pseudo, n) for n in range(5)]
    lam3 = sp.Matrix(u3_generators()[3])
    print('  V_aaaa =', sp.factor(constant_vertex_closed_form(hat, [], 4, B_0, f_pi, f_a)))
    print('  V_aaP3P3 =', sp.factor(constant_vertex_closed_form(hat, [lam3, lam3], 2, B_0, f_pi, f_a)))
    hat_ap = [diagonal_linear_spurion_dressed_mass(masses, axial, (0, 0, 0), n) for n in range(5)]
    waves, bounds, r = contact_unitarity_bounds(hat_ap, lam3, B_0, f_pi)
    print('  c_q^P = 0, c_q^A = c_q^{AP} = diag(c_u, c_d, c_s); bounds on f_pi/f_a from |Re a_0| <= 1/2:')
    for name in waves:
        print(f'    {name}: a_0 = {sp.factor(waves[name])};  f_pi/f_a <= {sp.simplify(bounds[name])}')

    print('\nIllustrative: c_u = 1 only, c_q^P = 0; inputs m_u = 2.16 MeV, m_d = 4.70 MeV,')
    print('  m_pi0 = 134.977 MeV, B_0 = m_pi0^2/(m_u+m_d), f_pi = 92.07 MeV:')
    numbers = {mu: sp.Rational(216, 100000), md: sp.Rational(470, 100000), ms: sp.Rational(935, 10000),
               cu: 1, cd: 0, cs: 0, f_pi: sp.Rational(9207, 100000)}
    numbers[B_0] = sp.Rational(134977, 1000000)**2 / (numbers[mu] + numbers[md])
    for name in bounds:
        print(f'    {name}: f_pi/f_a <= {sp.N(bounds[name].subs(numbers), 4)}')

    print('\nH5: external potential V_ext(a); V_tot = V_QCD + V_ext; a_0(aa->aa) at s >> m_a^2 = -V_tot\'\'\'\'(0)/(32 pi)')
    alp = sp.Symbol('a', real=True)
    theta = sp.symbols('Theta0:5', real=True)
    hat_general = [sp.diag(*sp.symbols(f'M{n}_u M{n}_d M{n}_s')) for n in range(5)]
    m0 = sp.Symbol('m_0', positive=True)
    v_qcd = qcd_alp_potential(hat_general, theta, B_0, f_pi, f_a, m0**2, alp)
    print('  -V_QCD\'\'\'\'(0) =', sp.factor(-potential_derivative_at_zero(v_qcd, alp, 4)))
    s_sym, m2_sym, quartic_sym, cubic_sym = sp.symbols('s m_a2 V4 V3', real=True)
    print('  a_0 with exchange =', alp_elastic_s_wave(quartic_sym, cubic_sym, m2_sym, s_sym))
    print('  phase model: -lambda|phi|^4 =', *phase_model_quartic_check())
    k_value = (16 * B_0 * sum(m * c**4 for m, c in zip(masses, axial)) / f_pi**2).subs(numbers)
    print('  c_u = 1, c_q^P = 0, inputs as above: allowed f_pi/f_a from |V_aaaa^QCD - V_ext\'\'\'\'| <= 16 pi')
    for external in (0, 8 * sp.pi, -8 * sp.pi, -16 * sp.pi, 100):
        interval = quartic_unitarity_interval(k_value, external)
        text = 'none' if interval is None else f'[{sp.N(interval[0], 4)}, {sp.N(interval[1], 4)}]'
        print(f'    V_ext\'\'\'\' = {external}: {text}')

    print('\nGeneral dressed-frame S waves (mass mixing neglected), P_i P_i identical pair:')
    print('  a_0(aa->aa) = f_pi^2 B_0 Re tr(hatM_4) / (32 pi f_a^4)')
    print('  a_0(aP_i->aP_i) = -B_0 Re tr(hatM_2 lambda_i^2) / (16 pi f_a^2)')
    print('  a_0(P_iP_i->aa) = -B_0 Re tr(hatM_2 lambda_i^2) / (32 pi f_a^2)')

    print('\nH3: Lambda_unit^2 = 16 pi f^2/N + 2 m^2/N^2 (|Re a_0|<=1/2), 32 pi f^2/N + 2 m^2/N^2 (|a_0|<=1);')
    print('  f_pi = 92.07 MeV; mass term for a common meson mass m (here m = m_pi0 = 134.977 MeV):')
    f_value = sp.Rational(9207, 100000)
    m_value = sp.Rational(134977, 1000000)
    for label, flavours in (('U(3)/SU(3) singlet, N=3', 3), ('pi pi I=0, N=2', 2)):
        for bound in (sp.Rational(1, 2), 1):
            chiral = sp.sqrt(unitarity_scale_squared(flavours, f_value, 0, bound))
            massive = sp.sqrt(unitarity_scale_squared(flavours, f_value, m_value**2, bound))
            print(f'    {label}, bound {bound}: Lambda = {sp.N(chiral * 1000, 4)} MeV (chiral limit), '
                  f'{sp.N(massive * 1000, 4)} MeV (with m)')

    print('\nH2 example (diagonal sources, f_pi/f_a = 1, s = 3 GeV^2, t = -3/5 GeV^2):')
    results = compare_frames(sp.Integer(1), [(sp.Integer(3), -sp.Rational(3, 5))])
    for key, value in results['original']['amplitudes'].items():
        print(f'    {key[0]}: original {sp.N(value, 15)}, dressed {sp.N(results["dressed"]["amplitudes"][key], 15)}')


if __name__ == '__main__':
    main()
