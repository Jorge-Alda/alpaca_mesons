#!/usr/bin/env python3
"""Is the LO lower bound on the bare ALP mass invariant under ALP-dependent chiral redefinitions?

Run:   .venv/bin/python scripts/derive_alp_lo_stability_bound_invariance.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_stability_bound_invariance.py

Question. scripts/derive_alp_lo_stability_bound.py derived, in the convention "quark
spurion linear in a, no independent a^2 term" (memory/alp/shared/finite-alp-meson-mixing.md),

  m_a0^2 >= eps^2 [2 tr(c^P chi_0 c^P)
                   - 2 (m_0^2/N_f)(tr c^P + c_G)^2 / (1 + (m_0^2/N_f) sum_q 1/chi_q)],   eps = f_pi/f_a.

Here we test, by exact computation, whether this statement survives ALP-dependent chiral
field redefinitions g_R(x), g_L(x), x = a/f_a, and write the frame-independent form.

Conventions (Alda:2026nrd; memory/alp/shared/alp-general-chiral-redefinition-invariants.md,
alp-chiral-source-transformations.md, alp-pole-alignment-chiral-generator.md,
alp-dressed-mass-first-order-invariants.md, alp-dressed-topological-invariant.md):
 * N_f = 3, nine Hermitian U(3) generators, tr(lambda_i lambda_j) = 2 delta_ij, ordered
   (lambda_0, ..., lambda_8); U = exp(i Phi/f_pi), Phi = P_i lambda_i.
 * LO Lagrangian (as in scripts/derive_alp_pole_alignment_chiral_generator.py):
     1/2 (da)^2 + f_pi^2/4 <u u> + f_pi^2/4 <chi U^dag + U chi^dag>
     - f_pi^2 m_0^2/(4 N_f) X_TC^2 - 1/2 m_a0^2 a^2,
   u = -dPhi/f_pi + 2 (dx) C_A(0) + (terms of cubic order in the fields),
   X_TC = -<Phi>/f_pi - 2 Theta(x).
 * A "frame" is the set of source functions: chi(x) = chi_0 + x chi_(1) + x^2 chi_(2)
   (power-series coefficients chi_(n), GeV^2), right/left connection coefficients
   C_R(x), C_L(x) (derivatives C_eta,n at x = 0), and the slope Theta_1 = dTheta/dx(0).
   chi = 2 B_0 M with the quark mass function M(x); chi_0 = 2 B_0 m_q, m_q diagonal, > 0.
 * Linear-spurion frame (finite-mixing record): chi_(1) = i[c^S, chi_0] + i{c^P, chi_0},
   chi_(2) = 0, C_R = c^V + c^A, C_L = c^V - c^A (constant), Theta_1 = c_G.
 * Constant-generator redefinition: g_R = exp(i x G_R), g_L = exp(i x G_L),
   G_R = G_V + G_A, G_L = G_V - G_A (Hermitian, noncommuting),
   chi' = g_R chi g_L^dag, C_eta' = g_eta C_eta g_eta^dag + G_eta,
   Theta' = Theta - (1/2i) log(det g_R/det g_L) = (Theta_1 - tr G_A) x.
 * Transport matrices dW_eta/dx = i C_eta W_eta, W_eta(0) = 1; dressed mass function
   Mhat(x) = W_R^dag M W_L; Mhat_n = d^n Mhat/dx^n(0) (derivative convention of the
   general record); Phat_1 = (Mhat_1 - Mhat_1^dag)/(2i) = -{c^AP, m_q};
   dressed topological slope Thetahat_1 = Theta_1 + tr C_A(0) = c_GA = c_G + tr c^A;
   vacuum-angle slope thetabar_1 = Theta_1 + (1/2) d/dx arg det M = c_G + tr c^P
   = c_GA - tr c^AP; c^AP = c^A - c^P, c^VS = c^V - c^S.
Units: f_pi, f_a, m_q in GeV; B_0 in GeV; chi, m_0^2, m_a0^2, M^2 in GeV^2; couplings
and generators dimensionless. All test values are exact rationals chosen for testing only.

Results established below (each one is computed, not assumed):
 1. A constant-generator redefinition acts on the quadratic terms as the congruence
    K' = T^T K T, M^2' = T^T M^2 T with the shift P = P' + w a,
    w_i = -(f_pi/f_a) tr(lambda_i G_A), PROVIDED the x^2 term chi'_(2) of the transformed
    spurion is kept in M^2'_aa. m_a0^2 enters only M^2_aa, with unit coefficient, in
    every frame. The Schur complement S = M^2_aa - M^2_aP (M^2_PP)^-1 M^2_Pa is invariant
    under the shift, and the inertia of M^2 is invariant under the congruence.
 2. In ANY frame, M^2 = (linear-spurion blocks built from c^P, Theta_1)
    - eps^2 Re tr chi_(2) e_a e_a^T, where c^P is read off from {c^P, chi_0} =
    (chi_(1) - chi_(1)^dag)/(2i); the Hermitian part of chi_(1) does not enter. Hence the
    frame-dependent statement of the bound is
      m_a0^2 >= eps^2 [Re tr chi_(2) + 2 tr(c^P chi_0 c^P)
                       - 2 (m_0^2/N_f)(tr c^P + Theta_1)^2/(1 + (m_0^2/N_f) sum_q 1/chi_q)].
 3. Under the redefinition c^P' = c^P + G_A, c^S' = c^S + G_V, c^A' = c^A + G_A,
    c^V' = c^V + G_V, c_G' = c_G - tr G_A; invariant: c^AP, c^VS, c_GA, tr c^P + c_G.
    Re tr chi'_(2) = 2 tr(c^P chi_0 c^P) - 2 tr(c^P' chi_0 c^P') (starting from the linear
    frame). So the previous bracket evaluated with c^P', c_G' differs from the original one
    by 2 tr(G_A chi_0 G_A) + 2 tr(G_A {c^P, chi_0}); the difference is exactly cancelled by
    the a^2 term of the transformed spurion. The topological part is invariant by itself.
 4. Frame-independent form (in terms of the dressed sources):
      m_a0^2 >= (f_pi^2/f_a^2) [ B_0 Re tr Mhat_2 + 2 tr(c^AP chi_0 c^AP)
                - 2 (m_0^2/N_f)(c_GA - tr c^AP)^2/(1 + (m_0^2/N_f) sum_q 1/chi_q) ]
             = (f_pi^2/f_a^2) [ B_0 d^2/dx^2 tr (Mhat^dag Mhat)^(1/2) |_0
                - 2 (m_0^2/N_f) thetabar_1^2/(1 + (m_0^2/N_f) sum_q 1/chi_q) ].
    The second line is the second derivative of the sum of the singular values of M(x)
    (equal to those of Mhat(x)), i.e. of the vacuum energy min_U V(U, x). With Phat_1:
    2 tr(c^AP chi_0 c^AP) = 2 B_0 sum_{q,q'} |Phat_1,qq'|^2/(m_q + m_q'),
    c_GA - tr c^AP = Thetahat_1 + tr(m_q^-1 Phat_1)/2. c^VS (Shat_1) drops out.
    The only second-order datum not fixed by the first-order invariants is the single
    real number Re tr Mhat_2. Theta_2 and the anti-Hermitian part of tr Mhat_2 do not enter.
 5. Linear-spurion convention (M_2 = 0, C_eta constant): Re tr Mhat_2 =
    4 tr(c^P m_q c^P) - 4 tr(c^AP m_q c^AP), independent of c^V, c^S and of C_eta,1;
    then item 4 reproduces the previous bracket exactly.
 6. Frame without kinetic mixing (G_R = -C_R, G_L = -C_L): C' = 0, K = 1, chi' = 2 B_0 Mhat,
    c^P' = -c^AP, Theta_1' = c_GA, chi'_(2) = B_0 Mhat_2.

Numerical checks: smallest eigenvalue of M^2 just above/below the invariant bound in three
frames with the same physics; the previous bracket with the new c^P', c_G' is not a
threshold; the minimum of the full nonlinear chiral potential over U(3) is the same function
of x in two frames, and its second derivative reproduces item 4.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from mpmath import mp
from scipy.linalg import expm
from scipy.optimize import minimize

from derive_alp_charged_meson_decoupling import symbolic_hermitian
from derive_alp_dressed_topological_invariant import (
    ALP_FIELD_RATIO,
    matrix_transport_series,
    taylor_coefficients,
    truncate_matrix,
)
from derive_alp_lo_neutral_blocks import ALL_GENERATOR_INDICES, NEUTRAL_GENERATOR_INDICES, closed_form_blocks
from derive_alp_lo_stability_bound import bound_coefficient
from derive_alp_pole_alignment_chiral_generator import quadratic_hessian
from derive_alp_single_pole_alignment import shift_matrix
from symbolic_utils import anticommutator, commutator, u3_generators

FLAVOR_COUNT = 3
DIAGONAL_GENERATOR_INDICES = (3, 8, 0)
x = ALP_FIELD_RATIO


# ------------------------------------------------------------------------- frames


def linear_spurion_frame(vector, axial, scalar, pseudoscalar, quark_mass, condensate, gluon):
    """Sources of the linear-spurion convention: chi = chi_0 + i x([c^S, chi_0] + {c^P, chi_0})."""
    chi_0 = 2 * condensate * quark_mass
    return {
        "mass_spurion": [chi_0, sp.I * commutator(scalar, chi_0) + sp.I * anticommutator(pseudoscalar, chi_0),
                         sp.zeros(3)],
        "right_current": [vector + axial, sp.zeros(3)],
        "left_current": [vector - axial, sp.zeros(3)],
        "topological_slope": gluon,
        "quark_mass": quark_mass,
        "condensate": condensate,
    }


def second_order_group_element(generator):
    """exp(i x G) through x^2 (enough for every quadratic-order quantity below)."""
    return sp.eye(3) + sp.I * x * generator - x**2 * generator * generator / 2


def series_coefficients(matrix, order):
    """Power-series coefficients of a matrix polynomial in x, degrees 0..order."""
    expanded = matrix.applyfunc(sp.expand)
    return [expanded.applyfunc(lambda entry: entry.coeff(x, degree)) for degree in range(order + 1)]


def redefine(frame, vector_generator, axial_generator):
    """Apply g_R = exp(i x (G_V + G_A)), g_L = exp(i x (G_V - G_A)) to the sources of a frame."""
    right_generator = vector_generator + axial_generator
    left_generator = vector_generator - axial_generator
    right_element = second_order_group_element(right_generator)
    left_element = second_order_group_element(left_generator)
    chi = sum((coefficient * x**degree for degree, coefficient in enumerate(frame["mass_spurion"])), sp.zeros(3))
    new_chi = series_coefficients(right_element * chi * left_element.H, 2)

    def new_current(element, generator, coefficients):
        current = coefficients[0] + x * coefficients[1]
        coefficients_new = series_coefficients(element * current * element.H + generator, 1)
        return [coefficients_new[0], coefficients_new[1]]  # derivatives at x = 0 (degree <= 1)

    return dict(frame,
                mass_spurion=new_chi,
                right_current=new_current(right_element, right_generator, frame["right_current"]),
                left_current=new_current(left_element, left_generator, frame["left_current"]),
                topological_slope=frame["topological_slope"] - axial_generator.trace())


def dressing_redefinition(frame):
    """Constant generators G_eta = -C_eta(0) that remove the connections of a frame with constant C."""
    right, left = frame["right_current"][0], frame["left_current"][0]
    return -(right + left) / 2, -(right - left) / 2


# ------------------------------------------------------------------------ Hessian


def lagrangian_hessian(frame, f_pi, f_a, anomaly_mass, bare_mass, generator_indices=NEUTRAL_GENERATOR_INDICES):
    """K and M^2 on (P_i for i in generator_indices, a), from the LO Lagrangian of the frame."""
    generators = [u3_generators()[index] for index in generator_indices]
    count = len(generators)
    mesons = sp.symbols(f"P_0:{count}", real=True)
    meson_derivatives = sp.symbols(f"dP_0:{count}", real=True)
    alp, alp_derivative = sp.symbols("a da", real=True)
    fields, derivatives = [*mesons, alp], [*meson_derivatives, alp_derivative]
    axial_current = (frame["right_current"][0] - frame["left_current"][0]) / 2

    derivative_matrix = sum((d * g for d, g in zip(meson_derivatives, generators)), sp.zeros(3))
    vielbein = -derivative_matrix / f_pi + 2 * alp_derivative / f_a * axial_current
    kinetic_lagrangian = alp_derivative**2 / 2 + f_pi**2 / 4 * (vielbein * vielbein).trace()

    meson_matrix = sum((p * g for p, g in zip(mesons, generators)), sp.zeros(3))
    goldstone = sp.eye(3) + sp.I * meson_matrix / f_pi - meson_matrix * meson_matrix / (2 * f_pi**2)
    ratio = alp / f_a
    chi = sum((c * ratio**n for n, c in enumerate(frame["mass_spurion"])), sp.zeros(3))
    topological = -meson_matrix.trace() / f_pi - 2 * frame["topological_slope"] * ratio
    mass_lagrangian = (f_pi**2 / 4 * (chi * goldstone.H + goldstone * chi.H).trace()
                       - f_pi**2 * anomaly_mass / (4 * FLAVOR_COUNT) * topological**2
                       - bare_mass * alp**2 / 2)
    zero = {field: 0 for field in fields}
    size = count + 1
    kinetic = sp.Matrix(size, size, lambda r, c: sp.diff(kinetic_lagrangian, derivatives[r], derivatives[c]))
    mass_lagrangian = sp.expand(mass_lagrangian)
    mass_squared = sp.Matrix(size, size, lambda r, c: -sp.diff(mass_lagrangian, fields[r], fields[c]).subs(zero))
    return kinetic.applyfunc(sp.expand), mass_squared.applyfunc(sp.expand)


def schur_complement(mass_squared):
    """M^2_aa - M^2_aP (M^2_PP)^-1 M^2_Pa (ALP is the last coordinate)."""
    meson, mixing = mass_squared[:-1, :-1], mass_squared[:-1, -1]
    return mass_squared[-1, -1] - (mixing.T * meson.LUsolve(mixing))[0]


def meson_shift_from_generator(axial_generator, f_pi, f_a, generator_indices=NEUTRAL_GENERATOR_INDICES):
    """w_i = -(f_pi/f_a) tr(lambda_i G_A) (alp-pole-alignment-chiral-generator.md, item 2)."""
    generators = u3_generators()
    return sp.Matrix([-f_pi / f_a * (generators[index] * axial_generator).trace() for index in generator_indices])


# ----------------------------------------------------------- reading off the sources


def spurion_diagonal(frame):
    return tuple(frame["mass_spurion"][0][q, q] for q in range(3))


def pseudoscalar_from_spurion(frame):
    """c^P with {c^P, chi_0} = (chi_(1) - chi_(1)^dag)/(2i), chi_0 diagonal."""
    chi_values = spurion_diagonal(frame)
    first = frame["mass_spurion"][1]
    antihermitian = (first - first.H) / (2 * sp.I)
    return sp.Matrix(3, 3, lambda r, c: sp.expand(antihermitian[r, c] / (chi_values[r] + chi_values[c])))


def scalar_from_spurion(frame):
    """Off-diagonal c^S with i[c^S, chi_0] = (chi_(1) + chi_(1)^dag)/2 (diagonal entries not fixed)."""
    chi_values = spurion_diagonal(frame)
    first = frame["mass_spurion"][1]
    hermitian = (first + first.H) / 2
    return sp.Matrix(3, 3, lambda r, c: 0 if r == c else
                     sp.expand(hermitian[r, c] / (sp.I * (chi_values[c] - chi_values[r]))))


def real_trace(matrix):
    return sp.expand((matrix.trace() + sp.conjugate(matrix.trace())) / 2)


def frame_bound_coefficient(frame, anomaly_mass):
    """Bracket of item 2 for the given frame: Re tr chi_(2) + previous bracket(c^P, Theta_1)."""
    return sp.expand(real_trace(frame["mass_spurion"][2])
                     + bound_coefficient(pseudoscalar_from_spurion(frame), spurion_diagonal(frame), anomaly_mass,
                                         frame["topological_slope"]))


def previous_bracket_in_frame(frame, anomaly_mass):
    """The bracket of derive_alp_lo_stability_bound.py evaluated with the frame's c^P and Theta_1 (no a^2 term)."""
    return bound_coefficient(pseudoscalar_from_spurion(frame), spurion_diagonal(frame), anomaly_mass,
                             frame["topological_slope"])


# ------------------------------------------------------------------ dressed sources


def dressed_mass_derivatives(frame):
    """[Mhat_0, Mhat_1, Mhat_2] (derivatives at x = 0) from the transport equation."""
    right_transport = matrix_transport_series(frame["right_current"], order=2)
    left_transport = matrix_transport_series(frame["left_current"], order=2)
    mass = sum((c * x**n for n, c in enumerate(frame["mass_spurion"])), sp.zeros(3)) / (2 * frame["condensate"])
    dressed = truncate_matrix(right_transport.H * mass * left_transport, 2)
    return [m.applyfunc(sp.expand) for m in taylor_coefficients(dressed, order=2)]


def dressed_topological_slope(frame):
    """Thetahat_1 = Theta_1 + tr C_A(0)."""
    return sp.expand(frame["topological_slope"]
                     + (frame["right_current"][0] - frame["left_current"][0]).trace() / 2)


def axial_pseudoscalar_from_dressed_mass(dressed_first, quark_mass):
    """c^AP from Phat_1 = -{c^AP, m_q} (m_q diagonal)."""
    phat = (dressed_first - dressed_first.H) / (2 * sp.I)
    return sp.Matrix(3, 3, lambda r, c: sp.expand(-phat[r, c] / (quark_mass[r, r] + quark_mass[c, c])))


def topological_term(anomaly_mass, chi_values, slope):
    singlet = anomaly_mass / FLAVOR_COUNT
    return 2 * singlet * slope**2 / (1 + singlet * sum(1 / value for value in chi_values))


def invariant_bound_coefficient(dressed, quark_mass, condensate, anomaly_mass, dressed_slope):
    """B_0 Re tr Mhat_2 + 2 tr(c^AP chi_0 c^AP) - 2 (m_0^2/N_f)(c_GA - tr c^AP)^2/(1 + (m_0^2/N_f) sum 1/chi_q)."""
    chi_0 = 2 * condensate * quark_mass
    axial_pseudoscalar = axial_pseudoscalar_from_dressed_mass(dressed[1], quark_mass)
    return sp.expand(condensate * real_trace(dressed[2])
                     + 2 * (axial_pseudoscalar * chi_0 * axial_pseudoscalar).trace()
                     - topological_term(anomaly_mass, chi_0.diagonal(), dressed_slope - axial_pseudoscalar.trace()))


def invariant_bound_from_dressed_mass_only(dressed, quark_mass, condensate, anomaly_mass, dressed_slope):
    """Same bracket written with Phat_1 only: 2 B_0 sum |Phat_1,qq'|^2/(m_q+m_q'), tr(m^-1 Phat_1)/2."""
    phat = (dressed[1] - dressed[1].H) / (2 * sp.I)
    masses = [quark_mass[q, q] for q in range(3)]
    chi_values = [2 * condensate * m for m in masses]
    absolute_square_sum = sum(phat[r, c] * sp.conjugate(phat[r, c]) / (masses[r] + masses[c])
                              for r in range(3) for c in range(3))
    slope = dressed_slope + sum(phat[q, q] / masses[q] for q in range(3)) / 2
    return sp.expand(condensate * real_trace(dressed[2]) + 2 * condensate * absolute_square_sum
                     - topological_term(anomaly_mass, chi_values, slope))


def singular_value_second_derivative(mass_derivatives):
    """d^2/dx^2 tr (M^dag M)^(1/2) at x = 0 for M_0 diagonal positive (Daleckii-Krein formula).

    With H = M^dag M = M_0^2 + x H_1 + x^2 H_2 and f = sqrt: the x^2 coefficient of tr f(H) is
    tr(f'(M_0^2) H_2) + (1/2) sum_ij |H_1,ij|^2 (f'(m_i^2) - f'(m_j^2))/(m_i^2 - m_j^2),
    the divided difference being -1/(2 m_i m_j (m_i + m_j)) (also for i = j).
    """
    mass_zero, mass_one, mass_two = mass_derivatives
    masses = [mass_zero[q, q] for q in range(3)]
    first = mass_zero * mass_one + mass_one.H * mass_zero
    second = mass_one.H * mass_one + (mass_zero * mass_two + mass_two.H * mass_zero) / 2
    coefficient = (sum(second[q, q] / (2 * masses[q]) for q in range(3))
                   - sum(first[r, c] * sp.conjugate(first[r, c]) / (4 * masses[r] * masses[c] * (masses[r] + masses[c]))
                         for r in range(3) for c in range(3)))
    return sp.expand(2 * coefficient)


def vacuum_angle_slope(frame):
    """thetabar_1 = Theta_1 + (1/2) Im tr(m^-1 M_1)."""
    mass_zero = frame["mass_spurion"][0]
    log_derivative = (mass_zero.inv() * frame["mass_spurion"][1]).trace()
    return sp.expand(frame["topological_slope"] + (log_derivative - sp.conjugate(log_derivative)) / (4 * sp.I))


def singular_value_bound_coefficient(frame, anomaly_mass):
    """B_0 d^2/dx^2 tr|M| - 2 (m_0^2/N_f) thetabar_1^2/(1 + ...), using the frame's own M(x)."""
    condensate = frame["condensate"]
    mass = [frame["mass_spurion"][0] / (2 * condensate), frame["mass_spurion"][1] / (2 * condensate),
            2 * frame["mass_spurion"][2] / (2 * condensate)]
    return sp.expand(condensate * singular_value_second_derivative(mass)
                     - topological_term(anomaly_mass, spurion_diagonal(frame), vacuum_angle_slope(frame)))


# ---------------------------------------------------------------- symbolic checks


def schur_complement_shift_invariance():
    """Generic symmetric 4x4 M^2 (3 mesons + a), symbolic shift w: S(T^T M^2 T) = S(M^2)."""
    entries = sp.symbols("M_0:10", real=True)
    index = iter(entries)
    matrix = sp.zeros(4)
    for r in range(4):
        for c in range(r, 4):
            matrix[r, c] = matrix[c, r] = next(index)
    conversion = shift_matrix(sp.Matrix(sp.symbols("w_0:3", real=True)))
    return sp.simplify(schur_complement(conversion.T * matrix * conversion) - schur_complement(matrix)) == 0


SYMBOLIC_PARAMETERS = {
    "quark_mass": sp.diag(*sp.symbols("m_u m_d m_s", positive=True)),
    "condensate": sp.Symbol("B_0", positive=True),
    "gluon": sp.Symbol("c_G", real=True),
    "f_pi": sp.Symbol("f_pi", positive=True),
    "f_a": sp.Symbol("f_a", positive=True),
    "anomaly_mass": sp.Symbol("m_0_sq", positive=True),
    "bare_mass": sp.Symbol("m_a0_sq", real=True),
}


def symbolic_linear_frame():
    p = SYMBOLIC_PARAMETERS
    sources = {name: symbolic_hermitian(name) for name in ("cV", "cA", "cS", "cP")}
    frame = linear_spurion_frame(sources["cV"], sources["cA"], sources["cS"], sources["cP"], p["quark_mass"],
                                 p["condensate"], p["gluon"])
    return frame, sources


def general_frame_reduction():
    """Item 2 (symbolic, neutral generators, complex d-s entries, generic complex chi_(2))."""
    p = SYMBOLIC_PARAMETERS
    frame, sources = symbolic_linear_frame()
    second = sp.Matrix(3, 3, lambda r, c: sp.Symbol(f"X_{r}{c}_re", real=True) + sp.I * sp.Symbol(f"X_{r}{c}_im", real=True))
    frame = dict(frame, mass_spurion=[frame["mass_spurion"][0], frame["mass_spurion"][1], second])
    kinetic, mass_squared = lagrangian_hessian(frame, p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"])
    closed_kinetic, closed_mass = closed_form_blocks(sources["cA"], sources["cP"], frame["mass_spurion"][0], p["gluon"],
                                                     p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"])
    epsilon_squared = p["f_pi"]**2 / p["f_a"]**2
    expected = closed_mass.copy()
    expected[-1, -1] -= epsilon_squared * real_trace(second)
    return {
        "K from the Lagrangian depends only on C_A(0) (closed form)":
            (kinetic - closed_kinetic).applyfunc(sp.simplify) == sp.zeros(6),
        "M^2 = linear-spurion blocks(c^P, Theta_1) - eps^2 Re tr chi_(2) e_a e_a^T; c^S does not enter":
            (mass_squared - expected).applyfunc(sp.simplify) == sp.zeros(6),
        "dM^2/dm_a0^2 = e_a e_a^T": mass_squared.diff(p["bare_mass"]) == sp.diag(0, 0, 0, 0, 0, 1),
        "c^P read off from chi_(1) equals the input c^P":
            (pseudoscalar_from_spurion(frame) - sources["cP"]).applyfunc(sp.simplify) == sp.zeros(3),
    }


def symbolic_redefinition():
    """Items 1 and 3 with symbolic sources and symbolic noncommuting G_V, G_A (neutral, d-s entries)."""
    p = SYMBOLIC_PARAMETERS
    frame, sources = symbolic_linear_frame()
    vector_generator, axial_generator = symbolic_hermitian("GV"), symbolic_hermitian("GA")
    new = redefine(frame, vector_generator, axial_generator)
    kinetic, mass_squared = lagrangian_hessian(frame, p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"])
    kinetic_new, mass_new = lagrangian_hessian(new, p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"])
    conversion = shift_matrix(meson_shift_from_generator(axial_generator, p["f_pi"], p["f_a"]))
    chi_0 = frame["mass_spurion"][0]
    pseudoscalar, pseudoscalar_new = sources["cP"], pseudoscalar_from_spurion(new)
    axial_new = (new["right_current"][0] - new["left_current"][0]) / 2
    vector_new = (new["right_current"][0] + new["left_current"][0]) / 2
    zero = sp.zeros(3)

    def vanishes(matrix):
        return matrix.applyfunc(sp.simplify) == sp.zeros(*matrix.shape)

    bracket_difference = sp.expand(previous_bracket_in_frame(new, p["anomaly_mass"])
                                   - previous_bracket_in_frame(frame, p["anomaly_mass"]))
    expected_difference = sp.expand(2 * (axial_generator * chi_0 * axial_generator).trace()
                                    + 2 * (axial_generator * anticommutator(pseudoscalar, chi_0)).trace())
    checks = {
        "K' = T^T K T (w_i = -(f_pi/f_a) tr(lambda_i G_A))": vanishes(kinetic_new - conversion.T * kinetic * conversion),
        "M^2' = T^T M^2 T including the a^2 term of chi'": vanishes(mass_new - conversion.T * mass_squared * conversion),
        "chi'_(2) != 0 (the redefinition generates an a^2 spurion term)": not vanishes(new["mass_spurion"][2]),
        "dM^2'/dm_a0^2 = e_a e_a^T (m_a0^2 untouched)": mass_new.diff(p["bare_mass"]) == sp.diag(0, 0, 0, 0, 0, 1),
        "c^P' = c^P + G_A": vanishes(pseudoscalar_new - pseudoscalar - axial_generator),
        "c^S' = c^S + G_V (off-diagonal part)": vanishes(scalar_from_spurion(new) - scalar_from_spurion(frame) - (
            vector_generator - sp.diag(*vector_generator.diagonal()))),
        "c^A' = c^A + G_A, c^V' = c^V + G_V": vanishes(axial_new - sources["cA"] - axial_generator)
                                              and vanishes(vector_new - sources["cV"] - vector_generator),
        "c_G' = c_G - tr G_A": sp.expand(new["topological_slope"] - p["gluon"] + axial_generator.trace()) == 0,
        "c^A' - c^P' = c^A - c^P (c^AP invariant)": vanishes(axial_new - pseudoscalar_new - sources["cA"] + pseudoscalar),
        "tr c^P' + c_G' = tr c^P + c_G": sp.expand(pseudoscalar_new.trace() + new["topological_slope"]
                                                  - pseudoscalar.trace() - p["gluon"]) == 0,
        "c_G' + tr c^A' = c_G + tr c^A (c_GA invariant)": sp.expand(new["topological_slope"] + axial_new.trace()
                                                                   - p["gluon"] - sources["cA"].trace()) == 0,
        "Re tr chi'_(2) = 2 tr(c^P chi_0 c^P) - 2 tr(c^P' chi_0 c^P')": sp.simplify(
            real_trace(new["mass_spurion"][2]) - 2 * (pseudoscalar * chi_0 * pseudoscalar).trace()
            + 2 * (pseudoscalar_new * chi_0 * pseudoscalar_new).trace()) == 0,
        "previous bracket(c^P', c_G') - previous bracket(c^P, c_G) = 2 tr(G_A chi_0 G_A) + 2 tr(G_A {c^P, chi_0})":
            sp.simplify(bracket_difference - expected_difference) == 0,
        "frame bracket (with Re tr chi'_(2)) is unchanged": sp.simplify(
            frame_bound_coefficient(new, p["anomaly_mass"]) - frame_bound_coefficient(frame, p["anomaly_mass"])) == 0,
        "G_V alone leaves the previous bracket unchanged": sp.simplify(
            previous_bracket_in_frame(redefine(frame, vector_generator, zero), p["anomaly_mass"])
            - previous_bracket_in_frame(frame, p["anomaly_mass"])) == 0,
    }
    return {"checks": checks, "bracket difference": expected_difference}


def dressed_mass_second_order_linear_convention():
    """Item 5: Mhat_2 in the linear convention, noncommuting symbols, with C_eta,1 kept.

    Expected: Mhat_2 = 2 (C_R Mt - Mt C_L) - (C_R^2 m + m C_L^2) + 2 C_R m C_L - i (C_R1 m - m C_L1),
    Mt = [c^S, m] + {c^P, m}, C_R = c^V + c^A, C_L = c^V - c^A, C_eta1 = dC_eta/dx(0).
    """
    from derive_alp_invariant_coefficient_functions import transport_series
    names = ("m", "cV", "cA", "cS", "cP", "CR1", "CL1")
    m, c_v, c_a, c_s, c_p, right_one, left_one = (sp.Symbol(n, commutative=False, hermitian=True) for n in names)
    right, left = c_v + c_a, c_v - c_a
    right_transport = transport_series([right, right_one], x, 2)
    left_transport = transport_series([left, left_one], x, 2)
    mass = m + sp.I * x * (c_s * m - m * c_s + c_p * m + m * c_p)
    dressed = sp.expand(sp.adjoint(right_transport) * mass * left_transport)
    second = sp.expand(sp.diff(dressed, x, 2).subs(x, 0))
    mixed = c_s * m - m * c_s + c_p * m + m * c_p
    expected = sp.expand(2 * (right * mixed - mixed * left) - (right * right * m + m * left * left)
                         + 2 * right * m * left - sp.I * (right_one * m - m * left_one))
    return sp.expand(second - expected) == 0


def dressed_mass_real_trace_linear_convention():
    """Item 5 with explicit symbolic Hermitian 3x3 matrices (all entries, incl. charged)."""
    m = SYMBOLIC_PARAMETERS["quark_mass"]
    matrices = {name: symbolic_hermitian(name, charge_conserving=False) for name in ("cV", "cA", "cS", "cP", "CR1", "CL1")}
    frame = linear_spurion_frame(matrices["cV"], matrices["cA"], matrices["cS"], matrices["cP"], m,
                                 SYMBOLIC_PARAMETERS["condensate"], SYMBOLIC_PARAMETERS["gluon"])
    frame["right_current"][1] = matrices["CR1"]
    frame["left_current"][1] = matrices["CL1"]
    dressed = dressed_mass_derivatives(frame)
    axial_pseudoscalar = matrices["cA"] - matrices["cP"]
    expected = 4 * (matrices["cP"] * m * matrices["cP"]).trace() - 4 * (axial_pseudoscalar * m * axial_pseudoscalar).trace()
    real = real_trace(dressed[2])
    return {
        "Re tr Mhat_2 = 4 tr(c^P m c^P) - 4 tr(c^AP m c^AP)": sp.simplify(real - expected) == 0,
        "Re tr Mhat_2 independent of c^V, c^S, C_eta,1": not (real.free_symbols & set().union(
            *(matrices[n].free_symbols for n in ("cV", "cS", "CR1", "CL1")))),
        "Mhat_1 = -i[c^VS, m] - i{c^AP, m}": (dressed[1] + sp.I * commutator(matrices["cV"] - matrices["cS"], m)
                                             + sp.I * anticommutator(axial_pseudoscalar, m)).applyfunc(sp.simplify)
                                            == sp.zeros(3),
    }


def invariant_form_checks():
    """Items 4-5 symbolically (neutral sources with d-s entries)."""
    p = SYMBOLIC_PARAMETERS
    frame, sources = symbolic_linear_frame()
    dressed = dressed_mass_derivatives(frame)
    slope = dressed_topological_slope(frame)
    invariant = invariant_bound_coefficient(dressed, p["quark_mass"], p["condensate"], p["anomaly_mass"], slope)
    previous = bound_coefficient(sources["cP"], spurion_diagonal(frame), p["anomaly_mass"], p["gluon"])
    only_dressed = invariant_bound_from_dressed_mass_only(dressed, p["quark_mass"], p["condensate"], p["anomaly_mass"], slope)
    singular = singular_value_bound_coefficient(frame, p["anomaly_mass"])

    # Generic dressed data: Mhat_1 from symbolic c^VS, c^AP, generic complex Mhat_2.
    m = p["quark_mass"]
    vector_scalar, axial_pseudoscalar = symbolic_hermitian("cVS"), symbolic_hermitian("cAP")
    generic_second = sp.Matrix(3, 3, lambda r, c: sp.Symbol(f"Y_{r}{c}_re", real=True) + sp.I * sp.Symbol(f"Y_{r}{c}_im", real=True))
    generic_first = -sp.I * commutator(vector_scalar, m) - sp.I * anticommutator(axial_pseudoscalar, m)
    generic = [m, generic_first, generic_second]
    gluon_axial = sp.Symbol("c_GA", real=True)
    generic_invariant = invariant_bound_coefficient(generic, m, p["condensate"], p["anomaly_mass"], gluon_axial)
    dressed_frame = {"mass_spurion": [2 * p["condensate"] * m, 2 * p["condensate"] * generic_first,
                                      p["condensate"] * generic_second],
                     "condensate": p["condensate"], "topological_slope": gluon_axial}
    generic_singular = singular_value_bound_coefficient(dressed_frame, p["anomaly_mass"])
    chi = 2 * p["condensate"] * m
    ds_entry = sp.Symbol("cAP_ds_re", real=True)**2 + sp.Symbol("cAP_ds_im", real=True)**2
    diagonal_part = axial_pseudoscalar.copy()
    diagonal_part[1, 2] = diagonal_part[2, 1] = 0
    return {
        "linear convention: invariant form = previous bracket": sp.simplify(invariant - previous) == 0,
        "linear convention: Thetahat_1 - tr c^AP = c_G + tr c^P": sp.expand(
            slope - (sources["cA"] - sources["cP"]).trace() - p["gluon"] - sources["cP"].trace()) == 0,
        "Phat_1-only form equals the c^AP form": sp.simplify(only_dressed - invariant) == 0,
        "singular-value form equals the c^AP form (linear frame)": sp.simplify(singular - invariant) == 0,
        "generic Mhat_1, Mhat_2: singular-value form = B_0 Re tr Mhat_2 + 2 tr(c^AP chi c^AP) - top. (c^VS drops out)":
            sp.simplify(generic_singular - generic_invariant) == 0,
        "generic: c^VS does not enter": not (generic_invariant.free_symbols & vector_scalar.free_symbols),
        "d-s entry of c^AP adds 2 (chi_d + chi_s)|c^AP_ds|^2": sp.simplify(
            2 * (axial_pseudoscalar * chi * axial_pseudoscalar).trace()
            - 2 * (diagonal_part * chi * diagonal_part).trace() - 2 * (chi[1, 1] + chi[2, 2]) * ds_entry) == 0,
    }


def exact_to_mpmath(entry):
    """Exact rational complex number to an mpmath number at the working precision."""
    real, imaginary = (sp.Rational(part) for part in sp.expand(entry).as_real_imag())
    return mp.mpc(mp.mpf(real.p) / real.q, mp.mpf(imaginary.p) / imaginary.q)


def singular_value_formula_numeric_check(precision=40, step=sp.Rational(1, 10**8)):
    """mpmath check of the Daleckii-Krein second derivative of tr (M^dag M)^(1/2) (complex d-s entries)."""
    mass_zero = sp.diag(sp.Rational(3, 1000), sp.Rational(1, 200), sp.Rational(9, 100))
    mass_one = sp.Matrix([[sp.Rational(1, 300) * sp.I, 0, 0], [0, sp.Rational(-1, 100), sp.Rational(1, 50) + sp.I / 70],
                          [0, sp.Rational(1, 90) - sp.I / 40, sp.Rational(1, 30) + sp.I / 20]])
    mass_two = sp.Matrix([[sp.Rational(1, 70), 0, sp.I / 90], [0, sp.Rational(-1, 30) * sp.I, sp.Rational(1, 80)],
                          [sp.Rational(1, 60), 0, sp.Rational(-1, 25) + sp.I / 30]])
    formula = singular_value_second_derivative([mass_zero, mass_one, mass_two])
    mp.dps = precision

    def trace_modulus(value):
        matrix = mass_zero + value * mass_one + value**2 * mass_two / 2
        numeric = mp.matrix([[exact_to_mpmath(entry) for entry in row] for row in matrix.tolist()])
        return mp.re(sum(mp.sqrtm(numeric.H * numeric)[q, q] for q in range(3)))

    h = mp.mpf(step.p) / step.q
    finite_difference = (trace_modulus(step) - 2 * trace_modulus(0) + trace_modulus(-step)) / h**2
    return float(sp.re(formula)), float(finite_difference)


# ------------------------------------------------------------------ exact rational frames


def rational_hermitian(entries):
    diagonal, ds_entry = entries
    matrix = sp.diag(*diagonal)
    matrix[1, 2], matrix[2, 1] = ds_entry, sp.conjugate(ds_entry)
    return matrix


RATIONAL_PARAMETERS = {
    "quark_mass": sp.diag(sp.Rational(1, 400), sp.Rational(1, 200), sp.Rational(1, 10)),
    "condensate": sp.Rational(5, 2),
    "gluon": sp.Rational(1, 3),
    "f_pi": sp.Rational(23, 250),
    "f_a": sp.Rational(3, 10),
    "anomaly_mass": sp.Rational(4, 5),
    "bare_mass": sp.Symbol("m_a0_sq", real=True),
}
RATIONAL_SOURCES = {
    "vector": rational_hermitian(((sp.Rational(1, 5), sp.Rational(-1, 3), sp.Rational(1, 2)), sp.Rational(1, 4) + sp.I / 3)),
    "axial": rational_hermitian(((sp.Rational(1, 2), sp.Rational(-1, 4), sp.Rational(2, 3)), sp.Rational(1, 3) + sp.I / 4)),
    "scalar": rational_hermitian(((sp.Rational(-1, 6), sp.Rational(1, 7), 0), sp.Rational(-1, 5) + sp.I / 6)),
    "pseudoscalar": rational_hermitian(((sp.Rational(-1, 3), sp.Rational(1, 2), sp.Rational(1, 7)), sp.Rational(1, 2) - sp.I / 3)),
}
RATIONAL_GENERATORS = {  # noncommuting with each other and with the sources
    "vector": rational_hermitian(((sp.Rational(2, 5), sp.Rational(-1, 2), sp.Rational(1, 9)), sp.Rational(-1, 3) + sp.I / 5)),
    "axial": rational_hermitian(((sp.Rational(-3, 4), sp.Rational(1, 6), sp.Rational(5, 4)), sp.Rational(2, 7) - sp.I / 2)),
}


def rational_frames():
    """Linear-spurion frame, a random-generator frame and the frame without kinetic mixing (same physics)."""
    p, s, g = RATIONAL_PARAMETERS, RATIONAL_SOURCES, RATIONAL_GENERATORS
    linear = linear_spurion_frame(s["vector"], s["axial"], s["scalar"], s["pseudoscalar"], p["quark_mass"],
                                  p["condensate"], p["gluon"])
    rotated = redefine(linear, g["vector"], g["axial"])
    dressed = redefine(linear, *dressing_redefinition(linear))
    return {"linear spurion": linear, "noncommuting G_V, G_A": rotated, "no kinetic mixing (dressed)": dressed}


def rational_frame_checks():
    """Items 1-6 at an exact rational point, on the five neutral mesons and on all nine."""
    p = RATIONAL_PARAMETERS
    frames = rational_frames()
    linear = frames["linear spurion"]
    epsilon_squared = p["f_pi"]**2 / p["f_a"]**2
    reference_dressed = dressed_mass_derivatives(linear)
    reference_slope = dressed_topological_slope(linear)
    reference_bound = invariant_bound_coefficient(reference_dressed, p["quark_mass"], p["condensate"],
                                                  p["anomaly_mass"], reference_slope)
    reference_previous = previous_bracket_in_frame(linear, p["anomaly_mass"])
    rows, checks = {}, {}
    for label, frame in frames.items():
        _, mass_squared = lagrangian_hessian(frame, p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"])
        _, mass_all = lagrangian_hessian(frame, p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"],
                                         ALL_GENERATOR_INDICES)
        schur_bound = sp.nsimplify(sp.expand(-(schur_complement(mass_squared) - p["bare_mass"]) / epsilon_squared))
        schur_bound_all = sp.expand(-(schur_complement(mass_all) - p["bare_mass"]) / epsilon_squared)
        dressed = dressed_mass_derivatives(frame)
        previous = previous_bracket_in_frame(frame, p["anomaly_mass"])
        rows[label] = {"Schur bound / eps^2": sp.simplify(schur_bound), "previous bracket with this frame's c^P, Theta_1":
                       sp.simplify(previous)}
        checks[f"[{label}] Mhat_0..2 equal to the linear frame"] = all(
            (d - r).applyfunc(sp.simplify) == sp.zeros(3) for d, r in zip(dressed, reference_dressed))
        checks[f"[{label}] Thetahat_1 equal to the linear frame"] = sp.simplify(dressed_topological_slope(frame)
                                                                                - reference_slope) == 0
        checks[f"[{label}] -(S - m_a0^2)/eps^2 = invariant bound (5 neutral mesons)"] = sp.simplify(
            schur_bound - reference_bound) == 0
        checks[f"[{label}] same with all nine mesons"] = sp.simplify(schur_bound_all - reference_bound) == 0
        checks[f"[{label}] frame bracket Re tr chi_(2) + bracket(c^P, Theta_1) = invariant bound"] = sp.simplify(
            frame_bound_coefficient(frame, p["anomaly_mass"]) - reference_bound) == 0
        checks[f"[{label}] singular-value form = invariant bound"] = sp.simplify(
            singular_value_bound_coefficient(frame, p["anomaly_mass"]) - reference_bound) == 0
    rotated, dressed_frame = frames["noncommuting G_V, G_A"], frames["no kinetic mixing (dressed)"]
    checks["linear frame: invariant bound = previous bracket"] = sp.simplify(reference_bound - reference_previous) == 0
    checks["G_V, G_A frame: previous bracket(c^P', c_G') != invariant bound"] = sp.simplify(
        previous_bracket_in_frame(rotated, p["anomaly_mass"]) - reference_bound) != 0
    kinetic_dressed, _ = lagrangian_hessian(dressed_frame, p["f_pi"], p["f_a"], p["anomaly_mass"], p["bare_mass"])
    checks["dressed frame: C_R' = C_L' = 0 (all orders here: C' = 0 at x^0 and x^1)"] = all(
        c == sp.zeros(3) for c in dressed_frame["right_current"] + dressed_frame["left_current"])
    checks["dressed frame: K = 1"] = kinetic_dressed == sp.eye(6)
    checks["dressed frame: chi' = 2 B_0 Mhat through x^2"] = all(
        (dressed_frame["mass_spurion"][n] - 2 * p["condensate"] * reference_dressed[n] / sp.factorial(n))
        .applyfunc(sp.simplify) == sp.zeros(3) for n in range(3))
    checks["dressed frame: c^P' = -c^AP, Theta_1' = c_GA"] = (
        (pseudoscalar_from_spurion(dressed_frame) + RATIONAL_SOURCES["axial"] - RATIONAL_SOURCES["pseudoscalar"])
        .applyfunc(sp.simplify) == sp.zeros(3)
        and sp.simplify(dressed_frame["topological_slope"] - p["gluon"] - RATIONAL_SOURCES["axial"].trace()) == 0)

    # Consistency of lagrangian_hessian with the 10x10 helper of the pole-alignment script.
    data = {"f_pi": p["f_pi"], "f_a": p["f_a"], "m_0_sq": p["anomaly_mass"], "m_a0_sq": sp.Rational(1, 4)}
    axial_rotated = (rotated["right_current"][0] - rotated["left_current"][0]) / 2
    helper = quadratic_hessian(data, axial_rotated, rotated["mass_spurion"], rotated["topological_slope"])
    own = lagrangian_hessian(rotated, p["f_pi"], p["f_a"], p["anomaly_mass"], sp.Rational(1, 4), ALL_GENERATOR_INDICES)
    checks["lagrangian_hessian equals quadratic_hessian of the pole-alignment script (10x10)"] = all(
        (a - b).applyfunc(sp.simplify) == sp.zeros(10) for a, b in zip(helper, own))
    return {"checks": checks, "rows": rows, "invariant bound": sp.simplify(reference_bound)}


# --------------------------------------------------------------- numerical checks


def numeric_mass_matrix(frame, bare_mass):
    p = RATIONAL_PARAMETERS
    _, mass_squared = lagrangian_hessian(frame, p["f_pi"], p["f_a"], p["anomaly_mass"], bare_mass)
    return np.array(mass_squared.evalf(30), dtype=complex).real


def eigenvalue_sign_checks(relative_offset=1e-6):
    """Smallest eigenvalue of M^2 at m_a0^2 = eps^2 (bound) (1 +- offset), in three frames; also the
    count of negative eigenvalues of M^2 equals that of the linear frame (congruence)."""
    p = RATIONAL_PARAMETERS
    epsilon_squared = p["f_pi"]**2 / p["f_a"]**2
    frames = rational_frames()
    linear = frames["linear spurion"]
    invariant = invariant_bound_coefficient(dressed_mass_derivatives(linear), p["quark_mass"], p["condensate"],
                                            p["anomaly_mass"], dressed_topological_slope(linear))
    threshold = epsilon_squared * invariant
    offset = relative_offset * (abs(threshold) + 1)
    results = {}
    for label, frame in frames.items():
        above = np.linalg.eigvalsh(numeric_mass_matrix(frame, threshold + offset))[0]
        below = np.linalg.eigvalsh(numeric_mass_matrix(frame, threshold - offset))[0]
        wrong = epsilon_squared * previous_bracket_in_frame(frame, p["anomaly_mass"])
        wrong_above = np.linalg.eigvalsh(numeric_mass_matrix(frame, wrong + offset))[0]
        wrong_below = np.linalg.eigvalsh(numeric_mass_matrix(frame, wrong - offset))[0]
        results[label] = {"threshold": float(threshold), "above": above, "below": below,
                          "previous-bracket threshold": float(wrong), "at previous +": wrong_above,
                          "at previous -": wrong_below}
    return results


def vacuum_energy(frame_numeric, ratio, f_pi, anomaly_mass, start=None):
    """min over the nine meson coordinates of V(U, x) = -(f_pi^2/4)<chi U^dag + U chi^dag>
    + f_pi^2 m_0^2/(4 N_f) (-<Phi>/f_pi - 2 Theta(x))^2, with exact exponentials in chi(x)."""
    chi_function, theta_function = frame_numeric
    chi = chi_function(ratio)
    theta = theta_function(ratio)
    generators = [np.array(g, dtype=complex) for g in u3_generators()]

    def potential(coordinates):
        meson = sum(c * g for c, g in zip(coordinates, generators))
        goldstone = expm(1j * meson / f_pi)
        topological = -np.trace(meson).real / f_pi - 2 * theta
        return (-(f_pi**2 / 4) * np.trace(chi @ goldstone.conj().T + goldstone @ chi.conj().T).real
                + f_pi**2 * anomaly_mass / (4 * FLAVOR_COUNT) * topological**2)

    result = minimize(potential, np.zeros(9) if start is None else start, method="BFGS",
                      options={"gtol": 1e-14, "maxiter": 10000})
    return result.fun, result.x


def nonlinear_frames():
    """Exact chi(x) and Theta(x) of the linear frame and of the G_V, G_A frame (no truncation in x)."""
    p, s, g = RATIONAL_PARAMETERS, RATIONAL_SOURCES, RATIONAL_GENERATORS
    to_array = lambda matrix: np.array(matrix.evalf(30), dtype=complex)
    chi_0 = to_array(2 * p["condensate"] * p["quark_mass"])
    first = to_array(sp.I * commutator(s["scalar"], 2 * p["condensate"] * p["quark_mass"])
                     + sp.I * anticommutator(s["pseudoscalar"], 2 * p["condensate"] * p["quark_mass"]))
    gluon = float(p["gluon"])
    right, left = to_array(g["vector"] + g["axial"]), to_array(g["vector"] - g["axial"])
    axial_trace = float(g["axial"].trace())
    linear = (lambda ratio: chi_0 + ratio * first, lambda ratio: gluon * ratio)
    rotated = (lambda ratio: expm(1j * ratio * right) @ (chi_0 + ratio * first) @ expm(1j * ratio * left).conj().T,
               lambda ratio: (gluon - axial_trace) * ratio)
    return linear, rotated


def nonlinear_vacuum_energy_checks(step=2e-2, ratios=(0.1, 0.3)):
    """E(x) = min_U V equal in the two frames at finite x; E''(0)/f_a^2 = -eps^2 (invariant bound)."""
    p = RATIONAL_PARAMETERS
    f_pi, f_a, anomaly_mass = float(p["f_pi"]), float(p["f_a"]), float(p["anomaly_mass"])
    linear, rotated = nonlinear_frames()
    finite = [(ratio, vacuum_energy(linear, ratio, f_pi, anomaly_mass)[0],
               vacuum_energy(rotated, ratio, f_pi, anomaly_mass)[0]) for ratio in ratios]
    energies = {}
    for label, frame in (("linear", linear), ("rotated", rotated)):
        values = [vacuum_energy(frame, ratio, f_pi, anomaly_mass)[0] for ratio in (-2 * step, -step, 0.0, step, 2 * step)]
        # fourth-order central difference for the second derivative in x
        second = (-values[0] + 16 * values[1] - 30 * values[2] + 16 * values[3] - values[4]) / (12 * step**2)
        energies[label] = second / f_a**2  # d^2 E / da^2
    frame = rational_frames()["linear spurion"]
    invariant = float(invariant_bound_coefficient(dressed_mass_derivatives(frame), p["quark_mass"], p["condensate"],
                                                  p["anomaly_mass"], dressed_topological_slope(frame)))
    expected = -float(p["f_pi"]**2 / p["f_a"]**2) * invariant
    return {"finite x": finite, "second derivative": energies, "expected second derivative": expected}


# ---------------------------------------------------------------------------- main


def report(label, passed):
    print(f"  [{'pass' if passed else 'FAIL'}] {label}")


def main():
    print("1. Schur complement and congruence")
    report("S(T^T M^2 T) = S(M^2) for the shift P = P' + w a (generic symbolic)", schur_complement_shift_invariance())
    print("2. Any frame: M^2 from the Lagrangian (symbolic, neutral sources, complex d-s entries)")
    for label, passed in general_frame_reduction().items():
        report(label, passed)
    print("3. Constant-generator redefinition, symbolic sources and noncommuting G_V, G_A")
    redefinition = symbolic_redefinition()
    for label, passed in redefinition["checks"].items():
        report(label, passed)
    print(f"  previous bracket(c^P', c_G') - previous bracket(c^P, c_G) = {redefinition['bracket difference']}")
    print("4-5. Dressed sources and the invariant form")
    report("Mhat_2 (linear convention) = 2(C_R Mt - Mt C_L) - (C_R^2 m + m C_L^2) + 2 C_R m C_L - i(C_R1 m - m C_L1)",
           dressed_mass_second_order_linear_convention())
    for label, passed in dressed_mass_real_trace_linear_convention().items():
        report(label, passed)
    for label, passed in invariant_form_checks().items():
        report(label, passed)
    formula, finite_difference = singular_value_formula_numeric_check()
    report(f"d^2 tr|M|/dx^2: formula {formula:.15g}, mpmath finite difference {finite_difference:.15g}",
           abs(formula - finite_difference) < 1e-9 * abs(formula))
    print("6. Exact rational frames with the same physics")
    frames = rational_frame_checks()
    for label, passed in frames["checks"].items():
        report(label, passed)
    print(f"  invariant bound / eps^2 = {frames['invariant bound']} GeV^2")
    for label, row in frames["rows"].items():
        print(f"  {label}: " + ", ".join(f"{key} = {value}" for key, value in row.items()))
    print("Numerical: smallest eigenvalue of M^2 (GeV^2) at m_a0^2 = threshold (1 +- 1e-6)")
    for label, row in eigenvalue_sign_checks().items():
        report(f"{label}: invariant threshold {row['threshold']:.8g}: above {row['above']:.3g}, below {row['below']:.3g}",
               row["above"] > 0 > row["below"])
        report(f"{label}: previous bracket in this frame {row['previous-bracket threshold']:.8g} is a threshold only if "
               f"it equals the invariant one: above {row['at previous +']:.3g}, below {row['at previous -']:.3g}",
               (row["at previous +"] > 0 > row["at previous -"])
               == np.isclose(row["previous-bracket threshold"], row["threshold"], rtol=1e-12, atol=0))
    vacuum = nonlinear_vacuum_energy_checks()
    print("Numerical: minimum of the full nonlinear chiral potential over U(3)")
    for ratio, linear, rotated in vacuum["finite x"]:
        report(f"x = {ratio}: E_linear = {linear:.15g}, E_rotated = {rotated:.15g}", abs(linear - rotated) < 1e-12)
    for label, value in vacuum["second derivative"].items():
        report(f"{label}: d^2E/da^2 = {value:.10g}, expected -eps^2 (bound) = {vacuum['expected second derivative']:.10g}",
               abs(value - vacuum["expected second derivative"]) < 1e-5 * abs(vacuum["expected second derivative"]))


if __name__ == "__main__":
    main()
