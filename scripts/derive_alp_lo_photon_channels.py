#!/usr/bin/env python3
"""LO photon channels a -> gamma gamma and a -> gamma P_1 P_2 at finite ALP-meson mixing.

Run:   .venv/bin/python scripts/derive_alp_lo_photon_channels.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_photon_channels.py

Conventions (Alda:2026nrd). Natural units, metric (+---), epsilon^{0123} = +1,
U = exp(i Phi/f_pi), Phi = sum_i P_i lambda_i with tr(lambda_i lambda_j) =
2 delta_ij (U(3) basis of scripts/symbolic_utils.py), N_f = N_c = 3,
Q = diag(2/3, -1/3, -1/3), epsilon = f_pi/f_a (variable
`decay_constant_ratio`, to avoid a clash with polarization vectors).
Units: f_pi, f_a, B_0, m_q in GeV; squared masses, s, m_0^2 in GeV^2; all
c coefficients and Taylor coefficients of the partonic functions are
dimensionless. All numerical inputs are illustrative rationals, not fits.

Mass operator on the neutral flavour space (Hermitian matrices commuting with Q):
    M_s(X) = {B_0 m_q, X} + (m_0^2/N_f) tr(X) 1 - s X,
and M_s^{-1}(X) is the matrix Y with M_s(Y) = X (memory/alp/shared/
alp-lo-exact-pole-equations.md). Its tadpole argument is written with the
first-order dressed invariants (memory/alp/shared/
alp-dressed-topological-invariant.md, item 3; alp-dressed-mass-first-order-
invariants.md):
    B_0 hatP_1 - (m_0^2/N_f) hatTheta_1 1,
    hatP_1 = anti-Hermitian part of hatM_1, hatM_1 = M_1 - i(C_{R,0} m_q - m_q C_{L,0}),
    hatTheta_1 = Theta_1 + tr C_{A,0} (= c_GA),
and it equals -{B_0 m_q, c_q^{AP}} - (m_0^2/N_f) c_GA 1 (checked). The photon
input is hatTheta_{gamma,1} = Theta_{gamma,1} + 2 N_c tr(Q^2 C_{A,0})
= c_gamma + 2 N_c tr(Q^2 c_q^A) (memory/alp/photon/alp-dressed-photon-invariant.md).

Final results (no expansion in epsilon; m_a, m_1, m_2 physical pole masses):

 a -> gamma gamma (Lagrangian (alpha_em/4pi)[c_gamma a/f_a + (N_c/f_pi) tr(Q^2 Phi)] F Ftilde,
 memory/alp/photon/alp-photon-invariant-from-mixing.md):
    M = -e^2/(4 pi^2 f_a) E
        [hatTheta_{gamma,1} + 2 N_c tr(Q^2 M_{m_a^2}^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1))]
        / sqrt(1 + 2 epsilon^2 tr[(M_{m_a^2}^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1))^2]),
    E = epsilon_{mu nu rho sigma} q_1^mu eps_1^nu q_2^rho eps_2^sigma.

 a -> gamma P_1 P_2 (leading WZW gamma Phi^3 vertex,
 memory/alp/photon-two-pseudoscalar/wzw-alp-photon-pp-amplitude.md):
    M = epsilon_{mu nu rho sigma} eps_gamma^mu p_1^nu p_2^rho k^sigma A,
    A = -orientation N_c e/(24 pi^2 f_pi^3) F(Q; X_a, X_1, X_2),  orientation = -1,
    F(Q;A,B,C) = tr(QABC) - tr(QACB) - tr(QBAC) + tr(QBCA) + tr(QCAB) - tr(QCBA).
 The flavour matrix of a leg at a pole of mass m with nonzero ALP component is
    2 epsilon M_{m^2}^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1)
      / sqrt(1 + 2 epsilon^2 tr[(M_{m^2}^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1))^2]),
 with m = m_a for the ALP and m = m_1, m_2 for neutral daughters. Charged
 daughters (and neutral poles with zero ALP component, e.g. the exact K^0
 pole at m_{K^0,0}^2) carry their outgoing generator lambda^dagger.

The common WZW orientation is -1: it is the output "Common WZW orientation
required by the signed neutral-pion comparison" of
scripts/derive_wzw_external_currents.py (memory/rcht/odd-parity/
levi-civita-convention.md, wzw-gamma3pi-contribution.md); the test file
re-runs that derivation.

Derivation (each step checked in this script):
 1. Canonical frame P_i = P_i' + epsilon tr(lambda_i c_q^A) a: K' = 1 exactly,
    the mixed mass entry is -epsilon tr(lambda_i [tadpole argument]), the
    WZW combination Phi - 2 epsilon c_q^A a becomes Phi', and the photon
    coefficient of a/f_a becomes hatTheta_{gamma,1}. Only first-order dressed
    invariants and m_q enter; the ALP mass entry contains the second-order
    invariants, which reach one-ALP amplitudes only through m_a.
 2. Pole vector at s = m^2 in the canonical frame: meson matrix
    sum_i u_i' lambda_i = 2 epsilon y M_s^{-1}(tadpole argument), ALP component y;
    K' = 1 normalization gives y = 1/sqrt(1 + 2 epsilon^2 tr[(M_s^{-1}(...))^2]).
    In the original fields the mixing vector is
    w_i = epsilon tr(lambda_i [c_q^A + M_{m_a^2}^{-1}(...)]); P_i = P_i'' + w_i a
    removes the ALP-meson mixing at the ALP pole.
 3. The photon vertex contracted with the pole vector, and the WZW legs
    sum_i u_i lambda_i - 2 epsilon c_q^A u_a, give the formulas above.
Identities (not assumptions): one-ALP amplitudes receive no contribution from
WZW terms with two ALP legs, because epsilon^{mu nu rho sigma} p_mu p_nu = 0
for the two ALP derivatives (checked); the Ward identities follow from the
same antisymmetry (checked).

Assumptions and scope: LO, tree level, stable simple poles; the ALP is
neutral and all sources commute with Q; chiral redefinitions are restricted to
[Q, g_{R,L}(x)] = 0 (U(1)_u x U(2)_ds), since redefinitions that do not
commute with Q turn the photon coupling into an x-dependent chiral charge,
not treated in the records (photon dressed record, Assumption 2). The
gamma-gamma Lagrangian has no direct c_q^A term, as in the photon-invariant
record. The vector current C_{V,0} does not enter the leading WZW
gamma P P a vertex used here (scope of the WZW record). No vector mesons.
"""

from __future__ import annotations

import contextlib
import io
import itertools
import random
from dataclasses import dataclass, replace

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS, NEUTRAL_GENERATOR_INDICES, closed_form_blocks
from derive_wzw_alp_photon_two_pseudoscalar import outgoing_generators, photon_three_meson_coefficient
from derive_wzw_alp_photon_two_pseudoscalar import orientation as record_orientation
from solve_alp_meson_spectrum import solve_spectrum_numerical
from symbolic_utils import levi_civita_contract, u3_generators

NUMBER_OF_COLOURS = 3
NUMBER_OF_FLAVOURS = 3
CHARGE_MATRIX = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
NEUTRAL_GENERATORS = tuple(sp.Matrix(u3_generators()[index]) for index in NEUTRAL_GENERATOR_INDICES)
WZW_ORIENTATION = -1  # scripts/derive_wzw_external_currents.py, signed neutral-pion matching
HIGH_PRECISION_DIGITS = 50
electric_charge = sp.Symbol("e", positive=True)


def anticommutator(first, second):
    return first * second + second * first


def six_term_flavour_trace(charge, first, second, third):
    """F(Q;A,B,C) = tr(QABC)-tr(QACB)-tr(QBAC)+tr(QBCA)+tr(QCAB)-tr(QCBA)."""
    def trace(*matrices):
        product = charge
        for matrix in matrices:
            product = product * matrix
        return product.trace()
    return sp.expand(trace(first, second, third) - trace(first, third, second) - trace(second, first, third)
                     + trace(second, third, first) + trace(third, first, second) - trace(third, second, first))


# ------------------------------------------------------------ theory inputs


@dataclass(frozen=True)
class PartonicTheory:
    """First Taylor coefficients of the partonic ALP functions and LO meson constants.

    quark_mass m_q = M(0) (Hermitian, here diagonal), mass_first_derivative M_1,
    right_current C_{R,0}, left_current C_{L,0}, gluon_coefficient Theta_1
    (= c_G), photon_coefficient Theta_{gamma,1} (= c_gamma). The ALP entry of
    the canonical mass matrix, which carries the second-order invariants, is
    set by bare_alp_mass_squared (the m_{a,0}^2 of closed_form_blocks).
    """

    quark_mass: sp.Matrix
    mass_first_derivative: sp.Matrix
    right_current: sp.Matrix
    left_current: sp.Matrix
    gluon_coefficient: sp.Expr
    photon_coefficient: sp.Expr
    condensate_parameter: sp.Expr
    anomaly_mass_squared: sp.Expr
    pion_decay_constant: sp.Expr
    decay_constant_ratio: sp.Expr
    bare_alp_mass_squared: sp.Expr

    @property
    def alp_decay_constant(self):
        return self.pion_decay_constant / self.decay_constant_ratio

    @property
    def axial_source(self):
        """c_q^A = C_{A,0} = (C_{R,0} - C_{L,0})/2."""
        return (self.right_current - self.left_current) / 2

    @property
    def pseudoscalar_source(self):
        """c_q^P from P_1 = {c_q^P, m_q} (Alda:2026nrd parametrization), diagonal m_q."""
        pseudoscalar_part = (self.mass_first_derivative - self.mass_first_derivative.H) / (2 * sp.I)
        return sp.Matrix(3, 3, lambda row, column: pseudoscalar_part[row, column]
                         / (self.quark_mass[row, row] + self.quark_mass[column, column]))


def theory_from_alda_couplings(vector, axial, scalar, pseudoscalar, gluon, photon, **constants):
    """Taylor coefficients from c_q^{V,A,S,P}, c_G, c_gamma: M_1 = i[c^S,m_q] + i{c^P,m_q}, C_{R,L} = c^V +- c^A."""
    quark_mass = constants["quark_mass"]
    first_derivative = sp.I * (scalar * quark_mass - quark_mass * scalar) + sp.I * anticommutator(pseudoscalar, quark_mass)
    return PartonicTheory(mass_first_derivative=sp.expand(first_derivative), right_current=vector + axial,
                          left_current=vector - axial, gluon_coefficient=gluon, photon_coefficient=photon,
                          **constants)


def illustrative_theory(decay_constant_ratio=sp.Rational(1, 2), bare_alp_mass_squared=sp.Rational(3, 5)):
    """Charge-conserving test theory with complex d-s entries (generic neutral sources).

    chi_0 = 2 B_0 m_q = diag(3/200, 1/40, 47/100) GeV^2 as in
    illustrative_lo_inputs("neutral") of derive_alp_lo_neutral_blocks.py.
    """
    rational = sp.Rational
    condensate_parameter = rational(27, 10)
    chiral_mass = sp.diag(rational(3, 200), rational(1, 40), rational(47, 100))
    axial = sp.Matrix([[1, 0, 0], [0, rational(1, 2), rational(1, 3) + sp.I / 4],
                       [0, rational(1, 3) - sp.I / 4, rational(-1, 3)]])
    pseudoscalar = sp.Matrix([[rational(1, 4), 0, 0], [0, rational(-1, 2), rational(1, 2) - sp.I / 3],
                              [0, rational(1, 2) + sp.I / 3, rational(2, 3)]])
    vector = sp.Matrix([[rational(2, 7), 0, 0], [0, rational(-1, 5), rational(1, 6) - sp.I / 5],
                        [0, rational(1, 6) + sp.I / 5, rational(3, 8)]])
    scalar = sp.Matrix([[rational(-1, 9), 0, 0], [0, rational(1, 3), rational(-1, 4) + sp.I / 7],
                        [0, rational(-1, 4) - sp.I / 7, rational(1, 5)]])
    return theory_from_alda_couplings(
        vector, axial, scalar, pseudoscalar, rational(1, 2), rational(1, 7),
        quark_mass=chiral_mass / (2 * condensate_parameter), condensate_parameter=condensate_parameter,
        anomaly_mass_squared=rational(4, 5), pion_decay_constant=rational(23, 250),
        decay_constant_ratio=sp.nsimplify(decay_constant_ratio), bare_alp_mass_squared=bare_alp_mass_squared)


# -------------------------------------------------------- dressed invariants


def dressed_first_order_invariants(theory):
    """hatP_1, hatS_1, hatTheta_1, hatTheta_{gamma,1} from the first Taylor coefficients.

    hatM_1 = M_1 - i(C_{R,0} m_q - m_q C_{L,0}) (first-order dressed-mass record);
    hatTheta_1 = Theta_1 + tr C_{A,0}; hatTheta_{gamma,1} = Theta_{gamma,1} + 2 N_c tr(Q^2 C_{A,0}).
    """
    dressed_mass = theory.mass_first_derivative - sp.I * (theory.right_current * theory.quark_mass
                                                          - theory.quark_mass * theory.left_current)
    axial = theory.axial_source
    return {
        "dressed_pseudoscalar_source": ((dressed_mass - dressed_mass.H) / (2 * sp.I)).applyfunc(sp.expand),
        "dressed_scalar_source": ((dressed_mass + dressed_mass.H) / 2).applyfunc(sp.expand),
        "dressed_gluon_coefficient": sp.expand(theory.gluon_coefficient + axial.trace()),
        "dressed_photon_coefficient": sp.expand(theory.photon_coefficient
                                                + 2 * NUMBER_OF_COLOURS * (CHARGE_MATRIX**2 * axial).trace()),
    }


def tadpole_argument(theory):
    """B_0 hatP_1 - (m_0^2/N_f) hatTheta_1 1 (the argument of M_s^{-1} in all results)."""
    invariants = dressed_first_order_invariants(theory)
    return (theory.condensate_parameter * invariants["dressed_pseudoscalar_source"]
            - theory.anomaly_mass_squared / NUMBER_OF_FLAVOURS * invariants["dressed_gluon_coefficient"]
            * sp.eye(3)).applyfunc(sp.expand)


# ------------------------------------------------------------ mass operator


def mass_operator(matrix, squared_momentum, theory):
    """M_s(X) = {B_0 m_q, X} + (m_0^2/N_f) tr(X) 1 - s X."""
    return (anticommutator(theory.condensate_parameter * theory.quark_mass, matrix)
            + theory.anomaly_mass_squared / NUMBER_OF_FLAVOURS * matrix.trace() * sp.eye(3)
            - squared_momentum * matrix)


def neutral_coordinates(matrix):
    """x_i = tr(lambda_i X)/2 for the neutral generators (lambda_3, lambda_8, lambda_0, lambda_6, lambda_7)."""
    return sp.Matrix([sp.expand((generator * matrix).trace() / 2) for generator in NEUTRAL_GENERATORS])


def from_neutral_coordinates(coordinates):
    return sum((coordinates[index] * generator for index, generator in enumerate(NEUTRAL_GENERATORS)),
               sp.zeros(3))


def inverse_mass_operator(matrix, squared_momentum, theory):
    """M_s^{-1}(X): the neutral Y with M_s(Y) = X (exact for rational or symbolic s, Floats otherwise)."""
    operator_matrix = sp.Matrix.hstack(*[neutral_coordinates(mass_operator(generator, squared_momentum, theory))
                                         for generator in NEUTRAL_GENERATORS])
    solution = operator_matrix.LUsolve(neutral_coordinates(matrix))
    return from_neutral_coordinates(solution.applyfunc(sp.cancel if squared_momentum.free_symbols else sp.expand))


# ------------------------------------------------------------ final formulas


def inverse_on_tadpole_argument(squared_momentum, theory):
    """M_s^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1)."""
    return inverse_mass_operator(tadpole_argument(theory), squared_momentum, theory)


def alp_component_squared_inverse(squared_momentum, theory):
    """1 + 2 epsilon^2 tr[(M_s^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1))^2] (inverse squared ALP component)."""
    inverse = inverse_on_tadpole_argument(squared_momentum, theory)
    return sp.expand(1 + 2 * theory.decay_constant_ratio**2 * (inverse * inverse).trace())


def two_photon_bracket(squared_momentum, theory):
    """hatTheta_{gamma,1} + 2 N_c tr(Q^2 M_s^{-1}(B_0 hatP_1 - m_0^2/N_f hatTheta_1 1))."""
    invariants = dressed_first_order_invariants(theory)
    inverse = inverse_on_tadpole_argument(squared_momentum, theory)
    return sp.expand(invariants["dressed_photon_coefficient"]
                     + 2 * NUMBER_OF_COLOURS * (CHARGE_MATRIX**2 * inverse).trace())


def two_photon_amplitude_factor(alp_squared_mass, theory):
    """Coefficient of E in M(a -> gamma gamma): -e^2/(4 pi^2 f_a) bracket / sqrt(normalization)."""
    return (-electric_charge**2 / (4 * sp.pi**2 * theory.alp_decay_constant)
            * two_photon_bracket(alp_squared_mass, theory)
            / sp.sqrt(alp_component_squared_inverse(alp_squared_mass, theory)))


def leg_flavour_matrix(pole_squared_mass, theory):
    """2 epsilon M_{m^2}^{-1}(...)/sqrt(1 + 2 epsilon^2 tr[(M_{m^2}^{-1}(...))^2]) for a pole with ALP component."""
    inverse = inverse_on_tadpole_argument(pole_squared_mass, theory)
    normalization = 1 + 2 * theory.decay_constant_ratio**2 * (inverse * inverse).trace()
    return 2 * theory.decay_constant_ratio * inverse / sp.sqrt(normalization)


def photon_two_pseudoscalar_factor(parent_leg, first_leg, second_leg, theory, orientation=WZW_ORIENTATION):
    """A = -orientation N_c e/(24 pi^2 f_pi^3) F(Q; X_a, X_1, X_2)."""
    return (-orientation * NUMBER_OF_COLOURS * electric_charge / (24 * sp.pi**2 * theory.pion_decay_constant**3)
            * six_term_flavour_trace(CHARGE_MATRIX, parent_leg, first_leg, second_leg))


# ----------------------------------------------------- Lorentz derivations


def two_photon_factor_from_lagrangian():
    """Amplitude of L = (alpha_em/4pi) G a F Ftilde, two outgoing photons, divided by G E.

    F_{mu nu} -> i(q_mu eps_nu - q_nu eps_mu) per outgoing photon; F Ftilde =
    (1/2) eps^{mu nu rho sigma} F F gives both photon assignments.
    Returns the number multiplying G E (expected -e^2/(4 pi^2)).
    """
    momenta = [sp.Matrix(sp.symbols(f"q{label}_0:4")) for label in (1, 2)]
    polarizations = [sp.Matrix(sp.symbols(f"eps{label}_0:4")) for label in (1, 2)]
    metric = sp.diag(1, -1, -1, -1)
    lowered = [(metric * vector) for vector in momenta + polarizations]
    strengths = [sp.Matrix(4, 4, lambda mu, nu: sp.I * (lowered[k][mu] * lowered[k + 2][nu]
                                                        - lowered[k][nu] * lowered[k + 2][mu])) for k in (0, 1)]
    cross_term = sum(sp.LeviCivita(mu, nu, rho, sigma) * strengths[0][mu, nu] * strengths[1][rho, sigma]
                     for mu, nu, rho, sigma in itertools.permutations(range(4)))
    lorentz_factor = levi_civita_contract(momenta[0], polarizations[0], momenta[1], polarizations[1])
    fine_structure = electric_charge**2 / (4 * sp.pi)
    return sp.simplify(fine_structure / (4 * sp.pi) * cross_term / lorentz_factor)


def wzw_vertex_sum(flavour_matrices, momenta, photon_polarization):
    """sum over assignments of legs to the three d Phi of tr(Q X X X) eps(eps_gamma, p, p, p)."""
    total = 0
    for order in itertools.permutations(range(3)):
        trace = (CHARGE_MATRIX * flavour_matrices[order[0]] * flavour_matrices[order[1]]
                 * flavour_matrices[order[2]]).trace()
        total += trace * levi_civita_contract(photon_polarization, momenta[order[0]], momenta[order[1]],
                                              momenta[order[2]])
    return sp.expand(total)


def wzw_amplitude_prefactor(seed=3):
    """Derive A/F from the gamma Phi^3 coefficient of the WZW record and the Fourier factors.

    L = c eps^{mu nu rho sigma} A_mu tr(Q d_nu Phi d_rho Phi d_sigma Phi), c from
    derive_wzw_alp_photon_two_pseudoscalar.py (keeps the symbol `orientation`).
    Incoming ALP: d -> -i p_a; outgoing mesons: d -> +i p. M = c * sum (iM = i L).
    Returns M / [F(Q; X_a, X_1, X_2) eps(eps_gamma, p_1, p_2, k)] with p_a = k + p_1 + p_2.
    """
    generator = random.Random(seed)
    legs = [sp.Matrix(3, 3, lambda r, c: sp.Rational(generator.randint(-5, 5), generator.randint(1, 4)))
            for _ in range(3)]
    photon_momentum, first, second, polarization = [sp.Matrix(sp.symbols(f"{name}_0:4"))
                                                    for name in ("k", "p1", "p2", "eps_gamma")]
    alp_momentum = photon_momentum + first + second
    amplitude = photon_three_meson_coefficient * wzw_vertex_sum(
        legs, [-sp.I * alp_momentum, sp.I * first, sp.I * second], polarization)
    reference = (six_term_flavour_trace(CHARGE_MATRIX, *legs)
                 * levi_civita_contract(polarization, first, second, photon_momentum))
    return sp.simplify(amplitude / reference)


def wzw_orientation_from_neutral_pion_matching():
    """Re-run scripts/derive_wzw_external_currents.py's signed neutral-pion comparison (returns -1)."""
    import derive_wzw_external_currents
    with contextlib.redirect_stdout(io.StringIO()):
        return derive_wzw_external_currents.derive_neutral_pion_two_photons()


# ----------------------------------------------------- original-field route


def lo_blocks(theory):
    """K, M^2 on (P_3, P_8, P_0, P_6, P_7, a) from closed_form_blocks with chi_0 = 2 B_0 m_q."""
    return closed_form_blocks(
        axial_source=theory.axial_source, pseudoscalar_source=theory.pseudoscalar_source,
        quark_spurion=2 * theory.condensate_parameter * theory.quark_mass,
        gluon_coupling=theory.gluon_coefficient, pion_decay_constant=theory.pion_decay_constant,
        alp_decay_constant=theory.alp_decay_constant, anomaly_mass_squared=theory.anomaly_mass_squared,
        bare_alp_mass_squared=theory.bare_alp_mass_squared)


def canonical_frame_matrix(theory):
    """T with Phi_original = T Phi': P_i = P_i' + epsilon tr(lambda_i c_q^A) a."""
    transform = sp.eye(6)
    for index, generator in enumerate(NEUTRAL_GENERATORS):
        transform[index, 5] = sp.expand(theory.decay_constant_ratio * (generator * theory.axial_source).trace())
    return transform


def canonical_alp_mass_entry(theory):
    """(T^T M^2 T)_{aa}: the canonical ALP mass entry (contains the second-order invariants)."""
    _, mass_squared = lo_blocks(theory)
    transform = canonical_frame_matrix(theory)
    return sp.expand((transform.T * mass_squared * transform)[5, 5])


def schur_mixing_vector(squared_momentum, theory):
    """w(s) = -D_PP(s)^{-1} D_Pa(s) from the original blocks (the pole direction (w, 1) at a pole)."""
    kinetic, mass_squared = lo_blocks(theory)
    inverse_propagator = squared_momentum * kinetic - mass_squared
    return -inverse_propagator[:5, :5].LUsolve(inverse_propagator[:5, 5])


def leg_matrix_from_pole_vector(pole_vector, theory):
    """WZW leg in the original fields: sum_i u_i lambda_i - 2 epsilon c_q^A u_a."""
    meson = sum((pole_vector[index] * generator for index, generator in enumerate(NEUTRAL_GENERATORS)), sp.zeros(3))
    return meson - 2 * theory.decay_constant_ratio * theory.axial_source * pole_vector[5]


def two_photon_bracket_from_pole_vector(pole_vector, theory):
    """c_gamma u_a + (N_c/epsilon) sum_i u_i tr(Q^2 lambda_i): coefficient of (alpha/4pi)(1/f_a) F Ftilde."""
    meson = sum((pole_vector[index] * (CHARGE_MATRIX**2 * generator).trace()
                 for index, generator in enumerate(NEUTRAL_GENERATORS)), 0)
    return theory.photon_coefficient * pole_vector[5] + NUMBER_OF_COLOURS / theory.decay_constant_ratio * meson


def solve_float_spectrum(theory):
    kinetic, mass_squared = lo_blocks(theory)
    return solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                    np.array(mass_squared.evalf(30), dtype=float), field_labels=NEUTRAL_FIELD_LABELS)


def unmixed_meson_squared_masses(theory):
    """Eigenvalues of M_0 on the neutral space (exact)."""
    operator_matrix = sp.Matrix.hstack(*[neutral_coordinates(mass_operator(generator, 0, theory))
                                         for generator in NEUTRAL_GENERATORS])
    return list(operator_matrix.eigenvals())


def high_precision_poles(theory, digits=HIGH_PRECISION_DIGITS):
    """Zeros of det(s K - M^2) to `digits` digits (Poly.nroots on the exact polynomial)."""
    kinetic, mass_squared = lo_blocks(theory)
    squared_momentum = sp.Symbol("s")
    polynomial = sp.Poly(sp.expand((squared_momentum * kinetic - mass_squared).det(method="berkowitz")),
                         squared_momentum)
    return sorted(sp.re(root) for root in polynomial.nroots(n=digits, maxsteps=500))


def poles_with_alp_component(theory, digits=HIGH_PRECISION_DIGITS):
    """High-precision poles away from the unmixed meson masses, with their inverse squared ALP component."""
    unmixed = [sp.N(value, digits) for value in unmixed_meson_squared_masses(theory)]
    result = []
    for pole in high_precision_poles(theory, digits):
        if min(abs(pole - value) for value in unmixed) < sp.Float(10)**(-digits // 2):
            continue
        result.append((pole, sp.re(alp_component_squared_inverse(sp.Float(pole, digits), theory))))
    return result


def alp_like_pole(theory, digits=HIGH_PRECISION_DIGITS):
    """Pole with the largest ALP component (smallest 1 + 2 epsilon^2 tr[(M^{-1})^2])."""
    return min(poles_with_alp_component(theory, digits), key=lambda item: item[1])[0]


def alp_secular_function(squared_momentum, theory):
    """D^eff_aa(s) = s - (T^T M^2 T)_aa + 2 epsilon^2 tr(arg M_s^{-1}(arg)), arg = B_0 hatP_1 - m_0^2/N_f hatTheta_1 1."""
    argument = tadpole_argument(theory)
    return sp.expand(squared_momentum - canonical_alp_mass_entry(theory) + 2 * theory.decay_constant_ratio**2
                     * (argument * inverse_mass_operator(argument, squared_momentum, theory)).trace())


def two_photon_weak_mixing_coefficients(theory):
    """Exact epsilon^0 and epsilon^2 coefficients of bracket(m_a^2)/sqrt(normalization(m_a^2)).

    Pole: m_a^2 = m_{a,0}^2 + epsilon^2 sigma + O(epsilon^4), sigma = d(T^T M^2 T)_aa/d(epsilon^2)
    - 2 tr(arg M_{m_{a,0}^2}^{-1}(arg)), since dD^eff_aa/ds = 1 at epsilon = 0. With
    d M_s^{-1}(X)/ds = M_s^{-1}(M_s^{-1}(X)), the epsilon^2 coefficient is
    2 N_c tr(Q^2 M^{-1}(M^{-1}(arg))) sigma - tr[(M^{-1}(arg))^2] bracket(m_{a,0}^2).
    """
    ratio_symbol = sp.Symbol("epsilon", positive=True)
    symbolic = replace(theory, decay_constant_ratio=ratio_symbol)
    bare = theory.bare_alp_mass_squared
    slope = sp.expand(canonical_alp_mass_entry(symbolic)).coeff(ratio_symbol, 2)
    argument = tadpole_argument(theory)
    inverse = inverse_mass_operator(argument, bare, theory)
    pole_slope = sp.expand(slope - 2 * (argument * inverse).trace())
    zeroth = two_photon_bracket(bare, theory)
    first = sp.expand(2 * NUMBER_OF_COLOURS * (CHARGE_MATRIX**2 * inverse_mass_operator(inverse, bare, theory)).trace()
                      * pole_slope - (inverse * inverse).trace() * zeroth)
    return zeroth, first, pole_slope


# ------------------------------------------------- chiral field redefinitions


def random_charge_commuting_hermitian(generator):
    """Random rational Hermitian matrix commuting with Q (u entry, d-s block)."""
    def rational():
        return sp.Rational(generator.randint(-6, 6), generator.randint(1, 5))
    off_diagonal = rational() + sp.I * rational()
    return sp.Matrix([[rational(), 0, 0], [0, rational(), off_diagonal], [0, sp.conjugate(off_diagonal), rational()]])


def redefinition_taylor_coefficients(first_generator, second_generator):
    """g(x) = exp(i [x G_1 + x^2 G_2]): Taylor coefficients g_0, g_1, g_2 (derivatives at x = 0)."""
    return [sp.eye(3), sp.I * first_generator, 2 * sp.I * second_generator - first_generator * first_generator]


def random_redefinition(seed):
    """Nonlinear g_R(x), g_L(x) with g(0) = 1 and [Q, g] = 0; returns their Taylor coefficients."""
    generator = random.Random(seed)
    return tuple(redefinition_taylor_coefficients(random_charge_commuting_hermitian(generator),
                                                  random_charge_commuting_hermitian(generator))
                 for _ in ("right", "left"))


def transform_theory(theory, right_coefficients, left_coefficients, second_mass_derivative=None):
    """Apply the Taylor laws of memory/alp/shared/alp-general-chiral-redefinition-invariants.md.

    C'_{eta,0} = C_{eta,0} + i g_{eta,1}^dagger,  M'_1 = g_{R,1} m_q + M_1 + m_q g_{L,1}^dagger,
    Theta'_1 = Theta_1 - (1/2i) tr(g_{R,1} - g_{L,1}),
    Theta'_{gamma,1} = Theta_{gamma,1} - (2 N_c/2i) tr(Q^2 (g_{R,1} - g_{L,1}))  ([Q, g] = 0).
    Also returns M'_2 and C'_{eta,1} (second-order coefficients) for the
    statement that they do not enter. The bare ALP entry is shifted so that the
    canonical ALP mass entry, and hence m_a, is unchanged (it contains the
    second-order invariants, which closed_form_blocks does not carry).
    """
    mass = [theory.quark_mass, theory.mass_first_derivative,
            sp.zeros(3) if second_mass_derivative is None else second_mass_derivative]
    second_mass = sp.zeros(3)
    for right_order, left_order in itertools.product(range(3), repeat=2):
        mass_order = 2 - right_order - left_order
        if mass_order >= 0:
            second_mass += (sp.factorial(2) / (sp.factorial(right_order) * sp.factorial(left_order)
                                               * sp.factorial(mass_order))
                            * right_coefficients[right_order] * mass[mass_order] * left_coefficients[left_order].H)
    currents = {}
    for name, current, coefficients in (("right", theory.right_current, right_coefficients),
                                        ("left", theory.left_current, left_coefficients)):
        currents[name] = (current + sp.I * coefficients[1].H).applyfunc(sp.expand)
        currents[name + "_first"] = (coefficients[1] * current + current * coefficients[1].H
                                     + sp.I * coefficients[2].H + sp.I * coefficients[1] * coefficients[1].H
                                     ).applyfunc(sp.expand)
    difference = right_coefficients[1] - left_coefficients[1]
    transformed = replace(
        theory,
        mass_first_derivative=(right_coefficients[1] * theory.quark_mass + theory.mass_first_derivative
                               + theory.quark_mass * left_coefficients[1].H).applyfunc(sp.expand),
        right_current=currents["right"], left_current=currents["left"],
        gluon_coefficient=sp.expand(theory.gluon_coefficient - difference.trace() / (2 * sp.I)),
        photon_coefficient=sp.expand(theory.photon_coefficient
                                     - 2 * NUMBER_OF_COLOURS * (CHARGE_MATRIX**2 * difference).trace() / (2 * sp.I)))
    shift = canonical_alp_mass_entry(theory) - canonical_alp_mass_entry(transformed)
    transformed = replace(transformed, bare_alp_mass_squared=sp.expand(theory.bare_alp_mass_squared + shift))
    second_order = {"second_mass_derivative": second_mass.applyfunc(sp.expand),
                    "right_current_first_derivative": currents["right_first"],
                    "left_current_first_derivative": currents["left_first"]}
    return transformed, second_order


# ------------------------------------------------------- numerical channels


def float_channel_comparison(theory):
    """Compare pole-vector amplitudes (original fields, float solver) with the invariant formulas."""
    spectrum = solve_float_spectrum(theory)
    vectors = spectrum.pole_vectors
    alp_index = int(np.argmax(np.abs(vectors[5, :])))
    result = {"poles": list(spectrum.pole_squared_masses), "alp_index": alp_index, "legs": {}}
    alp_pole = spectrum.pole_squared_masses[alp_index]
    alp_vector = vectors[:, alp_index]
    inverse = np.array(inverse_on_tadpole_argument(sp.Float(alp_pole, 30), theory).evalf(), dtype=complex)
    ratio = float(theory.decay_constant_ratio)
    mixing_vector = np.array([ratio * np.trace(np.array(generator, dtype=complex)
                                               @ (np.array(theory.axial_source, dtype=complex) + inverse)).real
                              for generator in NEUTRAL_GENERATORS])
    result["mixing_vector_difference"] = float(np.max(np.abs(alp_vector[:5] / alp_vector[5] - mixing_vector)))
    result["normalization_difference"] = abs(
        float(alp_component_squared_inverse(sp.Float(alp_pole, 30), theory)) - 1 / alp_vector[5]**2)
    bracket_vector = complex(two_photon_bracket_from_pole_vector(list(alp_vector), theory))
    bracket_invariant = complex(two_photon_bracket(sp.Float(alp_pole, 30), theory)
                                / sp.sqrt(alp_component_squared_inverse(sp.Float(alp_pole, 30), theory)))
    result["two_photon"] = (bracket_vector, bracket_invariant)
    for index, pole in enumerate(spectrum.pole_squared_masses):
        if abs(vectors[5, index]) < 1e-12:
            continue
        from_vector = np.array(leg_matrix_from_pole_vector(list(vectors[:, index]), theory).evalf(), dtype=complex)
        invariant = np.array(leg_flavour_matrix(sp.Float(pole, 30), theory).evalf(), dtype=complex)
        result["legs"][index] = (from_vector, invariant, float(np.max(np.abs(from_vector - invariant))))
    return result


def main():
    print("LO finite-mixing photon channels; natural units, epsilon^{0123}=+1, N_c=3, epsilon=f_pi/f_a.")
    print("Illustrative rational inputs (GeV units), not phenomenology.\n")
    print("WZW orientation from the neutral-pion matching:", wzw_orientation_from_neutral_pion_matching())
    print("gamma gamma: amplitude / (G E) from (alpha/4pi) G a F Ftilde:", two_photon_factor_from_lagrangian())
    print("gamma P P: A / F(Q;X_a,X_1,X_2) =", wzw_amplitude_prefactor(),
          "-> with orientation = -1:", wzw_amplitude_prefactor().subs(record_orientation, WZW_ORIENTATION), "\n")

    theory = illustrative_theory(sp.Rational(1, 2))
    invariants = dressed_first_order_invariants(theory)
    print("First-order dressed invariants: hatTheta_1 =", invariants["dressed_gluon_coefficient"],
          ", hatTheta_gamma1 =", invariants["dressed_photon_coefficient"])
    print("hatP_1 =", invariants["dressed_pseudoscalar_source"].tolist(), "\n")

    for ratio in (sp.Rational(1, 2), sp.Rational(1, 10)):
        theory = illustrative_theory(ratio)
        comparison = float_channel_comparison(theory)
        print(f"epsilon = {ratio}: poles (GeV^2) {np.round(comparison['poles'], 8).tolist()}")
        print(f"  ALP pole {comparison['poles'][comparison['alp_index']]:.12f}; max|u_i/u_a - w_i| ="
              f" {comparison['mixing_vector_difference']:.2e}; |normalization - 1/u_a^2| ="
              f" {comparison['normalization_difference']:.2e}")
        vector_value, invariant_value = comparison["two_photon"]
        print(f"  gamma gamma bracket: pole vector {vector_value.real:.14f}, invariant {invariant_value.real:.14f}")
        print("  leg matrices max difference per pole:",
              {index: f"{value[2]:.1e}" for index, value in comparison["legs"].items()})

    print("\nRandom nonlinear Q-commuting redefinitions (epsilon = 1/2):")
    theory = illustrative_theory(sp.Rational(1, 2))
    for seed in (11, 12):
        transformed, _ = transform_theory(theory, *random_redefinition(seed))
        before, after = float_channel_comparison(theory), float_channel_comparison(transformed)
        print(f"  seed {seed}: max pole shift {max(abs(a - b) for a, b in zip(before['poles'], after['poles'])):.1e},"
              f" gamma gamma bracket {before['two_photon'][0].real:.12f} -> {after['two_photon'][0].real:.12f}")


if __name__ == "__main__":
    main()
