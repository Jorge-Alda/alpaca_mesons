#!/usr/bin/env python3
"""Dimensionally consistent loop-size estimate for a -> pi+ pi- pi0 (up-only model), physical units.

NO LOOP IS COMPUTED.  No loop integral, Feynman parameter, angular average, tensor reduction
or scalar loop function appears.  Each diagram is replaced by a product of tree-level vertices,
pole data and a power of its loop scale fixed by dimensional analysis:

    M1_est = sum_d S_d (m_d^2)^(2 - r_d) / (16 pi^2) * A_d,      R = |M1_est| / |M0|,

    d     one-loop diagram of the classes 1-8 of evaluate_explicit_estimate.py (explicit list DIAGRAMS,
          labels "class.index", e.g. 6.7);
    r_d   number of loop propagators (= number of vertices on the loop);
    m_d   loop scale: m_a for loops with an ALP-like (A) line, m_sigma for loops with no A line
          (radial lines, possibly together with meson lines);
    A_d   flavour trace tr[V(B_1) W_1 V(B_2) W_2 ... V(B_r) W_r] of the TREE vertices in physical
          units, each internal line contributing its pole-projected weight W (pole residue, the
          line weights of evaluate_explicit_estimate.py).  The explicit eps = f_pi/f_a and N(m_a^2)
          factors are inside the vertices and pole vectors and are not multiplied again;
    S_d   symmetry factor 1/|Aut_d| (checked by enumerate_murayama_up_only_loop_diagrams.py).
Vertex set: the tree vertices of evaluate_explicit_estimate.py, including the sigma^3 vertex
V_sigma3 = -6 lambda f_a of the radial potential (enters class 6: diagrams 6.1, 6.3, 6.4 and the
sigma self-energy 6.7).
Every momentum product that contains the loop momentum is replaced by m_d^2; products of
external momenta (pole masses and channel invariants s_uv) are kept.

Justification of the power (assumption: a dimensional size rule, not a loop evaluation): the loop
measure d^4k/(2 pi)^4 gives m^4/(16 pi^2) and every propagator gives 1/m^2, so a loop with r
propagators is (m^2)^(2 - r)/(16 pi^2) times the product of its vertices in physical units.
A vertex with n legs has mass dimension 4 - n, the r vertices of the loop have 2r - 4 external
legs in total (a radial branch counts as one radial leg times its exchange V3/(m_sigma^2 - s),
of dimension -1), so every term of M1_est is dimensionless like M0.
External-leg residues are one-propagator loops (r = 1): (1/2) dSigma_ee/dp^2 M0 with
Sigma ~ m^2/(16 pi^2) V4: a closed A line (class 7, diagrams 7.1-7.4, loop scale m_a) on all four
legs a, pi+, pi-, pi0 (on the charged legs through the pi+ pi- A A vertex, which charge
conservation allows), and a closed sigma line (class 8, diagrams 8.1-8.2, loop scale m_sigma) on
the a and pi0 legs (the sigma sigma pi+ pi- vertex vanishes in the up-only model).  The pi+ slope has the closed form
d/dp^2 [(1/2) V(l_+, l_-, l_A, l_A)] = -eps^2 [(M^-1)_11 - (M^-1)_22]^2 / (3 N(m_a^2) f_pi^2),
M^-1 = M^{-1}_{m_a^2}(C).  M0 is the tree of evaluate_explicit_estimate.py in physical units.
Result record: memory/alp/shared/murayama-radial-up-only-estimate-result.md.

Physical units: every dimensionful quantity in GeV (f_pi, m_q, B0, m_1^2 in GeV^2), with the
inputs of memory/alp/shared/only-up-finite-mixing-inputs.md (scripts/implementation/only_up_inputs.py).
The quartic coupling lambda is dimensionless; f_a = sqrt(2) m_u / kappa_u, eps = f_pi/f_a.
Vertices: the vertex formulas of evaluate_explicit_estimate.py are written for f_pi = 1; in
physical units every term of an n-leg vertex carries f_pi^(2 - n) (the chiral terms
(f^2/4)<...> with U = exp(i Phi/f) give f^(2 - m) for m meson legs and every a/f_a, sigma/f_a gives
eps/f_pi; the singlet and radial kinetic terms have n = 2).  The pole data (M_s, F(s), N(s),
m_sigma^2 = B0 eps^2 tr(S_hat_10 - 2 S_hat_20) + 2 lambda f_a^2) contain no extra f_pi.
Sources enter only through the dressed invariants S_hat_rs, P_hat_rs of a frame (constant axial
redefinitions leave them unchanged).  Kinematic point: symmetric Dalitz point
s_{+-} = s_{+0} = s_{-0} = (m_a^2 + 2 m_pi+^2 + m_pi0^2)/3 (labelled choice).

Run:
    .venv/bin/python scripts/estimate_loop_size_dimensionless.py
Tests:
    PYTHONPATH=scripts .venv/bin/python -m unittest test_estimate_loop_size_dimensionless
"""
import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq

import evaluate_explicit_estimate as explicit
from evaluate_explicit_estimate import E_UP, RESIDUE_DIAGRAM_LABELS, Leg, diagram_groups, label_range
from only_up_inputs import B_ZERO, DOWN_MASS, M1_SQUARED, PION_DECAY_CONSTANT, STRANGE_MASS, UP_MASS

CLASS_ORDER = ('1', '2', '3', '4', '5', '6', '7', '8')
MAXIMUM_ALP_ORDER = 6
LOOP_FACTOR = 1 / (16 * math.pi ** 2)        # loop measure d^4k/(2 pi)^4 -> m^4/(16 pi^2)


# ---------------------------------------------------------------------------
# Frames and dressed invariants
#
# A frame is the data (M(x, y), C_R(x), C_L(x), Theta(x, y)) of the polar Lagrangian
# written in the field U of that frame, with x = a/f_a = eps a, y = sigma/f_a and the
# covariant derivative D U = dU - i dx (C_R U - U C_L).  The transports
# dW/dx = i C W, W(0) = 1, give the dressed source M_hat = W_R^dag M W_L, whose
# Hermitian parts are S_hat_rs and P_hat_rs, and Theta_hat = Theta + int_0^x tr C_A,
# C_A = (C_R - C_L)/2 (memory/alp/shared/alp-loop-estimate-frame-invariance.md, in
# the sign convention of this script: L contains <chi U^dag + U chi^dag>).
# Under U = xi U' xi, xi = exp(i kappa x):  M' = xi^dag M xi^dag,
# C_R' = xi^dag C_R xi - kappa,  C_L' = xi C_L xi^dag + kappa.  Only traceless kappa is
# used: then Theta (absent in this model) and tr C_A are unchanged.  The dressed sources are
# the same in every frame, and the Lagrangian of the field U_hat = W_R^dag U W_L has no
# connection and the sources M_hat.
# ---------------------------------------------------------------------------


class FrameData:
    """Sources {(r, s): 3x3} (coefficient of y^r x^s), connections {s: 3x3} and Theta {(r, s): float}."""

    def __init__(self, sources, right_connection, left_connection, theta, alp_order):
        self.sources, self.alp_order = sources, alp_order
        self.right_connection, self.left_connection, self.theta = right_connection, left_connection, theta


def series_multiply(first, second, alp_order):
    """Product of x-series {s: matrix} (or {(r, s): matrix} for the first factor), truncated in x."""
    result = {}
    for key, left in first.items():
        r, s = key if isinstance(key, tuple) else (None, key)
        for t, right in second.items():
            if s + t <= alp_order:
                new_key = (r, s + t) if r is not None else s + t
                result[new_key] = result.get(new_key, 0) + left @ right
    return result


def left_multiply(first, second, alp_order):
    """first(x) second(x, y) with first an x-series and second a double series."""
    result = {}
    for s, left in first.items():
        for (r, t), right in second.items():
            if s + t <= alp_order:
                result[(r, s + t)] = result.get((r, s + t), 0) + left @ right
    return result


def exponential_series(generator, sign, alp_order):
    """exp(sign i generator x) as an x-series."""
    return {s: np.linalg.matrix_power(sign * 1j * np.asarray(generator, dtype=complex), s) / math.factorial(s)
            for s in range(alp_order + 1)}


def transport_series(connection, alp_order):
    """W(x) with dW/dx = i C(x) W, W(0) = 1, order by order in x."""
    transport = {0: np.eye(3, dtype=complex)}
    for s in range(alp_order):
        value = sum((connection[t] @ transport[s - t] for t in connection if t <= s), np.zeros((3, 3), complex))
        transport[s + 1] = 1j * value / (s + 1)
    return transport


def axially_redefined_frame(frame, generator):
    """Frame data after U = xi U' xi, xi = exp(i kappa x), for a constant Hermitian kappa.

    Used only to produce transformed invariant-building data (sources and connections that
    enter dressed_invariants); no vertex or pole vector is ever built from these raw data.
    """
    order = frame.alp_order
    rotation, rotation_conjugate = exponential_series(generator, 1, order), exponential_series(generator, -1, order)
    sources = left_multiply(rotation_conjugate, frame.sources, order)
    multiplied = {}                                      # right multiplication by xi^dag
    for (r, s), value in sources.items():
        for t, right in rotation_conjugate.items():
            if s + t <= order:
                multiplied[(r, s + t)] = multiplied.get((r, s + t), 0) + value @ right
    generator = np.asarray(generator, dtype=complex)
    right = series_multiply(series_multiply(rotation_conjugate, frame.right_connection, order - 1), rotation, order - 1)
    left = series_multiply(series_multiply(rotation, frame.left_connection, order - 1), rotation_conjugate, order - 1)
    right[0] = right.get(0, 0) - generator
    left[0] = left.get(0, 0) + generator
    if abs(np.trace(generator)) > 1e-12:
        raise NotImplementedError('only traceless generators: then Theta and tr C_A are unchanged')
    return FrameData(multiplied, right, left, dict(frame.theta), order)


def dressed_invariants(frame):
    """Dressed sources {(r, s): S_hat_rs + i P_hat_rs} and {(r, s): Theta_hat_rs} of a frame."""
    order = frame.alp_order
    right_transport = transport_series(frame.right_connection, order)
    left_transport = transport_series(frame.left_connection, order)
    right_adjoint = {s: value.conj().T for s, value in right_transport.items()}
    dressed = left_multiply(right_adjoint, frame.sources, order)
    result = {}                                          # right multiplication by W_L
    for (r, s), value in dressed.items():
        for t, right in left_transport.items():
            if s + t <= order:
                result[(r, s + t)] = result.get((r, s + t), 0) + value @ right
    theta = dict(frame.theta)
    for s in set(frame.right_connection) | set(frame.left_connection):
        axial = 0.5 * (frame.right_connection.get(s, 0) - frame.left_connection.get(s, 0))
        theta[(0, s + 1)] = theta.get((0, s + 1), 0.) + np.trace(np.asarray(axial)).real / (s + 1)
    return result, theta


@dataclass(frozen=True)
class PhysicalInputs:
    """Dimensionful inputs in GeV (singlet_mass_squared in GeV^2); quartic_coupling dimensionless."""
    pion_decay_constant: float = PION_DECAY_CONSTANT
    up_mass: float = UP_MASS
    down_mass: float = DOWN_MASS
    strange_mass: float = STRANGE_MASS
    b_zero: float = B_ZERO
    singlet_mass_squared: float = M1_SQUARED
    quartic_coupling: float = explicit.DEFAULT_QUARTIC_COUPLING

    def rescaled(self, factor):
        """All dimensionful inputs multiplied by factor (mass dimension d -> factor^d)."""
        return replace(self, pion_decay_constant=factor * self.pion_decay_constant,
                       up_mass=factor * self.up_mass, down_mass=factor * self.down_mass,
                       strange_mass=factor * self.strange_mass, b_zero=factor * self.b_zero,
                       singlet_mass_squared=factor ** 2 * self.singlet_mass_squared)


DEFAULT_INPUTS = PhysicalInputs()


def polar_frame(inputs, maximum_alp_order=MAXIMUM_ALP_ORDER):
    """Polar frame of the up-only model in GeV: S + i P = m_u (1 + y) exp(-i x) e_1 + diag(0, m_d, m_s)."""
    sources = {}
    for radial_order in (0, 1):
        for alp_order in range(maximum_alp_order + 1):
            sources[(radial_order, alp_order)] = (inputs.up_mass * (-1j) ** alp_order
                                                  / math.factorial(alp_order) * E_UP)
    sources[(0, 0)] = sources[(0, 0)] + np.diag([0., inputs.down_mass, inputs.strange_mass])
    return FrameData(sources, {}, {}, {}, maximum_alp_order)


class DimensionfulEstimate(explicit.ExplicitEstimate):
    """Tree M0 and the estimate M1_est in physical units at one kappa_u.

    The vertex, line-weight, slot and cycle machinery is inherited from
    evaluate_explicit_estimate.ExplicitEstimate; only the units and the loop power change.
    """

    def __init__(self, kappa_u, inputs=DEFAULT_INPUTS, frame=None):
        self.kappa_u = kappa_u
        self.inputs = inputs
        self.pion_decay_constant = inputs.pion_decay_constant
        self.quartic_coupling = inputs.quartic_coupling
        self.epsilon = inputs.pion_decay_constant * kappa_u / (math.sqrt(2) * inputs.up_mass)
        self.alp_decay_constant = inputs.pion_decay_constant / self.epsilon
        self.b_zero = inputs.b_zero
        self.up_mass, self.down_mass, self.strange_mass = inputs.up_mass, inputs.down_mass, inputs.strange_mass
        self.singlet_mass_squared = inputs.singlet_mass_squared
        dressed_sources, dressed_theta = dressed_invariants(frame or polar_frame(inputs))
        assert all(abs(value) < 1e-12 for value in dressed_theta.values()), 'model assumes Theta_hat = 0'
        self.scalar_sources = {key: (value + value.conj().T) / 2 for key, value in dressed_sources.items()}
        self.pseudoscalar_sources = {key: (value - value.conj().T) / 2j for key, value in dressed_sources.items()}
        self.source_c = self.b_zero * self.pseudoscalar_sources[(0, 1)]
        self.radial_mass_squared = (self.b_zero * self.epsilon ** 2
                                    * np.trace(self.scalar_sources[(1, 0)]
                                               - 2 * self.scalar_sources.get((2, 0), 0 * E_UP)).real
                                    + 2 * inputs.quartic_coupling * self.alp_decay_constant ** 2)
        self.alp_mass_squared, self.neutral_pion_mass_squared = self.neutral_poles()
        self.charged_pion_mass_squared = self.b_zero * (self.up_mass + self.down_mass)
        self.legs = {
            'A': self.pole_leg(self.alp_mass_squared, 'A'),
            '+': Leg(math.sqrt(2) * explicit.unit_matrix(1, 0), 0j, 0j, '+'),
            '-': Leg(math.sqrt(2) * explicit.unit_matrix(0, 1), 0j, 0j, '-'),
            '0': self.pole_leg(self.neutral_pion_mass_squared, '0'),
        }
        self.masses_squared = {'A': self.alp_mass_squared, '+': self.charged_pion_mass_squared,
                               '-': self.charged_pion_mass_squared, '0': self.neutral_pion_mass_squared}
        self.pair_invariant = sum(self.masses_squared.values()) / 3
        self.external_dots = self.symmetric_point_dots()
        self.legs['beta'] = self.radial_branch('+', '-')
        self.legs['gamma'] = self.radial_branch('A', '0')
        self._matrix_cache = {}

    def neutral_poles(self):
        """As in the parent class, with brackets and tolerances proportional to f_pi^2."""
        unit_squared = self.pion_decay_constant ** 2
        unmixed = self.unmixed_neutral_masses()
        edges = [unmixed[0] - 1e3 * unit_squared] + list(unmixed) + [unmixed[-1] + 1e6 * unit_squared]
        roots = []
        for low, high in zip(edges[:-1], edges[1:]):
            width = high - low
            roots.append(brentq(self.pole_function, low + 1e-12 * width, high - 1e-12 * width,
                                xtol=1e-14 * unit_squared, rtol=1e-15, maxiter=500))
        alp_like = max(roots, key=lambda s: 1 / self.normalization(s))
        lambda_three = np.diag([1., -1., 0.])

        def isovector_overlap(s):
            meson = 2 * self.epsilon * self.pole_operator_inverse(s) / math.sqrt(self.normalization(s))
            return abs(np.trace(lambda_three @ meson)) / 2

        pion_like = max((s for s in roots if s != alp_like), key=isovector_overlap)
        return alp_like, pion_like

    def vertex(self, legs, dot):
        """Physical n-leg vertex (mass dimension 4 - n): f_pi^(2 - n) times the f_pi = 1 formula."""
        return self.pion_decay_constant ** (2 - len(legs)) * super().vertex(legs, dot)

    def neutral_leg_radial_slope(self, name):
        """d/dp^2 of (1/2) V(l_e(p), l_e(-p), l_sigma, l_sigma) = eps^2 xi_e^2 / f_pi^2."""
        return self.epsilon ** 2 * self.legs[name].alp ** 2 / self.pion_decay_constant ** 2

    def loop_weight(self, lines):
        """(m_d^2)^(2 - r_d) / (16 pi^2), r_d = number of lines, m_d = m_a with an A line, else m_sigma."""
        return LOOP_FACTOR * self.loop_scale_squared(lines) ** (2 - len(lines))

    def diagram_terms(self):
        """{(class, diagram labels): sum_d S_d (m_d^2)^(2 - r)/(16 pi^2) A_d} including the external-leg
        residues (r = 1): class 7 on the a, pi+, pi-, pi0 legs (closed A line, diagrams 7.1-7.4),
        class 8 on the a and pi0 legs (closed sigma line, diagrams 8.1-8.2).  The diagram labels are
        a range such as '4.7-4.14' (see evaluate_explicit_estimate.diagram_groups)."""
        groups, tree = self.group_terms()
        terms = {(name, labels): value for (name, labels, _, _), value in groups.items()}
        return terms, tree

    def group_terms(self):
        """{(class, diagram labels, r, kind): value}, kind = 'cycle' for the listed diagrams and
        'residue' for the external-leg residues; the entries of the printed table."""
        tree = self.tree()[0]
        groups = {}
        for name, labels, propagators, members in diagram_groups():
            value = 0j
            for _, symmetry_factor, blocks, lines in members:
                value += (symmetry_factor * self.loop_weight(lines)
                          * self.cycle(blocks, lines, self.loop_scale_squared(lines)))
            groups[(name, labels, propagators, 'cycle')] = value
        for name, (lines, value) in self.external_residues(tree).items():
            groups[(name, label_range(RESIDUE_DIAGRAM_LABELS[name]), 1, 'residue')] = self.loop_weight(lines) * value
        return groups, tree

    def contributions(self):
        """C_k summed over r for every class, and the tree M0."""
        terms, tree = self.diagram_terms()
        classes = {label: 0j for label in CLASS_ORDER}
        for (name, _), value in terms.items():
            classes[name] = classes.get(name, 0j) + value
        return classes, tree

    def ratio(self):
        classes, tree = self.contributions()
        return abs(sum(classes.values())) / abs(tree)


def channel_open(estimate):
    masses = estimate.masses_squared
    return math.sqrt(masses['A']) > 2 * math.sqrt(masses['+']) + math.sqrt(masses['0'])


def threshold_kappa(low=0.01, high=1.0):
    """kappa_u at which m_a = 2 m_pi+ + m_pi0 (pole data only)."""
    def gap(kappa_u):
        masses = DimensionfulEstimate(kappa_u).masses_squared
        return math.sqrt(masses['A']) - 2 * math.sqrt(masses['+']) - math.sqrt(masses['0'])
    return brentq(gap, low, high, xtol=1e-8)


def ratio_at(kappa_u):
    return DimensionfulEstimate(kappa_u).ratio()


def signed_ratio_at(kappa_u):
    """Re(M1_est / M0); M1_est and M0 are real at the symmetric Dalitz point."""
    estimate = DimensionfulEstimate(kappa_u)
    terms, tree = estimate.diagram_terms()
    return (sum(terms.values()) / tree).real


def crossings(grid, signed_values):
    """kappa_u with R = 1 between neighbouring grid points: sign changes of R - 1, and the pairs of
    roots around every zero of M1_est (a sign change of Re(M1_est/M0) between two points with R > 1)."""
    roots, zeros = [], []
    for low, high in zip(grid[:-1], grid[1:]):
        first, second = signed_values[low], signed_values[high]
        if (abs(first) - 1) * (abs(second) - 1) < 0:
            roots.append(brentq(lambda k: ratio_at(k) - 1, low, high, xtol=1e-9))
        elif first * second < 0 and abs(first) > 1:
            zero = brentq(signed_ratio_at, low, high, xtol=1e-10)
            zeros.append(zero)
            roots.append(brentq(lambda k: ratio_at(k) - 1, low, zero, xtol=1e-9))
            roots.append(brentq(lambda k: ratio_at(k) - 1, zero, high, xtol=1e-9))
    return roots, zeros


def main():
    print('a -> pi+ pi- pi0, up-only model, symmetric Dalitz point, lambda = 1e-3, physical units (GeV).')
    print('NO LOOP COMPUTED: S_d (m_d^2)^(2 - r_d)/(16 pi^2) x tree-vertex trace A_d.\n')
    threshold = threshold_kappa()
    print(f'Channel opens at kappa_u = {threshold:.6f}\n')
    grid = sorted({round(threshold * (1 + 1e-6), 8), 0.1, 0.15, 0.2, 0.205, 0.21, 0.22, 0.23, 0.25,
                   0.3, 0.4, 0.5, 0.7, 1.0})
    signed_values = {}
    for kappa_u in grid:
        estimate = DimensionfulEstimate(kappa_u)
        terms, tree = estimate.diagram_terms()
        signed_values[kappa_u] = (sum(terms.values()) / tree).real
        print(f'kappa_u = {kappa_u:.6f}: eps = {estimate.epsilon:.4f}, m_a = {math.sqrt(estimate.alp_mass_squared):.5f} GeV, '
              f'm_sigma = {math.sqrt(estimate.radial_mass_squared):.5f} GeV, open = {channel_open(estimate)}, '
              f'M0 = {tree.real:+.6e}, R = {abs(signed_values[kappa_u]):.6e}')
        if kappa_u in (0.3, 1.0):
            groups, _ = estimate.group_terms()
            for (name, labels, propagators, kind), value in groups.items():
                print(f'    class {name} diagrams {labels:9s} r = {propagators} {kind:7s}: '
                      f'C/|M0| = {value.real / abs(tree):+.6e} (imaginary {abs(value.imag) / abs(tree):.1e})')
    open_grid = [k for k in grid if k >= threshold]
    roots, zeros = crossings(open_grid, signed_values)
    print('\nzeros of M1_est (accidental cancellation between classes) in the open-channel range:',
          ', '.join(f'{zero:.6f}' for zero in zeros) or 'none')
    print('kappa* with R = 1 in the open-channel range:', ', '.join(f'{root:.6f}' for root in roots) or 'none')


if __name__ == '__main__':
    main()
