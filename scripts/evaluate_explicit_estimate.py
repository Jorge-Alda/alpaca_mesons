#!/usr/bin/env python3
"""Tree-level vertices, pole data, tree amplitude and diagram list for a -> pi+ pi- pi0 (up-only model).

NO LOOP IS COMPUTED.  No loop integral, Feynman parameter, angular average, tensor reduction
or scalar loop function appears.  This module provides, in units f_pi = 1, the ingredients of
the loop-size estimate: the tree-level vertices, the line weights built from the pole data, the
flavour traces A_d = tr[V(B_1) W_1 ... V(B_r) W_r] of the listed diagrams, the closed-line
self-energy slopes of the external-leg residues and the tree M0.  It assigns no loop weight:
the estimate M1_est = sum_d S_d (m_d^2)^(2 - r_d)/(16 pi^2) A_d in physical units is assembled
by scripts/estimate_loop_size_dimensionless.py (result record
memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md).
The diagram list DIAGRAMS is written out below and is the one authoritative list; it contains
only the diagrams allowed by the selection rules, which
scripts/enumerate_murayama_up_only_loop_diagrams.py checks.

Diagram classes (internal lines of the loop) and diagram labels "class.index":
    class 1  one A line, one vertex (r = 1)             diagrams 1.1-1.3
    class 2  one A line and one M line                  diagrams 2.1-2.5
    class 3  two A lines                                diagrams 3.1-3.3
    class 4  one S line and M lines                     diagrams 4.1-4.14 (bubbles 4.1-4.6, triangles 4.7-4.14)
    class 5  one S line, one A line and M lines         diagrams 5.1-5.12 (bubbles 5.1-5.4, triangles 5.5-5.12)
    class 6  two S lines and M lines                    diagrams 6.1-6.7 (bubble 6.1, triangles 6.2-6.4,
                                                        boxes 6.5-6.6, sigma self-energy 6.7)
    class 7  external-leg residues, closed A line       diagrams 7.1 (a), 7.2 (pi+), 7.3 (pi-), 7.4 (pi0)
    class 8  external-leg residues, closed S line       diagrams 8.1 (a), 8.2 (pi0)
    Within a class the index follows the order of DIAGRAMS (classes 1-6) and of
    RESIDUE_DIAGRAM_LABELS (classes 7, 8).

Decisions (2026-11-25):
    class-7 external residues (closed A line) on all four legs a, pi+, pi-, pi0; on the charged legs
    the closed A line needs only the pi+ pi- A A vertex, which charge conservation allows
    (decision of 2026-11-25, replacing the earlier choice without charged-pion residues);
    the sigma^3 vertex of the radial potential is included, V_sigma3 = -6 lambda f_a
    (multilinear coefficient of L = -V, no factor i); it enters class 6 (diagrams 6.1, 6.3, 6.4
    and the sigma self-energy 6.7);
    loop scale m_a for cycles with an A line, m_sigma otherwise.

Conventions (Alda:2026nrd, model of DiLuzio:2026poh, polar fields, connection-free frame):
    L = (f^2/4)<dU dU^dag> + (f^2/4)<chi U^dag + U chi^dag> - (m_1^2/12)<Phi>^2
        + (1/2)(d sigma)^2 + (1/2)(1 + sigma/f_a)^2 (d a)^2 - lambda f_a sigma^3,
    U = exp(i Phi/f),  Phi = sum_i P_i lambda_i,  <lambda_i lambda_j> = 2 delta_ij,
    chi = 2 B0 sum_{r,s} (S_hat_rs + i P_hat_rs) (sigma/f_a)^r (a/f_a)^s,
    up-only: S_hat_rs + i P_hat_rs = m_u (-i)^s/s! e_1 (r = 0, 1) + delta_{r0} delta_{s0} diag(0, m_d, m_s).
    f = f_pi = 1 in the evaluation (all masses in units of f_pi), eps = f_pi/f_a = f_pi kappa_u/(sqrt2 m_u).
    m_sigma^2 = B0 eps^2 (<S_hat_10> - 2<S_hat_20>) + 2 lambda f_a^2, lambda = 1e-3 (labelled default).
Vertices are multilinear coefficients of L (no factor i), all momenta incoming.  A leg is
(meson matrix Lambda, ALP component xi, radial component rho, momentum label).
Momentum rule (prescription, not a loop evaluation): every product of two momenta that contains
a loop slot equals +m_k^2; products of external momenta are the on-shell values at the symmetric
Dalitz point s_{+-} = s_{+0} = s_{-0} = (m_a^2 + 2 m_{pi+}^2 + m_{pi0}^2)/3.
Line weights on the slot basis (E_bc, bare a, pole A, sigma):
    A line: slots -> (l_A, l_A);  S line: (l_sigma, l_sigma);
    M line: 2 sum_{b,c} (E_bc, E_cb) + (l_a, l_a) - (l_A, l_A)
    (U(3) completeness sum_i (lambda_i)_{bc} (lambda_i)_{de} = 2 delta_be delta_cd).

Run:
    .venv/bin/python scripts/evaluate_explicit_estimate.py
Tests:
    PYTHONPATH=scripts .venv/bin/python -m unittest test_evaluate_explicit_estimate
"""
import itertools
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS / 'implementation') not in sys.path:
    sys.path.insert(0, str(_SCRIPTS / 'implementation'))

from only_up_inputs import (B_ZERO, DOWN_MASS, M1_SQUARED, PION_DECAY_CONSTANT,  # noqa: E402
                            STRANGE_MASS, UP_MASS)

DEFAULT_QUARTIC_COUPLING = 1e-3
LOOP = 'loop'


def unit_matrix(first_index, second_index):
    """3 x 3 matrix E_bc with a single entry 1 at (b, c) = (first_index, second_index)."""
    matrix = np.zeros((3, 3), dtype=complex)
    matrix[first_index, second_index] = 1.
    return matrix


E_UP = unit_matrix(0, 0)


@dataclass(frozen=True)
class Leg:
    """Meson matrix (coefficient of Phi/f), ALP component, radial component, momentum label."""
    meson: np.ndarray
    alp: complex
    radial: complex
    momentum: str


class ExplicitEstimate:
    """Pole data, legs, vertices and class sums of the explicit formulas at one kappa_u."""

    def __init__(self, kappa_u, quartic_coupling=DEFAULT_QUARTIC_COUPLING):
        self.kappa_u = kappa_u
        self.pion_decay_constant = 1.     # all masses in units of f_pi in this class
        self.quartic_coupling = quartic_coupling
        self.epsilon = PION_DECAY_CONSTANT * kappa_u / (math.sqrt(2) * UP_MASS)
        self.b_zero = B_ZERO / PION_DECAY_CONSTANT
        self.up_mass, self.down_mass, self.strange_mass = (m / PION_DECAY_CONSTANT for m in
                                                            (UP_MASS, DOWN_MASS, STRANGE_MASS))
        self.singlet_mass_squared = M1_SQUARED / PION_DECAY_CONSTANT ** 2
        self.scalar_sources, self.pseudoscalar_sources = self.up_only_sources()
        self.source_c = self.b_zero * self.pseudoscalar_sources[(0, 1)]
        alp_decay_constant = 1 / self.epsilon
        self.radial_mass_squared = (self.b_zero * self.epsilon ** 2
                                    * np.trace(self.scalar_sources[(1, 0)]
                                               - 2 * self.scalar_sources.get((2, 0), 0 * E_UP)).real
                                    + 2 * quartic_coupling * alp_decay_constant ** 2)
        self.alp_mass_squared, self.neutral_pion_mass_squared = self.neutral_poles()
        self.charged_pion_mass_squared = self.b_zero * (self.up_mass + self.down_mass)
        self.legs = {
            'A': self.pole_leg(self.alp_mass_squared, 'A'),
            '+': Leg(math.sqrt(2) * unit_matrix(1, 0), 0j, 0j, '+'),
            '-': Leg(math.sqrt(2) * unit_matrix(0, 1), 0j, 0j, '-'),
            '0': self.pole_leg(self.neutral_pion_mass_squared, '0'),
        }
        masses = {'A': self.alp_mass_squared, '+': self.charged_pion_mass_squared,
                  '-': self.charged_pion_mass_squared, '0': self.neutral_pion_mass_squared}
        self.masses_squared = masses
        self.pair_invariant = (masses['A'] + masses['+'] + masses['-'] + masses['0']) / 3
        self.external_dots = self.symmetric_point_dots()
        self.legs['beta'] = self.radial_branch('+', '-')
        self.legs['gamma'] = self.radial_branch('A', '0')
        self._matrix_cache = {}

    # ------------------------------------------------------------------ sources and poles
    def up_only_sources(self):
        """{(r, s): S_hat_rs}, {(r, s): P_hat_rs} of the up-only model, r <= 1, s <= 6."""
        scalar, pseudoscalar = {}, {}
        for r in (0, 1):
            for s in range(7):
                value = self.up_mass * (-1j) ** s / math.factorial(s)
                scalar[(r, s)] = value.real * E_UP
                pseudoscalar[(r, s)] = value.imag * E_UP
        scalar[(0, 0)] = scalar[(0, 0)] + np.diag([0., self.down_mass, self.strange_mass])
        return scalar, pseudoscalar

    def pole_operator(self, mass_squared, matrix):
        """M_s(X) = B0 {S_hat_00, X} + (m_1^2/3) <X> 1 - s X."""
        scalar = self.scalar_sources[(0, 0)]
        return (self.b_zero * (scalar @ matrix + matrix @ scalar)
                + self.singlet_mass_squared / 3 * np.trace(matrix) * np.eye(3) - mass_squared * matrix)

    def pole_operator_inverse(self, mass_squared):
        """M^{-1}_s(C), C = B0 P_hat_01, by solving the 9 x 9 linear system on 3 x 3 matrices."""
        columns = np.array([self.pole_operator(mass_squared, unit_matrix(b, c)).ravel()
                            for b in range(3) for c in range(3)]).T
        return np.linalg.solve(columns, self.source_c.ravel()).reshape(3, 3)

    def normalization(self, mass_squared):
        """N(s) = 1 + 2 eps^2 <M^{-1}_s(C)^2>."""
        inverse = self.pole_operator_inverse(mass_squared)
        return 1 + 2 * self.epsilon ** 2 * np.trace(inverse @ inverse).real

    def pole_function(self, mass_squared):
        """F(s) = s + 2 B0 eps^2 <S_hat_02> + 2 eps^2 <C M^{-1}_s(C)>; poles are the roots, N = F'."""
        inverse = self.pole_operator_inverse(mass_squared)
        return (mass_squared + 2 * self.b_zero * self.epsilon ** 2 * np.trace(self.scalar_sources[(0, 2)]).real
                + 2 * self.epsilon ** 2 * np.trace(self.source_c @ inverse).real)

    def unmixed_neutral_masses(self):
        """Eigenvalues of M_0 on diagonal matrices (unmixed pi0, eta, eta')."""
        basis = [unit_matrix(i, i) for i in range(3)]
        matrix = np.array([[self.pole_operator(0., b)[i, i].real for b in basis] for i in range(3)])
        return np.sort(np.linalg.eigvals(matrix).real)

    def neutral_poles(self):
        """(m_a^2, m_pi0^2): the root of F with the largest ALP component 1/N, and the root with
        the largest lambda_3 component among the others."""
        unmixed = self.unmixed_neutral_masses()
        edges = [min(-1., unmixed[0] - 1e3)] + list(unmixed) + [unmixed[-1] + 1e6]
        roots = []
        for low, high in zip(edges[:-1], edges[1:]):
            width = high - low
            roots.append(brentq(self.pole_function, low + 1e-12 * width, high - 1e-12 * width,
                                xtol=1e-14, rtol=1e-15, maxiter=500))
        alp_like = max(roots, key=lambda s: 1 / self.normalization(s))
        lambda_three = np.diag([1., -1., 0.])

        def isovector_overlap(s):
            meson = 2 * self.epsilon * self.pole_operator_inverse(s) / math.sqrt(self.normalization(s))
            return abs(np.trace(lambda_three @ meson)) / 2

        pion_like = max((s for s in roots if s != alp_like), key=isovector_overlap)
        return alp_like, pion_like

    def pole_leg(self, mass_squared, momentum):
        """Leg of a neutral pole: meson matrix 2 eps M^{-1}(C)/sqrt(N), ALP component 1/sqrt(N)."""
        norm = self.normalization(mass_squared)
        return Leg(2 * self.epsilon * self.pole_operator_inverse(mass_squared) / math.sqrt(norm),
                   1 / math.sqrt(norm) + 0j, 0j, momentum)

    # ------------------------------------------------------------------ momenta
    def symmetric_point_dots(self):
        """k_i.k_j of incoming external momenta at s_{+-} = s_{+0} = s_{-0} = s_0."""
        m, s0 = self.masses_squared, self.pair_invariant
        dots = {}
        for first, second in itertools.combinations(('A', '+', '-', '0'), 2):
            value = (s0 - m[first] - m[second]) / 2
            dots[(first, second)] = dots[(second, first)] = value
        return dots

    def dot_rule(self, scale_squared, extra=None):
        """Momentum products: slot products = scale_squared, external products at the symmetric point."""
        def dot(first, second):
            if first == LOOP or second == LOOP:
                return scale_squared
            if extra and (first, second) in extra:
                return extra[(first, second)]
            return self.external_dots[(first, second)]
        return dot

    # ------------------------------------------------------------------ vertex
    def vertex(self, legs, dot):
        """V(l_1..l_n) = potential + meson kinetic + ALP kinetic (+ singlet and radial kinetic at n = 2,
        + radial cubic sigma^3 at n = 3)."""
        return (self.potential_vertex(legs) + self.meson_kinetic_vertex(legs, dot)
                + self.alp_kinetic_vertex(legs, dot) + self.two_leg_terms(legs, dot)
                + self.radial_cubic_vertex(legs))

    def radial_cubic_vertex(self, legs):
        """sigma^3 vertex of the radial potential V = (m_phi^2/2)(f_a + sigma)^2 + (lambda/4)(f_a + sigma)^4.

        L contains -lambda f_a sigma^3; the vertex is the multilinear coefficient of L (no factor i),
        V_sigma3 = d^3 L/d sigma^3 = -6 lambda f_a, times rho_1 rho_2 rho_3.  It is written here as
        -6 lambda f_pi^2/eps = -6 lambda f_a f_pi, the mass-dimension-2 form of the f_pi = 1 formulas
        (DimensionfulEstimate.vertex multiplies by f_pi^(2 - 3) and gives -6 lambda f_a).
        """
        if len(legs) != 3:
            return 0j
        coupling = -6 * self.quartic_coupling * self.pion_decay_constant ** 2 / self.epsilon
        return coupling * legs[0].radial * legs[1].radial * legs[2].radial

    def potential_vertex(self, legs):
        """B0 sum_{(R,A,P)} r! s! eps^{r+s} prod rho prod xi (-1)^{floor(m/2)}/m! sum_pi <X_rs Lambda_pi>,
        X_rs = S_hat_rs (m even), P_hat_rs (m odd)."""
        roles_per_leg = []
        for leg in legs:
            roles = []
            if leg.radial != 0:
                roles.append('R')
            if leg.alp != 0:
                roles.append('A')
            if np.any(leg.meson != 0):
                roles.append('P')
            roles_per_leg.append(roles)
        total = 0j
        for assignment in itertools.product(*roles_per_leg):
            r, s = assignment.count('R'), assignment.count('A')
            if (r, s) not in self.scalar_sources:
                continue
            mesons = [leg.meson for leg, role in zip(legs, assignment) if role == 'P']
            m = len(mesons)
            source = self.scalar_sources[(r, s)] if m % 2 == 0 else self.pseudoscalar_sources[(r, s)]
            if not np.any(source):
                continue
            factor = (math.factorial(r) * math.factorial(s) * self.epsilon ** (r + s)
                      * (-1) ** (m // 2) / math.factorial(m))
            for leg, role in zip(legs, assignment):
                factor = factor * (leg.radial if role == 'R' else leg.alp if role == 'A' else 1)
            ordered = sum((np.trace(source @ np.linalg.multi_dot([np.eye(3)] + [mesons[i] for i in order] + [np.eye(3)]))
                           for order in itertools.permutations(range(m))), 0j)
            total += self.b_zero * factor * ordered
        return total

    @staticmethod
    def ordered_sum(matrices):
        """sum over orderings pi of Lambda_pi1 ... Lambda_pim (unit matrix for m = 0)."""
        if not matrices:
            return np.eye(3, dtype=complex)
        return sum((np.linalg.multi_dot([np.eye(3)] + [matrices[i] for i in order] + [np.eye(3)])
                    for order in itertools.permutations(range(len(matrices)))), np.zeros((3, 3), complex))

    def meson_kinetic_vertex(self, legs, dot):
        """(1/4) sum_{T, R} (-K_T.K_R) i^{|T|} (-i)^{|R|}/(|T|! |R|!) <(sum_pi Lambda_pi^T)(sum_pi Lambda_pi^R)>."""
        n = len(legs)
        if n < 2 or n % 2 or any(not np.any(leg.meson) for leg in legs):
            return 0j
        total = 0j
        for size in range(1, n):
            for first in itertools.combinations(range(n), size):
                second = [j for j in range(n) if j not in first]
                momentum_product = sum(dot(legs[i].momentum, legs[j].momentum) for i in first for j in second)
                coefficient = (-momentum_product * 1j ** size * (-1j) ** (n - size)
                               / (math.factorial(size) * math.factorial(n - size)))
                total += 0.25 * coefficient * np.trace(self.ordered_sum([legs[i].meson for i in first])
                                                       @ self.ordered_sum([legs[j].meson for j in second]))
        return total

    def alp_kinetic_vertex(self, legs, dot):
        """-sum_{i<j} xi_i xi_j k_i.k_j C(rest), C() = 1, C(m) = 2 eps rho_m, C(m, m') = 2 eps^2 rho_m rho_m'."""
        n = len(legs)
        if n < 2 or n > 4:
            return 0j
        total = 0j
        for i, j in itertools.combinations(range(n), 2):
            if legs[i].alp == 0 or legs[j].alp == 0:
                continue
            rest = [legs[k].radial for k in range(n) if k not in (i, j)]
            radial_factor = {0: 1., 1: 2 * self.epsilon * (rest[0] if rest else 0),
                             2: 2 * self.epsilon ** 2 * (rest[0] * rest[1] if len(rest) == 2 else 0)}[len(rest)]
            total += -legs[i].alp * legs[j].alp * dot(legs[i].momentum, legs[j].momentum) * radial_factor
        return total

    def two_leg_terms(self, legs, dot):
        """Singlet term -(m_1^2/6)<Lambda_1><Lambda_2> and radial kinetic -rho_1 rho_2 k_1.k_2 (n = 2 only)."""
        if len(legs) != 2:
            return 0j
        first, second = legs
        value = -self.singlet_mass_squared / 6 * np.trace(first.meson) * np.trace(second.meson)
        if first.radial != 0 and second.radial != 0:
            value -= first.radial * second.radial * dot(first.momentum, second.momentum)
        return value

    # ------------------------------------------------------------------ tree and branches
    def radial_branch(self, first, second):
        """Radial leg rho = V(l_sigma, l_u, l_v)/(m_sigma^2 - s_uv) entering the loop (momentum k_u + k_v)."""
        sigma = Leg(np.zeros((3, 3), complex), 0j, 1. + 0j, 'sigma')
        value = self.vertex([sigma, self.legs[first], self.legs[second]], self.dot_rule(0.))
        return Leg(np.zeros((3, 3), complex), 0j, value / (self.radial_mass_squared - self.pair_invariant),
                   'branch')

    def tree(self):
        """M0 = V(A,+,-,0) + V(sigma,A,0) V(sigma,+,-)/(m_sigma^2 - s_{+-}) (other channels vanish)."""
        dot = self.dot_rule(0.)
        legs = self.legs
        sigma = Leg(np.zeros((3, 3), complex), 0j, 1. + 0j, 'sigma')
        contact = self.vertex([legs['A'], legs['+'], legs['-'], legs['0']], dot)
        exchange = (self.vertex([sigma, legs['A'], legs['0']], dot) * self.vertex([sigma, legs['+'], legs['-']], dot)
                    / (self.radial_mass_squared - self.pair_invariant))
        return contact + exchange, contact, exchange

    # ------------------------------------------------------------------ slots and lines
    def slot_legs(self):
        """Slot basis: E_bc (9), bare ALP l_a, pole l_A, radial l_sigma; all carry the loop momentum."""
        zero = np.zeros((3, 3), complex)
        basis = [Leg(unit_matrix(b, c), 0j, 0j, LOOP) for b in range(3) for c in range(3)]
        alp_pole = self.legs['A']
        basis += [Leg(zero, 1. + 0j, 0j, LOOP), Leg(alp_pole.meson, alp_pole.alp, 0j, LOOP),
                  Leg(zero, 0j, 1. + 0j, LOOP)]
        return basis

    @staticmethod
    def line_weight(line):
        """Line weights W_L on the slot basis (index 9 = bare a, 10 = pole A, 11 = sigma)."""
        weight = np.zeros((12, 12))
        if line == 'A':
            weight[10, 10] = 1.
        elif line == 'S':
            weight[11, 11] = 1.
        elif line == 'M':
            for b in range(3):
                for c in range(3):
                    weight[3 * b + c, 3 * c + b] = 2.     # U(3) completeness
            weight[9, 9], weight[10, 10] = 1., -1.
        else:
            raise ValueError(line)
        return weight

    def vertex_matrix(self, block, scale_squared):
        """V(B; x, y) on the slot basis at the class scale."""
        key = (block, scale_squared)
        if key not in self._matrix_cache:
            fixed = [self.legs[name] for name in block]
            slots = self.slot_legs()
            dot = self.dot_rule(scale_squared)
            matrix = np.zeros((12, 12), dtype=complex)
            for x in range(12):
                for y in range(x, 12):
                    matrix[x, y] = matrix[y, x] = self.vertex(fixed + [slots[x], slots[y]], dot)
            self._matrix_cache[key] = matrix
        return self._matrix_cache[key]

    def cycle(self, blocks, lines, scale_squared):
        """tr[V(B_1) W_{L_1} V(B_2) W_{L_2} ... V(B_n) W_{L_n}] (line L_j joins B_j and B_{j+1})."""
        product = np.eye(12, dtype=complex)
        for block, line in zip(blocks, lines):
            product = product @ self.vertex_matrix(block, scale_squared) @ self.line_weight(line)
        return np.trace(product)

    # ------------------------------------------------------------------ self-energy slopes
    def neutral_leg_radial_slope(self, name):
        """d/dp^2 of (1/2) V(l_e(p), l_e(-p), l_sigma, l_sigma) = eps^2 xi_e^2 (closed form)."""
        return self.epsilon ** 2 * self.legs[name].alp ** 2

    def closed_alp_line_slope(self, first, second):
        """d/dp^2 of (1/2) V(l_first(p), l_second(-p), l_A, l_A): affine in k_first.k_second = -p^2.

        Closed-A-line self-energy slope of the class-7 external residues: (first, second) = (A, A)
        (diagram 7.1), (+, -) for the pi+ leg (7.2), (-, +) for the pi- leg (7.3), (0, 0) (7.4).
        """
        legs = [Leg(self.legs[first].meson, self.legs[first].alp, 0j, 'p'),
                Leg(self.legs[second].meson, self.legs[second].alp, 0j, 'q')]
        alp_pole = self.legs['A']
        slots = [Leg(alp_pole.meson, alp_pole.alp, 0j, LOOP)] * 2
        values = []
        for momentum_squared in (0., 1.):
            extra = {('p', 'q'): -momentum_squared, ('q', 'p'): -momentum_squared}
            values.append(0.5 * self.vertex(legs + slots, self.dot_rule(self.alp_mass_squared, extra)))
        return values[1] - values[0]

    def charged_leg_alp_slope(self):
        """d/dp^2 of (1/2) V(l_+(p), l_-(-p), l_A, l_A) = -eps^2 [(M^-1)_11 - (M^-1)_22]^2/(3 N f_pi^2),
        M^-1 = M^{-1}_{m_a^2}(C), N = N(m_a^2) (closed form, checked in the unit tests)."""
        return self.closed_alp_line_slope('+', '-')

    # ------------------------------------------------------------------ classes
    def loop_scale_squared(self, lines):
        """m_d^2: m_a^2 if the cycle has an A line, m_sigma^2 otherwise."""
        return self.alp_mass_squared if 'A' in lines else self.radial_mass_squared

    def external_residues(self, tree):
        """{class: (lines of the closed loop, (1/2) sum_e dSigma_ee/dp^2 M0)} without the loop weight.

        The diagrams of a class are summed (labels in RESIDUE_DIAGRAM_LABELS).
        Class 7: closed A line on every external leg, e = A, pi+, pi-, pi0 (diagrams 7.1-7.4; the
        pi+ and pi- legs use the pi+ pi- A A vertex).
        Class 8: closed sigma line on the sigma sigma a a vertex, e = A, pi0 (diagrams 8.1, 8.2; the
        sigma sigma pi+ pi- vertex vanishes in the up-only model, r <= 1 in the sources; see
        test_vanishing_vertices).
        """
        alp_line_slopes = sum(self.closed_alp_line_slope(first, second)
                              for first, second in (('A', 'A'), ('+', '-'), ('-', '+'), ('0', '0')))
        return {
            '7': (('A',), 0.5 * alp_line_slopes * tree),
            '8': (('S',), 0.5 * (self.neutral_leg_radial_slope('A') + self.neutral_leg_radial_slope('0')) * tree),
        }


# ---------------------------------------------------------------------- explicit diagram list
# The one authoritative list of the estimate: (diagram label "class.index", symmetry factor
# S_d = 1/|Aut_d|, blocks of external legs or radial branches at each vertex, lines).  Line j joins block j and block j+1
# (cyclically).  beta = radial branch of (pi+, pi-), gamma = radial branch of (A, pi0).
# Only the diagrams allowed by intrinsic parity, charge and the radial vertex content are listed;
# scripts/enumerate_murayama_up_only_loop_diagrams.py checks that every allowed cycle appears
# exactly once.  Diagrams marked sigma^3 contain the vertex V_sigma3 = -6 lambda f_a.
DIAGRAMS = (
    # class 1: closed A line at one vertex (r = 1), diagrams 1.1-1.3
    ('1.1', 0.5, (('A', '+', '-', '0'),), ('A',)),
    ('1.2', 0.5, (('beta', 'A', '0'),), ('A',)),
    ('1.3', 0.5, (('gamma', '+', '-'),), ('A',)),
    # class 2: one A line and one M line (r = 2), diagrams 2.1-2.5
    ('2.1', 1.0, (('A', '+'), ('-', '0')), ('A', 'M')),
    ('2.2', 1.0, (('A', '-'), ('+', '0')), ('A', 'M')),
    ('2.3', 1.0, (('A', '0'), ('+', '-')), ('A', 'M')),
    ('2.4', 1.0, (('beta',), ('A', '0')), ('A', 'M')),
    ('2.5', 1.0, (('gamma',), ('+', '-')), ('A', 'M')),
    # class 3: two A lines (r = 2), diagrams 3.1-3.3
    ('3.1', 0.5, (('A', '0'), ('+', '-')), ('A', 'A')),
    ('3.2', 0.5, (('beta',), ('A', '0')), ('A', 'A')),
    ('3.3', 0.5, (('gamma',), ('+', '-')), ('A', 'A')),
    # class 4, bubbles: one S line and one M line (r = 2), diagrams 4.1-4.6
    ('4.1', 1.0, (('A',), ('+', '-', '0')), ('S', 'M')),
    ('4.2', 1.0, (('+',), ('A', '-', '0')), ('S', 'M')),
    ('4.3', 1.0, (('-',), ('A', '+', '0')), ('S', 'M')),
    ('4.4', 1.0, (('0',), ('A', '+', '-')), ('S', 'M')),
    ('4.5', 1.0, (('A',), ('beta', '0')), ('S', 'M')),
    ('4.6', 1.0, (('0',), ('beta', 'A')), ('S', 'M')),
    # class 4, triangles: one S line and two M lines (r = 3), diagrams 4.7-4.14
    ('4.7', 1.0, (('A', '+'), ('-',), ('0',)), ('M', 'S', 'M')),
    ('4.8', 1.0, (('A', '-'), ('+',), ('0',)), ('M', 'S', 'M')),
    ('4.9', 1.0, (('A', '0'), ('+',), ('-',)), ('M', 'S', 'M')),
    ('4.10', 1.0, (('+', '-'), ('A',), ('0',)), ('M', 'S', 'M')),
    ('4.11', 1.0, (('+', '0'), ('A',), ('-',)), ('M', 'S', 'M')),
    ('4.12', 1.0, (('-', '0'), ('A',), ('+',)), ('M', 'S', 'M')),
    ('4.13', 1.0, (('beta',), ('A',), ('0',)), ('M', 'S', 'M')),
    ('4.14', 1.0, (('gamma',), ('+',), ('-',)), ('M', 'S', 'M')),
    # class 5, bubbles: one S line and one A line (r = 2), diagrams 5.1-5.4
    ('5.1', 1.0, (('A',), ('+', '-', '0')), ('S', 'A')),
    ('5.2', 1.0, (('0',), ('A', '+', '-')), ('S', 'A')),
    ('5.3', 1.0, (('A',), ('beta', '0')), ('S', 'A')),
    ('5.4', 1.0, (('0',), ('beta', 'A')), ('S', 'A')),
    # class 5, triangles: one S, one A and one M line (r = 3), diagrams 5.5-5.12
    ('5.5', 1.0, (('A', '+'), ('-',), ('0',)), ('M', 'S', 'A')),
    ('5.6', 1.0, (('A', '-'), ('+',), ('0',)), ('M', 'S', 'A')),
    ('5.7', 1.0, (('+', '-'), ('A',), ('0',)), ('A', 'S', 'M')),
    ('5.8', 1.0, (('+', '-'), ('A',), ('0',)), ('M', 'S', 'A')),
    ('5.9', 1.0, (('+', '0'), ('A',), ('-',)), ('A', 'S', 'M')),
    ('5.10', 1.0, (('-', '0'), ('A',), ('+',)), ('A', 'S', 'M')),
    ('5.11', 1.0, (('beta',), ('A',), ('0',)), ('A', 'S', 'M')),
    ('5.12', 1.0, (('beta',), ('A',), ('0',)), ('M', 'S', 'A')),
    # class 6, bubble: two S lines (r = 2), sigma^3 vertex at beta, diagram 6.1
    ('6.1', 0.5, (('beta',), ('A', '0')), ('S', 'S')),
    # class 6, triangles: two S lines and one M line (r = 3), diagrams 6.2-6.4; 6.3 and 6.4 with a sigma^3 vertex
    ('6.2', 1.0, (('A', '0'), ('+',), ('-',)), ('S', 'M', 'S')),
    ('6.3', 1.0, (('beta',), ('A',), ('0',)), ('S', 'M', 'S')),
    ('6.4', 1.0, (('gamma',), ('+',), ('-',)), ('S', 'M', 'S')),
    # class 6, boxes: two S lines and two M lines, alternating (r = 4), diagrams 6.5-6.6
    ('6.5', 1.0, (('A',), ('+',), ('-',), ('0',)), ('S', 'M', 'S', 'M')),
    ('6.6', 1.0, (('+',), ('-',), ('A',), ('0',)), ('M', 'S', 'M', 'S')),
    # class 6, sigma self-energy: self-energy of the internal sigma line of the tree exchange,
    # two sigma^3 vertices (r = 2), diagram 6.7
    ('6.7', 0.5, (('beta',), ('gamma',)), ('S', 'S')),
)
# Labels of the external-leg residues (one diagram per external leg, summed in external_residues).
RESIDUE_DIAGRAM_LABELS = {
    '7': ('7.1', '7.2', '7.3', '7.4'),      # closed A line on the a, pi+, pi-, pi0 legs
    '8': ('8.1', '8.2'),                    # closed S line on the a, pi0 legs
}


def diagram_class(label):
    """Class of a diagram label "class.index", e.g. '6' for '6.7'."""
    return label.split('.')[0]


def label_range(labels):
    """'first-last' for several consecutive diagram labels, the label itself for one."""
    return labels[0] if len(labels) == 1 else f'{labels[0]}-{labels[-1]}'


def diagram_groups():
    """[(class, label range, r, diagrams)]: consecutive entries of DIAGRAMS with the same class and
    the same number of lines r (bubbles, triangles, boxes of a class, and the sigma self-energy 6.7
    separately from the bubble 6.1)."""
    groups = []
    for diagram in DIAGRAMS:
        label, lines = diagram[0], diagram[3]
        key = (diagram_class(label), len(lines))
        if groups and groups[-1][0] == key:
            groups[-1][1].append(diagram)
        else:
            groups.append((key, [diagram]))
    return [(name, label_range([d[0] for d in members]), propagators, tuple(members))
            for (name, propagators), members in groups]


def main():
    print('a -> pi+ pi- pi0, up-only model, symmetric Dalitz point, lambda = 1e-3, units f_pi = 1.')
    print('NO LOOP COMPUTED: pole data, tree M0 and the flavour traces A_d of the listed diagrams.')
    print('The loop-size estimate is assembled by scripts/estimate_loop_size_dimensionless.py.\n')
    for kappa_u in (0.3, 1.0):
        estimate = ExplicitEstimate(kappa_u)
        tree, contact, exchange = estimate.tree()
        print(f'kappa_u = {kappa_u}: eps = {estimate.epsilon:.6f}, m_a = '
              f'{math.sqrt(estimate.alp_mass_squared) * PION_DECAY_CONSTANT:.6f} GeV, '
              f'm_pi0 = {math.sqrt(estimate.neutral_pion_mass_squared) * PION_DECAY_CONSTANT:.6f} GeV, '
              f'm_sigma = {math.sqrt(estimate.radial_mass_squared) * PION_DECAY_CONSTANT:.6f} GeV, '
              f'N(m_a^2) = {estimate.normalization(estimate.alp_mass_squared):.6f}')
        print(f'  M0 = {tree.real:+.10e} (contact {contact.real:+.6e}, radial exchange {exchange.real:+.6e})')
        for label, weight, blocks, lines in DIAGRAMS:
            trace = estimate.cycle(blocks, lines, estimate.loop_scale_squared(lines))
            print(f'  {label:5s} S_d = {weight:.1f} {lines} {blocks}: A_d/|M0| = {trace.real / abs(tree):+.6e}')
        print()


if __name__ == '__main__':
    main()
