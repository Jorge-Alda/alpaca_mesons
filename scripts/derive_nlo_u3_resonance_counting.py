#!/usr/bin/env python3
"""Assign delta orders to the U(3) NLO Lagrangian, resonance exchange and loops.

Purpose
    Power counting for the OpenSpec change `fix-u3-nlo-resonance-counting`
    (tasks 2.2, 3.1, 4.1 and 6.1).
    --orders: tabulates the order in the combined expansion parameter delta
    of Kaiser:1998ds and Kaiser:2000gs for every operator of
    memory/nlo/lagrangian.md, for the resonance couplings and masses of
    Ecker:1988te, and for the low-energy limit of single-resonance exchange.
    --loops: classifies the one-loop topologies (and a two-loop reference
    row) that contribute to Goldstone/ALP two-point functions and to ALP
    decay vertices, gives their delta order in the low-energy limit and
    their status at the declared order delta^2 ("included",
    "excluded: beyond declared order", "excluded: declared approximation",
    "excluded: first order in 1/f_a", or, for the resonance self-energy
    near its pole, "excluded from strict order; line-shape prescription
    near pole"). It also prints the lowest vector + pseudoscalar
    thresholds, where loops with an internal resonance line become
    non-analytic.
    --scheme (task 4.1): states the chiral subtraction scheme
    L_i = L_i^r(mu) + Gamma_i lambda(mu), tabulates the U(3) coefficients
    Gamma_i of Kaiser:1998ds next to the SU(3) ones, implements the running
    of L_i^r(mu) and checks that the bare L_i is mu-independent, lists the
    combinations of L_i with zero Gamma (independent of the chiral scale
    mu_chi). It gives the chiral renormalization of B,
    B = B^r(1 + 4 tau lambda/F^4), and asserts that the pole term is of
    order delta^2 (same as a one-loop G-LO graph); it asserts that the part
    -v_31/2 of Gamma_8 dropped by Kaiser:1998ds is delta^3. It tabulates
    separately the dependence on the QCD renormalization scale through
    Z_A, checks that L_A = (2 L_5^r + 3 L_18^r)/sqrt(1 + Lambda_1) is
    independent of both scales, and asserts that the failure of the
    multiplicative Z_A law for renormalized combinations with a pole
    (2 L_8 - 3 L_25, B) is (Z_A - 1) x pole = O(delta^3).
    --truncation [SQRT_S ...] (task 6.1): truncation uncertainty of an
    observable computed at the declared order delta^2, per unit |Delta|,
    where Delta is the relative NLO shift of the observable (delta^1 plus
    delta^2 terms), supplied by the downstream calculation. It prints the
    three terms of delta_est(s) = max(1/N_c, s/Lambda_chi^2,
    M_K^2/Lambda_chi^2) and the largest one, u_trunc/|Delta| = delta_est,
    the uncertainty from the excluded resonance loops (class R-int of
    --loops) u_Rloop/|Delta| = eps_R theta, and the quadrature total
    u_tot/|Delta| = sqrt(delta_est^2 + eps_R^2 theta), for both settings of
    inputs_from_data; and two flags, "in declared range" (sqrt(s) <= 1.0 GeV)
    and "delta_est < 1". It prints the crossover Lambda_chi/sqrt(N_c) where
    s/Lambda_chi^2 = 1/N_c, and asserts that delta_est is non-decreasing in
    s, that the M_K term never dominates, that 0.7 < delta_est(1 GeV) < 0.8,
    that delta_est < 1 iff sqrt(s) < Lambda_chi above the crossover, the
    threshold rule for u_Rloop and that the total is not smaller than each
    term. Default grid: sqrt(s) = 0.3, 0.5, 0.7, 0.9, 1.0, 1.2 GeV.
    Neither --orders nor --loops expands a propagator or evaluates a loop
    integral; they only classify orders. --scheme evaluates no loop either;
    it takes the Gamma_i as inputs.

Conventions
    - delta counting: p = O(sqrt(delta)), m_q = O(delta), 1/N_c = O(delta).
    - A term with N_c power k and chiral dimension 2n is of order
      delta^(n - k). The leading-order Lagrangian is delta^0.
    - Chiral dimension: each derivative, external field strength f_pm^{mu nu}
      or u_mu counts p^1; each chi_pm (proportional to m_q) counts p^2;
      X_TC = -(psi + theta), and the ALP through c_G a / f_a inside X_TC,
      count p^0. Resonance masses squared M_R^2 count p^0 N_c^0.
    - N_c power: read from the explicit prefactor of each term
      (F^2 ~ N_c, f_pi = F at this order). For the Goldstone Lagrangian it is
      cross-checked against the rule k = 2 - (number of flavour traces
      + number of X_TC factors), i.e. the pure-gluon vacuum energy is N_c^2,
      a single flavour trace is N_c, and each further flavour trace or
      factor of X_TC = -(psi + theta) costs 1/N_c (psi = -i<log U> is itself
      a flavour trace, and theta enters as theta / N_c).
    - sqrt(N_c) is represented by k = 1/2 (exact sympy Rationals).
    - The ALP field and its sources are kept at first order in 1/f_a; 1/f_a
      does not enter delta, and the ALP mass squared counts as p^2.
    - Spin-one resonances are antisymmetric tensor fields (Ecker:1988te),
      so their couplings to Goldstones start at p^2.
    - Loop rules (--loops), stated as assumptions of the counting:
      a loop whose internal lines are all Goldstones (eta' included, since
      m_0^2 = O(delta)) has loop momenta set by the external momenta and
      masses and multiplies the tree graph built from the same vertices by
      p^2/(16 pi^2 F^2), i.e. (k, 2n) += (-1, 2). A loop with an internal
      resonance line has loop momenta of order M_R and multiplies by
      N_c^{-1} without p^2, i.e. (k, 2n) += (-1, 0). Since the loop
      momentum ~ M_R supplies the momentum factors of the vertices (e.g.
      G_V^2 M_V^2/(16 pi^2 F^4) relative to the tree kinetic term), the
      lowest R-int contribution to a Goldstone two-point function is the
      generic N_c^{-1} loop factor times the chirally allowed minimum p^2
      (or m_q) of the F^2-normalized two-point function: N_c^0 p^2, i.e.
      delta^1. This is the O(p^2) divergence of the V pi and S pi loops in
      the pion self-energy (Rosell:2004mn, subsection "Pion self-energy",
      delta Z_pi ~ 6 G_V^2 M_V^2/(16 pi^2 F^4) plus a c_d^2 M_S^2 term),
      which multiplies <u.u> and <chi_+> (gamma_1, gamma_2 of
      Rosell:2005ai), so it renormalizes F^2 and F^2 B_0. The O(p^4) part
      (tree exchange N_c p^4 times N_c^{-1}) is delta^2. A resonance
      self-energy Sigma_R relative to M_R^2 counts N_c^{-1} p^4 in the
      low-energy limit and N_c^{-1} near the pole (Gamma_R/M_R ~ 1/N_c).
      A loop with an internal ALP line is O(1/f_a^2); the 1/f_a^2
      suppression is not part of delta, so no delta order is assigned to it.
    - Declared order: delta^2, first order in 1/f_a, Goldstone loops only
      (loops with internal resonance lines are a declared approximation).
    - Scheme (--scheme): d-dimensional regularization with
      lambda(mu) = mu^(d-4)/(16 pi^2) [1/(d-4) - (log 4 pi + Gamma'(1) + 1)/2]
      (Gasser:1983yg convention, as used by Kaiser:1998ds and
      Herrera-Siklody:1996tqr Eq. (3.20)), L_i = L_i^r(mu)
      + Gamma_i lambda(mu); hence mu dL_i^r/dmu = -Gamma_i/(16 pi^2).
      Kaiser:1998ds writes mu_chi for this chiral subtraction scale; it is
      independent of the QCD renormalization scale mu of the singlet axial
      current, mu dZ_A/dmu = gamma_A Z_A, Z_A = 1 + O(1/N_c).
    - In the --scheme counting lambda ~ 1/(16 pi^2) counts N_c^0 p^0; the
      1/N_c of a loop is carried by 1/F^2. Z_A - 1 counts N_c^-1 p^0.
    - Truncation (--truncation): the omitted terms are O(delta^3), one
      factor of delta beyond the delta^2 terms, so u_trunc = |Delta| delta_est.
      N_c = 3, Lambda_chi = 4 pi f_pi, f_pi = F (project convention,
      130.2 MeV / sqrt(2)). theta = 1 if the excluded resonance loops enter:
      for inputs_from_data=True (F, B_0, tau, Lambda_1, Lambda_2, L_i^r(mu),
      resonance masses and couplings fixed from data, short-distance
      constraints imposed; the condition stated in --loops) only the
      non-analytic part of R-int is missing, so theta = 1 only for sqrt(s)
      >= m_rho + m_pi = 0.91483 GeV (lowest R pi threshold); for
      inputs_from_data=False (large-N_c, saturation or LO short-distance
      estimates) theta = 1 at all sqrt(s).
      ASSUMPTION: eps_R = 0.30. Basis: the NLO-in-1/N_c shifts of the
      central values of L_10^r and C_87^r are about 17% and 28%
      (Pich:2008jm), and the second multiplets V', A' account for about
      30-40% of the correction there. The two uncertainties are added in
      quadrature (assumption: independent). Exact arithmetic in sympy
      (Rationals and pi); comparisons at 50 digits.
      Near a resonance pole (|M_R^2 - s| ~ M_R Gamma_R) the counting does not
      apply; see the line-shape prescription record under Inputs.

Units
    Natural units. The order tables use only powers of N_c and of the
    chiral dimension. Masses and thresholds in --loops are in GeV. In
    --scheme mu is in GeV (log mu means log(mu / 1 GeV); the reference
    scale drops out of every difference), and shifts of L_i^r are printed
    in units of 10^-3. In --truncation sqrt(s), f_pi, M_K, Lambda_chi and
    thresholds are in GeV; delta_est and all uncertainties are
    dimensionless, per unit |Delta|.

Inputs
    memory/nlo/lagrangian.md (L^(0), L^(1), L^(2) operators, Alda:2026nrd
    conventions), Kaiser:1998ds, Kaiser:2000gs (delta counting),
    Ecker:1988te (resonance couplings and their N_c scaling),
    Rosell:2004mn, Rosell:2005ai, Rosell:2006dt, Pich:2008jm, Pich:2010sm
    (one-loop resonance chiral theory: loops with internal resonances, their
    O(p^2) divergences, their absorption into resonance parameters through
    NLO short-distance constraints, and their R pi and R R' cuts),
    memory/rcht/foundations/vector-width-lineshape-prescriptions.md (line
    shape near the pole).
    --scheme (line numbers refer to the source TeX files, see
    memory/extracts/Kaiser:1998ds/u3-loop-renormalization-and-gamma-table.md,
    memory/extracts/Kaiser:1998ds/qcd-scale-dependence-za.md and
    memory/extracts/Herrera-Siklody:1996tqr/u3-one-loop-counterterms-and-gamma-table.md):
    Gasser:1983yg (convention for lambda(mu)); Kaiser:1998ds adehep.tex
    (lambda: lines 674-678; B renormalization: lines 679-686; no
    renormalization of F, tau, Lambda_1, Lambda_2: lines 687-689; U(3)
    Gamma_i: lines 692-696; Z_A laws: lines 583-594, 600, 621-634;
    invariant L_A: lines 777-781); Herrera-Siklody:1996tqr un.tex
    (lambda_infinity, Eq. (3.20); SU(3) Gamma_6, Gamma_8: lines 1507-1508;
    B_8 at next order in 1/N_c: lines 2044-2050; Gamma_47 = 0 by large-N_c
    truncation; --orders: O_47 at line 2516, theta^ = i theta at lines
    674-676, C_mu, M at lines 1208-1210, N_c rule (4.1) at lines 1904-1907); Gasser:1984gg (SU(3) Gamma_i). Gamma_25 = 0 is quoted
    by Kaiser:1998ds and not reproduced from Herrera-Siklody:1996tqr
    Appendix D. v_31 = -Lambda_2/3 from
    memory/extracts/Herrera-Siklody:1996tqr/u3-notation-dictionary-with-kaiser-1998ds.md.
    Masses from ParticleDataGroup:2024cfk, in GeV:
    M_rho = 0.77526, M_omega = 0.78266, M_K*(892)^+- = 0.89167,
    m_pi^+- = 0.13957, m_K^+- = 0.493677, m_eta = 0.547862.
    --truncation: f_pi = 0.0921 GeV (130.2 MeV / sqrt(2),
    ParticleDataGroup:2024cfk), M_K = m_K^+- above, M_rho and m_pi^+- above
    for the rho pi threshold; eps_R = 0.30 (assumption, from Pich:2008jm);
    exclusion condition of R-int from --loops; Delta from the downstream
    calculation (not an input of this script).

Invocation
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --orders
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --loops
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --scheme
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --truncation          # default grid
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --truncation 0.5 1.0  # sqrt(s) in GeV
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py      # all checks
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --test
"""

import argparse
import math
import sys
import unittest
from dataclasses import dataclass
from itertools import product

import sympy as sp

delta = sp.Symbol("delta", positive=True)
HALF = sp.Rational(1, 2)


@dataclass(frozen=True)
class CountedTerm:
    """One row of the order table.

    nc_power is the exponent k of N_c, chiral_dimension is 2n.
    flavor_traces and topological_factors (powers of X_TC) are None for rows
    that are not Goldstone Lagrangian operators (resonance couplings,
    exchange graphs, loops), for which the trace rule does not apply.
    """

    name: str
    description: str
    nc_power: sp.Rational
    chiral_dimension: int
    flavor_traces: object = None
    topological_factors: object = None
    expected_exponent: object = None


def delta_order(nc_power, chiral_dimension):
    """Return delta^(n - k) for N_c power k and chiral dimension 2n."""
    return delta ** (sp.Rational(chiral_dimension, 2) - sp.Rational(nc_power))


def delta_exponent(term):
    """Return the exponent of delta of a CountedTerm as an exact Rational."""
    return sp.Rational(term.chiral_dimension, 2) - sp.Rational(term.nc_power)


def nc_power_from_traces(flavor_traces, topological_factors):
    """N_c power of a Goldstone operator from its trace structure.

    k = 2 - (flavour traces + powers of X_TC); see the module docstring.
    """
    return 2 - flavor_traces - topological_factors


def combine(name, description, factors, expected_exponent=None):
    """Combine factors (vertices, propagator terms) by adding k and 2n."""
    return CountedTerm(
        name=name,
        description=description,
        nc_power=sum((sp.Rational(factor.nc_power) for factor in factors), sp.Integer(0)),
        chiral_dimension=sum(factor.chiral_dimension for factor in factors),
        expected_exponent=expected_exponent,
    )


def goldstone_lagrangian_rows():
    """Operators of L^(0), L^(1), L^(2) in memory/nlo/lagrangian.md."""
    return [
        # Leading order, delta^0.
        CountedTerm("F^2 kinetic", "F^2/4 <u_mu u^mu>", 1, 2, 1, 0, 0),
        CountedTerm("F^2 mass", "F^2/4 <chi_+>", 1, 2, 1, 0, 0),
        # tau = f_pi^2 m_0^2 / (4 N_f) is the topological susceptibility of
        # pure Yang-Mills, O(N_c^0); hence m_0^2 = O(1/N_c) = O(delta).
        CountedTerm("tau X_TC^2", "-(f_pi^2 m_0^2/(4 N_f)) X_TC^2", 0, 0, 0, 2, 0),
        # L^(1), delta^1.
        CountedTerm("L_5", "L_5 <u_mu u^mu chi_+>", 1, 4, 1, 0, 1),
        CountedTerm("L_8", "(L_8/2) <chi_+^2 + chi_-^2>", 1, 4, 1, 0, 1),
        # Lambda_1, Lambda_2 ~ 1/N_c, times F^2 ~ N_c: net N_c^0.
        CountedTerm("Lambda_1", "(f_pi^2/12) Lambda_1 <u_mu><u^mu>", 0, 2, 2, 0, 1),
        CountedTerm("Lambda_2", "i (f_pi^2/12) Lambda_2 X_TC <chi_->", 0, 2, 1, 1, 1),
        # L^(2), delta^2.
        CountedTerm("L_4", "L_4 <u_mu u^mu><chi_+>", 0, 4, 2, 0, 2),
        CountedTerm("L_6", "L_6 <chi_+>^2", 0, 4, 2, 0, 2),
        # L_7 is O(N_c^0) in U(3). In SU(3) ChPT eta' exchange is integrated
        # out, L_7 ~ F^2/m_0^2 ~ N_c^2, and the counting is different.
        CountedTerm("L_7", "L_7 <chi_->^2", 0, 4, 2, 0, 2),
        CountedTerm("L_18", "L_18 <u_mu><u^mu chi_+>", 0, 4, 2, 0, 2),
        # Single trace times one X_TC: L_25 ~ N_c^0 (not N_c).
        CountedTerm("L_25", "-i L_25 X_TC <chi_+ chi_->", 0, 4, 1, 1, 2),
        CountedTerm("N_c p^6", "generic single-trace p^6 term", 1, 6, 1, 0, 2),
    ]


def source_derivative_rows():
    """Source-derivative operator O_47 of Herrera-Siklody:1996tqr (un.tex line 2516).

    O_47 = D_mu theta^ <C^mu><M>, theta^ = i theta, C^mu = U^dagger D^mu U,
    M = U^dagger chi + chi^dagger U (un.tex lines 674-676, 1208-1210).
    Chiral dimension: D_mu theta^ (p^1) + C^mu (p^1) + M (p^2) = p^4.
    Flavour traces: <C^mu> and <M>, i.e. 2. One factor of theta^, which,
    like X_TC, costs 1/N_c (theta enters QCD as theta/N_c, un.tex line 365).
    Hence k = 2 - 2 - 1 = -1, which is also rule (4.1) of
    Herrera-Siklody:1996tqr, G = N_c^(2 - #traces - #theta^) (un.tex lines
    1904-1907), and the order is delta^(2 + 1) = delta^3.
    """
    return [
        CountedTerm("L_47", "L_47 D_mu theta^ <C^mu><M>", -1, 4, 2, 1, 3),
    ]


def resonance_rows():
    """Resonance masses and couplings (Ecker:1988te), and OZI-violating parts."""
    sqrt_nc = HALF
    rows = [
        CountedTerm("M_R^2", "resonance mass squared", 0, 0, expected_exponent=0),
        CountedTerm("c_d", "c_d <S u_mu u^mu>", sqrt_nc, 2),
        CountedTerm("c_m", "c_m <S chi_+>", sqrt_nc, 2),
        CountedTerm("d_m", "i d_m <P chi_->", sqrt_nc, 2),
        CountedTerm("F_V", "F_V/(2 sqrt2) <V_{mu nu} f_+^{mu nu}>", sqrt_nc, 2),
        CountedTerm("G_V", "i G_V/sqrt2 <V_{mu nu} u^mu u^nu>", sqrt_nc, 2),
        CountedTerm("F_A", "F_A/(2 sqrt2) <A_{mu nu} f_-^{mu nu}>", sqrt_nc, 2),
        # Singlet couplings equal the nonet values in the large-N_c limit.
        CountedTerm("c~_d", "singlet c~_d S_1 <u_mu u^mu>", sqrt_nc, 2),
        CountedTerm("c~_m", "singlet c~_m S_1 <chi_+>", sqrt_nc, 2),
        CountedTerm("d~_m", "singlet i d~_m P_1 <chi_->", sqrt_nc, 2),
        # OZI-violating difference from the nonet value: extra 1/N_c.
        CountedTerm("OZI diff.", "c~ - c/sqrt(3) (any S, P coupling)", sqrt_nc - 1, 2),
    ]
    for row in rows[1:]:
        object.__setattr__(row, "expected_exponent", 1 - row.nc_power)
    return rows


def exchange_rows():
    """Low-energy limit of single-resonance exchange.

    The propagator 1/(M_R^2 - s) = 1/M_R^2 + s/M_R^4 + ... is only used to
    classify orders. At the declared order the full propagator is kept.
    """
    vertex = CountedTerm("vertex", "", HALF, 2)
    ozi_vertex = CountedTerm("OZI vertex", "", HALF - 1, 2)
    leading_propagator = CountedTerm("1/M_R^2", "", 0, 0)
    next_propagator = CountedTerm("s/M_R^4", "", 0, 2)
    return [
        combine("exch. 1/M^2", "two vertices x 1/M_R^2 (N_c p^4)",
                [vertex, vertex, leading_propagator], 1),
        combine("exch. s/M^4", "two vertices x s/M_R^4 (N_c p^6)",
                [vertex, vertex, next_propagator], 2),
        combine("exch. OZI", "one OZI-violating vertex x 1/M_R^2",
                [vertex, ozi_vertex, leading_propagator], 2),
    ]


def pseudoscalar_mixing_rows():
    """P-Goldstone bilinear mixing induced by d_m <P chi_->.

    <chi_-> contains -i {Phi/F, chi_0}, so the bilinear coefficient is
    B ~ d_m chi_0 / F: (k, 2n) = (1/2 - 1/2, 2) = (0, 2).
    Integrating out P at tree level shifts the mass matrix of the
    canonically normalized Goldstone fields by B^2/M_R^2 ~ N_c^0 p^4. The
    Lagrangian of memory/nlo/lagrangian.md is written in terms of Phi/F
    (u_mu, chi_pm), so the corresponding coefficient there is
    F^2 B^2/M_R^2 ~ N_c p^4: an L_8-type term, order delta. The two
    estimates differ by the factor F^2 ~ N_c of the field normalization,
    and the delta order must be read from the F^2-normalized form.
    """
    mixing = CountedTerm("B (P-Phi mix)", "B ~ d_m chi_0 / F (canonical fields)",
                         HALF - HALF, 2, expected_exponent=1)
    canonical_shift = CountedTerm("B^2/M^2 canon.", "Goldstone mass shift, canonical fields",
                                  0, 4)
    lagrangian_shift = CountedTerm("F^2 B^2/M^2", "same, F^2-normalized coefficient",
                                   1, 4, expected_exponent=1)
    return [mixing, canonical_shift, lagrangian_shift]


def loop_rows():
    """One-loop Goldstone graph with LO vertices.

    Each loop adds p^2/(16 pi^2 F^2) relative to the LO N_c p^2 term,
    giving N_c^0 p^4.
    """
    loop_number = 1
    return [
        CountedTerm("1-loop Goldstone", "one loop, LO vertices (N_c^0 p^4)",
                    1 - loop_number, 2 + 2 * loop_number, expected_exponent=2),
    ]


def all_rows():
    """Every row of the order table, in printing order."""
    return (goldstone_lagrangian_rows() + source_derivative_rows()
            + resonance_rows() + exchange_rows()
            + pseudoscalar_mixing_rows() + loop_rows())


def format_table(rows):
    """Format the order table as text."""
    header = f"{'name':<17}{'operator':<42}{'k':>6}{'2n':>4}{'tr':>4}{'X':>3}  order"
    lines = [header, "-" * len(header)]
    for row in rows:
        traces = "-" if row.flavor_traces is None else str(row.flavor_traces)
        topological = "-" if row.topological_factors is None else str(row.topological_factors)
        lines.append(
            f"{row.name:<17}{row.description:<42}{str(row.nc_power):>6}"
            f"{row.chiral_dimension:>4}{traces:>4}{topological:>3}  "
            f"{sp.sstr(delta_order(row.nc_power, row.chiral_dimension))}"
        )
    return "\n".join(lines)


def check_orders():
    """Print the table and assert every expected order. Return True on success."""
    rows = all_rows()
    print("delta counting: p ~ sqrt(delta), m_q ~ delta, 1/N_c ~ delta;"
          " order = delta^(n-k)")
    print("k = N_c power, 2n = chiral dimension, tr = flavour traces,"
          " X = powers of X_TC\n")
    print(format_table(rows))

    for row in rows:
        if row.flavor_traces is not None:
            from_traces = nc_power_from_traces(row.flavor_traces, row.topological_factors)
            assert sp.Rational(row.nc_power) == from_traces, (
                f"{row.name}: explicit N_c power {row.nc_power} != trace rule {from_traces}")
        if row.expected_exponent is not None:
            assert delta_exponent(row) == sp.Rational(row.expected_exponent), (
                f"{row.name}: delta^{delta_exponent(row)} != delta^{row.expected_exponent}")

    by_name = {row.name: row for row in rows}
    for name in ["L_5", "L_8", "Lambda_1", "Lambda_2"]:
        assert delta_exponent(by_name[name]) == 1, name
    for name in ["L_4", "L_6", "L_7", "L_18", "L_25"]:
        assert delta_exponent(by_name[name]) == 2, name
    assert delta_exponent(by_name["exch. 1/M^2"]) == 1
    assert delta_exponent(by_name["exch. s/M^4"]) == delta_exponent(by_name["N_c p^6"]) == 2
    assert delta_exponent(by_name["F^2 B^2/M^2"]) == delta_exponent(by_name["exch. 1/M^2"])
    assert delta_exponent(by_name["F^2 B^2/M^2"]) == delta_exponent(by_name["L_8"])
    assert delta_exponent(by_name["1-loop Goldstone"]) == 2
    for name in ["L_4", "L_6", "L_7", "L_18", "L_25"]:
        assert delta_exponent(by_name[name]) == delta_exponent(by_name["1-loop Goldstone"])

    # L_47: the N_c power from rule (4.1) of Herrera-Siklody:1996tqr
    # (theta^ counted as a source) equals the one from the Alda:2026nrd form
    # D_mu theta^ = i(<u_mu> - d_mu X_TC): <u_mu><u^mu><chi_+> has 3 traces
    # and no X_TC, d_mu X_TC <u^mu><chi_+> has 2 traces and one X_TC.
    l47_nc_power = by_name["L_47"].nc_power
    assert l47_nc_power == nc_power_from_traces(2, 1) == -1
    assert nc_power_from_traces(3, 0) == nc_power_from_traces(2, 1) == l47_nc_power
    assert delta_exponent(by_name["L_47"]) == 3 > DECLARED_EXPONENT

    # m_0^2 = 4 N_f tau / F^2 with tau ~ N_c^0, F^2 ~ N_c and chiral dimension 0.
    m0_squared_exponent = 0 - (0 - 1)
    assert m0_squared_exponent == 1

    print("\nNotes:")
    print("- m_0^2 = 4 N_f tau / F^2 ~ N_c^0 / N_c: O(delta^1).")
    print("- L_7 is O(N_c^0) in U(3); in SU(3) eta' exchange makes it O(N_c^2).")
    print("- Exchange: the expansion of 1/(M_R^2 - s) only classifies orders."
          " At the declared order the full propagator is kept unexpanded,"
          " since s/M_R^2 is not small for sqrt(s) <~ 1 GeV.")
    print("- P-Phi mixing: B^2/M_R^2 ~ N_c^0 p^4 in canonical fields equals"
          " F^2 B^2/M_R^2 ~ N_c p^4 (L_8-type, delta^1) in the F^2-normalized"
          " Lagrangian. The two canonical-field rows are offset by one power"
          " of delta (canonical LO mass term N_c^0 p^2 -> delta); only"
          " their ratio to LO, delta^1, is meaningful.")
    print("- L_47 (O_47 = D_mu theta^ <C^mu><M>, Herrera-Siklody:1996tqr"
          " un.tex line 2516): 2n = 1 + 1 + 2 = 4, 2 traces, one theta^."
          " N_c power: module rule k = 2 - tr - X = -1; rule (4.1) of"
          " Herrera-Siklody:1996tqr N_c^(2 - #tf - #theta^) = N_c^-1 (lines"
          " 1904-1907); they agree. Herrera-Siklody:1996tqr does not quote"
          " the order of L_47 explicitly (its list, lines 1986-2008, covers"
          " only operators without derivatives of theta^). Order delta^3,"
          " beyond the declared delta^2. D_mu theta^ counts p^1, unlike"
          " X_TC (p^0).")
    print("- ALP at first order in 1/f_a: 1/f_a is not part of delta;"
          " c_G a/f_a in X_TC counts p^0, m_a^2 counts p^2.")
    print("\n--orders: all assertions passed")
    return True


class OrderTableTests(unittest.TestCase):
    """Unit tests of the counting functions and of the order table."""

    def test_delta_order_function(self):
        self.assertEqual(delta_order(1, 2), 1)
        self.assertEqual(delta_order(HALF, 2), delta ** HALF)
        self.assertEqual(delta_order(0, 4), delta ** 2)

    def test_trace_rule_matches_explicit_prefactors(self):
        for row in goldstone_lagrangian_rows():
            self.assertEqual(
                sp.Rational(row.nc_power),
                nc_power_from_traces(row.flavor_traces, row.topological_factors),
                row.name)

    def test_expected_orders(self):
        for row in all_rows():
            if row.expected_exponent is not None:
                self.assertEqual(delta_exponent(row), row.expected_exponent, row.name)

    def test_exchange_combination(self):
        vertex = CountedTerm("v", "", HALF, 2)
        combined = combine("x", "", [vertex, vertex, CountedTerm("p", "", 0, 0)])
        self.assertEqual((combined.nc_power, combined.chiral_dimension), (1, 4))

    def test_check_orders_runs(self):
        self.assertTrue(check_orders())

    def test_l47_row(self):
        (row,) = source_derivative_rows()
        self.assertEqual(row.name, "L_47")
        self.assertEqual((row.chiral_dimension, row.flavor_traces,
                          row.topological_factors), (1 + 1 + 2, 2, 1))
        self.assertEqual(row.nc_power, nc_power_from_traces(2, 1))
        self.assertEqual(row.nc_power, -1)
        # Alda:2026nrd form of D_mu theta^: <u_mu> term (3 traces, no X_TC).
        self.assertEqual(nc_power_from_traces(3, 0), row.nc_power)
        self.assertEqual(delta_exponent(row), 3)
        self.assertEqual(delta_order(row.nc_power, row.chiral_dimension), delta ** 3)
        self.assertGreater(delta_exponent(row), DECLARED_EXPONENT)
        self.assertIn(row, all_rows())


# --------------------------------------------------------------------------
# Loop topologies (--loops, task 3.1)
# --------------------------------------------------------------------------

DECLARED_EXPONENT = 2
INCLUDED = "included"
BEYOND_ORDER = "excluded: beyond declared order"
DECLARED_APPROXIMATION = "excluded: declared approximation"
FIRST_ORDER_IN_INVERSE_FA = "excluded: first order in 1/f_a"
NEAR_POLE_PRESCRIPTION = "excluded from strict order; line-shape prescription near pole"

# PDG 2024 masses in GeV (ParticleDataGroup:2024cfk).
MASSES_GEV = {
    "rho": 0.77526,
    "omega": 0.78266,
    "K*": 0.89167,
    "pi": 0.13957,
    "K": 0.493677,
    "eta": 0.547862,
}
VECTOR_PSEUDOSCALAR_PAIRS = [("rho", "pi"), ("omega", "pi"), ("K*", "K"),
                             ("rho", "eta"), ("omega", "eta")]
VALIDITY_LIMIT_GEV = 1.0


@dataclass(frozen=True)
class LoopClass:
    """One loop topology class.

    nc_power (k) and chiral_dimension (2n) give the low-energy order
    delta^(n - k); both are None when no delta order applies (internal ALP
    line, resonance self-energy near its pole).
    """

    name: str
    description: str
    nc_power: object
    chiral_dimension: object
    status: str
    regime: str = "low energy"

    def exponent(self):
        """Exponent of delta, or None if no delta order is assigned."""
        if self.nc_power is None:
            return None
        return sp.Rational(self.chiral_dimension, 2) - sp.Rational(self.nc_power)


def loop_factors():
    """Multiplicative (k, 2n) factors used to build the loop classes."""
    return {
        "LO vertex tree": CountedTerm("LO tree", "L^(0) tree graph", 1, 2),
        "N_c p^4 vertex tree": CountedTerm("L_5, L_8 tree", "L^(1) N_c p^4 vertex", 1, 4),
        "N_c^0 p^2 vertex tree": CountedTerm("Lambda tree", "L^(1) N_c^0 p^2 vertex", 0, 2),
        # Chirally allowed minimum of the F^2-normalized Goldstone two-point
        # function: F^2 p^2 or F^2 B_0 m_q, i.e. N_c p^2.
        "minimal two-point tree": CountedTerm("F^2 p^2 tree", "F^2 x minimal p^2 (or m_q)", 1, 2),
        "resonance exchange": CountedTerm("exchange", "two vertices x 1/M_R^2", 1, 4),
        "Goldstone loop": CountedTerm("G loop", "p^2/(16 pi^2 F^2)", -1, 2),
        "resonance loop": CountedTerm("R loop", "1/(16 pi^2 F^2) x M_R^2", -1, 0),
        "self-energy, low energy": CountedTerm("Sigma_R/M_R^2", "N_c^-1 p^4", -1, 4),
    }


def classify_status(exponent, declared_approximation=False):
    """Status of a class with delta exponent `exponent` at declared order."""
    if exponent > DECLARED_EXPONENT:
        return BEYOND_ORDER
    if declared_approximation:
        return DECLARED_APPROXIMATION
    return INCLUDED


def loop_class_from_factors(name, description, factor_names, declared_approximation=False):
    """Build a LoopClass by multiplying the named loop factors."""
    factors = loop_factors()
    combined = combine(name, description, [factors[label] for label in factor_names])
    exponent = delta_exponent(combined)
    return LoopClass(name, description, combined.nc_power, combined.chiral_dimension,
                     classify_status(exponent, declared_approximation))


def loop_classes():
    """All loop topology classes, in printing order."""
    return [
        loop_class_from_factors(
            "G-LO", "Goldstone-only 1 loop, L^(0) vertices (eta' incl.)",
            ["LO vertex tree", "Goldstone loop"]),
        # The ALP source vertex (c_q^A, c_q^V, c_q^S, c_q^P, c_G) is a vertex
        # of L^(0); 1/f_a is not part of delta.
        loop_class_from_factors(
            "G-ALP", "as G-LO with one L^(0) ALP-source vertex",
            ["LO vertex tree", "Goldstone loop"]),
        loop_class_from_factors(
            "G-NLO", "Goldstone-only 1 loop, one L_5 or L_8 vertex",
            ["N_c p^4 vertex tree", "Goldstone loop"]),
        loop_class_from_factors(
            "G-NLO", "Goldstone-only 1 loop, one Lambda_1/2 vertex",
            ["N_c^0 p^2 vertex tree", "Goldstone loop"]),
        # R-int, O(p^2) part: the loop momentum ~ M_R supplies the momentum
        # factors of the resonance vertices, so the loop gives the generic
        # N_c^{-1} factor (e.g. G_V^2 M_V^2/(16 pi^2 F^4)) times the minimal
        # p^2 (or m_q) of the F^2-normalized two-point function. These are
        # the O(p^2) divergences of the V pi and S pi loops in the pion
        # self-energy (Rosell:2004mn, subsection "Pion self-energy") and
        # gamma_1, gamma_2 of <u.u>, <chi_+> (Rosell:2005ai): they
        # renormalize F^2 and F^2 B_0.
        loop_class_from_factors(
            "R-int", "1 loop with internal resonance line(s): O(p^2) part"
            " (renormalizes F^2, F^2 B_0)",
            ["minimal two-point tree", "resonance loop"], declared_approximation=True),
        # R-int, O(p^4) part: tree exchange (N_c p^4) times N_c^{-1}.
        loop_class_from_factors(
            "R-int", "same, O(p^4) part",
            ["resonance exchange", "resonance loop"], declared_approximation=True),
        loop_class_from_factors(
            "R-vertex", "G loop on resonance vertex in tree exchange",
            ["resonance exchange", "Goldstone loop"]),
        loop_class_from_factors(
            "R-dress", "G loop in resonance propagator (Sigma_R)",
            ["resonance exchange", "self-energy, low energy"]),
        LoopClass("R-dress", "Sigma_R/M_R^2 ~ N_c^-1 (Gamma_R/M_R)", None, None,
                  NEAR_POLE_PRESCRIPTION, regime="s ~ M_R^2"),
        LoopClass("A-int", "loop with internal ALP line, O(1/f_a^2)", None, None,
                  FIRST_ORDER_IN_INVERSE_FA),
        loop_class_from_factors(
            "G-2loop", "Goldstone-only 2 loops, L^(0) vertices",
            ["LO vertex tree", "Goldstone loop", "Goldstone loop"]),
    ]


def vector_pseudoscalar_thresholds():
    """Sorted list of (label, threshold in GeV) for vector + pseudoscalar."""
    thresholds = [(f"{vector}{pseudoscalar}", MASSES_GEV[vector] + MASSES_GEV[pseudoscalar])
                  for vector, pseudoscalar in VECTOR_PSEUDOSCALAR_PAIRS]
    return sorted(thresholds, key=lambda item: item[1])


def format_loop_table(classes):
    """Format the loop classification as text."""
    width = max(52, max(len(loop_class.description) for loop_class in classes) + 2)
    header = (f"{'class':<9}{'description':<{width}}{'regime':<12}{'k':>4}{'2n':>4}"
              f"  {'order':<9}status")
    lines = [header, "-" * (len(header) + 40)]
    for loop_class in classes:
        if loop_class.nc_power is None:
            nc_text, dimension_text, order_text = "-", "-", "-"
        else:
            nc_text = str(loop_class.nc_power)
            dimension_text = str(loop_class.chiral_dimension)
            order_text = sp.sstr(delta_order(loop_class.nc_power, loop_class.chiral_dimension))
        lines.append(f"{loop_class.name:<9}{loop_class.description:<{width}}{loop_class.regime:<12}"
                     f"{nc_text:>4}{dimension_text:>4}  {order_text:<9}{loop_class.status}")
    return "\n".join(lines)


def check_loops():
    """Print the loop classification and thresholds; assert the rules."""
    classes = loop_classes()
    print("Loop topologies for Goldstone/ALP two-point functions and ALP decay"
          " vertices")
    print(f"declared order: delta^{DECLARED_EXPONENT}, first order in 1/f_a,"
          " Goldstone loops only")
    print("Goldstone loop: (k, 2n) += (-1, 2); resonance loop: (k, 2n) += (-1, 0)\n")
    print(format_loop_table(classes))

    included_names = {item.name for item in classes if item.status == INCLUDED}
    assert included_names == {"G-LO", "G-ALP"}, included_names
    for item in classes:
        if item.status == BEYOND_ORDER:
            assert item.exponent() > DECLARED_EXPONENT, item.name
        if item.status == INCLUDED:
            assert item.exponent() == DECLARED_EXPONENT, item.name
    by_name = {}
    for item in classes:
        by_name.setdefault(item.name, []).append(item)
    expected_low_energy = {"G-LO": 2, "G-ALP": 2, "G-NLO": 3,
                           "R-vertex": 3, "R-dress": 4, "G-2loop": 4}
    for name, exponent in expected_low_energy.items():
        rows = [item for item in by_name[name] if item.regime == "low energy"]
        assert rows and all(item.exponent() == exponent for item in rows), name
    assert len(by_name["G-NLO"]) == 2
    assert all(item.status == BEYOND_ORDER for item in by_name["G-NLO"])
    resonance_internal = by_name["R-int"]
    assert sorted(item.exponent() for item in resonance_internal) == [1, 2]
    assert min(item.exponent() for item in resonance_internal) < DECLARED_EXPONENT
    assert all(item.status == DECLARED_APPROXIMATION for item in resonance_internal)
    near_pole = [item for item in by_name["R-dress"] if item.regime != "low energy"]
    assert len(near_pole) == 1 and near_pole[0].status == NEAR_POLE_PRESCRIPTION
    (alp_internal,) = by_name["A-int"]
    assert alp_internal.status == FIRST_ORDER_IN_INVERSE_FA

    thresholds = vector_pseudoscalar_thresholds()
    print("\nVector + pseudoscalar thresholds (GeV, ParticleDataGroup:2024cfk):")
    for label, value in thresholds:
        print(f"  {label:<9}{value:.6f}")
    threshold_by_label = dict(thresholds)
    assert threshold_by_label["rhopi"] < VALIDITY_LIMIT_GEV
    assert [value for _, value in thresholds] == sorted(value for _, value in thresholds)

    print("\nNotes:")
    print("- R-int starts at delta^1: resonance loops renormalize the O(p^2)"
          " constants F^2 and F^2 B_0 at 1/N_c (Rosell:2004mn, Rosell:2005ai),"
          " and the O(p^4) constants at delta^2. It is excluded by declared"
          " choice (Goldstone loops only).")
    print("- Assumption: below the lowest resonance + Goldstone threshold R-int is"
          " analytic in s and m_q, so its effect has the form of local terms"
          " at O(p^2) and O(p^4). Sources absorb it into the Goldstone"
          " constants (F^2, F^2 B_0, pion wave function: Rosell:2004mn,"
          " Rosell:2005ai) and into the resonance parameters (F_V, F_A, M_V,"
          " M_A via NLO Weinberg sum rules, c_m, d_m, G_V: Pich:2008jm,"
          " Rosell:2006dt). The exclusion is consistent only if F, B_0, tau"
          " (m_0^2), Lambda_1, Lambda_2, L_i^r(mu) and the resonance masses"
          " and couplings are fixed from data, and if the short-distance"
          " constraints hold at NLO. Its non-analytic part starts at the"
          " R pi threshold, and its size enters the truncation uncertainty.")
    print("- Resonance-resonance cuts lie at (M_R+M_R')^2 >~ 2 GeV^2"
          " (Pich:2010sm), outside the validity range; only R pi (and R K,"
          " R eta) cuts can lie below 1 GeV. In the pseudoscalar channel"
          " V pi and S pi contribute (Rosell:2006dt).")
    print(f"- The declared validity range (<~ {VALIDITY_LIMIT_GEV:g} GeV) overlaps the"
          f" excluded R-int non-analytic region above ~{thresholds[0][1]:.2f} GeV"
          " (rho pi) in channels where a vector-pseudoscalar state is allowed."
          " Quantum-number selection rules are not applied here.")
    print("- R-dress near s ~ M_R^2: Sigma_R/M_R^2 ~ Gamma_R/M_R ~ 1/N_c with no p"
          " suppression, so it is not small near the pole. It enters only"
          " through the line-shape prescription of memory/rcht/foundations/"
          "vector-width-lineshape-prescriptions.md, as a resummation outside"
          " the strict counting.")
    print("- A-int: O(1/f_a^2), excluded at first order in 1/f_a; 1/f_a is not"
          " part of delta, so no delta order is assigned.")
    print("\n--loops: all assertions passed")
    return True


class LoopClassificationTests(unittest.TestCase):
    """Unit tests of the loop classification and of the thresholds."""

    def setUp(self):
        self.classes = loop_classes()

    def rows(self, name, regime="low energy"):
        return [item for item in self.classes if item.name == name and item.regime == regime]

    def test_included_set(self):
        included = {item.name for item in self.classes if item.status == INCLUDED}
        self.assertEqual(included, {"G-LO", "G-ALP"})

    def test_goldstone_one_loop_is_nc0_p4(self):
        for name in ["G-LO", "G-ALP"]:
            (row,) = self.rows(name)
            self.assertEqual((row.nc_power, row.chiral_dimension), (0, 4))

    def test_nlo_vertex_loops_beyond_order_for_both_vertex_types(self):
        rows = self.rows("G-NLO")
        self.assertEqual({(row.nc_power, row.chiral_dimension) for row in rows},
                         {(0, 6), (-1, 4)})
        for row in rows:
            self.assertEqual(row.exponent(), 3)
            self.assertEqual(row.status, BEYOND_ORDER)

    def test_beyond_order_strictly_above_declared(self):
        for item in self.classes:
            if item.status == BEYOND_ORDER:
                self.assertGreater(item.exponent(), DECLARED_EXPONENT, item.name)

    def test_resonance_internal_starts_below_declared_order(self):
        rows = self.rows("R-int")
        self.assertEqual({(row.nc_power, row.chiral_dimension) for row in rows},
                         {(0, 2), (0, 4)})
        lowest = min(row.exponent() for row in rows)
        self.assertEqual(lowest, 1)
        self.assertLess(lowest, DECLARED_EXPONENT)
        for row in rows:
            self.assertEqual(row.status, DECLARED_APPROXIMATION)

    def test_resonance_internal_rule_matches_goldstone_two_point_minimum(self):
        # N_c^-1 loop factor times the F^2-normalized minimum N_c p^2.
        factors = loop_factors()
        combined = combine("x", "", [factors["minimal two-point tree"],
                                     factors["resonance loop"]])
        self.assertEqual(delta_exponent(combined), 1)
        self.assertEqual(delta_exponent(combined),
                         delta_exponent(next(row for row in goldstone_lagrangian_rows()
                                             if row.name == "Lambda_1")))

    def test_resonance_vertex_and_dressing(self):
        (vertex_row,) = self.rows("R-vertex")
        self.assertEqual((vertex_row.nc_power, vertex_row.chiral_dimension, vertex_row.exponent()),
                         (0, 6, 3))
        (dress_row,) = self.rows("R-dress")
        self.assertEqual((dress_row.nc_power, dress_row.chiral_dimension, dress_row.exponent()),
                         (0, 8, 4))
        (near_pole,) = self.rows("R-dress", regime="s ~ M_R^2")
        self.assertEqual(near_pole.status, NEAR_POLE_PRESCRIPTION)
        self.assertIsNone(near_pole.exponent())

    def test_alp_internal_and_two_loop(self):
        (alp_row,) = self.rows("A-int")
        self.assertEqual(alp_row.status, FIRST_ORDER_IN_INVERSE_FA)
        self.assertIsNone(alp_row.exponent())
        (two_loop,) = self.rows("G-2loop")
        self.assertEqual((two_loop.nc_power, two_loop.chiral_dimension, two_loop.exponent()),
                         (-1, 6, 4))
        self.assertEqual(two_loop.status, BEYOND_ORDER)

    def test_classify_status(self):
        self.assertEqual(classify_status(2), INCLUDED)
        self.assertEqual(classify_status(3), BEYOND_ORDER)
        self.assertEqual(classify_status(2, declared_approximation=True), DECLARED_APPROXIMATION)
        self.assertEqual(classify_status(3, declared_approximation=True), BEYOND_ORDER)

    def test_thresholds(self):
        thresholds = vector_pseudoscalar_thresholds()
        labels = [label for label, _ in thresholds]
        self.assertEqual(labels, ["rhopi", "omegapi", "rhoeta", "omegaeta", "K*K"])
        values = dict(thresholds)
        self.assertAlmostEqual(values["rhopi"], 0.91483, places=10)
        self.assertAlmostEqual(values["K*K"], 1.385347, places=10)
        self.assertLess(values["rhopi"], VALIDITY_LIMIT_GEV)

    def test_check_loops_runs(self):
        self.assertTrue(check_loops())


# --------------------------------------------------------------------------
# Chiral subtraction scheme and running of L_i^r (--scheme, task 4.1)
# --------------------------------------------------------------------------

dimension = sp.Symbol("d")
scale = sp.Symbol("mu", positive=True)
LOOP_FACTOR = 16 * sp.pi ** 2
# Scales in GeV as exact rationals; default mu = M_rho = 0.77 GeV.
REFERENCE_SCALE_GEV = sp.Rational(77, 100)
SCALES_GEV = [sp.Rational(1, 2), REFERENCE_SCALE_GEV, sp.Integer(1)]
FLAVOUR_NUMBER = 3

# Coefficients Gamma_i of the U(3) scheme, Kaiser:1998ds Eq. (scaling),
# adehep.tex lines 692-696 (same values in Herrera-Siklody:1996tqr,
# un.tex lines 1467-1491, for i = 4, ..., 8, 18).
KAISER_1998DS_U3_GAMMA = {
    "L_4": sp.Rational(1, 8),
    "L_5": sp.Rational(3, 8),
    "L_6": sp.Rational(1, 16),
    "L_7": sp.Integer(0),
    "L_8": sp.Rational(3, 16),
    "L_18": sp.Rational(-1, 4),
    "L_25": sp.Integer(0),
}
GAMMA_SOURCE = "Kaiser:1998ds, adehep.tex lines 692-696"
# Remarks on individual U(3) Gamma_i.
U3_GAMMA_FLAGS = {
    "L_8": ("leading value; B_8 = n_l/16 - v_31/2 + O(N_c^-2)"
            " (Herrera-Siklody:1996tqr, un.tex lines 2044-2050); -v_31/2 is"
            " delta^3, beyond declared order"),
    "L_25": "quoted by Kaiser:1998ds; not reproduced from Herrera-Siklody:1996tqr Appendix D",
}
NO_CHIRAL_LOOP_RENORMALIZATION = ("no chiral-loop renormalization at this order"
                                  " (Kaiser:1998ds, adehep.tex lines 687-689)")
# Further U(3) constants with zero Gamma, with the reason.
U3_ZERO_GAMMA_FLAGS = {
    "F": NO_CHIRAL_LOOP_RENORMALIZATION,
    "tau": NO_CHIRAL_LOOP_RENORMALIZATION,
    "Lambda_1": NO_CHIRAL_LOOP_RENORMALIZATION,
    "Lambda_2": NO_CHIRAL_LOOP_RENORMALIZATION,
    "L_47": ("zero by truncation; fails at next 1/N_c order"
             " (B_47(0) = v_12 (v_31 - 3) = O(N_c^-2), Herrera-Siklody:1996tqr)"),
}
# SU(3) values for comparison: Gasser:1984gg; Gamma_6 and Gamma_8 also quoted
# in Herrera-Siklody:1996tqr, un.tex lines 1507-1508.
SU3_GAMMA = {
    "L_4": sp.Rational(1, 8),
    "L_5": sp.Rational(3, 8),
    "L_6": sp.Rational(11, 144),
    "L_7": sp.Integer(0),
    "L_8": sp.Rational(5, 48),
}
CONSTANT_NAMES = list(KAISER_1998DS_U3_GAMMA)
MAXIMAL_COEFFICIENT = 3

# Dependence on the QCD renormalization scale mu of the singlet axial
# current (Kaiser:1998ds, adehep.tex lines 621-630; B -> Z_M B at line 600).
# mu dZ_A/dmu = gamma_A Z_A, Z_A = 1 + O(1/N_c) (lines 583-594).
# Not part of the mu_chi running.
QCD_SCALE_LAWS = [
    ("tau", "Z_A^2 tau"),
    ("1 + Lambda_1", "Z_A^2 (1 + Lambda_1)"),
    ("1 + Lambda_2", "Z_A (1 + Lambda_2)"),
    ("2 L_5 + 3 L_18", "Z_A (2 L_5 + 3 L_18)"),
    ("2 L_8 - 3 L_25", "Z_A (2 L_8 - 3 L_25)"),
    ("L_4, ..., L_8", "QCD-scale independent"),
    ("F", "QCD-scale independent"),
    ("B", "Z_M B (Z_M: quark-mass factor, line 600)"),
]

# Symbols for the QCD-scale and 1/N_c checks.
colour_number = sp.Symbol("N_c", positive=True)
z_a_coefficient = sp.Symbol("z_A")
singlet_axial_factor = sp.Symbol("Z_A", positive=True)
pole = sp.Symbol("lambda")
tau_symbol, decay_constant = sp.symbols("tau F", positive=True)


def chiral_scale_lambda(mu=scale, spacetime_dimension=dimension):
    """lambda(mu) in the Gasser:1983yg convention, with Gamma'(1) = -EulerGamma.

    Same expression as Kaiser:1998ds (adehep.tex lines 674-675) and
    lambda_infinity of Herrera-Siklody:1996tqr Eq. (3.20).
    """
    return (mu ** (spacetime_dimension - 4) / LOOP_FACTOR
            * (1 / (spacetime_dimension - 4)
               - HALF * (sp.log(4 * sp.pi) + sp.polygamma(0, 1) + 1)))


def lambda_scale_dependence():
    """Expansion of mu^(d-4)/(d-4) about d = 4, up to O(d-4)."""
    epsilon = sp.Symbol("epsilon")
    expansion = sp.series(scale ** epsilon / epsilon, epsilon, 0, 1).removeO()
    return expansion.subs(epsilon, dimension - 4)


def running_constant(reference_value, gamma, mu, reference_mu=REFERENCE_SCALE_GEV):
    """L_i^r(mu) = L_i^r(mu_1) - Gamma_i/(16 pi^2) log(mu/mu_1)."""
    return reference_value - gamma / LOOP_FACTOR * sp.log(mu / reference_mu)


def bare_finite_part(renormalized_value, gamma, mu):
    """mu-dependent finite part of the bare L_i: L_i^r(mu) + Gamma_i log(mu)/(16 pi^2).

    The pole 1/(d-4) and the constant of lambda(mu) are mu-independent and
    are dropped.
    """
    return renormalized_value + gamma * sp.log(mu) / LOOP_FACTOR


def gamma_combination(coefficients, gammas=KAISER_1998DS_U3_GAMMA):
    """Sum_i c_i Gamma_i for a dict {constant name: integer coefficient}."""
    return sum((coefficient * gammas[name] for name, coefficient in coefficients.items()),
               sp.Integer(0))


def null_space_basis(gammas=KAISER_1998DS_U3_GAMMA):
    """Integer basis of {c : sum_i c_i Gamma_i = 0} (sympy nullspace).

    Each vector is made primitive and its first nonzero entry positive.
    """
    row = sp.Matrix([[gammas[name] for name in CONSTANT_NAMES]])
    basis = []
    for vector in row.nullspace():
        denominator = sp.ilcm(*[sp.fraction(entry)[1] for entry in vector])
        integer_vector = vector * denominator
        entries = [int(entry) for entry in integer_vector]
        common = math.gcd(*entries)
        sign = 1 if next(entry for entry in entries if entry) > 0 else -1
        basis.append([sign * entry // common for entry in entries])
    return basis


def small_integer_null_combinations(gammas=KAISER_1998DS_U3_GAMMA,
                                    maximal_coefficient=MAXIMAL_COEFFICIENT):
    """All nonzero integer vectors c with |c_i| <= maximal_coefficient and c.Gamma = 0."""
    common_denominator = sp.ilcm(*[sp.fraction(gammas[name])[1] for name in CONSTANT_NAMES])
    integer_gammas = [int(gammas[name] * common_denominator) for name in CONSTANT_NAMES]
    values = range(-maximal_coefficient, maximal_coefficient + 1)
    return [vector for vector in product(values, repeat=len(CONSTANT_NAMES))
            if any(vector) and sum(c * g for c, g in zip(vector, integer_gammas)) == 0]


def format_combination(vector):
    """Text form of sum_i c_i L_i."""
    terms = []
    for coefficient, name in zip(vector, CONSTANT_NAMES):
        if coefficient == 0:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = "" if abs(coefficient) == 1 else f"{abs(coefficient)} "
        terms.append(f"{sign} {magnitude}{name}")
    if not terms:
        return "0"
    text = " ".join(terms)
    if text.startswith("+ "):
        return text[2:]
    return "-" + text[2:]


def is_primitive_with_positive_lead(vector):
    """True if the gcd of the entries is 1 and the first nonzero entry is positive."""
    nonzero = [entry for entry in vector if entry]
    return math.gcd(*nonzero) == 1 and nonzero[0] > 0


def two_constant_combinations(combinations, gammas=KAISER_1998DS_U3_GAMMA):
    """Primitive vectors (sign fixed) with exactly two entries, both with Gamma != 0."""
    selected = []
    for vector in combinations:
        support = [name for entry, name in zip(vector, CONSTANT_NAMES) if entry]
        if (len(support) == 2 and all(gammas[name] != 0 for name in support)
                and is_primitive_with_positive_lead(vector)):
            selected.append(vector)
    return selected


# --- Chiral renormalization of B (Kaiser:1998ds, adehep.tex lines 679-686) ---

def b_correction_routes():
    """Three countings of the relative correction 4 tau lambda / F^4 to B.

    Returns a dict {route: CountedTerm}. lambda ~ 1/(16 pi^2) counts
    N_c^0 p^0; the N_c suppression of the loop is carried by 1/F^2.
    - "tau/F^4": tau ~ N_c^0 p^0 times 1/F^4 ~ N_c^-2.
    - "m_0^2/(16 pi^2 F^2)": with tau = F^2 m_0^2/(4 N_f), 4 tau/F^4 =
      m_0^2/(N_f F^2), i.e. the eta' tadpole: m_0^2 ~ N_c^-1 p^0 times
      1/(16 pi^2 F^2) ~ N_c^-1.
    - "Goldstone loop rule": the generic factor p^2/(16 pi^2 F^2),
      (k, 2n) = (-1, 2), with the loop momentum set by the internal mass.
    """
    tau_term = CountedTerm("tau", "", 0, 0)
    inverse_f4 = CountedTerm("1/F^4", "", -2, 0)
    eta_prime_mass_squared = CountedTerm("m_0^2", "", -1, 0)
    inverse_f2 = CountedTerm("1/(16 pi^2 F^2)", "", -1, 0)
    return {
        "tau/F^4": combine("4 tau lambda/F^4", "", [tau_term, inverse_f4]),
        "m_0^2/(16 pi^2 F^2)": combine("eta' tadpole", "",
                                       [eta_prime_mass_squared, inverse_f2]),
        "Goldstone loop rule": loop_factors()["Goldstone loop"],
    }


def b_correction_term():
    """Counting of F^2 B^r <chi> x 4 tau lambda/F^4 in the Lagrangian."""
    return combine("B pole term", "F^2/4 <chi_+> x 4 tau lambda/F^4",
                   [goldstone_lagrangian_rows()[1], b_correction_routes()["tau/F^4"]])


def b_running_check():
    """Return mu d(bare B)/dmu at first order in tau for the running B^r.

    B^r(mu) = B^r(mu_1) exp(-tau/(4 pi^2 F^4) log(mu/mu_1)) (Kaiser:1998ds,
    lines 683-686) and bare B = B^r(mu)(1 + 4 tau lambda(mu)/F^4); only the
    log mu of lambda is kept. The result must vanish at O(tau).
    """
    reference_b = sp.Symbol("B^r(0.77)", positive=True)
    running_b = reference_b * sp.exp(-tau_symbol / (4 * sp.pi ** 2 * decay_constant ** 4)
                                     * sp.log(scale / REFERENCE_SCALE_GEV))
    bare_b = running_b * (1 + 4 * tau_symbol / decay_constant ** 4 * sp.log(scale) / LOOP_FACTOR)
    derivative = scale * sp.diff(bare_b, scale)
    return sp.series(derivative, tau_symbol, 0, 2).removeO()


def tau_in_terms_of_eta_prime_mass():
    """Match -tau X^2/2 to -f_pi^2 m_0^2 X^2/(4 N_f)."""
    eta_prime_mass_squared = sp.Symbol("m_0^2", positive=True)
    tau_value = sp.solve(-tau_symbol / 2 + decay_constant ** 2 * eta_prime_mass_squared
                         / (4 * FLAVOUR_NUMBER), tau_symbol)[0]
    return sp.simplify(4 * tau_value / decay_constant ** 4
                       - 2 * eta_prime_mass_squared / (FLAVOUR_NUMBER * decay_constant ** 2))


# --- Gamma_8 beyond leading order (Herrera-Siklody:1996tqr lines 2044-2050) ---

def gamma_8_correction_term():
    """Counting of the dropped -v_31/2 lambda term in the L_8 coefficient.

    v_31 = -Lambda_2/3 ~ N_c^-1 relative to Gamma_8 ~ N_c^0, times the
    one-loop pole term Gamma_8 lambda <chi_+^2 + chi_-^2> ~ N_c^0 p^4
    (the G-LO order).
    """
    loop_pole_term = CountedTerm("Gamma_8 lambda", "N_c^0 p^4", 0, 4)
    relative_v31 = CountedTerm("v_31/Gamma_8", "N_c^-1", -1, 0)
    return combine("-v_31/2 lambda", "dropped part of B_8", [loop_pole_term, relative_v31])


# --- QCD-scale dependence (Kaiser:1998ds lines 621-630, 777-781) ---

def singlet_axial_invariant(l5, l18, lambda1):
    """L_A = (2 L_5^r + 3 L_18^r)/sqrt(1 + Lambda_1) (Kaiser:1998ds lines 777-781)."""
    return (2 * l5 + 3 * l18) / sp.sqrt(1 + lambda1)


def singlet_axial_invariant_chiral_scale_derivative():
    """mu_chi dL_A/dmu_chi with L_5^r, L_18^r running with their Gamma_i."""
    lambda1 = sp.Symbol("Lambda_1", positive=True)
    l5 = running_constant(sp.Symbol("L_5^r(0.77)"), KAISER_1998DS_U3_GAMMA["L_5"], scale)
    l18 = running_constant(sp.Symbol("L_18^r(0.77)"), KAISER_1998DS_U3_GAMMA["L_18"], scale)
    return sp.simplify(scale * sp.diff(singlet_axial_invariant(l5, l18, lambda1), scale))


def singlet_axial_invariant_qcd_scale_change():
    """L_A after the Z_A transformation minus L_A before; must vanish."""
    numerator, one_plus_lambda1 = sp.symbols("N_A s_1", positive=True)
    before = numerator / sp.sqrt(one_plus_lambda1)
    after = (singlet_axial_factor * numerator
             / sp.sqrt(singlet_axial_factor ** 2 * one_plus_lambda1))
    return sp.simplify(after - before)


def su3_l7_qcd_scale_change():
    """Change of F^4 (1 + Lambda_2)^2/(288 tau) under the Z_A laws; must vanish.

    Kaiser:1998ds checks the laws of tau and Lambda_2 against the scale
    independence of L_7^SU(3) = -F^4 (1 + Lambda_2)^2/(288 tau) + L_7
    (lines 631-634).
    """
    one_plus_lambda2 = sp.Symbol("s_2", positive=True)
    before = decay_constant ** 4 * one_plus_lambda2 ** 2 / (288 * tau_symbol)
    after = (decay_constant ** 4 * (singlet_axial_factor * one_plus_lambda2) ** 2
             / (288 * singlet_axial_factor ** 2 * tau_symbol))
    return sp.simplify(after - before)


def renormalized_scale_mismatch(gamma):
    """Mismatch of the Z_A law for a renormalized combination with pole gamma lambda.

    The bare combination C = C^r + gamma lambda transforms as C -> Z_A C.
    Then C^r -> Z_A C - gamma lambda = Z_A C^r + (Z_A - 1) gamma lambda.
    Returns (Z_A - 1) gamma lambda with Z_A = 1 + z_A / N_c.
    """
    z_a = 1 + z_a_coefficient / colour_number
    renormalized = sp.Symbol("C^r")
    transformed = z_a * (renormalized + gamma * pole) - gamma * pole
    return sp.expand(transformed - z_a * renormalized)


def b_renormalized_scale_mismatch():
    """Relative mismatch of B^r -> Z_M B^r, at first order in tau.

    B = B^r (1 + 4 tau lambda/F^4) -> Z_M B and tau -> Z_A^2 tau give
    B^r -> Z_M B^r (1 + 4 tau lambda/F^4)/(1 + 4 Z_A^2 tau lambda/F^4).
    Returns the O(tau) part of the ratio minus one, with Z_A = 1 + z_A/N_c.
    """
    z_a = 1 + z_a_coefficient / colour_number
    correction = 4 * tau_symbol * pole / decay_constant ** 4
    ratio = (1 + correction) / (1 + z_a ** 2 * correction)
    return sp.expand(sp.series(ratio - 1, tau_symbol, 0, 2).removeO())


def scale_mismatch_term(pole_term):
    """Counting of (Z_A - 1) x pole term, with Z_A - 1 ~ N_c^-1."""
    z_a_minus_one = CountedTerm("Z_A - 1", "N_c^-1", -1, 0)
    return combine("(Z_A - 1) x pole", "", [pole_term, z_a_minus_one])


def check_scheme():
    """Print scheme, Gamma table, running and scale-independent combinations; assert."""
    print("Chiral subtraction scheme (Gasser:1983yg convention, as used by"
          " Kaiser:1998ds and Herrera-Siklody:1996tqr Eq. (3.20)):")
    print("  lambda(mu) = mu^(d-4)/(16 pi^2) [1/(d-4) - (log 4 pi + Gamma'(1) + 1)/2]")
    print("  L_i = L_i^r(mu) + Gamma_i lambda(mu)")
    print("  mu here is the chiral subtraction scale; Kaiser:1998ds writes mu_chi"
          " to distinguish it from the QCD scale (adehep.tex lines 676-678)")
    print(f"  sympy: lambda(mu) = {sp.sstr(chiral_scale_lambda())}")
    expansion = lambda_scale_dependence()
    assert sp.simplify(expansion - (1 / (dimension - 4) + sp.log(scale))) == 0
    print(f"  mu^(d-4)/(d-4) = {sp.sstr(expansion)} + O(d-4)")

    print(f"\nGamma_i: U(3), {GAMMA_SOURCE}, vs SU(3), Gasser:1984gg"
          " (Gamma_6, Gamma_8 also Herrera-Siklody:1996tqr lines 1507-1508)")
    print(f"  {'constant':<10}{'U(3)':>8}{'SU(3)':>9}")
    for name in CONSTANT_NAMES:
        su3_text = str(SU3_GAMMA[name]) if name in SU3_GAMMA else "-"
        flag = f"   [{U3_GAMMA_FLAGS[name]}]" if name in U3_GAMMA_FLAGS else ""
        print(f"  {name:<10}{str(KAISER_1998DS_U3_GAMMA[name]):>8}{su3_text:>9}{flag}")
    for name, flag in U3_ZERO_GAMMA_FLAGS.items():
        print(f"  {name:<10}{'0':>8}{'-':>9}   [{flag}]")
    differing = [name for name in SU3_GAMMA if SU3_GAMMA[name] != KAISER_1998DS_U3_GAMMA[name]]
    assert differing == ["L_6", "L_8"], differing
    print(f"  U(3) and SU(3) differ for: {', '.join(differing)}")
    gamma_8_correction = gamma_8_correction_term()
    assert delta_exponent(gamma_8_correction) == 3 > DECLARED_EXPONENT
    print("  Gamma_8: dropped -v_31/2 (v_31 = -Lambda_2/3 ~ 1/N_c) is delta^1 relative"
          " to Gamma_8; times the one-loop pole term (delta^2) it is"
          f" {sp.sstr(delta_order(gamma_8_correction.nc_power, gamma_8_correction.chiral_dimension))},"
          " beyond the declared order")

    print("\nChiral renormalization of B (Kaiser:1998ds, adehep.tex lines 679-686):")
    print("  B = B^r(mu) (1 + 4 tau lambda(mu)/F^4),"
          "  mu dB^r/dmu = -tau B^r/(4 pi^2 F^4)")
    routes = b_correction_routes()
    for route, term in routes.items():
        print(f"  relative correction via {route:<22} (k, 2n) = ({term.nc_power}, "
              f"{term.chiral_dimension}): "
              f"{sp.sstr(delta_order(term.nc_power, term.chiral_dimension))}")
    route_exponents = {delta_exponent(term) for term in routes.values()}
    assert route_exponents == {2}, route_exponents
    assert tau_in_terms_of_eta_prime_mass() == 0
    correction_term = b_correction_term()
    g_lo = next(item for item in loop_classes() if item.name == "G-LO")
    assert delta_exponent(correction_term) == g_lo.exponent() == DECLARED_EXPONENT
    assert sp.simplify(b_running_check()) == 0
    print(f"  pole term F^2 B^r <chi_+> x 4 tau lambda/F^4: (k, 2n) ="
          f" ({correction_term.nc_power}, {correction_term.chiral_dimension}),"
          f" {sp.sstr(delta_order(correction_term.nc_power, correction_term.chiral_dimension))},"
          " same order as G-LO: inside the declared order, so B^r depends on mu_chi")
    print("  checked: 4 tau/f_pi^4 = 2 m_0^2/(N_f f_pi^2) for tau = f_pi^2 m_0^2/(2 N_f);"
          " bare B is mu-independent at O(tau)")

    print("\nRunning: L_i^r(mu_2) = L_i^r(mu_1) - Gamma_i/(16 pi^2) log(mu_2/mu_1),"
          " i.e. mu dL_i^r/dmu = -Gamma_i/(16 pi^2)")
    for name in CONSTANT_NAMES:
        gamma = KAISER_1998DS_U3_GAMMA[name]
        reference_value = sp.Symbol(f"{name}^r(0.77)")
        running = running_constant(reference_value, gamma, scale)
        assert sp.simplify(scale * sp.diff(running, scale) + gamma / LOOP_FACTOR) == 0, name
        assert sp.simplify(scale * sp.diff(bare_finite_part(running, gamma, scale), scale)) == 0
        reference_bare = bare_finite_part(reference_value, gamma, REFERENCE_SCALE_GEV)
        for mu_value in SCALES_GEV:
            difference = (bare_finite_part(running.subs(scale, mu_value), gamma, mu_value)
                          - reference_bare)
            assert sp.simplify(sp.expand_log(difference, force=True)) == 0, (name, mu_value)
            assert abs(float(difference.subs(reference_value, 0))) < 1e-15, (name, mu_value)
    print("  checked: d/dlog mu [L_i^r(mu) + Gamma_i log(mu)/(16 pi^2)] = 0 for every i;")
    print("  bare combination equal at mu = 0.5, 0.77, 1.0 GeV for symbolic L_i^r(0.77)")

    print("\nShifts L_i^r(mu) - L_i^r(0.77 GeV) in units of 10^-3:")
    print(f"  {'constant':<10}" + "".join(f"{f'mu={float(mu_value):g}':>12}"
                                         for mu_value in SCALES_GEV))
    for name in CONSTANT_NAMES:
        shifts = [1000 * float(running_constant(0, KAISER_1998DS_U3_GAMMA[name], mu_value))
                  for mu_value in SCALES_GEV]
        print(f"  {name:<10}" + "".join(f"{shift:>12.4f}" for shift in shifts))

    print("\nCombinations sum_i c_i L_i with sum_i c_i Gamma_i = 0 (independent of mu_chi),")
    print(f"  among {', '.join(CONSTANT_NAMES)}")
    basis = null_space_basis()
    assert len(basis) == len(CONSTANT_NAMES) - 1
    print("  null-space basis (sympy nullspace, rescaled to integers):")
    for vector in basis:
        assert gamma_combination(dict(zip(CONSTANT_NAMES, vector))) == 0
        print(f"    {format_combination(vector)}")
    combinations = small_integer_null_combinations()
    print(f"  integer combinations with |c_i| <= {MAXIMAL_COEFFICIENT}: {len(combinations)}"
          " nonzero vectors")
    zero_gamma_names = [name for name in CONSTANT_NAMES if KAISER_1998DS_U3_GAMMA[name] == 0]
    print(f"  constants with Gamma = 0 (any multiple can be added): {', '.join(zero_gamma_names)}")
    short = two_constant_combinations(combinations)
    print("  combinations of two constants with Gamma != 0 (primitive, overall sign fixed):")
    for vector in short:
        print(f"    {format_combination(vector)}")
    checks = [({"L_5": 2, "L_18": 3}, 0), ({"L_7": 1}, 0), ({"L_25": 1}, 0),
              ({"L_8": 2, "L_25": -3}, sp.Rational(3, 8))]
    for coefficients, expected in checks:
        value = gamma_combination(coefficients)
        assert value == expected, (coefficients, value)
        text = format_combination([coefficients.get(name, 0) for name in CONSTANT_NAMES])
        verdict = "mu_chi-independent" if value == 0 else "NOT mu_chi-independent"
        print(f"  {text:<16} Gamma = {str(value):<5} {verdict}")

    print("\nQCD-scale dependence (singlet axial current), Kaiser:1998ds,"
          " adehep.tex lines 621-630; not part of the mu_chi running:")
    print("  mu dZ_A/dmu = gamma_A Z_A, Z_A = 1 + O(1/N_c)")
    for quantity, transformation in QCD_SCALE_LAWS:
        print(f"  {quantity:<16}-> {transformation}")
    assert dict(QCD_SCALE_LAWS)["tau"] == "Z_A^2 tau"
    assert singlet_axial_invariant_chiral_scale_derivative() == 0
    assert singlet_axial_invariant_qcd_scale_change() == 0
    assert su3_l7_qcd_scale_change() == 0
    print("  invariant L_A = (2 L_5^r + 3 L_18^r)/sqrt(1 + Lambda_1) (lines 777-781):")
    print("    checked: mu_chi dL_A/dmu_chi = 0 (2 Gamma_5 + 3 Gamma_18 = 0, Lambda_1"
          " has no pole) and L_A unchanged under the Z_A laws")
    print("  checked: F^4 (1 + Lambda_2)^2/tau (L_7 of SU(3), lines 631-634)"
          " unchanged under the Z_A laws")

    print("\nZ_A law for renormalized combinations with a pole:")
    combination_gamma = gamma_combination({"L_8": 2, "L_25": -3})
    mismatch = renormalized_scale_mismatch(combination_gamma)
    assert sp.simplify(mismatch - z_a_coefficient / colour_number * combination_gamma * pole) == 0
    assert sp.limit(mismatch * colour_number, colour_number, sp.oo) != 0
    assert sp.limit(mismatch, colour_number, sp.oo) == 0
    pole_term = CountedTerm("Gamma lambda", "one-loop pole term, N_c^0 p^4", 0, 4)
    mismatch_term = scale_mismatch_term(pole_term)
    assert delta_exponent(pole_term) == DECLARED_EXPONENT
    assert delta_exponent(mismatch_term) == 3 > DECLARED_EXPONENT
    print(f"  2 L_8^r - 3 L_25^r -> Z_A (2 L_8^r - 3 L_25^r) + {sp.sstr(mismatch)}"
          " with Z_A = 1 + z_A/N_c")
    print("    mismatch (Z_A - 1)(3/8) lambda = O(1/N_c) x loop pole:"
          f" {sp.sstr(delta_order(mismatch_term.nc_power, mismatch_term.chiral_dimension))}"
          " relative to LO, beyond the declared order")
    b_mismatch = b_renormalized_scale_mismatch()
    assert sp.limit(b_mismatch, colour_number, sp.oo) == 0
    assert sp.limit(b_mismatch * colour_number, colour_number, sp.oo) == sp.simplify(
        -8 * z_a_coefficient * tau_symbol * pole / decay_constant ** 4)
    b_mismatch_term = scale_mismatch_term(correction_term)
    assert delta_exponent(b_mismatch_term) == 3 > DECLARED_EXPONENT
    print(f"  B^r -> Z_M B^r [1 + ({sp.sstr(sp.factor(b_mismatch))}) + O(tau^2)]")
    print(f"    same structure through tau -> Z_A^2 tau: O(1/N_c) x loop pole,"
          f" {sp.sstr(delta_order(b_mismatch_term.nc_power, b_mismatch_term.chiral_dimension))},"
          " beyond the declared order")

    print("\nNotes:")
    print("- mu_chi (chiral subtraction scale) and the QCD renormalization scale"
          " are independent. A quantity can be mu_chi-independent but QCD-scale"
          " dependent, e.g. 2 L_5 + 3 L_18 (Gamma = 0, transforms as Z_A); having"
          " no pole (Gamma_25 = 0, Lambda_1, Lambda_2, tau) does not remove the"
          " Z_A dependence.")
    print("- The Z_A laws of Kaiser:1998ds hold for the bare couplings. A combination"
          " with a nonzero pole (2 L_8 - 3 L_25, and B through tau) is then not"
          " multiplicative after subtraction: the renormalized one picks up"
          " (Z_A - 1) x pole. This is O(1/N_c) times a one-loop pole, i.e. delta^3,"
          " the same order as the dropped -v_31/2 in Gamma_8 and the unknown"
          " O(1/N_c) parts of Gamma_25; consistency at delta^3 would need the"
          " Gamma_i at next order in 1/N_c. At the declared order delta^2 the"
          " laws hold for the renormalized couplings.")
    print("- Gamma_47 = 0 holds only under the large-N_c truncation of"
          " Herrera-Siklody:1996tqr and fails at the next order in 1/N_c.")
    print("\n--scheme: all assertions passed")
    return True


class SchemeTests(unittest.TestCase):
    """Unit tests of the subtraction scheme, running and Gamma combinations."""

    def test_lambda_pole_and_scale_dependence(self):
        expansion = lambda_scale_dependence()
        self.assertEqual(sp.simplify(expansion - 1 / (dimension - 4) - sp.log(scale)), 0)
        residue = sp.limit(chiral_scale_lambda() * (dimension - 4), dimension, 4)
        self.assertEqual(sp.simplify(residue - 1 / LOOP_FACTOR), 0)

    def test_lambda_finite_part_constant(self):
        epsilon = sp.Symbol("epsilon")
        finite = sp.series(chiral_scale_lambda(1, 4 + epsilon), epsilon, 0, 1).removeO()
        finite = finite.coeff(epsilon, 0)
        expected = -(sp.log(4 * sp.pi) - sp.EulerGamma + 1) / (2 * LOOP_FACTOR)
        self.assertEqual(sp.simplify(finite - expected), 0)

    def test_u3_gamma_values(self):
        self.assertEqual(KAISER_1998DS_U3_GAMMA, {
            "L_4": sp.Rational(1, 8), "L_5": sp.Rational(3, 8), "L_6": sp.Rational(1, 16),
            "L_7": 0, "L_8": sp.Rational(3, 16), "L_18": sp.Rational(-1, 4), "L_25": 0})
        self.assertIn("Kaiser:1998ds", GAMMA_SOURCE)

    def test_u3_su3_difference(self):
        differing = {name for name in SU3_GAMMA
                     if SU3_GAMMA[name] != KAISER_1998DS_U3_GAMMA[name]}
        self.assertEqual(differing, {"L_6", "L_8"})

    def test_herrera_siklody_su3_relations(self):
        # Herrera-Siklody:1996tqr lines 1520-1525 and general-n_l forms.
        self.assertEqual(3 * KAISER_1998DS_U3_GAMMA["L_6"] + KAISER_1998DS_U3_GAMMA["L_8"],
                         sp.Rational(9, 8) * (3 * SU3_GAMMA["L_6"] + SU3_GAMMA["L_8"]))
        flavours = sp.Integer(3)
        self.assertEqual(SU3_GAMMA["L_6"], (2 + flavours ** 2) / (16 * flavours ** 2))
        self.assertEqual(SU3_GAMMA["L_8"], (flavours ** 2 - 4) / (16 * flavours))

    def test_flags(self):
        self.assertIn("v_31", U3_GAMMA_FLAGS["L_8"])
        self.assertIn("not reproduced", U3_GAMMA_FLAGS["L_25"])
        for name in ["F", "tau", "Lambda_1", "Lambda_2"]:
            self.assertEqual(U3_ZERO_GAMMA_FLAGS[name], NO_CHIRAL_LOOP_RENORMALIZATION)

    def test_running_and_bare_independence(self):
        reference_value = sp.Symbol("L_ref")
        gamma = sp.Rational(3, 8)
        running = running_constant(reference_value, gamma, scale)
        self.assertEqual(running.subs(scale, REFERENCE_SCALE_GEV), reference_value)
        self.assertEqual(sp.simplify(scale * sp.diff(running, scale)), -gamma / LOOP_FACTOR)
        bare = bare_finite_part(running, gamma, scale)
        self.assertEqual(sp.simplify(sp.diff(bare, scale)), 0)

    def test_running_shift_value(self):
        shift = float(running_constant(0, sp.Rational(3, 8), sp.Rational(1, 2)))
        expected = -3 / 8 / (16 * math.pi ** 2) * math.log(0.5 / 0.77)
        self.assertAlmostEqual(shift, expected, places=15)
        self.assertGreater(shift, 0)

    def test_running_composes(self):
        gamma = sp.Rational(1, 8)
        reference_value = sp.Symbol("L_ref")
        at_one = running_constant(reference_value, gamma, 1)
        back = running_constant(at_one, gamma, sp.Rational(1, 2), reference_mu=1)
        direct = running_constant(reference_value, gamma, sp.Rational(1, 2))
        self.assertEqual(sp.simplify(sp.expand_log(back - direct, force=True)), 0)

    def test_named_combinations(self):
        self.assertEqual(gamma_combination({"L_5": 2, "L_18": 3}), 0)
        self.assertEqual(gamma_combination({"L_7": 1}), 0)
        self.assertEqual(gamma_combination({"L_25": 1}), 0)
        self.assertEqual(gamma_combination({"L_8": 2, "L_25": -3}), sp.Rational(3, 8))

    def test_null_space_basis(self):
        basis = null_space_basis()
        self.assertEqual(len(basis), 6)
        self.assertEqual(sp.Matrix(basis).rank(), 6)
        for vector in basis:
            self.assertEqual(gamma_combination(dict(zip(CONSTANT_NAMES, vector))), 0)

    def test_small_integer_combinations(self):
        combinations = small_integer_null_combinations()
        index = {name: position for position, name in enumerate(CONSTANT_NAMES)}
        target = [0] * len(CONSTANT_NAMES)
        target[index["L_5"]], target[index["L_18"]] = 2, 3
        self.assertIn(tuple(target), combinations)
        excluded = [0] * len(CONSTANT_NAMES)
        excluded[index["L_8"]], excluded[index["L_25"]] = 2, -3
        self.assertNotIn(tuple(excluded), combinations)
        # Orthogonality to the Gamma row is membership of the null space.
        for vector in combinations:
            self.assertEqual(gamma_combination(dict(zip(CONSTANT_NAMES, vector))), 0)
        self.assertEqual(len(combinations), len(set(combinations)))
        self.assertTrue(all(tuple(-entry for entry in vector) in combinations
                            for vector in combinations[:50]))

    def test_small_enumeration_toy(self):
        toy = {name: sp.Integer(0) for name in CONSTANT_NAMES}
        toy["L_4"], toy["L_5"] = sp.Integer(1), sp.Integer(-1)
        # Gamma_4 = 1, Gamma_5 = -1, rest 0: c_4 = c_5 and five free entries.
        combinations = small_integer_null_combinations(toy, maximal_coefficient=1)
        self.assertEqual(len(combinations), 3 * 3 ** 5 - 1)

    def test_two_constant_combinations(self):
        texts = {format_combination(vector)
                 for vector in two_constant_combinations(small_integer_null_combinations())}
        self.assertEqual(texts, {"3 L_4 - L_5", "L_4 - 2 L_6", "3 L_4 - 2 L_8", "2 L_4 + L_18",
                                 "L_5 - 2 L_8", "2 L_5 + 3 L_18", "3 L_6 - L_8"})

    def test_format_combination(self):
        self.assertEqual(format_combination([0, 2, 0, 0, 0, 3, 0]), "2 L_5 + 3 L_18")
        self.assertEqual(format_combination([0, 0, 0, 0, 2, 0, -3]), "2 L_8 - 3 L_25")
        self.assertEqual(format_combination([-1, 0, 0, 0, 0, 0, 0]), "-L_4")

    def test_b_renormalization_order(self):
        for route, term in b_correction_routes().items():
            self.assertEqual(delta_exponent(term), 2, route)
        self.assertEqual((b_correction_routes()["tau/F^4"].nc_power,
                          b_correction_routes()["tau/F^4"].chiral_dimension), (-2, 0))
        term = b_correction_term()
        self.assertEqual((term.nc_power, term.chiral_dimension), (-1, 2))
        self.assertEqual(delta_exponent(term), DECLARED_EXPONENT)
        self.assertEqual(tau_in_terms_of_eta_prime_mass(), 0)

    def test_b_running(self):
        self.assertEqual(sp.simplify(b_running_check()), 0)

    def test_gamma_8_correction_beyond_order(self):
        term = gamma_8_correction_term()
        self.assertEqual((term.nc_power, term.chiral_dimension), (-1, 4))
        self.assertEqual(delta_exponent(term), 3)
        self.assertGreater(delta_exponent(term), DECLARED_EXPONENT)

    def test_qcd_scale_laws(self):
        laws = dict(QCD_SCALE_LAWS)
        self.assertEqual(laws["tau"], "Z_A^2 tau")
        self.assertEqual(laws["1 + Lambda_1"], "Z_A^2 (1 + Lambda_1)")
        self.assertEqual(laws["1 + Lambda_2"], "Z_A (1 + Lambda_2)")
        self.assertEqual(laws["2 L_5 + 3 L_18"], "Z_A (2 L_5 + 3 L_18)")
        self.assertEqual(laws["2 L_8 - 3 L_25"], "Z_A (2 L_8 - 3 L_25)")
        self.assertEqual(laws["F"], "QCD-scale independent")
        self.assertNotIn("none", laws.values())

    def test_singlet_axial_invariant(self):
        self.assertEqual(singlet_axial_invariant_chiral_scale_derivative(), 0)
        self.assertEqual(singlet_axial_invariant_qcd_scale_change(), 0)
        self.assertEqual(su3_l7_qcd_scale_change(), 0)

    def test_renormalized_scale_mismatch(self):
        mismatch = renormalized_scale_mismatch(sp.Rational(3, 8))
        self.assertEqual(sp.simplify(
            mismatch - sp.Rational(3, 8) * z_a_coefficient / colour_number * pole), 0)
        self.assertEqual(renormalized_scale_mismatch(0), 0)
        b_mismatch = b_renormalized_scale_mismatch()
        self.assertEqual(sp.limit(b_mismatch, colour_number, sp.oo), 0)
        self.assertEqual(sp.simplify(sp.limit(b_mismatch * colour_number, colour_number, sp.oo)
                                     + 8 * z_a_coefficient * tau_symbol * pole
                                     / decay_constant ** 4), 0)
        pole_term = CountedTerm("pole", "", 0, 4)
        self.assertEqual(delta_exponent(scale_mismatch_term(pole_term)), 3)
        self.assertEqual(delta_exponent(scale_mismatch_term(b_correction_term())), 3)

    def test_check_scheme_runs(self):
        self.assertTrue(check_scheme())


# --------------------------------------------------------------------------
# Truncation uncertainty per unit |Delta| (--truncation, task 6.1)
# --------------------------------------------------------------------------

COLOUR_NUMBER_VALUE = 3
# f_pi = F convention of the project: 130.2 MeV / sqrt(2) (ParticleDataGroup:2024cfk).
PION_DECAY_CONSTANT_GEV = sp.Rational("0.0921")
KAON_MASS_GEV = sp.Rational(str(MASSES_GEV["K"]))
CHIRAL_SCALE_GEV = 4 * sp.pi * PION_DECAY_CONSTANT_GEV
# Assumption: relative size of the excluded resonance loops (see docstring).
RESONANCE_LOOP_FRACTION = sp.Rational(3, 10)
LOWEST_RESONANCE_PION_THRESHOLD_GEV = (sp.Rational(str(MASSES_GEV["rho"]))
                                       + sp.Rational(str(MASSES_GEV["pi"])))
DECLARED_RANGE_LIMIT_GEV = sp.Rational(str(VALIDITY_LIMIT_GEV))
DEFAULT_SQRT_S_GRID_GEV = ["0.3", "0.5", "0.7", "0.9", "1.0", "1.2"]
INVERSE_COLOUR_TERM = "1/N_c"
ENERGY_TERM = "s/Lambda_chi^2"
KAON_MASS_TERM = "M_K^2/Lambda_chi^2"
NUMERICAL_DIGITS = 50
LINE_SHAPE_NOTE = (
    "Near a resonance pole, where |M_R^2 - s| is of order M_R Gamma_R, the"
    " narrow-width tree propagator fails and the counting is not applicable;"
    " the line-shape prescription of"
    " memory/rcht/foundations/vector-width-lineshape-prescriptions.md replaces"
    " it as a resummation outside the strict counting. The windows are set by"
    " the downstream calculation with the PDG widths.")
ASSUMPTION_NOTE = (
    "eps_R is an assumption (Pich:2008jm: 17% and 28% shifts of L_10^r and"
    " C_87^r; V', A' about 30-40% of the correction). Delta must be supplied"
    " by the downstream calculation.")


def numerical_value(expression):
    """High-precision float of an exact sympy expression (50 digits)."""
    return sp.N(expression, NUMERICAL_DIGITS)


def expansion_parameter_terms(sqrt_s):
    """The three candidate sizes of delta at sqrt(s) (GeV), exact sympy values."""
    return {
        INVERSE_COLOUR_TERM: sp.Rational(1, COLOUR_NUMBER_VALUE),
        ENERGY_TERM: sqrt_s ** 2 / CHIRAL_SCALE_GEV ** 2,
        KAON_MASS_TERM: KAON_MASS_GEV ** 2 / CHIRAL_SCALE_GEV ** 2,
    }


def estimated_expansion_parameter(sqrt_s):
    """delta_est(s) = max(1/N_c, s/Lambda_chi^2, M_K^2/Lambda_chi^2).

    Returns (value, name of the largest term).
    """
    terms = expansion_parameter_terms(sqrt_s)
    name = max(terms, key=lambda label: numerical_value(terms[label]))
    return terms[name], name


def crossover_sqrt_s():
    """sqrt(s) where s/Lambda_chi^2 = 1/N_c, i.e. Lambda_chi/sqrt(N_c), in GeV."""
    return CHIRAL_SCALE_GEV / sp.sqrt(COLOUR_NUMBER_VALUE)


def resonance_loop_indicator(sqrt_s, inputs_from_data):
    """theta = 1 if the excluded resonance loops enter the uncertainty, else 0.

    With inputs fixed from data only the non-analytic part, which starts at
    the rho pi threshold, is missing; otherwise the whole R-int class is.
    """
    if not inputs_from_data:
        return 1
    return 1 if sqrt_s >= LOWEST_RESONANCE_PION_THRESHOLD_GEV else 0


def truncation_uncertainty(sqrt_s):
    """u_trunc / |Delta| = delta_est(s)."""
    return estimated_expansion_parameter(sqrt_s)[0]


def resonance_loop_uncertainty(sqrt_s, inputs_from_data):
    """u_Rloop / |Delta| = eps_R theta."""
    return RESONANCE_LOOP_FRACTION * resonance_loop_indicator(sqrt_s, inputs_from_data)


def total_uncertainty(sqrt_s, inputs_from_data):
    """u_tot / |Delta| = sqrt(delta_est^2 + eps_R^2 theta) (quadrature)."""
    return sp.sqrt(truncation_uncertainty(sqrt_s) ** 2
                   + resonance_loop_uncertainty(sqrt_s, inputs_from_data) ** 2)


def in_declared_range(sqrt_s):
    """True if sqrt(s) <= 1.0 GeV, the declared validity range."""
    return sqrt_s <= DECLARED_RANGE_LIMIT_GEV


def expansion_parameter_below_one(sqrt_s):
    """True if delta_est(s) < 1."""
    return numerical_value(truncation_uncertainty(sqrt_s)) < 1


def yes_no(flag):
    """'yes' or 'no'."""
    return "yes" if flag else "no"


def format_number(expression, digits=6):
    """Fixed-point text of an exact expression."""
    return f"{float(numerical_value(expression)):.{digits}f}"


def check_truncation(sqrt_s_values=None):
    """Print delta_est and the uncertainties per unit |Delta|; assert the rules.

    sqrt_s_values: list of exact sympy Rationals in GeV (default grid if None
    or empty).
    """
    if not sqrt_s_values:
        sqrt_s_values = [sp.Rational(value) for value in DEFAULT_SQRT_S_GRID_GEV]
    crossover = crossover_sqrt_s()
    kaon_term = expansion_parameter_terms(0)[KAON_MASS_TERM]

    print("Truncation uncertainty at declared order delta^2 (omitted terms"
          " O(delta^3), one more factor of delta)")
    print("delta_est(s) = max(1/N_c, s/Lambda_chi^2, M_K^2/Lambda_chi^2);"
          " all uncertainties per unit |Delta|")
    print(f"N_c = {COLOUR_NUMBER_VALUE}, f_pi = {format_number(PION_DECAY_CONSTANT_GEV, 4)} GeV,"
          f" M_K = {format_number(KAON_MASS_GEV)} GeV (ParticleDataGroup:2024cfk)")
    print(f"Lambda_chi = 4 pi f_pi = {format_number(CHIRAL_SCALE_GEV)} GeV")
    print(f"crossover s/Lambda_chi^2 = 1/N_c at sqrt(s) = Lambda_chi/sqrt(3)"
          f" = {format_number(crossover)} GeV")
    print(f"M_K^2/Lambda_chi^2 = {format_number(kaon_term)} < 1/N_c ="
          f" {format_number(sp.Rational(1, COLOUR_NUMBER_VALUE))}: the M_K term never dominates")
    print(f"eps_R = {format_number(RESONANCE_LOOP_FRACTION, 2)} (assumption);"
          f" lowest R pi threshold m_rho + m_pi ="
          f" {format_number(LOWEST_RESONANCE_PION_THRESHOLD_GEV, 5)} GeV")
    print(f"declared validity range: sqrt(s) <= {format_number(DECLARED_RANGE_LIMIT_GEV, 1)} GeV")

    for sqrt_s in sqrt_s_values:
        terms = expansion_parameter_terms(sqrt_s)
        value, dominant = estimated_expansion_parameter(sqrt_s)
        print(f"\nsqrt(s) = {format_number(sqrt_s, 3)} GeV")
        for label, term in terms.items():
            print(f"  {label:<20}{format_number(term)}")
        print(f"  delta_est           {format_number(value)}  (maximum: {dominant})")
        print(f"  u_trunc/|Delta|     {format_number(truncation_uncertainty(sqrt_s))}")
        for inputs_from_data in (True, False):
            theta = resonance_loop_indicator(sqrt_s, inputs_from_data)
            print(f"  inputs_from_data={str(inputs_from_data):<5}  theta = {theta},"
                  f" u_Rloop/|Delta| = "
                  f"{format_number(resonance_loop_uncertainty(sqrt_s, inputs_from_data))},"
                  f" u_tot/|Delta| = {format_number(total_uncertainty(sqrt_s, inputs_from_data))}")
        print(f"  in declared range: {yes_no(in_declared_range(sqrt_s))};"
              f" delta_est < 1: {yes_no(expansion_parameter_below_one(sqrt_s))}")

    # Assertions on the printed values and on a fine grid.
    assert numerical_value(kaon_term) < sp.Rational(1, COLOUR_NUMBER_VALUE)
    fine_grid = [sp.Rational(step, 1000) for step in range(50, 1501, 5)]
    check_points = sorted(set(fine_grid) | set(sqrt_s_values)
                          | {LOWEST_RESONANCE_PION_THRESHOLD_GEV})
    previous = None
    for sqrt_s in check_points:
        value, dominant = estimated_expansion_parameter(sqrt_s)
        assert dominant != KAON_MASS_TERM, sqrt_s
        if sqrt_s < crossover:
            assert dominant == INVERSE_COLOUR_TERM, sqrt_s
        elif sqrt_s <= CHIRAL_SCALE_GEV:
            assert dominant == ENERGY_TERM, sqrt_s
        if previous is not None:
            assert numerical_value(value - previous) >= 0, sqrt_s
        previous = value
        if sqrt_s > crossover:
            assert expansion_parameter_below_one(sqrt_s) == (sqrt_s < CHIRAL_SCALE_GEV), sqrt_s
        assert resonance_loop_indicator(sqrt_s, False) == 1
        if sqrt_s < LOWEST_RESONANCE_PION_THRESHOLD_GEV:
            assert resonance_loop_uncertainty(sqrt_s, True) == 0, sqrt_s
        else:
            assert resonance_loop_uncertainty(sqrt_s, True) > 0, sqrt_s
        for inputs_from_data in (True, False):
            total = numerical_value(total_uncertainty(sqrt_s, inputs_from_data))
            assert total >= numerical_value(truncation_uncertainty(sqrt_s))
            assert total >= numerical_value(resonance_loop_uncertainty(sqrt_s, inputs_from_data))
    at_one_gev = numerical_value(truncation_uncertainty(sp.Integer(1)))
    assert 0.7 < at_one_gev < 0.8, at_one_gev
    assert LOWEST_RESONANCE_PION_THRESHOLD_GEV == sp.Rational("0.91483")
    assert vector_pseudoscalar_thresholds()[0][0] == "rhopi"

    print("\nNotes:")
    print(f"- {LINE_SHAPE_NOTE}")
    print(f"- {ASSUMPTION_NOTE}")
    print("- u_Rloop applies for inputs_from_data=True only at sqrt(s) >= m_rho + m_pi"
          " (non-analytic part of R-int); for inputs_from_data=False at all sqrt(s).")
    print("\n--truncation: all assertions passed")
    return True


class TruncationTests(unittest.TestCase):
    """Unit tests of delta_est and of the uncertainties per unit |Delta|."""

    def test_constants(self):
        self.assertAlmostEqual(float(CHIRAL_SCALE_GEV), 4 * math.pi * 0.0921, places=14)
        self.assertAlmostEqual(float(crossover_sqrt_s()), 4 * math.pi * 0.0921 / math.sqrt(3),
                               places=14)
        self.assertAlmostEqual(float(crossover_sqrt_s()), 0.668204, places=6)
        self.assertEqual(LOWEST_RESONANCE_PION_THRESHOLD_GEV, sp.Rational("0.91483"))
        self.assertAlmostEqual(float(PION_DECAY_CONSTANT_GEV) * math.sqrt(2), 0.1302, places=3)

    def test_crossover_is_exact(self):
        terms = expansion_parameter_terms(crossover_sqrt_s())
        self.assertEqual(sp.simplify(terms[ENERGY_TERM] - terms[INVERSE_COLOUR_TERM]), 0)

    def test_kaon_term_never_dominates(self):
        kaon_term = float(expansion_parameter_terms(0)[KAON_MASS_TERM])
        self.assertAlmostEqual(kaon_term, 0.493677 ** 2 / (4 * math.pi * 0.0921) ** 2, places=12)
        self.assertLess(kaon_term, 1 / 3)
        self.assertGreater(kaon_term, 0.18)

    def test_dominant_term(self):
        self.assertEqual(estimated_expansion_parameter(sp.Rational(1, 2))[1], INVERSE_COLOUR_TERM)
        self.assertEqual(estimated_expansion_parameter(sp.Rational(66, 100))[1],
                         INVERSE_COLOUR_TERM)
        self.assertEqual(estimated_expansion_parameter(sp.Rational(67, 100))[1], ENERGY_TERM)
        self.assertEqual(estimated_expansion_parameter(sp.Integer(1))[1], ENERGY_TERM)

    def test_values(self):
        self.assertEqual(truncation_uncertainty(sp.Rational(3, 10)), sp.Rational(1, 3))
        at_one_gev = float(truncation_uncertainty(sp.Integer(1)))
        self.assertAlmostEqual(at_one_gev, 1 / (4 * math.pi * 0.0921) ** 2, places=12)
        self.assertTrue(0.7 < at_one_gev < 0.8)

    def test_monotonic(self):
        grid = [sp.Rational(step, 100) for step in range(5, 151)]
        values = [float(truncation_uncertainty(sqrt_s)) for sqrt_s in grid]
        self.assertEqual(values, sorted(values))

    def test_below_one_iff_below_chiral_scale(self):
        self.assertTrue(expansion_parameter_below_one(sp.Rational(115, 100)))
        self.assertFalse(expansion_parameter_below_one(sp.Rational(116, 100)))
        self.assertFalse(expansion_parameter_below_one(CHIRAL_SCALE_GEV))

    def test_resonance_loop_rule(self):
        threshold = LOWEST_RESONANCE_PION_THRESHOLD_GEV
        just_below = threshold - sp.Rational(1, 10 ** 6)
        self.assertEqual(resonance_loop_uncertainty(just_below, True), 0)
        self.assertEqual(resonance_loop_uncertainty(threshold, True), RESONANCE_LOOP_FRACTION)
        self.assertEqual(resonance_loop_uncertainty(sp.Integer(1), True), sp.Rational(3, 10))
        for sqrt_s in (sp.Rational(3, 10), just_below, threshold, sp.Integer(1)):
            self.assertEqual(resonance_loop_uncertainty(sqrt_s, False), RESONANCE_LOOP_FRACTION)

    def test_quadrature(self):
        sqrt_s = sp.Integer(1)
        total = total_uncertainty(sqrt_s, True)
        expected = math.sqrt(float(truncation_uncertainty(sqrt_s)) ** 2 + 0.3 ** 2)
        self.assertAlmostEqual(float(total), expected, places=14)
        self.assertEqual(total_uncertainty(sp.Rational(1, 2), True), sp.Rational(1, 3))
        for value in (sp.Rational(1, 2), sp.Integer(1)):
            for inputs_from_data in (True, False):
                total = float(total_uncertainty(value, inputs_from_data))
                self.assertGreaterEqual(total, float(truncation_uncertainty(value)))
                self.assertGreaterEqual(
                    total, float(resonance_loop_uncertainty(value, inputs_from_data)))

    def test_declared_range(self):
        self.assertTrue(in_declared_range(sp.Integer(1)))
        self.assertFalse(in_declared_range(sp.Rational(6, 5)))

    def test_check_truncation_runs(self):
        self.assertTrue(check_truncation([sp.Rational(1, 2), sp.Integer(1)]))
        self.assertTrue(check_truncation())


def main():
    """Parse the command line and run the requested checks."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--orders", action="store_true",
                        help="print the delta order table and assert the orders")
    parser.add_argument("--loops", action="store_true",
                        help="classify loop topologies and print thresholds")
    parser.add_argument("--scheme", action="store_true",
                        help="subtraction scheme, Gamma_i, running of L_i^r and"
                             " mu_chi-independent combinations")
    parser.add_argument("--truncation", nargs="*", metavar="SQRT_S", default=None,
                        help="truncation and excluded resonance-loop uncertainty per"
                             " unit |Delta| at the given sqrt(s) values in GeV"
                             " (default grid: " + " ".join(DEFAULT_SQRT_S_GRID_GEV) + ")")
    parser.add_argument("--test", action="store_true", help="run the unit tests")
    arguments = parser.parse_args()

    if arguments.test:
        result = unittest.main(argv=[sys.argv[0]], exit=False).result
        return 0 if result.wasSuccessful() else 1

    selected = {"orders": arguments.orders, "loops": arguments.loops,
                "scheme": arguments.scheme, "truncation": arguments.truncation is not None}
    run_all = not any(selected.values())
    sqrt_s_values = [sp.Rational(value) for value in (arguments.truncation or [])]
    checks = {"orders": check_orders, "loops": check_loops, "scheme": check_scheme,
              "truncation": lambda: check_truncation(sqrt_s_values)}
    for index, (name, check) in enumerate(checks.items()):
        if run_all or selected[name]:
            if run_all and index > 0:
                print()
            check()
    return 0


if __name__ == "__main__":
    sys.exit(main())
