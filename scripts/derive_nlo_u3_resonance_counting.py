#!/usr/bin/env python3
"""Assign delta orders to the U(3) NLO Lagrangian, resonance exchange and loops.

Purpose
    Power counting for the OpenSpec change `fix-u3-nlo-resonance-counting`
    (tasks 2.2 and 3.1).
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
    Neither mode expands a propagator or evaluates a loop integral; they
    only classify orders.

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

Units
    Natural units. The order tables use only powers of N_c and of the
    chiral dimension. Masses and thresholds in --loops are in GeV.

Inputs
    memory/nlo/lagrangian.md (L^(0), L^(1), L^(2) operators, Alda:2026nrd
    conventions), Kaiser:1998ds, Kaiser:2000gs (delta counting),
    Ecker:1988te (resonance couplings and their N_c scaling),
    Rosell:2004mn, Rosell:2005ai, Rosell:2006dt, Pich:2010sm (one-loop
    resonance chiral theory: loops with internal resonances, their O(p^2)
    divergences and their R pi and R R' cuts), memory/rcht/foundations/
    vector-width-lineshape-prescriptions.md (line shape near the pole).
    Masses from ParticleDataGroup:2024cfk, in GeV:
    M_rho = 0.77526, M_omega = 0.78266, M_K*(892)^+- = 0.89167,
    m_pi^+- = 0.13957, m_K^+- = 0.493677, m_eta = 0.547862.

Invocation
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --orders
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --loops
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py      # all checks
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --test
"""

import argparse
import sys
import unittest
from dataclasses import dataclass

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
    return (goldstone_lagrangian_rows() + resonance_rows() + exchange_rows()
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
          " analytic in s and m_q. Its effect then has the form of local terms"
          " at O(p^2) and O(p^4), so it is absorbed only if F, B_0, tau"
          " (m_0^2), Lambda_1, Lambda_2 and L_i^r(mu) are fixed from data"
          " (physical or fitted values), not from large-N_c or"
          " resonance-saturation estimates. Its non-analytic part starts at"
          " the R pi threshold, and its size enters the truncation"
          " uncertainty.")
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


def main():
    """Parse the command line and run the requested checks."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--orders", action="store_true",
                        help="print the delta order table and assert the orders")
    parser.add_argument("--loops", action="store_true",
                        help="classify loop topologies and print thresholds")
    parser.add_argument("--test", action="store_true", help="run the unit tests")
    arguments = parser.parse_args()

    if arguments.test:
        result = unittest.main(argv=[sys.argv[0]], exit=False).result
        return 0 if result.wasSuccessful() else 1

    selected = {"orders": arguments.orders, "loops": arguments.loops}
    run_all = not any(selected.values())
    checks = {"orders": check_orders, "loops": check_loops}
    for index, (name, check) in enumerate(checks.items()):
        if run_all or selected[name]:
            if run_all and index > 0:
                print()
            check()
    return 0


if __name__ == "__main__":
    sys.exit(main())
