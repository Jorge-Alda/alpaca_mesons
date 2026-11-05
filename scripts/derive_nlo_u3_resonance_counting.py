#!/usr/bin/env python3
"""Assign delta orders to the U(3) NLO Lagrangian and to resonance exchange.

Purpose
    Power counting for the OpenSpec change `fix-u3-nlo-resonance-counting`
    (task 2.2). The script tabulates the order in the combined expansion
    parameter delta of Kaiser:1998ds and Kaiser:2000gs for every operator of
    memory/nlo/lagrangian.md, for the resonance couplings and masses of
    Ecker:1988te, and for the low-energy limit of single-resonance exchange.
    It only classifies orders; no propagator is expanded in the calculation.

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

Units
    Natural units. Only powers of N_c and of the chiral dimension are
    tabulated; no dimensionful number is used.

Inputs
    memory/nlo/lagrangian.md (L^(0), L^(1), L^(2) operators, Alda:2026nrd
    conventions), Kaiser:1998ds, Kaiser:2000gs (delta counting),
    Ecker:1988te (resonance couplings and their N_c scaling).

Invocation
    .venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --orders
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


def main():
    """Parse the command line and run the requested checks."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--orders", action="store_true",
                        help="print the delta order table and assert the orders")
    parser.add_argument("--test", action="store_true", help="run the unit tests")
    arguments = parser.parse_args()

    if arguments.test:
        result = unittest.main(argv=[sys.argv[0]], exit=False).result
        return 0 if result.wasSuccessful() else 1

    selected = {"orders": arguments.orders}
    run_all = not any(selected.values())
    checks = {"orders": check_orders}
    for name, check in checks.items():
        if run_all or selected[name]:
            check()
    return 0


if __name__ == "__main__":
    sys.exit(main())
