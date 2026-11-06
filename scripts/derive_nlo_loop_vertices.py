#!/usr/bin/env python3
"""Cubic and quartic vertices of the leading U(3) chiral Lagrangian with ALP sources.

Run:
    .venv/bin/python scripts/derive_nlo_loop_vertices.py            # all checks
    .venv/bin/python scripts/derive_nlo_loop_vertices.py --cubic    # Phi^3, a Phi^2
    .venv/bin/python scripts/derive_nlo_loop_vertices.py --quartic  # Phi^4, a Phi^3
    .venv/bin/python scripts/derive_nlo_loop_vertices.py --test     # unit tests only

Conventions (Alda:2026nrd, Appendix A and Sec. III; memory records
memory/extracts/Alda:2026nrd/appendix-a-expansion.md and
memory/extracts/Alda:2026nrd/lo-topological-action.md):

* u = exp(i phi), u^2 = U, phi = Phi/(2 f_pi), Phi = sum_P lambda_P P with
  tr(lambda_i lambda_j) = 2 delta_ij.  U = exp(i sqrt2 Phi/f) is NOT used.
* Commutator and anticommutator maps C_phi[A] = [phi, A], A_phi[A] = {phi, A}.
* Vielbein at first order in 1/f_a (Eq. A.13):
      u_mu = -2 sin(C)/C [d phi] + 2 (d a/f_a) cos(C)[c_q^A]
             - 2 i (d a/f_a) sin(C)[c_q^V].
* chi_+ = 2 cos(A)[chi_S] + 2 sin(A)[chi_P],
  chi_- = 2 i cos(A)[chi_P] - 2 i sin(A)[chi_S]  (chi_- only for the NLO
  cross-check), chi_S = chi_0 + i (a/f_a)[c_q^S, chi_0],
  chi_P = (a/f_a){c_q^P, chi_0}, chi_0 = 2 B_0 m_q.
* X_TC = 2 i (tr log u + i c_G a/f_a), with log u = i phi exactly.
* L^(0) = f_pi^2/4 <u_mu u^mu> + f_pi^2/4 <chi_+> - f_pi^2 m_0^2/(4 N_f) X_TC^2.

Units: natural units; Phi, a, f_pi, f_a in GeV, chi_0 in GeV^2.
The flavour matrices Phi, d Phi, chi_0, c_q^{A,V,S,P} are noncommuting
symbols; no commutation is assumed.  c_G is a real number.  Traces are
made cyclic by rotating each word to a canonical order.  The symbols
partial_mu_Phi and partial_mu_a denote d_mu Phi and d_mu a; repeated
Lorentz indices in a product are contracted.

Assumptions: the ALP enters at first order in 1/f_a (terms a^2 dropped);
no CP-violating source and no constant theta term are included.  The
all-orders map series are truncated at the total field order 4 needed
here; the truncation is exact for the coefficients that are extracted.
"""

from __future__ import annotations

import argparse
import functools
import sys
import unittest

import sympy as sp

from derive_nlo_quadratic_chipt import derive as derive_nlo_quadratic
from derive_nlo_quadratic_chipt import trace as cyclic_trace
from symbolic_utils import anticommutator, commutator, u3_generators

MESON_COUNTER, ALP_COUNTER = sp.symbols("t_meson t_alp", real=True)
F_PI, F_A, ALP_FIELD, ALP_DERIVATIVE, C_G = sp.symbols(
    "f_pi f_a a partial_mu_a c_G", real=True)
TOPOLOGICAL_MASS, NUMBER_OF_FLAVOURS = sp.symbols("m_0 N_f", positive=True)
MAXIMUM_TOTAL_ORDER = 4

NONCOMMUTING = dict(zip(
    ("Phi", "partial_mu_Phi", "chi_0", "c_q^A", "c_q^V", "c_q^S", "c_q^P"),
    sp.symbols("Phi partial_mu_Phi chi_0 c_q^A c_q^V c_q^S c_q^P",
               commutative=False)))
SOURCE_NAMES = ("c_q^A", "c_q^V", "c_q^S", "c_q^P", "c_G")
VERTEX_CLASSES = {"Phi^3": (3, 0), "a Phi^2": (2, 1),
                  "Phi^4": (4, 0), "a Phi^3": (3, 1)}


# ---------------------------------------------------------------- truncation
def term_orders(term):
    """Return (meson order, ALP order) of one expanded monomial."""
    commuting_part = sp.Mul(*term.args_cnc()[0])
    powers = commuting_part.as_powers_dict()
    return int(powers.get(MESON_COUNTER, 0)), int(powers.get(ALP_COUNTER, 0))


def truncate(expression, maximum_total_order):
    """Drop terms with ALP order > 1 or total field order above the maximum."""
    if isinstance(expression, sp.MatrixBase):
        return expression.applyfunc(lambda entry: truncate(entry, maximum_total_order))
    kept = []
    for term in sp.Add.make_args(sp.expand(expression)):
        meson_order, alp_order = term_orders(term)
        if alp_order <= 1 and meson_order + alp_order <= maximum_total_order:
            kept.append(term)
    return sp.Add(*kept)


def coefficient_of_orders(expression, meson_order, alp_order):
    """Coefficient of t_meson^meson_order t_alp^alp_order (counters removed)."""
    selected = [term for term in sp.Add.make_args(sp.expand(expression))
                if term_orders(term) == (meson_order, alp_order)]
    return sp.Add(*selected).subs({MESON_COUNTER: 1, ALP_COUNTER: 1})


# ------------------------------------------------------- all-orders map series
def map_function_series(function_of_x, linear_map, argument, maximum_order):
    """Apply f(map)[argument] = sum_n f_n map^n[argument] up to map^maximum_order.

    The Taylor coefficients f_n are taken from SymPy's series of f(x).  Each
    map application raises the meson order by one, so the truncation is exact
    for all terms up to meson order maximum_order.
    """
    x = sp.Symbol("x")
    taylor = sp.series(function_of_x(x), x, 0, maximum_order + 1).removeO()
    result = 0 * argument
    power_of_map = argument
    for order in range(maximum_order + 1):
        coefficient = taylor.coeff(x, order)
        if coefficient != 0:
            result = result + coefficient * power_of_map
        power_of_map = linear_map(power_of_map)
    return result


def building_blocks(fields, maximum_total_order=MAXIMUM_TOTAL_ORDER):
    """Return u_mu, chi_+, chi_-, phi from the Appendix A all-orders formulas.

    ``fields`` maps the names in NONCOMMUTING to either noncommuting symbols
    or explicit 3x3 matrices; the same algebra applies to both.
    """
    phi = MESON_COUNTER * fields["Phi"] / (2 * F_PI)
    derivative_phi = MESON_COUNTER * fields["partial_mu_Phi"] / (2 * F_PI)
    alp_over_fa = ALP_COUNTER * ALP_FIELD / F_A
    alp_derivative_over_fa = ALP_COUNTER * ALP_DERIVATIVE / F_A
    chi_0 = fields["chi_0"]

    def commutator_map(argument):
        return truncate(commutator(phi, argument), maximum_total_order)

    def anticommutator_map(argument):
        return truncate(anticommutator(phi, argument), maximum_total_order)

    order = maximum_total_order
    u_mu = (-2 * map_function_series(lambda x: sp.sin(x) / x, commutator_map,
                                      derivative_phi, order)
            + 2 * alp_derivative_over_fa
            * map_function_series(sp.cos, commutator_map, fields["c_q^A"], order)
            - 2 * sp.I * alp_derivative_over_fa
            * map_function_series(sp.sin, commutator_map, fields["c_q^V"], order))
    chi_scalar = chi_0 + sp.I * alp_over_fa * commutator(fields["c_q^S"], chi_0)
    chi_pseudoscalar = alp_over_fa * anticommutator(fields["c_q^P"], chi_0)
    chi_plus = (2 * map_function_series(sp.cos, anticommutator_map, chi_scalar, order)
                + 2 * map_function_series(sp.sin, anticommutator_map,
                                          chi_pseudoscalar, order))
    chi_minus = (2 * sp.I * map_function_series(sp.cos, anticommutator_map,
                                                chi_pseudoscalar, order)
                 - 2 * sp.I * map_function_series(sp.sin, anticommutator_map,
                                                  chi_scalar, order))
    return {"phi": phi, "u_mu": truncate(u_mu, order),
            "chi_plus": truncate(chi_plus, order),
            "chi_minus": truncate(chi_minus, order)}


def leading_lagrangian(fields, trace, maximum_total_order=MAXIMUM_TOTAL_ORDER):
    """Return L^(0) truncated at the given total field order, first order in 1/f_a."""
    blocks = building_blocks(fields, maximum_total_order)
    trace_log_u = trace(sp.I * blocks["phi"])
    topological = 2 * sp.I * (trace_log_u + sp.I * C_G * ALP_COUNTER * ALP_FIELD / F_A)
    kinetic_product = truncate(blocks["u_mu"] * blocks["u_mu"], maximum_total_order)
    lagrangian = (F_PI**2 / 4 * trace(kinetic_product)
                  + F_PI**2 / 4 * trace(blocks["chi_plus"])
                  - F_PI**2 * TOPOLOGICAL_MASS**2 / (4 * NUMBER_OF_FLAVOURS)
                  * topological**2)
    return truncate(lagrangian, maximum_total_order)


# --------------------------------------------------------- symbolic vertices
def symbolic_vertices():
    """Vertices of L^(0) with noncommuting flavour symbols, by (meson, ALP) order."""
    lagrangian = leading_lagrangian(NONCOMMUTING, cyclic_trace)
    orders = [(2, 0), (1, 1)] + list(VERTEX_CLASSES.values())
    return {pair: coefficient_of_orders(lagrangian, *pair) for pair in orders}


def split_by_source(expression):
    """Split a term linear in the ALP into the five source contributions."""
    contributions = {name: sp.Integer(0) for name in SOURCE_NAMES}
    for term in sp.Add.make_args(sp.expand(expression)):
        names = {str(matrix) for matrix in term.atoms(sp.MatrixSymbol)}
        if term.has(C_G):
            names.add("c_G")
        present = names.intersection(SOURCE_NAMES)
        if len(present) != 1:
            raise AssertionError(f"term not linear in a single source: {term}")
        contributions[present.pop()] += term
    return contributions


# ------------------------------------------------- explicit 3x3 matrix checks
def generic_matrix(label):
    """3x3 matrix with nine independent symbolic entries.

    No Hermiticity is imposed: the selection rules tested with it follow
    from trace cyclicity alone, so a generic matrix is a stronger test and
    keeps the explicit expansion fast (no sqrt or i in the entries).
    """
    return sp.Matrix(3, 3, lambda row, column: sp.Symbol(f"{label}_{row}{column}"))


def explicit_fields(quark_masses, sources, nonet_generators=False):
    """Explicit 3x3 fields: generic matrices, or the nonet sum_P lambda_P P."""
    if nonet_generators:
        generators = u3_generators()
        meson_fields = sp.symbols("P_0:9", real=True)
        meson_derivatives = sp.symbols("dP_0:9", real=True)
        meson_matrix = sum((field * generator for field, generator
                            in zip(meson_fields, generators)), sp.zeros(3))
        derivative_matrix = sum((field * generator for field, generator
                                 in zip(meson_derivatives, generators)), sp.zeros(3))
    else:
        meson_matrix = generic_matrix("Phi")
        derivative_matrix = generic_matrix("dPhi")
    fields = {
        "Phi": meson_matrix,
        "partial_mu_Phi": derivative_matrix,
        "chi_0": 2 * sp.Symbol("B_0", positive=True) * sp.diag(*quark_masses),
    }
    for name in ("c_q^A", "c_q^V", "c_q^S", "c_q^P"):
        fields[name] = sources.get(name, sp.zeros(3))
    return fields


def explicit_cubic_alp_vertex(quark_masses, sources, with_topological_source=False):
    """a Phi^2 coefficient of L^(0) for explicit matrices (one source at a time)."""
    fields = explicit_fields(quark_masses, sources)
    lagrangian = leading_lagrangian(fields, lambda matrix: sp.Matrix(matrix).trace(),
                                    maximum_total_order=3)
    if not with_topological_source:
        lagrangian = lagrangian.subs(C_G, 0)
    return sp.expand(coefficient_of_orders(lagrangian, 2, 1))


GENERIC_MASSES = sp.symbols("m_u m_d m_s", positive=True)
DEGENERATE_LIGHT_MASSES = (sp.Symbol("m_l", positive=True),) * 2 + (
    sp.Symbol("m_s", positive=True),)


def commuting_scalar_source():
    """c_q^S block-diagonal in (u,d)+(s): commutes with diag(m_l,m_l,m_s)."""
    generic = generic_matrix("cS")
    for light in range(2):
        generic[light, 2] = 0
        generic[2, light] = 0
    return generic


@functools.cache
def explicit_vanishing_checks():
    """Return {label: a Phi^2 coefficient} for the cases that must vanish."""
    return {
        "c_q^A generic": explicit_cubic_alp_vertex(
            GENERIC_MASSES, {"c_q^A": generic_matrix("cA")}),
        "c_q^P generic": explicit_cubic_alp_vertex(
            GENERIC_MASSES, {"c_q^P": generic_matrix("cP")}),
        "c_G": explicit_cubic_alp_vertex(
            GENERIC_MASSES, {}, with_topological_source=True),
        "c_q^V proportional to 1": explicit_cubic_alp_vertex(
            GENERIC_MASSES, {"c_q^V": sp.Symbol("c_V", real=True) * sp.eye(3)}),
        "c_q^S with [c_q^S, m_q]=0": explicit_cubic_alp_vertex(
            DEGENERATE_LIGHT_MASSES, {"c_q^S": commuting_scalar_source()}),
    }


@functools.cache
def explicit_nonvanishing_checks():
    """Return {label: a Phi^2 coefficient} for generic noncommuting c_q^V, c_q^S."""
    return {
        "c_q^V generic": explicit_cubic_alp_vertex(
            GENERIC_MASSES, {"c_q^V": generic_matrix("cV")}),
        "c_q^S generic": explicit_cubic_alp_vertex(
            GENERIC_MASSES, {"c_q^S": generic_matrix("cS")}),
    }


# ------------------------------------------------- quadratic cross-checks
def quadratic_record_blocks():
    """Building-block expansions quoted in memory/nlo/quadratic-lagrangian.md."""
    phi = NONCOMMUTING["Phi"] / (2 * F_PI)
    chi_0 = NONCOMMUTING["chi_0"]
    c_P, c_S = NONCOMMUTING["c_q^P"], NONCOMMUTING["c_q^S"]
    alp_over_fa = ALP_FIELD / F_A
    return {
        ("u_mu", 1, 0): -NONCOMMUTING["partial_mu_Phi"] / F_PI,
        ("u_mu", 0, 1): 2 * ALP_DERIVATIVE / F_A * NONCOMMUTING["c_q^A"],
        ("chi_plus", 0, 0): 2 * chi_0,
        ("chi_plus", 2, 0): -anticommutator(phi, anticommutator(phi, chi_0)),
        ("chi_plus", 1, 1): 2 * alp_over_fa * anticommutator(phi, anticommutator(c_P, chi_0)),
        ("chi_minus", 0, 1): 2 * sp.I * alp_over_fa * anticommutator(c_P, chi_0),
        ("chi_minus", 1, 0): -2 * sp.I * anticommutator(phi, chi_0),
        ("chi_minus", 1, 1): 2 * alp_over_fa * anticommutator(phi, commutator(c_S, chi_0)),
    }


def nlo_quadratic_from_blocks():
    """Recompute the NLO Phi^2 and a Phi coefficients from the all-orders blocks."""
    blocks = building_blocks(NONCOMMUTING, maximum_total_order=2)
    u_mu, chi_plus, chi_minus = blocks["u_mu"], blocks["chi_plus"], blocks["chi_minus"]
    topological = 2 * sp.I * (cyclic_trace(sp.I * blocks["phi"])
                              + sp.I * C_G * ALP_COUNTER * ALP_FIELD / F_A)
    trace = cyclic_trace
    operators = {
        "L_5": trace(truncate(u_mu * u_mu * chi_plus, 2)),
        "L_8": trace(truncate(chi_plus**2 + chi_minus**2, 2)) / 2,
        "Lambda_1": F_PI**2 * trace(u_mu)**2 / 12,
        "Lambda_2": sp.I * F_PI**2 * topological * trace(chi_minus) / 12,
        "L_4": trace(truncate(u_mu * u_mu, 2)) * trace(chi_plus),
        "L_6": trace(chi_plus)**2,
        "L_7": trace(chi_minus)**2,
        "L_18": trace(u_mu) * trace(truncate(u_mu * chi_plus, 2)),
        "L_25": -sp.I * topological * trace(truncate(chi_plus * chi_minus, 2)),
    }
    return {name: {"PhiPhi": coefficient_of_orders(expression, 2, 0),
                   "aPhi": coefficient_of_orders(expression, 1, 1)}
            for name, expression in operators.items()}


# ------------------------------------------------------------------ direct definitions
def exponential(sign, field, maximum_order):
    """exp(sign i field) as a finite series (field carries the meson counter)."""
    return sum(((sign * sp.I * field)**power / sp.factorial(power)
                for power in range(maximum_order + 1)), sp.Integer(0))


def differentiated_exponential(sign, field, derivative_field, maximum_order):
    """d_mu exp(sign i field) by the product rule on each noncommuting power."""
    total = sp.Integer(0)
    for power in range(1, maximum_order + 1):
        prefactor = (sign * sp.I)**power / sp.factorial(power)
        total += prefactor * sum((field**position * derivative_field
                                  * field**(power - position - 1)
                                  for position in range(power)), sp.Integer(0))
    return total


# ----------------------------------------------------------------- unit tests
class LoopVertexTests(unittest.TestCase):
    """Checks of the all-orders expansion and of the vertex selection rules."""

    @classmethod
    def setUpClass(cls):
        cls.vertices = symbolic_vertices()

    def test_series_match_direct_definitions(self):
        """sin/cos map series equal u^dagger (.) u products through order 4.

        Assumption used only here: right/left sources r = (d a/f_a)(c_V + c_A),
        l = (d a/f_a)(c_V - c_A), u_mu = i[u^dag(d - i r)u - u(d - i l)u^dag];
        chi_+ = u^dag chi u^dag + u chi^dag u, chi_- = u^dag chi u^dag - u chi^dag u.
        """
        order = MAXIMUM_TOTAL_ORDER
        blocks = building_blocks(NONCOMMUTING, order)
        phi = blocks["phi"]
        derivative_phi = MESON_COUNTER * NONCOMMUTING["partial_mu_Phi"] / (2 * F_PI)
        u_matrix = exponential(1, phi, order)
        u_dagger = exponential(-1, phi, order)
        derivative_u = differentiated_exponential(1, phi, derivative_phi, order)
        derivative_u_dagger = differentiated_exponential(-1, phi, derivative_phi, order)
        alp_derivative = ALP_COUNTER * ALP_DERIVATIVE / F_A
        right_source = alp_derivative * (NONCOMMUTING["c_q^V"] + NONCOMMUTING["c_q^A"])
        left_source = alp_derivative * (NONCOMMUTING["c_q^V"] - NONCOMMUTING["c_q^A"])
        direct_u_mu = sp.I * (u_dagger * derivative_u - sp.I * u_dagger * right_source * u_matrix
                              - u_matrix * derivative_u_dagger
                              + sp.I * u_matrix * left_source * u_dagger)
        self.assertEqual(sp.expand(truncate(direct_u_mu, order) - blocks["u_mu"]), 0)
        chi_0 = NONCOMMUTING["chi_0"]
        alp_over_fa = ALP_COUNTER * ALP_FIELD / F_A
        chi_scalar = chi_0 + sp.I * alp_over_fa * commutator(NONCOMMUTING["c_q^S"], chi_0)
        chi_pseudoscalar = alp_over_fa * anticommutator(NONCOMMUTING["c_q^P"], chi_0)
        chi = chi_scalar + sp.I * chi_pseudoscalar
        chi_dagger = chi_scalar - sp.I * chi_pseudoscalar
        direct_plus = u_dagger * chi * u_dagger + u_matrix * chi_dagger * u_matrix
        direct_minus = u_dagger * chi * u_dagger - u_matrix * chi_dagger * u_matrix
        self.assertEqual(sp.expand(truncate(direct_plus, order) - blocks["chi_plus"]), 0)
        self.assertEqual(sp.expand(truncate(direct_minus, order) - blocks["chi_minus"]), 0)

    def test_blocks_reproduce_quadratic_record(self):
        blocks = building_blocks(NONCOMMUTING, 2)
        for (name, meson_order, alp_order), expected in quadratic_record_blocks().items():
            obtained = coefficient_of_orders(blocks[name], meson_order, alp_order)
            self.assertEqual(sp.expand(obtained - expected), 0, (name, meson_order, alp_order))

    def test_nlo_quadratic_terms_reproduce_record_script(self):
        reference = derive_nlo_quadratic()["results"]
        obtained = nlo_quadratic_from_blocks()
        for name, sectors in reference.items():
            for sector, expression in sectors.items():
                self.assertEqual(sp.simplify(sp.expand(obtained[name][sector] - expression)),
                                 0, (name, sector))

    def test_leading_quadratic_terms(self):
        """Canonical kinetic term and LO masses for the neutral fields."""
        fields = explicit_fields(GENERIC_MASSES, {}, nonet_generators=True)
        lagrangian = leading_lagrangian(fields, lambda matrix: matrix.trace(), 2)
        quadratic = sp.expand(coefficient_of_orders(lagrangian, 2, 0))
        P0, P3 = sp.Symbol("P_0", real=True), sp.Symbol("P_3", real=True)
        dP3 = sp.Symbol("dP_3", real=True)
        b_zero = sp.Symbol("B_0", positive=True)
        m_u, m_d, m_s = GENERIC_MASSES
        self.assertEqual(sp.diff(quadratic, dP3, 2), 1)
        self.assertEqual(sp.simplify(-sp.diff(quadratic, P3, 2) - b_zero * (m_u + m_d)), 0)
        singlet_mass = -sp.diff(quadratic, P0, 2)
        expected = 2 * b_zero * (m_u + m_d + m_s) / 3 + 6 * TOPOLOGICAL_MASS**2 / (
            2 * NUMBER_OF_FLAVOURS)
        self.assertEqual(sp.simplify(singlet_mass - expected), 0)

    def test_no_three_meson_vertex(self):
        self.assertEqual(self.vertices[(3, 0)], 0)

    def test_cubic_alp_vertex_symbolic_selection(self):
        by_source = split_by_source(self.vertices[(2, 1)])
        for name in ("c_q^A", "c_q^P", "c_G"):
            self.assertEqual(by_source[name], 0, name)
        for name in ("c_q^V", "c_q^S"):
            self.assertNotEqual(by_source[name], 0, name)

    def test_cubic_alp_vertex_explicit_vanishing(self):
        for label, expression in explicit_vanishing_checks().items():
            self.assertEqual(expression, 0, label)

    def test_cubic_alp_vertex_explicit_nonvanishing(self):
        for label, expression in explicit_nonvanishing_checks().items():
            self.assertNotEqual(expression, 0, label)


# ------------------------------------------------------------------- printing
LATEX_PLACEHOLDERS = {"Phi": ("FIELDPHI", r"\Phi"),
                      "partial_mu_Phi": ("DERIVATIVEPHI", r"\partial_\mu\Phi"),
                      "chi_0": ("CHIZERO", r"\chi_0"),
                      "c_q^A": ("SOURCEA", r"c_q^A"), "c_q^V": ("SOURCEV", r"c_q^V"),
                      "c_q^S": ("SOURCES", r"c_q^S"), "c_q^P": ("SOURCEP", r"c_q^P")}


def to_latex(expression):
    r"""LaTeX with readable flavour-matrix names; repeated Lorentz indices contracted.

    Placeholder names avoid SymPy's sub/superscript parsing of the names;
    powers of a matrix are printed with parentheses, e.g. (\partial_\mu\Phi)^{2}.
    """
    replacements = {}
    for matrix in expression.atoms(sp.MatrixSymbol):
        placeholder = LATEX_PLACEHOLDERS[str(matrix)][0]
        replacements[matrix] = sp.MatrixSymbol(placeholder, *matrix.shape)
    replacements[ALP_DERIVATIVE] = sp.Symbol("DERIVATIVEALP", real=True)
    text = sp.latex(sp.collect(expression.xreplace(replacements), [F_PI, F_A], evaluate=True))
    for placeholder, latex_name in LATEX_PLACEHOLDERS.values():
        text = text.replace(placeholder + "^{", "(" + latex_name + ")^{")
        text = text.replace(placeholder, latex_name + " ")
    return text.replace("DERIVATIVEALP", r"\partial^\mu a ")


def print_class(vertices, label):
    meson_order, alp_order = VERTEX_CLASSES[label]
    expression = vertices[(meson_order, alp_order)]
    if alp_order == 0:
        print(f"[{label}] (source independent):\n{to_latex(expression) if expression != 0 else '0'}\n")
        return
    for source, contribution in split_by_source(expression).items():
        text = to_latex(contribution) if contribution != 0 else "0"
        print(f"[{label}] source {source}:\n{text}\n")


def summary_table(vertices, explicit_zero, explicit_nonzero):
    rows = []
    for label, (meson_order, alp_order) in VERTEX_CLASSES.items():
        expression = vertices[(meson_order, alp_order)]
        if alp_order == 0:
            rows.append(("none (a^0)", label, expression == 0))
            continue
        for source, contribution in split_by_source(expression).items():
            rows.append((source + " (generic)", label, contribution == 0))
    for label, value in explicit_zero.items():
        rows.append((label + " [3x3]", "a Phi^2", value == 0))
    for label, value in explicit_nonzero.items():
        rows.append((label + " [3x3]", "a Phi^2", value == 0))
    print(f"{'source':<32} {'vertex class':<12} result")
    for source, label, is_zero in rows:
        print(f"{source:<32} {label:<12} {'zero' if is_zero else 'non-zero'}")


def run_all(show_cubic, show_quartic, run_tests):
    vertices = symbolic_vertices()
    assert vertices[(3, 0)] == 0, "Phi^3 terms do not vanish"
    cubic_by_source = split_by_source(vertices[(2, 1)])
    for name in ("c_q^A", "c_q^P", "c_G"):
        assert cubic_by_source[name] == 0, f"a Phi^2 term from {name} does not vanish"
    for name in ("c_q^V", "c_q^S"):
        assert cubic_by_source[name] != 0, f"a Phi^2 term from {name} vanishes"
    explicit_zero = explicit_vanishing_checks()
    explicit_nonzero = explicit_nonvanishing_checks()
    for label, value in explicit_zero.items():
        assert value == 0, f"a Phi^2 term does not vanish for {label}"
    for label, value in explicit_nonzero.items():
        assert value != 0, f"a Phi^2 term vanishes for {label}"
    print("Conventions: u=exp(i Phi/(2 f_pi)); natural units (GeV); first order in 1/f_a;")
    print("partial_mu fields with repeated Lorentz index contracted; tr = flavour trace.\n")
    if show_cubic:
        print("[Phi^2] (cross-check, LO quadratic):\n" + to_latex(vertices[(2, 0)]) + "\n")
        print("[a Phi] (cross-check, LO quadratic):\n" + to_latex(vertices[(1, 1)]) + "\n")
        print_class(vertices, "Phi^3")
        print_class(vertices, "a Phi^2")
    if show_quartic:
        print_class(vertices, "Phi^4")
        print_class(vertices, "a Phi^3")
    summary_table(vertices, explicit_zero, explicit_nonzero)
    if run_tests:
        result = unittest.TextTestRunner(verbosity=1).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(LoopVertexTests))
        if not result.wasSuccessful():
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cubic", action="store_true", help="print Phi^3 and a Phi^2")
    parser.add_argument("--quartic", action="store_true", help="print Phi^4 and a Phi^3")
    parser.add_argument("--test", action="store_true", help="run the unit tests only")
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=[sys.argv[0]], defaultTest="LoopVertexTests")
        return
    no_option = not (arguments.cubic or arguments.quartic)
    run_all(show_cubic=arguments.cubic or no_option,
            show_quartic=arguments.quartic or no_option,
            run_tests=no_option)


if __name__ == "__main__":
    main()
