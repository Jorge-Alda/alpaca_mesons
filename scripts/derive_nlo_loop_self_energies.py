#!/usr/bin/env python3
"""One-loop meson self-energy Sigma_ij(s) of the U(3) theory and its pole cancellation.

OpenSpec change `calculate-nlo-loop-self-energies`, tasks 2.1 and 2.3.

Run:
    .venv/bin/python scripts/derive_nlo_loop_self_energies.py           # all checks
    .venv/bin/python scripts/derive_nlo_loop_self_energies.py --meson   # task 2.1
    .venv/bin/python scripts/derive_nlo_loop_self_energies.py --cancel  # task 2.3
    .venv/bin/python scripts/derive_nlo_loop_self_energies.py --test    # unit tests only

Conventions (Alda:2026nrd; memory/nlo/loops/vertices-and-topologies.md):

* u = exp(i Phi/(2 f_pi)), Phi = sum_a lambda_a P_a, tr(lambda_a lambda_b) =
  2 delta_ab, generators ordered (lambda_0, ..., lambda_8) as in
  symbolic_utils.u3_generators; chi_0 = 2 B_0 diag(m_u, m_d, m_s).
* L^(0) = f_pi^2/4 <u_mu u^mu> + f_pi^2/4 <chi_+> - f_pi^2 m_0^2/(4 N_f) X_TC^2.
  The LO kinetic matrix is the unit matrix; the LO mass matrix is
  M^2_ab = 1/4 tr({lambda_a, lambda_b} chi_0) + m_0^2 delta_a0 delta_b0.
* Quartic vertex (scripts/derive_nlo_loop_vertices.py):
  L_Phi^4 = (1/f_pi^2)[-1/24 <Phi^2 dPhi dPhi> + 1/24 <Phi dPhi Phi dPhi>
                       + 1/48 <Phi^4 chi_0>].
  The X_TC^2 term is exactly quadratic in Phi and gives no vertex.
* Inverse propagator Q(s) = s K - M^2 + Sigma(s), s = p^2
  (memory/nlo/loop-corrections-to-alp-meson-lsz-mixing.md). The one-loop
  quadratic Lagrangian is L_loop = 1/2 Sigma1_ab dP_a dP_b + 1/2 Sigma0_ab P_a P_b,
  so Sigma(s) = Sigma0 + s Sigma1, and the pole of one field is at
  s = (M^2 - Sigma0)/(1 + Sigma1).
* Tadpole (Gasser:1983yg): for a mass eigenstate k,
  <phi_k phi_k> = 2 m_k^2 lambda(mu) + m_k^2/(16 pi^2) log(m_k^2/mu^2),
  <d_mu phi_k d^mu phi_k> = m_k^2 <phi_k phi_k>, <phi_k d_mu phi_k> = 0
  (dimensional regularization). In the generator basis the tadpole matrix is
  <phi_c phi_d> = sum_k O_ck O_dk <phi_k phi_k>, with M^2 = O diag(m_k^2) O^T,
  i.e. the matrix function of M^2; <d phi_c d phi_d> = (M^2 <phi phi>)_cd.
* lambda(mu) is the pole of Gasser:1983yg; L_i = L_i^r(mu) + Gamma_i lambda(mu),
  B = B^r(mu)(1 + 4 tau lambda/F^4) (Kaiser:1998ds); F is identified with f_pi.
  The finite part (MS-bar - 1 scheme) is the tadpole with lambda set to zero.

Units: GeV for masses, f_pi, B_0, mu; GeV^2 for s, chi_0 and Sigma0; Sigma1
is dimensionless.

Assumptions: no CP-violating source, no theta term, no ALP (order a^0); the
eta' is a light field with m_0^2 of the order of the quark-mass terms
(memory/nlo/counting/power-counting.md); one-loop order, so products of a
loop term with a running of B^r or L_i^r are dropped.

--meson (task 2.1): prints the tadpole-contracted loop Lagrangian as a trace
formula, Sigma0 and Sigma1 for the full nonet in the (3,8,0) basis and for
the charged mesons, split into pole and finite parts, and their values at the
test point for mu in {0.5, 0.77, 1.0} GeV, also in the LO mass eigenbasis.
Asserts: symmetry; vanishing in the chiral limit (also with m_0 != 0);
in the SU(3) octet theory (singlet removed) the degenerate-mass pole mass
M^2 (1 + 2/3 mu_P) and F_pi = F (1 - 3 mu_P) of Gasser:1984gg, and the
isospin-limit masses and decay constants of Gasser:1984gg.
--cancel (task 2.3): collects the coefficient of lambda of every entry of
Q(s) = s K - M^2 + Sigma with L_i -> L_i^r + Gamma_i lambda,
B -> B^r(1 + k_B m_0^2 lambda/F^2) and Lambda_j -> Lambda_j^r + k_Lambdaj
m_0^2 lambda/F^2, plus the change A^T Q0 + Q0 A of the LO Q0 = s - M^2 under
a field redefinition P -> (1 + lambda A) P, A built from the chirally
covariant maps of Phi linear in chi_0 or m_0^2 (field_redefinition_maps).
The off-shell Q(s) of the field P is finite only up to such a redefinition;
it does not move the pole positions. It solves for all coefficients,
reports which Gamma-combinations Sigma_ij fixes (pole_coefficient_report,
the structured comparison with the U(3) table of alp-gii.1 imported from
derive_nlo_u3_resonance_counting.py), asserts that the SU(3) Gamma_i of
Gasser:1984gg cancel every pole of the octet theory, that the U(3) table
with k_B = 2/3 cancels every pole of the nonet theory (entry by entry, for K
and M^2), and that the renormalized Q(s) changes with mu only by a field
redefinition and the one-loop pole positions not at all, for mu in
{0.5, 0.77, 1.0} GeV. Differences with the table or the records are printed
as lines starting with "FINDING:"; they do not change the exit code.
"""

from __future__ import annotations

import argparse
import functools
import random
import sys
import unittest

import mpmath
import sympy as sp

from derive_nlo_u3_resonance_counting import KAISER_1998DS_U3_GAMMA, SU3_GAMMA
from symbolic_utils import u3_generators

GENERATORS = u3_generators()
NONET_LABELS = tuple(range(9))
OCTET_LABELS = tuple(range(1, 9))
CHARGED_PAIRS = {(1, 2): r"\pi^+", (4, 5): r"K^+", (6, 7): r"K^0"}
CHARGED_TADPOLE_NAMES = {(1, 2): "<pi+pi->", (4, 5): "<K+K->", (6, 7): "<K0K0bar>"}
CHARGED_TADPOLE_LATEX = {"<pi+pi->": r"\langle\pi^+\pi^-\rangle",
                         "<K+K->": r"\langle K^+K^-\rangle",
                         "<K0K0bar>": r"\langle K^0\bar K^0\rangle"}

MASS_UP, MASS_DOWN, MASS_STRANGE = sp.symbols("m_u m_d m_s", positive=True)
B_ZERO, F_PI, SINGLET_MASS = sp.symbols("B_0 f_pi m_0", positive=True)
MANDELSTAM_S = sp.Symbol("s", real=True)
QUARK_MASSES = (MASS_UP, MASS_DOWN, MASS_STRANGE)
MESON_FIELDS = sp.symbols("P_0:9", real=True)
MESON_FIELD_DERIVATIVES = sp.symbols("dP_0:9", real=True)

# Local NLO constants of memory/nlo/quadratic-lagrangian.md.
LOCAL_CONSTANTS = ("L_4", "L_5", "L_6", "L_7", "L_8", "L_18", "L_25",
                   "Lambda_1", "Lambda_2")

# Test point (GeV). Illustrative inputs, not a fit: quark masses near the
# MS-bar values at 2 GeV, B_0 such that B_0 (m_u + m_d) ~ (0.135 GeV)^2,
# f_pi = F = 0.0924 GeV, m_0^2 = 0.7 GeV^2 so that the LO eta' is near 1 GeV.
TEST_POINT = {MASS_UP: sp.Rational(216, 100000), MASS_DOWN: sp.Rational(467, 100000),
              MASS_STRANGE: sp.Rational(934, 10000), B_ZERO: sp.Rational(267, 100),
              F_PI: sp.Rational(924, 10000), SINGLET_MASS: sp.sqrt(sp.Rational(7, 10))}
REFERENCE_SCALE = sp.Rational(77, 100)
SCALES = (sp.Rational(1, 2), REFERENCE_SCALE, sp.Integer(1))
# Illustrative L_i^r(0.77 GeV) used only in the mu-independence check.
ILLUSTRATIVE_CONSTANTS = {"L_4": sp.Rational(-3, 10000), "L_5": sp.Rational(12, 10000),
                          "L_6": sp.Rational(1, 10000), "L_7": sp.Rational(-3, 10000),
                          "L_8": sp.Rational(5, 10000), "L_18": sp.Rational(2, 10000),
                          "L_25": sp.Rational(1, 10000), "Lambda_1": sp.Integer(0),
                          "Lambda_2": sp.Integer(0)}
mpmath.mp.dps = 40

# k_B of B = B^r(1 + 4 tau lambda/F^4) = B^r(1 + k_B m_0^2 lambda/F^2) (Kaiser:1998ds).
# Record memory/nlo/counting/subtraction-scheme-and-scales.md (and
# tau_in_terms_of_eta_prime_mass of derive_nlo_u3_resonance_counting.py) uses
# tau = F^2 m_0^2/(4 N_f), so 4 tau/F^4 = m_0^2/(3 F^2): k_B = 1/3.
# memory/extracts/Kaiser:1998ds/qcd-scale-dependence-za.md quotes
# M_eta'^2 = 6 tau/F^2 at tree level, i.e. tau = F^2 m_0^2/6 and k_B = 2/3.
B_POLE_FROM_RECORD = sp.Rational(1, 3)
B_POLE_FROM_KAISER_ETA_PRIME_MASS = sp.Rational(2, 3)


# ------------------------------------------------------------ LO quantities
def chi_zero(masses=QUARK_MASSES):
    """chi_0 = 2 B_0 diag(m_u, m_d, m_s) in GeV^2."""
    return 2 * B_ZERO * sp.diag(*masses)


def leading_mass_matrix(labels, masses=QUARK_MASSES, singlet_mass=SINGLET_MASS):
    """LO M^2_ab = 1/4 tr({lambda_a, lambda_b} chi_0) + m_0^2 delta_a0 delta_b0."""
    chi = chi_zero(masses)
    size = len(labels)
    matrix = sp.zeros(size)
    for row, first in enumerate(labels):
        for column, second in enumerate(labels):
            left, right = GENERATORS[first], GENERATORS[second]
            matrix[row, column] = sp.expand(((left * right + right * left) * chi).trace() / 4)
            if first == 0 and second == 0:
                matrix[row, column] += singlet_mass**2
    return matrix


def background(labels):
    """Background Phi = sum_a lambda_a P_a and d_mu Phi for the given generators."""
    field = sum((MESON_FIELDS[label] * GENERATORS[label] for label in labels), sp.zeros(3))
    derivative = sum((MESON_FIELD_DERIVATIVES[label] * GENERATORS[label] for label in labels),
                     sp.zeros(3))
    return sp.Matrix(field), sp.Matrix(derivative)


# ------------------------------------------------------------ vertex
def quartic_vertex(field, derivative_field, chi, f_pi=F_PI):
    """L_Phi^4 of memory/nlo/loops/vertices-and-topologies.md for explicit matrices.

    Products of two derivative matrices stand for the contraction over mu.
    """
    return ((-(field * field * derivative_field * derivative_field).trace() / 24
             + (field * derivative_field * field * derivative_field).trace() / 24
             + (field**4 * chi).trace() / 48) / f_pi**2)


# ------------------------------------------------------------ tadpoles
def tadpole_symbol_matrix(labels):
    """Symmetric matrix of <phi_c phi_d> symbols in the generator basis.

    Charged pairs (1,2), (4,5), (6,7) are diagonal and degenerate; the neutral
    block (3,8,0) or (3,8) is a general symmetric matrix (pi0-eta-eta' mixing).
    """
    position = {label: index for index, label in enumerate(labels)}
    matrix = sp.zeros(len(labels))
    for pair, name in CHARGED_TADPOLE_NAMES.items():
        symbol = sp.Symbol(name, real=True)
        for label in pair:
            matrix[position[label], position[label]] = symbol
    neutral = [label for label in labels if label in (3, 8, 0)]
    for first_index, first in enumerate(neutral):
        for second in neutral[first_index:]:
            symbol = sp.Symbol(f"<phi{first}phi{second}>", real=True)
            matrix[position[first], position[second]] = symbol
            matrix[position[second], position[first]] = symbol
    return matrix


TRACE_FORMULA = """  L_loop = 1/f_pi^2 sum_cd { <d phi_c d phi_d>/24 [<Phi l_c Phi l_d> - <Phi^2 l_c l_d>]
           + <phi_c phi_d>/24 [<l_c dPhi l_d dPhi> - <l_c l_d dPhi dPhi>]
           + <phi_c phi_d>/48 [<l_c l_d Phi^2 chi_0> + <Phi l_c l_d Phi chi_0>
             + <Phi^2 l_c l_d chi_0> + <l_c Phi l_d Phi chi_0>
             + <Phi l_c Phi l_d chi_0> + <l_c Phi^2 l_d chi_0>] }"""


def loop_lagrangian(labels, tadpole, derivative_tadpole, chi, f_pi=F_PI):
    """Contract two fields of L_Phi^4 with the tadpole matrices.

    Returns (kinetic part, non-derivative part) of L_loop, quadratic in the
    background P_a, dP_a. Trace formula (contractions with one derivative on
    a loop line and one on the background vanish):
      L_loop = 1/f_pi^2 sum_cd { <d phi_c d phi_d>/24 [<Phi l_c Phi l_d> - <Phi^2 l_c l_d>]
               + <phi_c phi_d>/24 [<l_c dPhi l_d dPhi> - <l_c l_d dPhi dPhi>]
               + <phi_c phi_d>/48 [<l_c l_d Phi^2 chi_0> + <Phi l_c l_d Phi chi_0>
                 + <Phi^2 l_c l_d chi_0> + <l_c Phi l_d Phi chi_0>
                 + <Phi l_c Phi l_d chi_0> + <l_c Phi^2 l_d chi_0>] }.
    """
    field, derivative_field = background(labels)
    kinetic = sp.Integer(0)
    non_derivative = sp.Integer(0)
    for row, first in enumerate(labels):
        for column, second in enumerate(labels):
            weight = tadpole[row, column]
            derivative_weight = derivative_tadpole[row, column]
            if weight == 0 and derivative_weight == 0:
                continue
            left, right = GENERATORS[first], GENERATORS[second]
            non_derivative += derivative_weight * (
                (field * left * field * right).trace()
                - (field * field * left * right).trace()) / 24
            kinetic += weight * (
                (left * derivative_field * right * derivative_field).trace()
                - (left * right * derivative_field * derivative_field).trace()) / 24
            words = (left * right * field * field, field * left * right * field,
                     field * field * left * right, left * field * right * field,
                     field * left * field * right, left * field * field * right)
            non_derivative += weight * sum((word * chi).trace() for word in words) / 48
    return sp.expand(kinetic / f_pi**2), sp.expand(non_derivative / f_pi**2)


@functools.lru_cache(maxsize=None)
def self_energy(field_content="nonet"):
    """Sigma0 and Sigma1 in the generator basis with symbolic tadpoles.

    field_content: "nonet" (U(3), with m_0^2) or "octet" (SU(3), traceless Phi,
    singlet removed from the field content).
    """
    labels = NONET_LABELS if field_content == "nonet" else OCTET_LABELS
    mass_matrix = leading_mass_matrix(labels)
    if field_content == "octet":
        mass_matrix = mass_matrix.subs(SINGLET_MASS, 0)
    tadpole = tadpole_symbol_matrix(labels)
    derivative_tadpole = ((mass_matrix * tadpole + tadpole * mass_matrix) / 2).applyfunc(sp.expand)
    kinetic, non_derivative = loop_lagrangian(labels, tadpole, derivative_tadpole, chi_zero())
    fields = [MESON_FIELDS[label] for label in labels]
    derivatives = [MESON_FIELD_DERIVATIVES[label] for label in labels]
    sigma_zero = sp.hessian(non_derivative, fields).applyfunc(sp.expand)
    sigma_one = sp.hessian(kinetic, derivatives).applyfunc(sp.expand)
    return {"labels": labels, "mass_matrix": mass_matrix, "tadpole": tadpole,
            "sigma_zero": sigma_zero, "sigma_one": sigma_one}


def tadpole_substitution(tadpole, values):
    """Map tadpole symbols to the entries of a matrix of values."""
    substitution = {}
    for row in range(tadpole.rows):
        for column in range(tadpole.cols):
            entry = tadpole[row, column]
            if isinstance(entry, sp.Symbol):
                substitution[entry] = values[row, column]
    return substitution


@functools.lru_cache(maxsize=None)
def pole_part(field_content="nonet"):
    """Coefficient of lambda in Sigma0, Sigma1: <phi phi> -> 2 lambda M^2, divided by lambda."""
    data = self_energy(field_content)
    substitution = tadpole_substitution(data["tadpole"], 2 * data["mass_matrix"])
    return (data["sigma_zero"].subs(substitution).applyfunc(sp.expand),
            data["sigma_one"].subs(substitution).applyfunc(sp.expand))


# ------------------------------------------------------------ local terms
def local_operators(field, derivative_field, chi, f_pi=F_PI):
    """Phi^2 terms of the NLO Lagrangian per constant (memory/nlo/quadratic-lagrangian.md)."""
    trace = lambda matrix: matrix.trace()  # noqa: E731
    kinetic = {
        "L_5": 2 * trace(chi * derivative_field**2) / f_pi**2,
        "Lambda_1": trace(derivative_field)**2 / 12,
        "L_4": 2 * trace(derivative_field**2) * trace(chi) / f_pi**2,
        "L_18": 2 * trace(chi * derivative_field) * trace(derivative_field) / f_pi**2,
    }
    non_derivative = {
        "L_8": -2 * (trace(field**2 * chi**2) + trace(field * chi * field * chi)) / f_pi**2,
        "Lambda_2": -trace(field * chi) * trace(field) / 6,
        "L_6": -4 * trace(field**2 * chi) * trace(chi) / f_pi**2,
        "L_7": -4 * trace(field * chi)**2 / f_pi**2,
        "L_25": 4 * trace(field * chi**2) * trace(field) / f_pi**2,
    }
    return kinetic, non_derivative


@functools.lru_cache(maxsize=None)
def local_matrices(field_content="nonet"):
    """delta K and delta M^2 per unit constant, in the generator basis."""
    labels = NONET_LABELS if field_content == "nonet" else OCTET_LABELS
    field, derivative_field = background(labels)
    kinetic, non_derivative = local_operators(field, derivative_field, chi_zero())
    fields = [MESON_FIELDS[label] for label in labels]
    derivatives = [MESON_FIELD_DERIVATIVES[label] for label in labels]
    zero = sp.zeros(len(labels))
    result = {}
    for name in LOCAL_CONSTANTS:
        kinetic_matrix = (sp.hessian(sp.expand(kinetic[name]), derivatives)
                          if name in kinetic else zero)
        mass_matrix = (-sp.hessian(sp.expand(non_derivative[name]), fields)
                       if name in non_derivative else zero)
        result[name] = (kinetic_matrix.applyfunc(sp.expand), mass_matrix.applyfunc(sp.expand))
    return result


# ------------------------------------------------------------ pole cancellation
def pole_coefficient_unknowns(field_content="nonet"):
    """Unknown pole coefficients: Gamma_i, and k_B, k_Lambda1, k_Lambda2 in units of m_0^2/F^2.

    B = B^r (1 + k_B m_0^2 lambda/F^2); Lambda_j = Lambda_j^r + k_Lambdaj m_0^2 lambda/F^2.
    In the octet theory (tr Phi = 0) only L_4, ..., L_8 and B can appear.
    """
    names = LOCAL_CONSTANTS[:7] if field_content == "nonet" else LOCAL_CONSTANTS[:5]
    unknowns = {name: sp.Symbol(f"Gamma[{name}]") for name in names} | {"B": sp.Symbol("k_B")}
    if field_content == "nonet":
        unknowns |= {"Lambda_1": sp.Symbol("k_Lambda1"), "Lambda_2": sp.Symbol("k_Lambda2")}
    return unknowns


def field_redefinition_maps(field_content="nonet"):
    """Linear field redefinitions Phi -> Phi + lambda R(Phi) allowed by chiral symmetry.

    R is linear in Phi and contains one chi_0 or one m_0^2, divided by f_pi^2.
    Returns {name: matrix S} with S_ab = 1/2 tr(lambda_a R(lambda_b)), so that
    P -> (1 + lambda sum_k c_k S_k) P. Maps that vanish for the field content
    are dropped.
    """
    labels = NONET_LABELS if field_content == "nonet" else OCTET_LABELS
    chi = chi_zero()
    maps = {"{chi_0,Phi}": lambda matrix: chi * matrix + matrix * chi,
            "<chi_0>Phi": lambda matrix: chi.trace() * matrix,
            "<chi_0 Phi>1": lambda matrix: (chi * matrix).trace() * sp.eye(3),
            "<Phi>chi_0": lambda matrix: matrix.trace() * chi,
            "<chi_0><Phi>1": lambda matrix: chi.trace() * matrix.trace() * sp.eye(3)}
    if field_content == "nonet":
        maps |= {"m_0^2 Phi": lambda matrix: SINGLET_MASS**2 * matrix,
                 "m_0^2<Phi>1": lambda matrix: SINGLET_MASS**2 * matrix.trace() * sp.eye(3)}
    result = {}
    size = len(labels)
    for name, linear_map in maps.items():
        matrix = sp.Matrix(size, size, lambda row, column: sp.expand(
            (GENERATORS[labels[row]] * linear_map(GENERATORS[labels[column]])).trace()
            / (2 * F_PI**2)))
        if not matrix.is_zero_matrix:
            result[name] = matrix
    return result


def field_redefinition_term(field_content, coefficients):
    """Change A^T Q0 + Q0 A of the LO Q0 = s - M^2 under P -> (1 + lambda A) P."""
    data = self_energy(field_content)
    size = len(data["labels"])
    leading = MANDELSTAM_S * sp.eye(size) - data["mass_matrix"]
    redefinition = sp.zeros(size)
    for name, matrix in field_redefinition_maps(field_content).items():
        redefinition += coefficients.get("field:" + name, 0) * matrix
    return redefinition.T * leading + leading * redefinition


def total_pole(field_content, coefficients):
    """Coefficient of lambda in Q(s) = s K - M^2 + Sigma, entry by entry.

    coefficients: dict name -> value for the Gamma_i of L_4..L_25, the
    k_B, k_Lambda1, k_Lambda2 above, and "field:<map>" for the coefficients
    c_k of the field redefinition (missing names count as zero).
    """
    sigma_zero_pole, sigma_one_pole = pole_part(field_content)
    local = local_matrices(field_content)
    data = self_energy(field_content)
    chiral_mass = data["mass_matrix"].subs(SINGLET_MASS, 0)
    singlet_factor = SINGLET_MASS**2 / F_PI**2
    total = sigma_zero_pole + MANDELSTAM_S * sigma_one_pole
    for name in LOCAL_CONSTANTS[:7]:
        kinetic, mass = local[name]
        total += coefficients.get(name, 0) * (MANDELSTAM_S * kinetic - mass)
    for name in ("Lambda_1", "Lambda_2"):
        kinetic, mass = local[name]
        total += coefficients.get(name, 0) * singlet_factor * (MANDELSTAM_S * kinetic - mass)
    total -= coefficients.get("B", 0) * singlet_factor * chiral_mass
    total += field_redefinition_term(field_content, coefficients)
    return total.applyfunc(sp.expand)


def pole_equations(total):
    """Coefficients of all monomials in (m_q, m_0, s, B_0, f_pi) of every entry."""
    variables = (*QUARK_MASSES, SINGLET_MASS, MANDELSTAM_S, B_ZERO, F_PI)
    equations = set()
    for entry in total:
        if entry != 0:
            equations.update(sp.expand(coefficient) for coefficient in
                             sp.Poly(sp.expand(entry * F_PI**2), *variables).coeffs())
    return equations


def solve_pole_coefficients(field_content="nonet", fixed=None, with_field_redefinition=True):
    """Solve total_pole = 0 entry by entry.

    fixed: dict of coefficients held at given values. Returns (values,
    free_symbols, consistent); values maps every unknown name to its
    solution, which may depend on free_symbols.
    """
    fixed = dict(fixed or {})
    unknowns = {name: symbol for name, symbol in pole_coefficient_unknowns(field_content).items()
                if name not in fixed}
    if with_field_redefinition:
        unknowns |= {"field:" + name: sp.Symbol(f"c[{name}]")
                     for name in field_redefinition_maps(field_content)}
    total = total_pole(field_content, fixed | unknowns)
    if not unknowns:
        return {}, set(), all(entry == 0 for entry in total)
    solutions = sp.linsolve(list(pole_equations(total)), list(unknowns.values()))
    if solutions == sp.EmptySet:
        return {}, set(), False
    (solution,) = solutions
    values = dict(zip(unknowns, solution))
    free = set().union(*(value.free_symbols for value in values.values()))
    return values, free, True


# Gamma-combinations that the self-energy fixes, i.e. that do not change under
# the field redefinitions; their number equals (number of Gamma_i) minus the
# number of free field-redefinition parameters (asserted).
INVARIANT_COMBINATIONS = {
    "octet": {"Gamma_6 - Gamma_4/2": {"L_6": 1, "L_4": sp.Rational(-1, 2)},
              "Gamma_8 - Gamma_5/2": {"L_8": 1, "L_5": sp.Rational(-1, 2)},
              "Gamma_7 + Gamma_5/6": {"L_7": 1, "L_5": sp.Rational(1, 6)}},
    "nonet": {"Gamma_6 - Gamma_4/2": {"L_6": 1, "L_4": sp.Rational(-1, 2)},
              "Gamma_8 - Gamma_5/2": {"L_8": 1, "L_5": sp.Rational(-1, 2)},
              "Gamma_25 - Gamma_4/2": {"L_25": 1, "L_4": sp.Rational(-1, 2)},
              "Gamma_18 + Gamma_4 - 2 Gamma_7": {"L_18": 1, "L_4": 1, "L_7": -2}},
}


def gamma_table(field_content="nonet", b_pole=B_POLE_FROM_KAISER_ETA_PRIME_MASS):
    """Table values: U(3) of Kaiser:1998ds (alp-gii.1) or SU(3) of Gasser:1984gg."""
    if field_content == "octet":
        return dict(SU3_GAMMA) | {"B": sp.Integer(0)}
    return dict(KAISER_1998DS_U3_GAMMA) | {"B": b_pole, "Lambda_1": sp.Integer(0),
                                           "Lambda_2": sp.Integer(0)}


@functools.lru_cache(maxsize=None)
def pole_coefficient_report(field_content="nonet"):
    """Derived pole coefficients against the table, in a quotable structure.

    Returns a dict with
      "without_field_redefinition": {name: value} of the Gamma_i and k that
          make Q(s) finite entry by entry with no field redefinition (A = 0),
      "table_without_field_redefinition": True if the table does that,
      "parametric": {name: derived value, possibly depending on the free
          field-redefinition parameters c[...]},
      "free_parameters": sorted names of the free parameters,
      "constants": {name: {"derived", "table", "fixed_by_self_energy",
          "agree"}} for every Gamma_i and k ("agree" for a constant not fixed
          by the self-energy means: the table value is reached by some field
          redefinition consistent with the whole table),
      "invariants": {label: {"derived", "table", "agree"}},
      "table_consistent": True if the full table (with k_B of
          B_POLE_FROM_KAISER_ETA_PRIME_MASS for the nonet) cancels every pole for
          some field redefinition,
      "table_field_redefinition": {c: value} for that redefinition,
      "record_b_consistent": same check with k_B = B_POLE_FROM_RECORD (nonet).
    """
    values, free, consistent = solve_pole_coefficients(field_content)
    if not consistent:
        raise AssertionError(f"poles of the {field_content} theory cannot be absorbed")
    table = gamma_table(field_content)
    table_values, table_free, table_consistent = solve_pole_coefficients(field_content, fixed=table)
    _, _, table_without_redefinition = solve_pole_coefficients(
        field_content, fixed=table, with_field_redefinition=False)
    without_redefinition = {name: value.subs({symbol: 0 for symbol in free})
                            for name, value in values.items() if not name.startswith("field:")}
    record_consistent = True
    if field_content == "nonet":
        _, _, record_consistent = solve_pole_coefficients(
            field_content, fixed=gamma_table(field_content, B_POLE_FROM_RECORD))
    constants = {}
    for name, value in values.items():
        if name.startswith("field:"):
            continue
        fixed_by_self_energy = not value.free_symbols
        agree = (sp.simplify(value - table[name]) == 0) if fixed_by_self_energy else table_consistent
        constants[name] = {"derived": value, "table": table[name],
                           "fixed_by_self_energy": fixed_by_self_energy, "agree": agree}
    invariants = {}
    for label, weights in INVARIANT_COMBINATIONS[field_content].items():
        derived = sp.expand(sum(weight * values[name] for name, weight in weights.items()))
        reference = sum(weight * table[name] for name, weight in weights.items())
        if derived.free_symbols:
            raise AssertionError(f"{label} is not fixed by the self-energy: {derived}")
        invariants[label] = {"derived": derived, "table": reference,
                             "agree": sp.simplify(derived - reference) == 0}
    gamma_names = [name for name in values if name.startswith("L_")]
    gamma_parameters = set().union(*(values[name].free_symbols for name in gamma_names))
    if len(invariants) != len(gamma_names) - len(gamma_parameters):
        raise AssertionError("invariant combinations do not span the fixed Gamma-space")
    return {"without_field_redefinition": without_redefinition,
            "table_without_field_redefinition": table_without_redefinition,
            "parametric": values, "free_parameters": sorted(map(str, free)),
            "constants": constants, "invariants": invariants,
            "table_consistent": table_consistent and not table_free,
            "table_field_redefinition": {name: value for name, value in table_values.items()
                                         if name.startswith("field:")},
            "record_b_consistent": record_consistent}


# ------------------------------------------------------------ numerical evaluation
def numerical_point(point=TEST_POINT):
    """mpmath values of the test-point symbols."""
    return {symbol: mpmath.mpf(sp.N(value, 50)) for symbol, value in point.items()}


def numerical_matrix(matrix, point):
    return mpmath.matrix([[mpmath.mpf(sp.N(matrix[row, column].subs(point), 50))
                           for column in range(matrix.cols)] for row in range(matrix.rows)])


def finite_tadpoles(mass_matrix, scale):
    """<phi phi> and <d phi d phi> with lambda = 0: O diag(m^2 log(m^2/mu^2)/(16 pi^2)) O^T."""
    eigenvalues, rotation = mpmath.eigsy(mass_matrix)
    size = mass_matrix.rows
    tadpole = mpmath.zeros(size)
    derivative_tadpole = mpmath.zeros(size)
    for k in range(size):
        mass_squared = eigenvalues[k]
        value = (mass_squared * mpmath.log(mass_squared / scale**2) / (16 * mpmath.pi**2)
                 if mass_squared > 0 else mpmath.mpf(0))
        for row in range(size):
            for column in range(size):
                product = rotation[row, k] * rotation[column, k]
                tadpole[row, column] += product * value
                derivative_tadpole[row, column] += product * mass_squared * value
    return tadpole, derivative_tadpole


@functools.lru_cache(maxsize=None)
def linear_tadpole_coefficients(field_content="nonet"):
    """Sigma0, Sigma1 as linear maps of the tadpole symbols (for numerical evaluation)."""
    data = self_energy(field_content)
    symbols = sorted(data["tadpole"].free_symbols, key=str)
    maps = []
    for matrix in (data["sigma_zero"], data["sigma_one"]):
        maps.append({symbol: matrix.applyfunc(lambda entry: sp.expand(entry).coeff(symbol))
                     for symbol in symbols})
    return data, maps


def finite_self_energy(field_content, scale, point=TEST_POINT):
    """Finite (lambda = 0) Sigma0, Sigma1 at the given scale, generator basis, mpmath."""
    data, (zero_map, one_map) = linear_tadpole_coefficients(field_content)
    numbers = numerical_point(point)
    mass_matrix = numerical_matrix(data["mass_matrix"], numbers)
    tadpole_values, _ = finite_tadpoles(mass_matrix, mpmath.mpf(sp.N(scale, 50)))
    symbol_values = {}
    for row in range(data["tadpole"].rows):
        for column in range(data["tadpole"].cols):
            entry = data["tadpole"][row, column]
            if isinstance(entry, sp.Symbol):
                symbol_values[entry] = tadpole_values[row, column]
    size = mass_matrix.rows
    results = []
    for linear_map in (zero_map, one_map):
        total = mpmath.zeros(size)
        for symbol, coefficient_matrix in linear_map.items():
            total += numerical_matrix(coefficient_matrix, numbers) * symbol_values[symbol]
        results.append(total)
    return results[0], results[1], mass_matrix


def neutral_block(matrix, labels, neutral=(3, 8, 0)):
    positions = [labels.index(label) for label in neutral]
    return mpmath.matrix([[matrix[row, column] for column in positions] for row in positions])


def renormalized_inverse_propagator(scale, coefficients, point=TEST_POINT):
    """Constant part and s-coefficient of the finite (MS-bar - 1) Q(s) of the nonet.

    L_i^r(mu) = L_i^r(0.77) - Gamma_i log(mu/0.77)/(16 pi^2) with the
    ILLUSTRATIVE_CONSTANTS at 0.77 GeV; Lambda_j^r likewise with k_Lambdaj m_0^2/F^2;
    B^r(mu) = B^r(0.77)(1 - k_B m_0^2 log(mu/0.77)/(16 pi^2 F^2)) in the LO mass
    matrix only (one-loop order: B inside Sigma and inside the NLO local terms
    stays at B^r(0.77)).
    """
    numbers = numerical_point(point)
    running_log = scale_logarithm(scale)
    singlet_factor = numbers[SINGLET_MASS]**2 / numbers[F_PI]**2
    sigma_zero, sigma_one, mass_matrix = finite_self_energy("nonet", scale, point)
    local = local_matrices("nonet")
    data = self_energy("nonet")
    chiral_mass = numerical_matrix(data["mass_matrix"].subs(SINGLET_MASS, 0), numbers)
    kinetic = mpmath.eye(9) + sigma_one
    b_pole = mpmath.mpf(sp.N(coefficients["B"], 50))
    mass = mass_matrix - chiral_mass * b_pole * singlet_factor * running_log - sigma_zero
    for name in LOCAL_CONSTANTS:
        gamma = mpmath.mpf(sp.N(coefficients[name], 50))
        if name.startswith("Lambda"):
            gamma *= singlet_factor
        value = mpmath.mpf(sp.N(ILLUSTRATIVE_CONSTANTS[name], 50)) - gamma * running_log
        kinetic += numerical_matrix(local[name][0], numbers) * value
        mass += numerical_matrix(local[name][1], numbers) * value
    return -mass, kinetic


def scale_logarithm(scale):
    """log(mu/0.77 GeV)/(16 pi^2) = lambda(mu) - lambda(0.77 GeV)."""
    return mpmath.log(mpmath.mpf(sp.N(scale, 50)) / mpmath.mpf(sp.N(REFERENCE_SCALE, 50))) \
        / (16 * mpmath.pi**2)


def one_loop_pole_masses(scale, coefficients, point=TEST_POINT):
    """First-order pole positions s_k = m_k^2 - v_k^T Q(m_k^2) v_k for the nine LO eigenstates.

    v_k, m_k^2 are the LO eigenvectors and masses (generator basis). These
    do not change under a field redefinition (Q0 v_k = 0 at s = m_k^2).
    """
    numbers = numerical_point(point)
    leading = numerical_matrix(self_energy("nonet")["mass_matrix"], numbers)
    eigenvalues, rotation = mpmath.eigsy(leading)
    constant, slope = renormalized_inverse_propagator(scale, coefficients, point)
    masses = []
    for k in range(9):
        vector = rotation[:, k]
        value = eigenvalues[k]
        masses.append(value - (vector.T * (constant + value * slope) * vector)[0, 0])
    return masses


def scale_dependence(coefficients, point=TEST_POINT):
    """mu dependence of the renormalized Q(s) at one loop, per scale.

    Returns {mu: (largest |Delta Q| of the constant part, of the s part,
    largest |Delta Q - log(mu/0.77)/(16 pi^2) (A^T Q0 + Q0 A)| over both parts,
    largest |Delta s_k| of the pole positions)}, Delta = value at mu minus
    value at 0.77 GeV, A the field redefinition in coefficients ("field:" keys).
    """
    numbers = numerical_point(point)
    redefinition = field_redefinition_term("nonet", coefficients)
    redefinition_constant = numerical_matrix(redefinition.subs(MANDELSTAM_S, 0), numbers)
    redefinition_slope = numerical_matrix(redefinition.applyfunc(
        lambda entry: sp.expand(entry).coeff(MANDELSTAM_S)), numbers)
    reference = renormalized_inverse_propagator(REFERENCE_SCALE, coefficients, point)
    reference_masses = one_loop_pole_masses(REFERENCE_SCALE, coefficients, point)
    result = {}
    for scale in SCALES:
        current = renormalized_inverse_propagator(scale, coefficients, point)
        differences = [current[part] - reference[part] for part in range(2)]
        logarithm = scale_logarithm(scale)
        remainder = max(max(abs(x) for x in differences[0] - logarithm * redefinition_constant),
                        max(abs(x) for x in differences[1] - logarithm * redefinition_slope))
        masses = one_loop_pole_masses(scale, coefficients, point)
        result[scale] = (max(abs(x) for x in differences[0]), max(abs(x) for x in differences[1]),
                         remainder, max(abs(a - b) for a, b in zip(masses, reference_masses)))
    return result


# ------------------------------------------------------------ Gasser:1984gg checks
def axial_current_loop_ratio(field_content, label, tadpole):
    """Relative one-loop correction delta_A of <0|tr(lambda_label J_mu)|P_label>.

    J_mu = (f_pi^2/2) i (d U U^dag + U^dag d U), U = exp(i Phi/f_pi), expanded to
    third order; the two undifferentiated fields are contracted with the tadpole
    matrix (<phi d phi> = 0). F_P/F = 1 + delta_A - Sigma1_PP/2.
    """
    labels = NONET_LABELS if field_content == "nonet" else OCTET_LABELS
    loop_fields = sp.symbols(f"x_0:{len(labels)}", real=True)
    counter = sp.Symbol("t_order")
    fluctuation = sum((counter * field * GENERATORS[index] for field, index in zip(loop_fields, labels)),
                      sp.zeros(3))
    derivative = GENERATORS[label] * 1
    unit = sp.I / F_PI

    def power(matrix, exponent):
        return sp.eye(3) if exponent == 0 else matrix**exponent

    def exponential(sign):
        return sum(((sign * unit)**n / sp.factorial(n) * power(fluctuation, n) for n in range(4)),
                   sp.zeros(3))

    def derivative_exponential(sign):
        total = sp.zeros(3)
        for n in range(1, 4):
            total += (sign * unit)**n / sp.factorial(n) * sum(
                (power(fluctuation, k) * derivative * power(fluctuation, n - 1 - k) for k in range(n)),
                sp.zeros(3))
        return total

    current = F_PI**2 / 2 * sp.I * (derivative_exponential(1) * exponential(-1)
                                     + exponential(-1) * derivative_exponential(1))
    projected = sp.expand((GENERATORS[label] * current).trace())
    tree = projected.coeff(counter, 0)
    quadratic = projected.coeff(counter, 2)
    position = {index: number for number, index in enumerate(labels)}
    loop = sp.Integer(0)
    for first, field_first in zip(labels, loop_fields):
        for second, field_second in zip(labels, loop_fields):
            coefficient = sp.diff(quadratic, field_first, field_second) / 2
            loop += coefficient * tadpole[position[first], position[second]]
    return sp.simplify(loop / tree)


def octet_isospin_tadpoles():
    """Symbolic isospin-limit tadpoles t_pi, t_K, t_eta for the octet theory."""
    t_pion, t_kaon, t_eta = sp.symbols("t_pi t_K t_eta", real=True)
    data = self_energy("octet")
    values = sp.zeros(8)
    for label, value in zip(OCTET_LABELS, (t_pion, t_pion, t_pion, t_kaon, t_kaon, t_kaon, t_kaon, t_eta)):
        values[label - 1, label - 1] = value
    return (t_pion, t_kaon, t_eta), tadpole_substitution(data["tadpole"], values), values


def gasser_leutwyler_checks():
    """Return dict of residuals that must vanish (Gasser:1984gg, loop parts, L_i^r = 0).

    mu_P = t_P/(2 F^2), with t_P = M_P^2/(16 pi^2) log(M_P^2/mu^2) the finite tadpole.
    Isospin limit m_u = m_d = m_hat; the pole mass at one loop is
    M^2 - Sigma0 - M^2 Sigma1.
    """
    data = self_energy("octet")
    light = sp.Symbol("m_hat", positive=True)
    isospin = {MASS_UP: light, MASS_DOWN: light}
    (t_pion, t_kaon, t_eta), substitution, tadpole_values = octet_isospin_tadpoles()
    sigma_zero = data["sigma_zero"].subs(substitution).subs(isospin)
    sigma_one = data["sigma_one"].subs(substitution).subs(isospin)
    mass = data["mass_matrix"].subs(isospin)
    mu_pion, mu_kaon, mu_eta = (t / (2 * F_PI**2) for t in (t_pion, t_kaon, t_eta))
    pion_mass, kaon_mass, eta_mass = mass[0, 0], mass[3, 3], mass[7, 7]

    def pole_mass(position):
        return mass[position, position] - sigma_zero[position, position] \
            - mass[position, position] * sigma_one[position, position]

    def decay_constant_ratio(label):
        position = label - 1
        delta = axial_current_loop_ratio("octet", label, tadpole_values).subs(isospin)
        return delta - sigma_one[position, position] / 2

    residuals = {
        "M_pi^2": pole_mass(0) - pion_mass * (1 + mu_pion - mu_eta / 3),
        "M_K^2": pole_mass(3) - kaon_mass * (1 + 2 * mu_eta / 3),
        "M_eta^2": pole_mass(7) - (eta_mass * (1 + 2 * mu_kaon - 4 * mu_eta / 3)
                                   + pion_mass * (-mu_pion + 2 * mu_kaon / 3 + mu_eta / 3)),
        "pi3-eta8 mixing": sigma_zero[2, 7] + sigma_one[2, 7],
        "F_pi/F": decay_constant_ratio(1) - (-2 * mu_pion - mu_kaon),
        "F_K/F": decay_constant_ratio(4) - (-3 * mu_pion / 4 - 3 * mu_kaon / 2 - 3 * mu_eta / 4),
        "F_eta/F": decay_constant_ratio(8) - (-3 * mu_kaon),
    }
    common = sp.Symbol("t", real=True)
    degenerate = {t_pion: common, t_kaon: common, t_eta: common, MASS_STRANGE: light}
    mu_common = common / (2 * F_PI**2)
    residuals["degenerate M^2"] = (pole_mass(0) - pion_mass * (1 + 2 * mu_common / 3)).subs(degenerate)
    residuals["degenerate F/F"] = (decay_constant_ratio(1) + 3 * mu_common).subs(degenerate)
    return {name: sp.simplify(sp.expand(value)) for name, value in residuals.items()}


def chiral_limit_self_energy():
    """Sigma0, Sigma1 of the nonet for m_q = 0: octet tadpoles vanish, eta' has mass m_0."""
    data = self_energy("nonet")
    eta_prime = sp.Symbol("t_eta'", real=True)
    values = sp.zeros(9)
    values[0, 0] = eta_prime
    substitution = tadpole_substitution(data["tadpole"], values)
    chiral = {mass: 0 for mass in QUARK_MASSES}
    return (data["sigma_zero"].subs(substitution).subs(chiral).applyfunc(sp.expand),
            data["sigma_one"].subs(substitution).subs(chiral).applyfunc(sp.expand))


# ------------------------------------------------------------ printing
def entry_name(labels, row, column):
    return f"({labels[row]},{labels[column]})"


def tadpole_latex(symbol):
    """LaTeX name of a tadpole symbol, e.g. <phi3phi8> -> \\langle\\phi_3\\phi_8\\rangle."""
    name = str(symbol)
    if name in CHARGED_TADPOLE_LATEX:
        return CHARGED_TADPOLE_LATEX[name]
    return rf"\langle\phi_{name[4]}\phi_{name[8]}\rangle"


def print_entry(label, expression):
    expression = sp.expand(expression)
    tadpoles = sorted((symbol for symbol in expression.free_symbols if str(symbol).startswith("<")),
                      key=str)
    collected = sp.collect(expression, tadpoles or sorted(expression.free_symbols, key=str))
    names = {symbol: tadpole_latex(symbol) for symbol in tadpoles}
    print(f"  {label}: {sp.latex(collected, symbol_names=names)}")


def run_meson():
    print("Conventions: Q(s) = s K - M^2 + Sigma(s), Sigma = Sigma0 + s Sigma1; GeV units;")
    print("Phi = sum lambda_a P_a, tr(lambda_a lambda_b) = 2 delta_ab; chi_0 = 2 B_0 m_q.")
    print("Tadpole: <phi_k phi_k> = 2 m_k^2 lambda(mu) + m_k^2/(16 pi^2) log(m_k^2/mu^2);")
    print("<d phi_c d phi_d> = (M^2 <phi phi>)_cd; <phi_c phi_d> = sum_k O_ck O_dk <phi_k phi_k>.\n")
    print("Loop Lagrangian (trace formula, l_c = lambda_c; Sigma0, Sigma1 are its Hessians):")
    print(TRACE_FORMULA)
    data = self_energy("nonet")
    labels = data["labels"]
    sigma_zero, sigma_one = data["sigma_zero"], data["sigma_one"]
    assert (sigma_zero - sigma_zero.T).is_zero_matrix, "Sigma0 not symmetric"
    assert (sigma_one - sigma_one.T).is_zero_matrix, "Sigma1 not symmetric"
    print("\n[Sigma, generator basis, tadpoles kept symbolic] (charged: entry (a,a) of the pair)")
    for pair, name in CHARGED_PAIRS.items():
        position = labels.index(pair[0])
        assert sp.expand(sigma_zero[position, position] - sigma_zero[position + 1, position + 1]) == 0
        print_entry(f"Sigma0[{name}]", sigma_zero[position, position])
        print_entry(f"Sigma1[{name}]", sigma_one[position, position])
    positions = [labels.index(label) for label in (3, 8, 0)]
    for first_index, row in enumerate(positions):
        for column in positions[first_index:]:
            print_entry("Sigma0" + entry_name(labels, row, column), sigma_zero[row, column])
            print_entry("Sigma1" + entry_name(labels, row, column), sigma_one[row, column])
    sigma_zero_pole, sigma_one_pole = pole_part("nonet")
    print("\n[Pole part: coefficient of lambda(mu), <phi phi> -> 2 lambda M^2]")
    for pair, name in CHARGED_PAIRS.items():
        position = labels.index(pair[0])
        print_entry(f"Sigma0[{name}]", sp.factor(sigma_zero_pole[position, position]))
        print_entry(f"Sigma1[{name}]", sp.factor(sigma_one_pole[position, position]))
    for first_index, row in enumerate(positions):
        for column in positions[first_index:]:
            print_entry("Sigma0" + entry_name(labels, row, column), sp.factor(sigma_zero_pole[row, column]))
            print_entry("Sigma1" + entry_name(labels, row, column), sp.factor(sigma_one_pole[row, column]))
    print("\n[Finite part]: same expressions with <phi_c phi_d> -> sum_k O_ck O_dk "
          "m_k^2/(16 pi^2) log(m_k^2/mu^2).")
    inputs = [f"{symbol}={sp.N(value, 6)}" for symbol, value in TEST_POINT.items()
              if symbol != SINGLET_MASS]
    inputs.append(f"m_0^2={sp.N(TEST_POINT[SINGLET_MASS]**2, 6)} GeV^2")
    print("Test point (GeV): " + ", ".join(inputs))
    for scale in SCALES:
        sigma_zero_value, sigma_one_value, mass_matrix = finite_self_energy("nonet", scale)
        print(f"\n mu = {sp.N(scale, 3)} GeV, finite part, neutral block (3,8,0):")
        print("  Sigma0 [GeV^2] =", mpmath.nstr(neutral_block(sigma_zero_value, labels), 6))
        print("  Sigma1 =", mpmath.nstr(neutral_block(sigma_one_value, labels), 6))
        for pair, name in CHARGED_PAIRS.items():
            position = labels.index(pair[0])
            print(f"  {name}: Sigma0 = {mpmath.nstr(sigma_zero_value[position, position], 6)} GeV^2,"
                  f" Sigma1 = {mpmath.nstr(sigma_one_value[position, position], 6)}")
        neutral_mass = neutral_block(mass_matrix, labels)
        eigenvalues, rotation = mpmath.eigsy(neutral_mass)
        print("  LO neutral masses [GeV]:", mpmath.nstr(mpmath.matrix([mpmath.sqrt(v) for v in eigenvalues]), 6))
        print("  Sigma0 in LO mass eigenbasis [GeV^2] =",
              mpmath.nstr(rotation.T * neutral_block(sigma_zero_value, labels) * rotation, 6))
        print("  Sigma1 in LO mass eigenbasis =",
              mpmath.nstr(rotation.T * neutral_block(sigma_one_value, labels) * rotation, 6))
    sigma_zero_chiral, sigma_one_chiral = chiral_limit_self_energy()
    assert sigma_zero_chiral.is_zero_matrix and sigma_one_chiral.is_zero_matrix, "chiral limit"
    print("\nChiral limit (m_q = 0, m_0 != 0): Sigma0 = Sigma1 = 0 (asserted).")
    residuals = gasser_leutwyler_checks()
    print("SU(3) octet checks against Gasser:1984gg (loop parts, L_i^r = 0):")
    for name, residual in residuals.items():
        assert residual == 0, f"{name}: residual {residual}"
        print(f"  {name}: residual 0")


def print_report(label, report):
    """Print a pole_coefficient_report."""
    print(f"{label}:")
    print("  coefficients that make Q finite with A = 0:", ", ".join(
        f"{name} = {value}" for name, value in report["without_field_redefinition"].items()))
    print("  table finite with A = 0:", report["table_without_field_redefinition"])
    print("  free field-redefinition parameters:", ", ".join(report["free_parameters"]) or "none")
    for name, row in report["constants"].items():
        status = ("fixed by Sigma" if row["fixed_by_self_energy"]
                  else "not fixed by Sigma alone")
        print(f"  {name}: derived {row['derived']}  | table {row['table']}  | {status},"
              f" {'agree' if row['agree'] else 'DIFFER'}")
    for name, row in report["invariants"].items():
        print(f"  {name}: derived {row['derived']}, table {row['table']},"
              f" {'agree' if row['agree'] else 'DIFFER'}")
    print("  whole table cancels every pole:", report["table_consistent"])
    if report["table_consistent"]:
        print("  field redefinition for the table:", ", ".join(
            f"{name[6:]} = {value}" for name, value in report["table_field_redefinition"].items()))


def print_scale_dependence(label, coefficients, assert_independent):
    tolerance = mpmath.mpf(10)**-30
    for scale, (constant, slope, remainder, pole_shift) in scale_dependence(coefficients).items():
        print(f"  {label}, mu = {sp.N(scale, 3)} GeV: |dQ| constant {mpmath.nstr(constant, 3)} GeV^2,"
              f" s part {mpmath.nstr(slope, 3)}; after field redefinition {mpmath.nstr(remainder, 3)};"
              f" pole positions {mpmath.nstr(pole_shift, 3)} GeV^2")
        if assert_independent:
            assert remainder < tolerance and pole_shift < tolerance, f"mu dependence at {scale}"
        elif remainder > tolerance or pole_shift > tolerance:
            print(f"FINDING: {label}: pole positions depend on mu at one loop"
                  f" (|Delta s| = {mpmath.nstr(pole_shift, 3)} GeV^2 at mu = {sp.N(scale, 3)} GeV)")


def run_cancel():
    print("Pole cancellation: coefficient of lambda in Q(s), entry by entry, with")
    print("L_i -> L_i^r + Gamma_i lambda, B -> B^r (1 + k_B m_0^2 lambda/F^2),")
    print("Lambda_j -> Lambda_j^r + k_Lambdaj m_0^2 lambda/F^2, and the field redefinition")
    print("P -> (1 + lambda A) P, A = sum_k c[k] S_k (Q -> Q + lambda (A^T Q0 + Q0 A)).\n")
    octet = pole_coefficient_report("octet")
    print_report("SU(3) octet theory, table of Gasser:1984gg", octet)
    assert octet["table_consistent"], "octet: Gasser:1984gg Gamma_i do not cancel the poles"
    assert all(row["agree"] for row in octet["invariants"].values()), "octet invariants"
    print()
    nonet = pole_coefficient_report("nonet")
    print_report("U(3) nonet theory, table of alp-gii.1 (Kaiser:1998ds); table k_B from"
                 " 4 tau/F^4 with tau = F^2 m_0^2/6", nonet)
    if not nonet["table_without_field_redefinition"]:
        print("FINDING: with the table Gamma_i (as with the SU(3) Gamma_i of Gasser:1984gg),"
              " Q(s) of the field P is finite entry by entry only after the divergent field"
              " redefinition P -> (1 + lambda A) P listed above (it changes Q off shell but"
              " not the pole positions). Sigma_ij fixes only the four combinations above,"
              " k_B and k_Lambda1; Gamma_4, Gamma_5, Gamma_7 (and the pole of Lambda_2, which"
              " moves with A) need another Green function, e.g. the axial current, i.e. the"
              " c_q^A mixed function of task 3.1.")
    for name, row in nonet["invariants"].items():
        if not row["agree"]:
            print(f"FINDING: {name}: derived {row['derived']} vs table {row['table']}")
    derived_b = nonet["constants"]["B"]["derived"]
    if derived_b != B_POLE_FROM_RECORD:
        print(f"FINDING: B pole: derived B = B^r(1 + {derived_b} m_0^2 lambda/F^2), i.e."
              f" 4 tau/F^4 = {derived_b} m_0^2/F^2 with tau = F^2 m_0^2/6 (M_eta'^2 = 6 tau/F^2,"
              f" Kaiser:1998ds). The record subtraction-scheme-and-scales.md and"
              f" derive_nlo_u3_resonance_counting.tau_in_terms_of_eta_prime_mass use"
              f" tau = F^2 m_0^2/(4 N_f), giving {B_POLE_FROM_RECORD}: a factor 2 in the record,"
              f" not in Kaiser:1998ds. Record value cancels the poles: {nonet['record_b_consistent']}.")
    assert nonet["table_consistent"] or not all(
        row["agree"] for row in nonet["invariants"].values()), "inconsistent report"
    if not nonet["table_consistent"]:
        print("FINDING: the U(3) table does not cancel the poles for any field redefinition")
    table = gamma_table("nonet") | nonet["table_field_redefinition"]
    residual = total_pole("nonet", table)
    if nonet["table_consistent"]:
        assert residual.is_zero_matrix, "table plus field redefinition leaves a pole"
        print("\nFiniteness: with the table Gamma_i, k_B = 2/3 and the field redefinition above,"
              " every entry of K and M^2 is finite (asserted).")
    print("\nmu independence (Delta = value at mu minus value at 0.77 GeV):")
    print_scale_dependence("table, k_B = 2/3", table, assert_independent=nonet["table_consistent"])
    record = gamma_table("nonet", B_POLE_FROM_RECORD) | nonet["table_field_redefinition"]
    print_scale_dependence("table, k_B = 1/3 (record)", record, assert_independent=False)
    return nonet


# ------------------------------------------------------------ tests
class LoopSelfEnergyTests(unittest.TestCase):
    """Checks of the vertex, the contraction, the local terms and the limits."""

    @staticmethod
    def random_matrix(generator, hermitian=False):
        entries = [[sp.Rational(generator.randint(-9, 9), generator.randint(1, 5)) for _ in range(3)]
                   for _ in range(3)]
        matrix = sp.Matrix(entries)
        return (matrix + matrix.T) if hermitian else matrix

    def test_vertex_matches_vertex_script(self):
        from derive_nlo_loop_vertices import coefficient_of_orders, leading_lagrangian
        from derive_nlo_loop_vertices import F_PI as vertex_f_pi
        generator = random.Random(1)
        field, derivative = self.random_matrix(generator), self.random_matrix(generator)
        chi = sp.diag(*[sp.Rational(generator.randint(1, 9), 7) for _ in range(3)])
        fields = {"Phi": field, "partial_mu_Phi": derivative, "chi_0": chi}
        fields |= {name: sp.zeros(3) for name in ("c_q^A", "c_q^V", "c_q^S", "c_q^P")}
        lagrangian = leading_lagrangian(fields, lambda matrix: sp.Matrix(matrix).trace())
        reference = coefficient_of_orders(lagrangian, 4, 0).subs(vertex_f_pi, F_PI)
        self.assertEqual(sp.simplify(reference - quartic_vertex(field, derivative, chi)), 0)

    def test_contraction_matches_wick_expansion(self):
        """Trace formula equals the explicit Wick contraction of the vertex."""
        generator = random.Random(2)
        for field_content in ("nonet", "octet"):
            labels = NONET_LABELS if field_content == "nonet" else OCTET_LABELS
            size = len(labels)
            tadpole = sp.Matrix(size, size, lambda r, c: sp.Symbol(f"T_{min(r, c)}{max(r, c)}"))
            derivative_tadpole = sp.Matrix(size, size, lambda r, c: sp.Symbol(f"U_{min(r, c)}{max(r, c)}"))
            loop_x = sp.symbols(f"x_0:{size}")
            loop_y = sp.symbols(f"y_0:{size}")
            numbers = [sp.Rational(generator.randint(-5, 5), 3) for _ in range(2 * size)]
            field = sum((n * GENERATORS[a] for n, a in zip(numbers[:size], labels)), sp.zeros(3))
            derivative = sum((n * GENERATORS[a] for n, a in zip(numbers[size:], labels)), sp.zeros(3))
            fluctuation = sum((x * GENERATORS[a] for x, a in zip(loop_x, labels)), sp.zeros(3))
            fluctuation_derivative = sum((y * GENERATORS[a] for y, a in zip(loop_y, labels)), sp.zeros(3))
            chi = chi_zero()
            vertex = sp.expand(quartic_vertex(field + fluctuation, derivative + fluctuation_derivative, chi))
            polynomial = sp.Poly(vertex, *loop_x, *loop_y)
            contracted = sp.Integer(0)
            for monomial, coefficient in polynomial.terms():
                if sum(monomial) != 2:
                    continue
                powers_x, powers_y = monomial[:size], monomial[size:]
                indices_x = [i for i, p in enumerate(powers_x) for _ in range(p)]
                indices_y = [i for i, p in enumerate(powers_y) for _ in range(p)]
                if len(indices_x) == 2:
                    contracted += coefficient * tadpole[indices_x[0], indices_x[1]]
                elif len(indices_y) == 2:
                    contracted += coefficient * derivative_tadpole[indices_y[0], indices_y[1]]
            kinetic, non_derivative = loop_lagrangian(labels, tadpole, derivative_tadpole, chi)
            substitution = {MESON_FIELDS[a]: n for n, a in zip(numbers[:size], labels)}
            substitution |= {MESON_FIELD_DERIVATIVES[a]: n for n, a in zip(numbers[size:], labels)}
            self.assertEqual(sp.expand(kinetic.subs(substitution) + non_derivative.subs(substitution)
                                       - contracted), 0, field_content)

    def test_local_neutral_block_matches_record_script(self):
        from derive_nlo_neutral_matrices import derive_matrices
        reference = derive_matrices()
        replace = {reference["B_0"]: B_ZERO, reference["f_pi"]: F_PI, reference["m_u"]: MASS_UP,
                   reference["m_d"]: MASS_DOWN, reference["m_s"]: MASS_STRANGE}
        local = local_matrices("nonet")
        positions = [NONET_LABELS.index(label) for label in (3, 8, 0)]
        for name, (kinetic, mass) in local.items():
            for sector, matrix in (("kinetic", kinetic), ("mass_squared", mass)):
                expected = reference[sector].get(name)
                block = sp.Matrix(3, 3, lambda r, c: matrix[positions[r], positions[c]])
                if expected is None:
                    self.assertTrue(block.is_zero_matrix, (name, sector))
                else:
                    self.assertEqual((block - expected.subs(replace)).applyfunc(sp.simplify),
                                     sp.zeros(3), (name, sector))

    def test_leading_mass_matrix(self):
        mass = leading_mass_matrix(NONET_LABELS)
        chi = [2 * B_ZERO * m for m in QUARK_MASSES]
        self.assertEqual(sp.simplify(mass[1, 1] - (chi[0] + chi[1]) / 2), 0)
        self.assertEqual(sp.simplify(mass[4, 4] - (chi[0] + chi[2]) / 2), 0)
        self.assertEqual(sp.simplify(mass[6, 6] - (chi[1] + chi[2]) / 2), 0)
        self.assertEqual(sp.simplify(mass[0, 0] - sum(chi) / 3 - SINGLET_MASS**2), 0)
        self.assertEqual(mass[1, 2], 0)

    def test_symmetry_and_isospin(self):
        for field_content in ("nonet", "octet"):
            data = self_energy(field_content)
            for name in ("sigma_zero", "sigma_one"):
                matrix = data[name]
                self.assertTrue((matrix - matrix.T).is_zero_matrix)
        sigma_zero_pole, sigma_one_pole = pole_part("nonet")
        isospin = {MASS_DOWN: MASS_UP}
        for matrix in (sigma_zero_pole, sigma_one_pole):
            self.assertEqual(sp.expand(matrix[3, 8].subs(isospin)), 0)
            self.assertEqual(sp.expand(matrix[3, 0].subs(isospin)), 0)

    def test_chiral_limit(self):
        sigma_zero, sigma_one = chiral_limit_self_energy()
        self.assertTrue(sigma_zero.is_zero_matrix)
        self.assertTrue(sigma_one.is_zero_matrix)

    def test_gasser_leutwyler_limits(self):
        for name, residual in gasser_leutwyler_checks().items():
            self.assertEqual(residual, 0, name)

    def test_su3_gamma(self):
        report = pole_coefficient_report("octet")
        self.assertTrue(report["table_consistent"])
        self.assertTrue(all(row["agree"] for row in report["invariants"].values()))
        self.assertEqual(report["constants"]["B"]["derived"], 0)

    def test_nonet_finiteness_and_scale_independence(self):
        report = pole_coefficient_report("nonet")
        self.assertTrue(report["table_consistent"])
        table = gamma_table("nonet") | report["table_field_redefinition"]
        self.assertTrue(total_pole("nonet", table).is_zero_matrix)
        for _, _, remainder, pole_shift in scale_dependence(table).values():
            self.assertLess(remainder, mpmath.mpf(10)**-30)
            self.assertLess(pole_shift, mpmath.mpf(10)**-30)

    def test_field_redefinition_leaves_pole_positions(self):
        """A pure field redefinition changes Q but not the first-order pole positions."""
        coefficients = {name: sp.Integer(0) for name in LOCAL_CONSTANTS} | {"B": sp.Integer(0)}
        shifted = coefficients | {"field:{chi_0,Phi}": sp.Rational(1, 3)}
        point = numerical_point()
        term = field_redefinition_term("nonet", shifted)
        self.assertFalse(term.is_zero_matrix)
        leading = numerical_matrix(self_energy("nonet")["mass_matrix"], point)
        eigenvalues, rotation = mpmath.eigsy(leading)
        for k in range(9):
            value = numerical_matrix(term.subs(MANDELSTAM_S, sp.Float(str(eigenvalues[k]), 40)), point)
            vector = rotation[:, k]
            self.assertLess(abs((vector.T * value * vector)[0, 0]), mpmath.mpf(10)**-25)

    def test_finite_tadpole_matrix_function(self):
        """Eigen route equals the closed form for a diagonal mass matrix."""
        mass = mpmath.diag([mpmath.mpf("0.02"), mpmath.mpf("0.25")])
        tadpole, derivative_tadpole = finite_tadpoles(mass, mpmath.mpf("0.77"))
        for k, value in enumerate((mpmath.mpf("0.02"), mpmath.mpf("0.25"))):
            expected = value * mpmath.log(value / mpmath.mpf("0.77")**2) / (16 * mpmath.pi**2)
            self.assertLess(abs(tadpole[k, k] - expected), mpmath.mpf(10)**-35)
            self.assertLess(abs(derivative_tadpole[k, k] - value * expected), mpmath.mpf(10)**-35)


def run_tests():
    result = unittest.TextTestRunner(verbosity=1).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(LoopSelfEnergyTests))
    return result.wasSuccessful()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--meson", action="store_true", help="task 2.1: Sigma_ij of the nonet")
    parser.add_argument("--cancel", action="store_true", help="task 2.3: pole cancellation")
    parser.add_argument("--test", action="store_true", help="run the unit tests only")
    arguments = parser.parse_args()
    if arguments.test:
        sys.exit(0 if run_tests() else 1)
    no_option = not (arguments.meson or arguments.cancel)
    if arguments.meson or no_option:
        run_meson()
    if arguments.cancel or no_option:
        print()
        run_cancel()
    if no_option and not run_tests():
        sys.exit(1)


if __name__ == "__main__":
    main()
