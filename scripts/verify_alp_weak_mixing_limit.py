#!/usr/bin/env python3
"""Weak-mixing limit f_pi/f_a -> 0 of the exact LO poles and pole vectors (task 2.1).

Run:   .venv/bin/python scripts/verify_alp_weak_mixing_limit.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_verify_alp_weak_mixing_limit.py

Conventions (Alda:2026nrd; memory/alp/shared/finite-alp-meson-mixing.md and
scripts/solve_alp_meson_spectrum.py): fields Phi = (P_i, a) with the ALP last,
L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi, D(s) = s K - M^2, pole vectors
v^(n) with v^(n)T K v^(n) = 1. The expansion parameter is
    epsilon = f_pi / f_a,
varied at fixed f_pi, c_q^A, c_q^P, c_G, chi_0 = 2 B_0 m_q, m_0^2, m_a0^2.
At LO (closed forms of derive_alp_lo_neutral_blocks.py)
    K_PP = 1,          K_Pa = epsilon k,     K_aa = 1 + epsilon^2 k.k,
    M^2_PP (no epsilon), M^2_Pa = epsilon mu, M^2_aa = m_a0^2 + epsilon^2 gamma,
with k_i = -tr(c_q^A lambda_i), mu_i = -tr(lambda_i {c_q^P, chi_0})/2
+ m_0^2 c_G tr(lambda_i)/N_f and gamma = 2 m_0^2 c_G^2/N_f. A constant
orthogonal rotation of the mesons (allowed field change, congruence record)
makes M^2_PP = diag(m_j^2) and keeps this structure, with k, mu rotated.

Part A (symbolic, generic): in that meson mass basis, with symbols m_j^2
(distinct), k_j, mu_j, m_a0^2, gamma, the exact pole conditions are expanded
in epsilon with sympy:
  ALP pole:   S(s) = s K_aa - M^2_aa - sum_j epsilon^2 c_j(s)^2/(s - m_j^2) = 0,
  meson n:    (s - m_n^2) S(s) = 0 near s = m_n^2,
with c_j(s) = s k_j - mu_j. The pole vectors are built from the exact null
vectors of D(s) and normalized with K. Results compared with:
  First order (Alda:2026nrd, no-mixing prescription; see
  memory/alp/baryon/finite-mixing-baryon-couplings.md and
  memory/extracts/Alda:2026nrd/appendix-a-expansion.md):
    m_a^2 = m_a0^2 + O(epsilon^2), m_n^2 = m_n0^2 + O(epsilon^2),
    ALP-pole meson components v_j = epsilon c_j(m_a0^2)/(m_j^2 - m_a0^2) + O(epsilon^3),
    i.e. sum_i w_i lambda_i = 2 epsilon [c_q^A + M_a^{-1}(C)] at s = m_a0^2,
    v_a = 1 + O(epsilon^2) (LSZ factor 1 at first order).
  Second order (derived here; not given in Alda:2026nrd, where mass shifts
  and normalization corrections are beyond first order):
    m_a^2 = m_a0^2 + epsilon^2 [gamma - m_a0^2 k.k + sum_j c_j(m_a0^2)^2/(m_a0^2 - m_j^2)] + O(epsilon^4),
      and the epsilon^4 coefficient is F F' at m_a0^2 with
      F(s) = s k.k - gamma - sum_j c_j(s)^2/(s - m_j^2) (so the epsilon^2 term is -F),
    m_n^2 = m_n0^2 + epsilon^2 c_n(m_n^2)^2/(m_n^2 - m_a0^2) + O(epsilon^4),
    meson-pole ALP component v^(n)_a = -epsilon c_n(m_n^2)/(m_n^2 - m_a0^2) + O(epsilon^3).
  Poles are even in epsilon (a -> -a is a congruence flipping epsilon), so
  the residual after the epsilon^2 term is O(epsilon^4); vector components
  odd in epsilon have residual O(epsilon^3).
Part B (symbolic, all nine U(3) generators, symbolic Hermitian sources):
  the first-order pole equation (M^2_PP - s) w = epsilon (s k - mu) with
  w_i = epsilon tr(lambda_i (c_q^A + Y)) is equivalent to the Alda:2026nrd
  equation M_s(Y) = C, M_s(X) = {chi_0, X}/2 + m_0^2 tr(X) 1/N_f - s X,
  C = -{chi_0, c_q^AP}/2 - m_0^2 c_GA 1/N_f, c_GA = c_G + tr c_q^A,
  c_q^AP = c_q^A - c_q^P.
Part C (numerical): the solver of solve_alp_meson_spectrum.py on the LO 6x6
  blocks (P_3, P_8, P_0, P_6, P_7, a) at decreasing epsilon. The first-order
  ALP-pole components are computed independently from the flavour matrices
  (solving M_{m_a0^2}(Y) = C as a 9x9 linear system). Degenerate meson poles
  at epsilon = 0 (the K0 doublet) are treated as a group: the ALP couples
  to the single combination U^T c, the other combination keeps m_K0^2.
  Residual scaling exponents are fitted on a log-log grid.

Units: f_pi, f_a in GeV; chi_0, m_0^2, m_a0^2, M^2, s in GeV^2; epsilon, K,
c_q^{A,P}, c_G and pole vectors dimensionless. Inputs are illustrative
rationals (derive_alp_lo_neutral_blocks.illustrative_lo_inputs), not fits.
Assumptions: LO chiral Lagrangian plus a canonical bare ALP term, meson
masses at epsilon = 0 distinct from m_a0^2 (no degeneracy with the ALP).
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from derive_alp_charged_meson_decoupling import symbolic_hermitian
from derive_alp_lo_neutral_blocks import (
    FLAVOR_COUNT,
    NEUTRAL_FIELD_LABELS,
    NEUTRAL_GENERATOR_INDICES,
    closed_form_blocks,
    illustrative_lo_inputs,
)
from solve_alp_meson_spectrum import solve_spectrum_numerical
from symbolic_utils import anticommutator, u3_generators

EXPANSION_PARAMETER = sp.Symbol("epsilon", positive=True)
SQUARED_MOMENTUM = sp.Symbol("s")


# ------------------------------------------------------- part A: generic model


def generic_mass_basis_model(meson_count=5):
    """Symbols of the LO structure in the meson mass basis; returns a dict.

    K = [[1, epsilon k], [epsilon k^T, 1 + epsilon^2 k.k]],
    M^2 = [[diag(m_j^2), epsilon mu], [epsilon mu^T, m_a0^2 + epsilon^2 gamma]].
    """
    meson_masses = sp.symbols(f"m_1:{meson_count + 1}_sq", positive=True)
    kinetic_mixing = sp.symbols(f"k_1:{meson_count + 1}", real=True)
    mass_mixing = sp.symbols(f"mu_1:{meson_count + 1}", real=True)
    bare_alp_mass = sp.Symbol("m_a0_sq", positive=True)
    gamma = sp.Symbol("gamma", real=True)
    epsilon = EXPANSION_PARAMETER
    size = meson_count + 1
    kinetic = sp.eye(size)
    mass_squared = sp.diag(*meson_masses, bare_alp_mass + epsilon**2 * gamma)
    for index in range(meson_count):
        kinetic[index, meson_count] = kinetic[meson_count, index] = epsilon * kinetic_mixing[index]
        mass_squared[index, meson_count] = mass_squared[meson_count, index] = epsilon * mass_mixing[index]
    kinetic[meson_count, meson_count] = 1 + epsilon**2 * sum(value**2 for value in kinetic_mixing)
    # d_j stands for 1/(m_a0^2 - m_j^2); used only in the ALP-pole expansion.
    inverse_gaps = sp.symbols(f"d_1:{meson_count + 1}", real=True)
    return {
        "inverse_gaps": inverse_gaps,
        "meson_masses": meson_masses, "kinetic_mixing": kinetic_mixing, "mass_mixing": mass_mixing,
        "bare_alp_mass": bare_alp_mass, "gamma": gamma, "kinetic": kinetic, "mass_squared": mass_squared,
    }


def _coupling(model, index, squared_momentum):
    """c_j(s) = s k_j - mu_j for a single value of s."""
    return squared_momentum * model["kinetic_mixing"][index] - model["mass_mixing"][index]


# Truncated power series in epsilon: lists [a_0, a_1, ..., a_order] of sympy expressions.
# Plain arithmetic on these lists replaces sympy.series, which is too slow here.


def _simplified(value):
    return sp.cancel(value)


def _series_constant(value, order):
    return [value] + [sp.Integer(0)] * order


def _series_add(*terms):
    return [sp.Add(*values) for values in zip(*terms)]


def _series_scale(series, factor):
    return [factor * value for value in series]


def _series_shift(series, power):
    """Multiply by epsilon^power and truncate."""
    order = len(series) - 1
    return [sp.Integer(0)] * power + series[:order + 1 - power]


def _series_multiply(left, right):
    order = len(left) - 1
    return [_simplified(sum(left[index] * right[power - index] for index in range(power + 1))) for power in range(order + 1)]


def _series_inverse(series, inverse_constant=None):
    """1/series; `inverse_constant` may supply a symbol standing for 1/a_0 (keeps polynomials)."""
    order = len(series) - 1
    inverse = [1 / series[0] if inverse_constant is None else inverse_constant]
    for power in range(1, order + 1):
        inverse.append(_simplified(-inverse[0] * sum(series[index] * inverse[power - index] for index in range(1, power + 1))))
    return inverse


def _series_square_root(series):
    """Square root with positive constant term sqrt(a_0)."""
    order = len(series) - 1
    root = [sp.sqrt(series[0])]
    for power in range(1, order + 1):
        cross = sum(root[index] * root[power - index] for index in range(1, power))
        root.append(_simplified((series[power] - cross) / (2 * root[0])))
    return root


def _coupling_series(model, index, momentum_series):
    """c_j(s) = s k_j - mu_j as a series (D_ja = epsilon c_j(s))."""
    return _series_add(_series_scale(momentum_series, model["kinetic_mixing"][index]),
                       _series_constant(-model["mass_mixing"][index], len(momentum_series) - 1))


def _alp_diagonal_series(model, momentum_series):
    """D_aa(s) = s (1 + epsilon^2 k.k) - m_a0^2 - epsilon^2 gamma."""
    order = len(momentum_series) - 1
    kinetic_norm = sum(value**2 for value in model["kinetic_mixing"])
    kinetic_aa = _series_add(_series_constant(1, order), _series_shift(_series_constant(kinetic_norm, order), 2))
    return _series_add(_series_multiply(momentum_series, kinetic_aa),
                       _series_constant(-model["bare_alp_mass"], order),
                       _series_shift(_series_constant(-model["gamma"], order), 2))


def _meson_propagator_series(model, index, momentum_series, near_alp=False):
    """1/(s - m_j^2) as a series; requires s(0) != m_j^2.

    With near_alp=True (s(0) = m_a0^2) the constant 1/(m_a0^2 - m_j^2) is
    written as the symbol d_j, so that all ALP-pole coefficients are
    polynomials (much faster than rational simplification).
    """
    order = len(momentum_series) - 1
    inverse_constant = model["inverse_gaps"][index] if near_alp else None
    return _series_inverse(_series_add(momentum_series, _series_constant(-model["meson_masses"][index], order)),
                           inverse_constant)


def _series_root(condition, start, order):
    """Coefficients [start, r_1, ..., r_order] of the root s(epsilon) of condition(s) = 0.

    `condition` maps a momentum series to a series. The coefficient of
    epsilon^p is linear in r_p once r_1..r_{p-1} are fixed; d condition/ds
    must be nonzero at (start, 0).
    """
    coefficients = [start]
    for power in range(1, order + 1):
        unknown = sp.Symbol("r_unknown")
        trial = coefficients + [unknown] + [sp.Integer(0)] * (order - power)
        value = sp.expand(sp.cancel(condition(trial[:power + 1])[power]))
        coefficients.append(_simplified(sp.solve(value, unknown)[0]))
    return coefficients


def _normalized_vector_series(model, entries):
    """Divide the series entries of a null vector by sqrt(v^T K v)."""
    order = len(entries[0]) - 1
    size = len(entries)
    norm_squared = _series_constant(0, order)
    for row in range(size):
        for column in range(size):
            kinetic_entry = model["kinetic"][row, column]
            if kinetic_entry == 0:
                continue
            kinetic_series = [_simplified(kinetic_entry.diff(EXPANSION_PARAMETER, power).subs(EXPANSION_PARAMETER, 0)
                                          / sp.factorial(power)) for power in range(order + 1)]
            norm_squared = _series_add(norm_squared, _series_multiply(kinetic_series, _series_multiply(entries[row], entries[column])))
    inverse_norm = _series_inverse(_series_square_root(norm_squared))
    return [[_simplified(value) for value in _series_multiply(entry, inverse_norm)] for entry in entries]


def alp_pole_expansion(model, order=4):
    """Series coefficients of the ALP-like pole (to epsilon^order) and of its K-normalized pole vector."""
    meson_count = len(model["meson_masses"])

    def schur_complement(momentum_series):
        value = _alp_diagonal_series(model, momentum_series)
        for index in range(meson_count):
            coupling = _coupling_series(model, index, momentum_series)
            term = _series_multiply(_series_multiply(coupling, coupling), _meson_propagator_series(model, index, momentum_series, near_alp=True))
            value = _series_add(value, _series_scale(_series_shift(term, 2), -1))
        return value

    pole = _series_root(schur_complement, model["bare_alp_mass"], order)
    # Null vector (w, 1): w_j = epsilon c_j(s)/(m_j^2 - s), to epsilon^(order - 1).
    truncated_pole = pole[:order]
    entries = [_series_scale(_series_shift(_series_multiply(_coupling_series(model, index, truncated_pole),
                                                            _meson_propagator_series(model, index, truncated_pole, near_alp=True)), 1), -1)
               for index in range(meson_count)]
    entries.append(_series_constant(1, order - 1))
    vector = _normalized_vector_series(model, entries)
    return {"pole": pole, "vector": vector}


def meson_pole_expansion(model, meson, order=4):
    """Series of the pole that tends to m_meson^2 and of its K-normalized pole vector (sign: v_meson > 0)."""
    meson_count = len(model["meson_masses"])
    others = [index for index in range(meson_count) if index != meson]

    def remaining_alp_diagonal(momentum_series):
        value = _alp_diagonal_series(model, momentum_series)
        for index in others:
            coupling = _coupling_series(model, index, momentum_series)
            term = _series_multiply(_series_multiply(coupling, coupling), _meson_propagator_series(model, index, momentum_series))
            value = _series_add(value, _series_scale(_series_shift(term, 2), -1))
        return value

    def reduced_condition(momentum_series):
        order_here = len(momentum_series) - 1
        own = _series_add(momentum_series, _series_constant(-model["meson_masses"][meson], order_here))
        coupling = _coupling_series(model, meson, momentum_series)
        return _series_add(_series_multiply(own, remaining_alp_diagonal(momentum_series)),
                           _series_scale(_series_shift(_series_multiply(coupling, coupling), 2), -1))

    pole = _series_root(reduced_condition, model["meson_masses"][meson], order)
    truncated_pole = pole[:order]
    # a row: epsilon c_n + y [D_aa - epsilon^2 sum_{j != n} c_j^2/(s - m_j^2)] = 0;
    # meson rows j != n: x_j = -epsilon c_j y/(s - m_j^2); x_n = 1.
    alp_component = _series_scale(_series_shift(_series_multiply(
        _coupling_series(model, meson, truncated_pole), _series_inverse(remaining_alp_diagonal(truncated_pole))), 1), -1)
    entries = []
    for index in range(meson_count):
        if index == meson:
            entries.append(_series_constant(1, order - 1))
        else:
            entries.append(_series_scale(_series_shift(_series_multiply(
                _series_multiply(_coupling_series(model, index, truncated_pole), alp_component),
                _meson_propagator_series(model, index, truncated_pole)), 1), -1))
    entries.append(alp_component)
    vector = _normalized_vector_series(model, entries)
    return {"pole": pole, "vector": vector}


def expected_low_orders(model):
    """First-order (Alda:2026nrd) and second-order (derived) closed forms in the mass basis.

    ALP-pole formulas use d_j = 1/(m_a0^2 - m_j^2).
    """
    masses, bare = model["meson_masses"], model["bare_alp_mass"]
    meson_count = len(masses)
    kinetic_norm = sum(value**2 for value in model["kinetic_mixing"])
    alp_shift = (model["gamma"] - bare * kinetic_norm
                 + sum(_coupling(model, index, bare)**2 * model["inverse_gaps"][index] for index in range(meson_count)))
    # S(s) = (s - m_a0^2) + epsilon^2 F(s), F(s) = s k.k - gamma - sum_j c_j(s)^2/(s - m_j^2);
    # the root is s = m_a0^2 - epsilon^2 F + epsilon^4 F F' + O(epsilon^6) at s = m_a0^2.
    derivative = kinetic_norm - sum(
        2 * _coupling(model, index, bare) * model["kinetic_mixing"][index] * model["inverse_gaps"][index]
        - _coupling(model, index, bare)**2 * model["inverse_gaps"][index]**2 for index in range(meson_count))
    return {
        "alp_pole_second_order": alp_shift,
        "alp_pole_fourth_order": -alp_shift * derivative,
        "alp_vector_first_order": [-_coupling(model, index, bare) * model["inverse_gaps"][index] for index in range(meson_count)],
        "meson_pole_second_order": [_coupling(model, index, masses[index])**2 / (masses[index] - bare)
                                    for index in range(meson_count)],
        "meson_alp_component_first_order": [-_coupling(model, index, masses[index]) / (masses[index] - bare)
                                            for index in range(meson_count)],
    }


def _rational_test_point(model):
    """Exact rational values of all symbols (distinct meson masses, none equal to m_a0^2)."""
    point = {model["bare_alp_mass"]: sp.Rational(3, 10), model["gamma"]: sp.Rational(2, 15)}
    for index, symbol in enumerate(model["meson_masses"]):
        point[symbol] = sp.Rational(1, 50) + sp.Rational(index * 37 + 11, 97)
    for index, symbol in enumerate(model["kinetic_mixing"]):
        point[symbol] = sp.Rational((-1)**index * (index + 2), 7)
    for index, symbol in enumerate(model["mass_mixing"]):
        point[symbol] = sp.Rational(index + 3, 41)
    for index, symbol in enumerate(model["inverse_gaps"]):
        point[symbol] = 1 / (point[model["bare_alp_mass"]] - point[model["meson_masses"][index]])
    return point


def symbolic_weak_mixing_checks(meson_count=3, order=4):
    """Part A: compare the sympy series with the first- and second-order closed forms."""
    model = generic_mass_basis_model(meson_count)
    expected = expected_low_orders(model)
    alp = alp_pole_expansion(model, order)
    test_point = _rational_test_point(model)

    def zero(expression):
        return sp.cancel(sp.together(expression)) == 0

    def nonzero(expression):
        """Nonzero as a function: a nonzero exact value at one rational point suffices."""
        return sp.nsimplify(expression).subs(test_point) != 0
    checks = {
        "ALP pole: no O(epsilon) term (first order of Alda:2026nrd: m_a^2 = m_a0^2)": zero(alp["pole"][1]),
        "ALP pole: epsilon^2 coefficient equals the second-order formula": zero(alp["pole"][2] - expected["alp_pole_second_order"]),
        "ALP pole: no O(epsilon^3) term (poles even in epsilon)": zero(alp["pole"][3]),
        "ALP pole: epsilon^4 coefficient is nonzero (residual scales as epsilon^4)": nonzero(alp["pole"][4]),
        "ALP pole: epsilon^4 coefficient = F(m_a0^2) F'(m_a0^2)": zero(alp["pole"][4] - expected["alp_pole_fourth_order"]),
        "ALP vector: meson components = epsilon c_j(m_a0^2)/(m_j^2 - m_a0^2) + O(epsilon^3)": all(
            zero(alp["vector"][index][0]) and zero(alp["vector"][index][1] - expected["alp_vector_first_order"][index])
            and zero(alp["vector"][index][2]) for index in range(meson_count)),
        "ALP vector: epsilon^3 meson coefficient nonzero (residual scales as epsilon^3)": nonzero(alp["vector"][0][3]),
        "ALP vector: v_a = 1 + O(epsilon^2), epsilon^2 coefficient nonzero": (
            alp["vector"][meson_count][0] == 1 and zero(alp["vector"][meson_count][1])
            and nonzero(alp["vector"][meson_count][2])),
    }
    meson_results = []
    for meson in range(meson_count):
        result = meson_pole_expansion(model, meson, order)
        meson_results.append(result)
        others = [index for index in range(meson_count) if index != meson]
        checks[f"meson {meson + 1}: no O(epsilon) and O(epsilon^3) pole terms"] = zero(result["pole"][1]) and zero(result["pole"][3])
        checks[f"meson {meson + 1}: epsilon^2 pole coefficient equals the second-order formula"] = zero(
            result["pole"][2] - expected["meson_pole_second_order"][meson])
        checks[f"meson {meson + 1}: epsilon^4 pole coefficient nonzero"] = nonzero(result["pole"][4])
        checks[f"meson {meson + 1}: ALP component = -epsilon c_n/(m_n^2 - m_a0^2) + O(epsilon^3)"] = (
            zero(result["vector"][meson_count][0])
            and zero(result["vector"][meson_count][1] - expected["meson_alp_component_first_order"][meson])
            and zero(result["vector"][meson_count][2]) and nonzero(result["vector"][meson_count][3]))
        checks[f"meson {meson + 1}: other meson components start at O(epsilon^2)"] = all(
            zero(result["vector"][index][0]) and zero(result["vector"][index][1]) and nonzero(result["vector"][index][2])
            for index in others)
    return {"checks": checks, "model": model, "alp": alp, "mesons": meson_results, "expected": expected}


# ---------------------------------------------- part B: Alda:2026nrd identity


def alda_first_order_operator(argument, quark_spurion, anomaly_mass_squared, squared_momentum):
    """M_s(X) = {chi_0, X}/2 + m_0^2 tr(X) 1/N_f - s X (Alda:2026nrd's M_a at s = m_a^2)."""
    return (anticommutator(quark_spurion, argument) / 2
            + anomaly_mass_squared * sp.trace(argument) * sp.eye(FLAVOR_COUNT) / FLAVOR_COUNT
            - squared_momentum * argument)


def alda_tadpole_source(axial_source, pseudoscalar_source, quark_spurion, anomaly_mass_squared, gluon_coupling):
    """C = -{chi_0, c_q^AP}/2 - m_0^2 c_GA 1/N_f with c_q^AP = c_q^A - c_q^P, c_GA = c_G + tr c_q^A."""
    axial_minus_pseudoscalar = axial_source - pseudoscalar_source
    return (-anticommutator(quark_spurion, axial_minus_pseudoscalar) / 2
            - anomaly_mass_squared * (gluon_coupling + sp.trace(axial_source)) * sp.eye(FLAVOR_COUNT) / FLAVOR_COUNT)


def symbolic_alda_identity_checks():
    """Part B: (M^2_PP - s) w - epsilon (s k - mu) = epsilon tr(lambda_i (M_s(Y) - C)) for all nine i."""
    generators = u3_generators()
    inputs = {
        "axial_source": symbolic_hermitian("cA", charge_conserving=False),
        "pseudoscalar_source": symbolic_hermitian("cP", charge_conserving=False),
        "quark_spurion": symbolic_hermitian("chi", charge_conserving=False),
        "gluon_coupling": sp.Symbol("c_G", real=True),
        "pion_decay_constant": sp.Symbol("f_pi", positive=True),
        "alp_decay_constant": sp.Symbol("f_a", positive=True),
        "anomaly_mass_squared": sp.Symbol("m_0_sq", real=True),
        "bare_alp_mass_squared": sp.Symbol("m_a0_sq", real=True),
    }
    kinetic, mass_squared = closed_form_blocks(**inputs, generator_indices=tuple(range(9)))
    epsilon = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
    squared_momentum = SQUARED_MOMENTUM
    unknown = symbolic_hermitian("Y", charge_conserving=False)
    overlap = sp.Matrix([epsilon * sp.trace(generator * (inputs["axial_source"] + unknown)) for generator in generators])
    pole_equation = (mass_squared[:9, :9] - squared_momentum * sp.eye(9)) * overlap - (
        squared_momentum * kinetic[:9, 9] - mass_squared[:9, 9])
    operator_value = alda_first_order_operator(unknown, inputs["quark_spurion"], inputs["anomaly_mass_squared"], squared_momentum)
    source = alda_tadpole_source(inputs["axial_source"], inputs["pseudoscalar_source"], inputs["quark_spurion"],
                                 inputs["anomaly_mass_squared"], inputs["gluon_coupling"])
    projected = sp.Matrix([epsilon * sp.trace(generator * (operator_value - source)) for generator in generators])
    difference = (pole_equation - projected).applyfunc(lambda entry: sp.simplify(sp.expand(entry)))
    return {"checks": {
        "(M^2_PP - s) w - (s K_Pa - M^2_Pa) = epsilon tr(lambda_i (M_s(Y) - C)), all nine generators":
            difference == sp.zeros(9, 1),
    }}


# --------------------------------------------------------- part C: numerical


def lo_blocks_at_ratio(ratio, source_structure="neutral", inputs=None):
    """LO 6x6 blocks (floats) with f_a = f_pi / ratio; other inputs fixed."""
    inputs = dict(illustrative_lo_inputs(source_structure) if inputs is None else inputs)
    inputs["alp_decay_constant"] = inputs["pion_decay_constant"] / sp.nsimplify(ratio)
    kinetic, mass_squared = closed_form_blocks(**inputs)
    return np.array(kinetic.evalf(30), dtype=float), np.array(mass_squared.evalf(30), dtype=float)


def alda_first_order_components(inputs, generator_indices=NEUTRAL_GENERATOR_INDICES):
    """tr(lambda_i (c_q^A + Y)) with M_{m_a0^2}(Y) = C, solved as a linear system (complex floats).

    The ALP-pole meson components at first order are epsilon times these.
    Also returns the components on the charged generators (zero for Q-commuting sources).
    """
    generators = u3_generators()
    as_array = lambda matrix: np.array(sp.Matrix(matrix).evalf(30), dtype=complex)  # noqa: E731
    spurion = as_array(inputs["quark_spurion"])
    anomaly = float(inputs["anomaly_mass_squared"])
    bare = float(inputs["bare_alp_mass_squared"])
    source = as_array(alda_tadpole_source(inputs["axial_source"], inputs["pseudoscalar_source"], inputs["quark_spurion"],
                                          inputs["anomaly_mass_squared"], inputs["gluon_coupling"]))
    operator = np.zeros((9, 9), dtype=complex)
    for column in range(9):
        unit = np.zeros(9, dtype=complex)
        unit[column] = 1
        argument = unit.reshape(3, 3)
        image = (spurion @ argument + argument @ spurion) / 2 + anomaly * np.trace(argument) * np.eye(3) / FLAVOR_COUNT - bare * argument
        operator[:, column] = image.reshape(9)
    solution = np.linalg.solve(operator, source.reshape(9)).reshape(3, 3)
    total = as_array(inputs["axial_source"]) + solution
    components = np.array([np.trace(as_array(generator) @ total) for generator in generators])
    if np.max(np.abs(components.imag)) > 1e-12:
        raise ValueError("first-order components are not real")
    components = components.real
    charged = [index for index in range(9) if index not in generator_indices]
    return components[list(generator_indices)], components[charged]


def mass_basis_groups(meson_mass_squared, relative_tolerance=1e-9):
    """Eigenvalues and eigenvectors of M^2_PP grouped by (near-)equal eigenvalues."""
    values, vectors = np.linalg.eigh(meson_mass_squared)
    scale = max(1.0, float(np.max(np.abs(values))))
    groups = [[0]]
    for index in range(1, len(values)):
        if abs(values[index] - values[groups[-1][-1]]) <= relative_tolerance * scale:
            groups[-1].append(index)
        else:
            groups.append([index])
    return [(float(np.mean(values[group])), vectors[:, group]) for group in groups]


def perturbative_prediction(kinetic, mass_squared, ratio, bare_alp_mass_squared):
    """First- and second-order predictions from the LO blocks at f_pi/f_a = ratio (floats).

    Returns the predicted meson poles (each degenerate group expanded with the
    rank-one shift), the ALP pole, and the size of the ALP component summed
    over each meson group.
    """
    meson = slice(0, kinetic.shape[0] - 1)
    alp = kinetic.shape[0] - 1
    kinetic_mixing = kinetic[meson, alp] / ratio
    mass_mixing = mass_squared[meson, alp] / ratio
    gamma = (mass_squared[alp, alp] - bare_alp_mass_squared) / ratio**2
    groups = mass_basis_groups(mass_squared[meson, meson])
    alp_shift = gamma - bare_alp_mass_squared * kinetic_mixing @ kinetic_mixing
    meson_entries = []
    for mass, basis in groups:
        alp_coupling = basis.T @ (bare_alp_mass_squared * kinetic_mixing - mass_mixing)
        alp_shift += float(alp_coupling @ alp_coupling) / (bare_alp_mass_squared - mass)
        own_coupling = basis.T @ (mass * kinetic_mixing - mass_mixing)
        coupling_size = float(np.linalg.norm(own_coupling))
        shifts = [coupling_size**2 / (mass - bare_alp_mass_squared)] + [0.0] * (basis.shape[1] - 1)
        meson_entries.append({
            "unperturbed": mass,
            "first_order_poles": [mass] * basis.shape[1],
            "second_order_poles": [mass + ratio**2 * shift for shift in shifts],
            "alp_component_size": ratio * coupling_size / abs(mass - bare_alp_mass_squared),
            "multiplicity": basis.shape[1],
        })
    return {
        "alp_first_order_pole": bare_alp_mass_squared,
        "alp_second_order_pole": bare_alp_mass_squared + ratio**2 * alp_shift,
        "mesons": meson_entries,
    }


def numerical_residuals(ratio, source_structure="neutral"):
    """Differences between the exact solver and the first/second-order predictions at f_pi/f_a = ratio."""
    inputs = illustrative_lo_inputs(source_structure)
    bare = float(inputs["bare_alp_mass_squared"])
    kinetic, mass_squared = lo_blocks_at_ratio(ratio, source_structure)
    spectrum = solve_spectrum_numerical(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
    poles = np.array(spectrum.pole_squared_masses)
    vectors = spectrum.pole_vectors
    alp = spectrum.size - 1
    prediction = perturbative_prediction(kinetic, mass_squared, ratio, bare)
    alp_pole_index = int(np.argmin(np.abs(poles - bare)))
    first_order_components, _ = alda_first_order_components(inputs)
    alp_vector = vectors[:, alp_pole_index]
    residuals = {
        "ALP pole minus first order": abs(poles[alp_pole_index] - prediction["alp_first_order_pole"]),
        "ALP pole minus second order": abs(poles[alp_pole_index] - prediction["alp_second_order_pole"]),
        "ALP-pole meson components minus first order (Alda:2026nrd)":
            float(np.max(np.abs(alp_vector[:alp] - ratio * first_order_components))),
        "ALP-pole v_a minus 1": abs(alp_vector[alp] - 1),
    }
    remaining = [index for index in range(spectrum.size) if index != alp_pole_index]
    first_order, second_order, component = [], [], []
    for entry in prediction["mesons"]:
        multiplicity = entry["multiplicity"]
        chosen = sorted(remaining, key=lambda index: abs(poles[index] - entry["unperturbed"]))[:multiplicity]
        remaining = [index for index in remaining if index not in chosen]
        exact_values = np.sort(poles[chosen])
        first_order.append(np.max(np.abs(exact_values - np.sort(entry["first_order_poles"]))))
        second_order.append(np.max(np.abs(exact_values - np.sort(entry["second_order_poles"]))))
        component.append(abs(np.sqrt(np.sum(vectors[alp, chosen]**2)) - entry["alp_component_size"]))
    residuals["meson poles minus first order"] = float(np.max(first_order))
    residuals["meson poles minus second order"] = float(np.max(second_order))
    residuals["meson-pole ALP components minus first order"] = float(np.max(component))
    return residuals


EXPECTED_EXPONENTS = {
    "ALP pole minus first order": 2,
    "ALP pole minus second order": 4,
    "ALP-pole meson components minus first order (Alda:2026nrd)": 3,
    "ALP-pole v_a minus 1": 2,
    "meson poles minus first order": 2,
    "meson poles minus second order": 4,
    "meson-pole ALP components minus first order": 3,
}

DEFAULT_RATIOS = tuple(0.02 / 2**power for power in range(5))


def residual_exponents(source_structure="neutral", ratios=DEFAULT_RATIOS):
    """Fitted exponents p of residual ~ epsilon^p (least squares in log-log), and the raw residuals."""
    table = [numerical_residuals(ratio, source_structure) for ratio in ratios]
    exponents = {}
    for label in EXPECTED_EXPONENTS:
        values = np.array([row[label] for row in table])
        exponents[label] = float(np.polyfit(np.log(ratios), np.log(values), 1)[0]) if np.all(values > 0) else None
    return {"ratios": ratios, "residuals": table, "exponents": exponents}


def numerical_weak_mixing_checks(source_structure="neutral", exponent_tolerance=0.15):
    """Part C: fitted exponents agree with EXPECTED_EXPONENTS; charged first-order components vanish."""
    result = residual_exponents(source_structure)
    _, charged = alda_first_order_components(illustrative_lo_inputs(source_structure))
    checks = {f"{label}: exponent {expected}": result["exponents"][label] is not None
              and abs(result["exponents"][label] - expected) <= exponent_tolerance
              for label, expected in EXPECTED_EXPONENTS.items()}
    checks["first-order components on charged generators vanish"] = bool(np.max(np.abs(charged)) < 1e-12)
    return {"checks": checks, **result}


def derive():
    """Run parts A, B and C; returns a dict of named groups."""
    groups = {
        "A symbolic series, 3 mesons": symbolic_weak_mixing_checks(meson_count=3),
        "B Alda:2026nrd first-order identity": symbolic_alda_identity_checks(),
    }
    for source_structure in ("flavor_diagonal", "neutral"):
        groups[f"C numerical, {source_structure} sources"] = numerical_weak_mixing_checks(source_structure)
    return groups


def main():
    groups = derive()
    for name, group in groups.items():
        print(name)
        for label, passed in group["checks"].items():
            print(f"  [{'pass' if passed else 'FAIL'}] {label}")
        if "exponents" in group:
            print(f"  epsilon values: {[float(f'{value:.6g}') for value in group['ratios']]}")
            for label, exponent in group["exponents"].items():
                values = [row[label] for row in group["residuals"]]
                print(f"  {label}: fitted exponent {exponent:.4f}; residuals {[float(f'{value:.3g}') for value in values]}")
    print("ALP pole: m_a^2 = m_a0^2 - epsilon^2 F(m_a0^2) + epsilon^4 F(m_a0^2) F'(m_a0^2) + O(epsilon^6),")
    print("  F(s) = s k.k - gamma - sum_j c_j(s)^2/(s - m_j^2), c_j(s) = s k_j - mu_j (mass basis).")


if __name__ == "__main__":
    main()
