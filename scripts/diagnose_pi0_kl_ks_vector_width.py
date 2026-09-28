#!/usr/bin/env python3
"""Diagnose the pi0 K_L K_S vector-width cancellation slowdown.

Invocation:
    .venv/bin/python scripts/diagnose_pi0_kl_ks_vector_width.py

This diagnostic deliberately imports the current implementation modules as
top-level modules because they currently use top-level imports internally.  It
does not call the unbounded nested production width integral.
"""

from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
import sympy as sp


# Preserve the current implementation's top-level import arrangement without
# importing scripts.implementation (whose package initializer is unrelated to
# this diagnostic).
IMPLEMENTATION = Path(__file__).resolve().parent / "implementation"
sys.path.insert(0, str(IMPLEMENTATION))

from three_meson_width import (  # noqa: E402
    EPSILON_K,
    _kaon_mixing_coefficients,
    _s23_bounds,
    build_amplitude,
)
from vector_meson_exchange import (  # noqa: E402
    build_vector_exchange_amplitude,
    build_vector_matching_amplitude,
)


def show(label: str, expression: sp.Basic) -> None:
    """Print one symbolic derivation step using SymPy's LaTeX printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def one_vector_term(model, contribution, s12: float, s23: float) -> complex:
    """Evaluate one retained term through the implementation's own callable."""
    return replace(model, contributions=(contribution,))(s12, s23)


def print_vector_terms(label: str, model, s12: float, s23: float) -> float:
    """Print retained metadata and terms; return their absolute-value sum."""
    print(f"{label}: {len(model.contributions)} retained contributions")
    terms = []
    for index, contribution in enumerate(model.contributions, start=1):
        term = one_vector_term(model, contribution, s12, s23)
        terms.append(term)
        print(
            f"  {index}: exchange={contribution.exchange_name}, "
            f"spectator={contribution.spectator_index + 1}, "
            f"pair=({contribution.pair_indices[0] + 1},"
            f"{contribution.pair_indices[1] + 1}), "
            f"flavor_factor={contribution.flavor_factor:.17e}, "
            f"term={term:.17e}"
        )
    term_sum = sum(terms, 0.0j)
    absolute_sum = sum(abs(term) for term in terms)
    print(f"  component sum       = {term_sum:.17e}")
    print(f"  absolute-term sum   = {absolute_sum:.17e}")
    print(f"  callable residual   = {(model(s12, s23) - term_sum):.17e}")
    return absolute_sum


def build_components(epsilon_k: complex, channel: tuple[str, str, str]):
    """Build current ChPT, vector-pole, and vector-matching amplitudes."""
    common = (M_A, C_Q_A, C_Q_P, *channel, F_A)
    return (
        build_amplitude(*common, epsilon_k=epsilon_k),
        build_vector_exchange_amplitude(*common, epsilon_k=epsilon_k),
        build_vector_matching_amplitude(*common, epsilon_k=epsilon_k),
    )


def bounded_inner_quad(
    models, s12: float, epsabs: float, limit: int = 24
) -> tuple[float, float, int, int, bool, str]:
    """Run one bounded inner Dalitz integral and expose QUADPACK status."""
    lower, upper = _s23_bounds(
        s12, models[0].alp_mass**2, *models[0].masses
    )

    def integrand(s23: float) -> float:
        return abs(sum(model(s12, s23) for model in models)) ** 2

    result = quad(
        integrand,
        lower,
        upper,
        epsabs=epsabs,
        epsrel=1.0e-6,
        limit=limit,
        full_output=1,
    )
    value, error, info = result[:3]
    converged = len(result) == 3
    message = "converged" if converged else str(result[3]).splitlines()[0]
    return value, error, info["neval"], info["last"], converged, message


# Representative assumptions, not a phenomenological benchmark.  Natural
# units are used, so masses, f_a, and integration invariants are in GeV,
# GeV, and GeV^2 respectively; invariant amplitudes are dimensionless.
M_A = 1.6
F_A = 100.0
C_Q_A = np.diag([0.8, -0.3, 0.15])
C_Q_P = np.zeros((3, 3))
PHYSICAL_CHANNEL = ("pi0", "KL", "KS")
FLAVOR_CHANNEL = ("pi0", "K0", "K0bar")
REVERSED_FLAVOR_CHANNEL = ("pi0", "K0bar", "K0")


print("Conventions and assumptions:")
print("- Natural units; amplitudes are dimensionless and Dalitz invariants are GeV^2.")
print("- K_S=p K0+q K0bar and K_L=p K0-q K0bar; final-state bras use p* and q*.")
print("- Only the strangeness-zero K0 K0bar and K0bar K0 flavor orderings are retained.")
print("- The flavor amplitude is Bose/order symmetric under interchange of the two equal-mass kaon legs.")
print("- Representative inputs: m_a=1.6 GeV, f_a=100 GeV, c_q^A=diag(0.8,-0.3,0.15), c_q^P=0.")
print("- Current default pole data, constant widths, ideal vector mixing, and G_V=f_pi/sqrt(2) are used.")
print("- Quadrature evidence is one bounded inner integral (limit=24), not the nested production width.\n")


# ---------------------------------------------------------------------------
# 1. Exact SymPy projection cancellation
# ---------------------------------------------------------------------------
p_bar, q_bar = sp.symbols(r"p^{*} q^{*}")
A_0bar, A_bar0 = sp.symbols(
    r"\mathcal{A}_{0\bar{0}} \mathcal{A}_{\bar{0}0}"
)
A_sym = sp.symbols(r"\mathcal{A}_{\mathrm{sym}}")

kl_bra = sp.Matrix([[p_bar, -q_bar]])
ks_bra = sp.Matrix([[p_bar, q_bar]])
ordered_coefficients = sp.Matrix(
    [kl_bra[0, 0] * ks_bra[0, 1], kl_bra[0, 1] * ks_bra[0, 0]]
)
ordered_amplitudes = sp.Matrix([A_0bar, A_bar0])
projected = sp.expand((ordered_coefficients.T * ordered_amplitudes)[0])
symmetry_substitution = {A_0bar: A_sym, A_bar0: A_sym}
projected_symmetric = sp.simplify(projected.subs(symmetry_substitution))

show("K_L final-state bra coefficients in the (K0,K0bar) basis", kl_bra)
show("K_S final-state bra coefficients in the (K0,K0bar) basis", ks_bra)
show("Coefficients of the two allowed ordered flavor amplitudes", ordered_coefficients)
show("Projected pi0 K_L K_S flavor amplitude", projected)
show("Bose/order-symmetry substitution", sp.Eq(A_0bar, A_bar0))
show("Exact projected amplitude after the symmetry substitution", projected_symmetric)
assert projected_symmetric == 0


# ---------------------------------------------------------------------------
# 2. Representative corrected-amplitude evaluation
# ---------------------------------------------------------------------------
default_models = build_components(EPSILON_K, PHYSICAL_CHANNEL)
flavor_models = build_components(EPSILON_K, FLAVOR_CHANNEL)
reversed_flavor_models = build_components(EPSILON_K, REVERSED_FLAVOR_CHANNEL)
assert default_models[1].contributions == ()
assert default_models[2].contributions == ()
assert default_models[1](0.0, 0.0) == 0.0j
assert default_models[2](0.0, 0.0) == 0.0j

s12_lower = (default_models[0].masses[0] + default_models[0].masses[1]) ** 2
s12_upper = (M_A - default_models[0].masses[2]) ** 2
S12 = 0.4 * s12_lower + 0.6 * s12_upper
s23_lower, s23_upper = _s23_bounds(
    S12, M_A**2, *default_models[0].masses
)
S23 = 0.3 * s23_lower + 0.7 * s23_upper

print("Representative Dalitz point:")
print(f"  s12={S12:.17e} GeV^2")
print(f"  s23={S23:.17e} GeV^2\n")

p_value, q_value = _kaon_mixing_coefficients(EPSILON_K)
projection_coefficient = np.conjugate(p_value) * np.conjugate(q_value)

print("Retained flavor-projection contributions (default complex epsilon_K):")
for label, model_0bar, model_bar0 in zip(
    ("ChPT", "vector pole", "vector matching"),
    flavor_models,
    reversed_flavor_models,
):
    terms = (
        projection_coefficient * model_0bar(S12, S23),
        -projection_coefficient * model_bar0(S12, S23),
    )
    print(f"  {label} K0 K0bar term  = {terms[0]:.17e}")
    print(f"  {label} K0bar K0 term  = {terms[1]:.17e}")
    print(f"  {label} projected sum  = {sum(terms):.17e}")
    print(f"  {label} absolute sum   = {sum(abs(term) for term in terms):.17e}")
print()

physical_absolute_sums = [
    2.0
    * abs(projection_coefficient)
    * abs(flavor_models[0](S12, S23))
]
physical_absolute_sums.append(
    print_vector_terms("pi0 KL KS vector pole", default_models[1], S12, S23)
)
physical_absolute_sums.append(
    print_vector_terms("pi0 KL KS vector matching", default_models[2], S12, S23)
)
print()

flavor_absolute_sums = [abs(flavor_models[0](S12, S23))]
flavor_absolute_sums.append(
    print_vector_terms("pi0 K0 K0bar vector pole", flavor_models[1], S12, S23)
)
flavor_absolute_sums.append(
    print_vector_terms("pi0 K0 K0bar vector matching", flavor_models[2], S12, S23)
)
print()

print("Representative-point sector sums:")
for channel_label, models, absolute_sums in (
    ("pi0 KL KS", default_models, physical_absolute_sums),
    ("pi0 K0 K0bar", flavor_models, flavor_absolute_sums),
):
    values = tuple(model(S12, S23) for model in models)
    total = sum(values)
    absolute_total = sum(absolute_sums)
    print(f"  {channel_label} ChPT            = {values[0]:.17e}")
    print(f"  {channel_label} vector pole     = {values[1]:.17e}")
    print(f"  {channel_label} vector matching = {values[2]:.17e}")
    print(f"  {channel_label} component sum    = {total:.17e}")
    print(f"  {channel_label} absolute-term sum= {absolute_total:.17e}")
    print(f"  {channel_label} residual ratio   = {abs(total) / absolute_total:.17e}")
    if channel_label == "pi0 KL KS":
        assert abs(total) <= 512.0 * np.finfo(float).eps * absolute_total
    else:
        assert absolute_total > 0.0 and abs(total) > 1.0e-6
print()


# ---------------------------------------------------------------------------
# 3. Bounded evidence for epsilon_K and epsabs behavior
# ---------------------------------------------------------------------------
cp_models = build_components(0.0j, PHYSICAL_CHANNEL)
cp_values = tuple(model(S12, S23) for model in cp_models)
assert all(value == 0.0j for value in cp_values)

inner_span = s23_upper - s23_lower
flavor_scale = inner_span * sum(
    abs(model(S12, 0.5 * (s23_lower + s23_upper))) for model in flavor_models
) ** 2
# Diagnostic tolerance assumption: resolving 1e-12 of the corresponding
# noncancelled flavor-channel integral scale is more than sufficient to
# classify an analytically zero physical channel.
scale_aware_epsabs = 1.0e-12 * flavor_scale

print("Bounded inner-quadrature diagnostic at fixed s12:")
print(f"  noncancelled flavor integral scale = {flavor_scale:.17e} GeV^2")
print(f"  scale-aware epsabs                 = {scale_aware_epsabs:.17e} GeV^2")
for epsilon_label, models in (
    ("default complex epsilon_K", default_models),
    ("epsilon_K=0", cp_models),
):
    for epsabs in (0.0, scale_aware_epsabs):
        value, error, neval, last, converged, message = bounded_inner_quad(
            models, S12, epsabs
        )
        print(
            f"  {epsilon_label}, epsabs={epsabs:.3e}: "
            f"integral={value:.17e}, error={error:.3e}, neval={neval}, "
            f"subintervals={last}/24, converged={converged}, status={message}"
        )

relative_only = bounded_inner_quad(default_models, S12, 0.0)
absolute_bounded = bounded_inner_quad(default_models, S12, scale_aware_epsabs)
cp_relative_only = bounded_inner_quad(cp_models, S12, 0.0)
assert relative_only[4] and relative_only[3] == 1
assert absolute_bounded[4] and absolute_bounded[3] == 1
assert cp_relative_only[4] and cp_relative_only[0] == 0.0

print("\nAssertions passed: exact symbolic cancellation, pre-kinematic vector pruning, retained nonzero flavor amplitudes, and one-panel bounded quadrature.")
