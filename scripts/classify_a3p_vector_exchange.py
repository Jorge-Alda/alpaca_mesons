#!/usr/bin/env python3
"""Classify LO RChT vector-exchange flavor factors in a -> 3P.

Invocation:
    .venv/bin/python scripts/classify_a3p_vector_exchange.py

Natural units are used.  The component convention is
U=exp(i*sqrt(2)*Phi/f_pi), f_pi=93/1000 GeV.  For each allowed three-body
multiset, this script evaluates

    Tr(T_V [T_spectator,C_A]) Tr(T_Vbar [T_pair1,T_pair2])

with explicit U(3) matrices.  Pair ordering follows MESON_ORDER; reversing a
pair reverses the displayed product but cannot change the classification.
"""

from dataclasses import dataclass
from itertools import combinations_with_replacement

import sympy as sp
try:
    from scripts.symbolic_utils import matrix_unit
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import matrix_unit


def commutator(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    """Return the explicit SymPy matrix commutator."""
    return left * right - right * left


def simplified(expression: sp.Expr) -> sp.Expr:
    """Put trace products in a deterministic, condition-visible form."""
    return sp.factor(sp.trigsimp(sp.expand(expression)))


def show(label: str, expression: sp.Basic) -> None:
    """Print one analytical input with SymPy's LaTeX renderer."""
    print(f"{label}: ${sp.latex(expression)}$")


def E(row: int, column: int) -> sp.Matrix:
    """Return the three-flavor matrix unit E_(row,column)."""
    return sp.Matrix(matrix_unit(row, column))


@dataclass(frozen=True)
class State:
    generator: sp.Matrix
    charge: int
    strangeness: int
    tex: str


@dataclass(frozen=True)
class Contribution:
    spectator: str
    spectator_multiplicity: int
    pair: tuple[str, str]
    vector: str
    first_trace: sp.Expr
    second_trace: sp.Expr
    product: sp.Expr


# ---------------------------------------------------------------------------
# 1. RChT component conventions and exact symbolic inputs
# ---------------------------------------------------------------------------
GeV = sp.Symbol("GeV", positive=True)
f_pi = sp.Rational(93, 1000) * GeV
theta_P = sp.Symbol("theta_P", real=True)
c_u, c_d, c_s, c_sd = sp.symbols("c_u c_d c_s c_sd")

pi0 = sp.diag(1, -1, 0) / sp.sqrt(2)
eta8 = sp.diag(1, 1, -2) / sp.sqrt(6)
eta0 = sp.eye(3) / sp.sqrt(3)
eta = sp.cos(theta_P) * eta8 - sp.sin(theta_P) * eta0
eta_prime = sp.sin(theta_P) * eta8 + sp.cos(theta_P) * eta0

PSEUDOSCALARS = {
    "pi+": State(E(0, 1), +1, 0, r"\pi^+"),
    "pi-": State(E(1, 0), -1, 0, r"\pi^-"),
    "pi0": State(pi0, 0, 0, r"\pi^0"),
    "K+": State(E(0, 2), +1, +1, r"K^+"),
    "K-": State(E(2, 0), -1, -1, r"K^-"),
    "K0": State(E(1, 2), 0, +1, r"K^0"),
    "K0bar": State(E(2, 1), 0, -1, r"\bar{K}^0"),
    "eta": State(eta, 0, 0, r"\eta"),
    "eta_prime": State(eta_prime, 0, 0, r"\eta'"),
}
MESON_ORDER = tuple(PSEUDOSCALARS)
MESON_INDEX = {name: index for index, name in enumerate(MESON_ORDER)}

rho0 = pi0
omega_ideal = sp.diag(1, 1, 0) / sp.sqrt(2)
phi_ideal = sp.diag(0, 0, 1)
VECTORS = {
    "rho+": State(E(0, 1), +1, 0, r"\rho^+"),
    "rho-": State(E(1, 0), -1, 0, r"\rho^-"),
    "rho0": State(rho0, 0, 0, r"\rho^0"),
    "K*+": State(E(0, 2), +1, +1, r"K^{*+}"),
    "K*-": State(E(2, 0), -1, -1, r"K^{*-}"),
    "K*0": State(E(1, 2), 0, +1, r"K^{*0}"),
    "K*0bar": State(E(2, 1), 0, -1, r"\bar{K}^{*0}"),
    "omega": State(omega_ideal, 0, 0, r"\omega"),
    "phi": State(phi_ideal, 0, 0, r"\phi"),
}
VECTOR_CONJUGATE = {
    "rho+": "rho-",
    "rho-": "rho+",
    "rho0": "rho0",
    "K*+": "K*-",
    "K*-": "K*+",
    "K*0": "K*0bar",
    "K*0bar": "K*0",
    "omega": "omega",
    "phi": "phi",
}

C_DIAGONAL = sp.diag(c_u, c_d, c_s)
C_SD = c_sd * E(2, 1)  # E_sd has S=-1 and selects final-state S=+1.
SCENARIOS = {
    "diag C_A": (C_DIAGONAL, 0, 0),
    "C_A=c_sd E_sd": (C_SD, 0, -1),
}

Q_FLAVOR = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
S_FLAVOR = sp.diag(0, 0, -1)

show("f_pi [GeV]", f_pi / GeV)
show("eta generator", eta)
show("eta-prime generator", eta_prime)
show("generic diagonal C_A", C_DIAGONAL)
show("Delta-S spurion C_A=c_sd E_sd", C_SD)
print("Vector convention: ideal omega=diag(1,1,0)/sqrt(2), phi=diag(0,0,1).")
print("Pair convention: pair1 precedes pair2 in MESON_ORDER.\n")


# ---------------------------------------------------------------------------
# 2. Enumerate all charge-neutral S=0,+1 multisets and all trace products
# ---------------------------------------------------------------------------
CHANNELS = tuple(
    channel
    for channel in combinations_with_replacement(MESON_ORDER, 3)
    if sum(PSEUDOSCALARS[name].charge for name in channel) == 0
    and sum(PSEUDOSCALARS[name].strangeness for name in channel) in (0, +1)
)


def pair_after_spectator(channel: tuple[str, str, str], spectator: str) -> tuple[str, str]:
    """Remove one occurrence of spectator and return the ordered residual pair."""
    residual = list(channel)
    residual.remove(spectator)
    return tuple(sorted(residual, key=MESON_INDEX.__getitem__))


def classify(channel: tuple[str, str, str], C_A: sp.Matrix) -> tuple[Contribution, ...]:
    """Calculate every spectator/vector product and retain symbolic nonzeros."""
    nonzero = []
    for spectator in dict.fromkeys(channel):
        pair = pair_after_spectator(channel, spectator)
        T_spectator = PSEUDOSCALARS[spectator].generator
        T_pair1 = PSEUDOSCALARS[pair[0]].generator
        T_pair2 = PSEUDOSCALARS[pair[1]].generator
        for vector, vector_state in VECTORS.items():
            conjugate = VECTOR_CONJUGATE[vector]
            T_V = vector_state.generator
            T_Vbar = VECTORS[conjugate].generator
            first_trace = simplified(sp.trace(T_V * commutator(T_spectator, C_A)))
            second_trace = simplified(
                sp.trace(T_Vbar * commutator(T_pair1, T_pair2))
            )
            product = simplified(first_trace * second_trace)

            # Trace-product assertion: the classified expression is obtained
            # directly from the two requested vertices, not from a lookup rule.
            raw_product = sp.trace(T_V * commutator(T_spectator, C_A)) * sp.trace(
                T_Vbar * commutator(T_pair1, T_pair2)
            )
            assert simplified(product - raw_product) == 0
            if product != 0:
                nonzero.append(
                    Contribution(
                        spectator,
                        channel.count(spectator),
                        pair,
                        vector,
                        first_trace,
                        second_trace,
                        product,
                    )
                )
    return tuple(nonzero)


all_results: dict[tuple[str, tuple[str, str, str]], tuple[Contribution, ...]] = {}
for scenario, (C_A, spurion_charge, spurion_strangeness) in SCENARIOS.items():
    for channel in CHANNELS:
        channel_charge = sum(PSEUDOSCALARS[name].charge for name in channel)
        channel_strangeness = sum(PSEUDOSCALARS[name].strangeness for name in channel)
        assert channel_charge == 0
        assert channel_strangeness in (0, +1)

        contributions = classify(channel, C_A)
        all_results[(scenario, channel)] = contributions
        for contribution in contributions:
            spectator_state = PSEUDOSCALARS[contribution.spectator]
            pair_states = [PSEUDOSCALARS[name] for name in contribution.pair]
            vector_state = VECTORS[contribution.vector]

            # Flavor-eigenvalue checks for both traces.  The conjugate vector
            # carries the opposite charge and strangeness at the VPP vertex.
            assert vector_state.charge + spectator_state.charge + spurion_charge == 0
            assert (
                vector_state.strangeness
                + spectator_state.strangeness
                + spurion_strangeness
                == 0
            )
            assert -vector_state.charge + sum(state.charge for state in pair_states) == 0
            assert (
                -vector_state.strangeness
                + sum(state.strangeness for state in pair_states)
                == 0
            )


# ---------------------------------------------------------------------------
# 3. Concise complete channel table (zeros included as '--')
# ---------------------------------------------------------------------------
print(f"Enumerated channels: {len(CHANNELS)} (charge Q=0, total S=0 or +1).")
for scenario, (_, _, spurion_strangeness) in SCENARIOS.items():
    print(f"\n[{scenario}]  required nonzero final S={-spurion_strangeness}")
    print("channel | nonzero spectator -> Vbar-pair : trace product\n-----------|------------")
    for channel in CHANNELS:
        channel_strangeness = sum(PSEUDOSCALARS[name].strangeness for name in channel)
        if channel_strangeness != -spurion_strangeness:
            continue
        contributions = all_results[(scenario, channel)]
        entries = []
        for item in contributions:
            pair_text = "\\,".join(PSEUDOSCALARS[name].tex for name in item.pair)
            entries.append(
                f"$a \\to {PSEUDOSCALARS[item.spectator].tex} "
                f"{'(\\times' + str(item.spectator_multiplicity) + ')' if item.spectator_multiplicity > 1 else ''}"
                f"{VECTORS[item.vector].tex}[\\to {pair_text}]$ : "
                f"${sp.latex(item.product)}$"
            )
        print(f"${' '.join([PSEUDOSCALARS[name].tex for name in channel])}$ | " + ("; ".join(entries) or "--"))


# Global selection and conjugation checks.
for (scenario, channel), contributions in all_results.items():
    final_strangeness = sum(PSEUDOSCALARS[name].strangeness for name in channel)
    required_strangeness = -SCENARIOS[scenario][2]
    assert not contributions or final_strangeness == required_strangeness

for vector, conjugate in VECTOR_CONJUGATE.items():
    assert VECTORS[vector].generator.T == VECTORS[conjugate].generator
    assert VECTOR_CONJUGATE[conjugate] == vector

print("\nAssertions passed: channel Q/S, both-vertex flavor flow, conjugate-vector maps, and all trace products.")
print("Approximation: ideal omega-phi mixing. Exact symbolic inputs: cu, cd, cs, c_sd and theta_P.")
print("S=-1 classification is obtained by Hermitian conjugating fields and E_sd -> E_ds.")
